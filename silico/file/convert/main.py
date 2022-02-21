# General imports.
import yaml
from pathlib import Path
from openbabel import pybel

# Silico imports.
from silico.exception.base import Silico_exception
from silico.file.convert.gaussian import Gaussian_input_parser
from silico.file.convert.babel import Openbabel_converter
import silico.logging

# Custom formats to allow literal strings in yaml output.
# Adapted from https://stackoverflow.com/questions/6432605/any-yaml-libraries-in-python-that-support-dumping-of-long-strings-as-block-liter

# Custom 'str' class that will be dumped literally.
class literal_str(str): pass
# The dumper (which recognises the literal_str class).
def literal_str_representer(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')
# Add to yaml.
yaml.add_representer(literal_str, literal_str_representer)


class Silico_coords():
    """
    Class that represents an input file in the silico input (.si) format.
     
    The .si format is yaml based, and stores the input geometry in a modified xyz format along with charge and multiplicity.
    """
    
    def __init__(self, geometry, *, charge = None, multiplicity = None, name = None, file_name = None):
        """
        Constructor for .si files.
         
        :param geometry: The molecular geometry in .si format. Use from_xyz() instead if your format is in xyz.
        :param charge: The molecular charge.
        :param multiplicity: The molecular multiplicity (as an integer).
        :param name: Name of the system/molecule.
        :param file_name: The name of the file which was loaded. This can be used as a back-up file name.
        """
        self.geometry = geometry
        self.charge = charge
        self.multiplicity = multiplicity
        self.file_name = file_name
        self.name = name
        
    
    @property
    def formula(self):
        """
        """
        molecule = pybel.readstring("xyz", self.to_format("xyz"))
        return molecule.formula
    
    @property
    def elements(self):
        """
        A unique list of the elements in this input file.
        """
        return list(set([atom.atomicnum for atom in pybel.readstring("xyz", self.to_format("xyz")).atoms]))
    
    @property
    def auto_name(self):
        """
        A more intelligent name for this molecule, taking into account our old file name if necessary.
        """
        # If a real name wasn't given, but a file name was, use it.
        if self.name is None and self.file_name is not None:
            name = Path(self.file_name).with_suffix("").name
            
        elif self.name is not None:
            name = self.name
        
        else:
            name = "molecule"
            
        return name

    @classmethod
    def from_xyz(self, geometry, **kwargs):
        """
        Create a Silico_coords object from a molecule in xyz format.
        
        :param geometry: The input geometry in xyz format.
        :param charge: The molecular charge.
        :param multiplicity: The molecular multiplicity (as an integer).
        :param name: Name of the system/molecule
        """
        # Split the xyz format on newlines.
        split_geom = geometry.split("\n")
        
        # Remove the first line (the number of atoms) and the second (optional comment line)
        split_geom = split_geom[2:]
        
        # Call our main constructor.
        return self("\n".join(split_geom), **kwargs)
    
    @classmethod
    def from_com(self, geometry, *, charge = None, multiplicity = None, **kwargs):
        """
        Create a Silico_coords object from a molecule in gaussian input format (.com, .gjf etc).
        
        :param geometry: The input geometry in gaussian format.
        :param charge: The molecular charge.
        :param multiplicity: The molecular multiplicity (as an integer).
        :param name: Name of the system/molecule
        """
        # Get a parser for our gaussian input file.
        parser = Gaussian_input_parser(geometry)
        
        # If we've not been given a charge and/or multi, use the values from the input file.
        if charge is None:
            charge = parser.charge
        if multiplicity is None:
            multiplicity = parser.multiplicity
            
        # Continue with other constructors.
        return self.from_xyz(parser.xyz, charge = charge, multiplicity = multiplicity, **kwargs)
        
        
    @classmethod
    def from_file(self, file_name, file_type = None, *, gen3D = None, **kwargs):
        """
        Create a Silico_coords object from a file in arbitrary format.
        
        :param file_name: Name/path of the input file to read from.
        :param file_type: The format of the file; a string recognised by openbabel. If not given, an attempt will be made to guess from the file name (see Openbabel_converter.type_from_file_name()).
        :param charge: The molecular charge.
        :param multiplicity: The molecular multiplicity (as an integer).
        :param name: Name of the system/molecule.
        """
        file_name = Path(file_name)
        silico.logging.get_logger().info("Parsing coordinate file '{}'".format(file_name))
        
        # Get the file format.
        if file_type is None:
            file_type = Openbabel_converter.type_from_file_name(file_name)
                
        # Certain formats we support natively; others we convert to an intermediate format.
        if file_type in ["com", "gau", "gjc", "gjf"]:
            # Gaussian input format.
            with open(file_name, "rt") as com_file:
                return self.from_com(com_file.read(), file_name = file_name, **kwargs)

        elif file_type == "si":
            # Silico input format.
            with open(file_name, "rt") as si_file:
                return self.from_si(si_file.read(), file_name = file_name, **kwargs)
            
        else:
            # Generic input format.
            
            # We convert all formats to gaussian input formats (because this format contains charge and multiplicity, which we can extract).
            com_file = Openbabel_converter.from_file(file_name, file_type).convert("com", gen3D = gen3D)         
        
            # Continue with other constructors.
            return self.from_com(com_file, file_name = file_name, **kwargs)
    
    @classmethod
    def from_yaml(self, yaml_dict, file_name = None, **kwargs):
        """
        Create a Silico_coords object from a loaded/parsed .si file.
        """
        yaml_dict = dict(yaml_dict)
        # Overwrite dictionary values if explicit ones have been given.
        for kwarg in kwargs:
            if kwargs[kwarg] is not None:
                yaml_dict[kwarg] = kwargs[kwarg]
            
        # Continue constructing.
        try:
            return self(**yaml_dict, file_name = file_name)
        except TypeError:
            raise Silico_exception("Failed to load .si file '{}'; is the file formatted correctly?".format(file_name))
    
    @classmethod
    def from_si(self, si_file, **kwargs):
        """
        Create a Silico_coords object from a raw .si file.
        """
        return self.from_yaml(yaml.safe_load(si_file), **kwargs)
        
    def to_file(self, file):
        """
        Write this silico input file to an open file object.
        
        :param file: A file opened with open().
        """
        file.write(self.yaml)
        
    def to_format(self, file_type, file = None):
        """
        Get this input file in an arbitrary format.
        
        :param file_type: The format of the file; a string recognised by openbabel.
        :param file: An optional file to write to, if not given the converted file is returned as a string.
        """
        if file_type.lower() == "si":
            # Convert to yaml
            return self.yaml
        else:
            # Convert.
            charge = int(self.charge) if self.charge is not None else None
            multiplicity = int(self.multiplicity) if self.multiplicity is not None else None
            return Openbabel_converter.get_cls("xyz")(input_file = self.xyz, input_file_path = self.auto_name, input_file_type = "xyz").convert(file_type, file, charge = charge, multiplicity = multiplicity)

    @classmethod
    def input_formats(self):
        """
        A dictionary of available input formats that this file can be created from.
        
        Each key is the short-code of the format (eg, si, com, xyz etc) while the value is a longer description.
        """
        formats = {
            "si": "Silico input format",
            "com": "Gaussian Input",
            "gau": "Gaussian Input",
            "gjc": "Gaussian Input",
            "gjf": "Gaussian Input"
        }
        formats.update(pybel.informats)
        return formats

    @classmethod
    def output_formats(self):
        """
        A dictionary of available output formats that this file can be converted to.
        
        Each key is the short-code of the format (eg, si, com, xyz etc) while the value is a longer description.
        """
        formats = {"si": "Silico input format"}
        formats.update(pybel.outformats)
        return formats

    @property
    def implicit_charge(self):
        """
        The charge of the molecule/system, accounting for cases where no explicit charge is set.
        """ 
        if self.charge is None:
            return 0
        else:
            return self.charge
        
    @property
    def implicit_multiplicity(self):
        """
        The multiplicity (as an integer) of the molecule/system, accounting for cases where no explicit multiplicity is set.
        """ 
        if self.multiplicity is None:
            return 1
        else:
            return self.multiplicity
        
    @property
    def dict(self):
        """
        Get this input file as a dict.
        """
        return {
            'name': self.name,
            'charge': float(self.charge),
            'multiplicity': float(self.multiplicity),
            'geometry': self.geometry,
        }
        
    @property
    def yaml(self):
        """
        Get this input file in yaml format.
        """
        # Get in dict format.
        intermediate = self.dict
        
        # Wrap geometry in literal format (so it will appear line by line).
        intermediate['geometry'] = literal_str(intermediate['geometry'])
        
        # Convert.
        return yaml.dump(intermediate, sort_keys=False)
        
        
    @property
    def xyz(self):
        """
        Get the geometry of this input file in XYZ format.
        """
        return "{}\n\n{}".format(len(self.geometry.strip().split("\n")), self.geometry) 
    