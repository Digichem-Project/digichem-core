# General imports.
import yaml
from pathlib import Path
from openbabel import pybel
import packaging.version

# Silico imports.
from silico.exception.base import Silico_exception
from silico.file.convert.gaussian import Gaussian_input_parser
from silico.file.convert.babel import Openbabel_converter
import silico.log
from silico.file.input import Input_file

# Custom formats to allow literal strings in yaml output.
# Adapted from https://stackoverflow.com/questions/6432605/any-yaml-libraries-in-python-that-support-dumping-of-long-strings-as-block-liter

# Custom 'str' class that will be dumped literally.
class literal_str(str): pass
# The dumper (which recognises the literal_str class).
def literal_str_representer(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')
# Add to yaml.
yaml.add_representer(literal_str, literal_str_representer)

class flow_mapping(dict): pass
def flow_mapping_representer(dumper, data):
    return dumper.represent_mapping( u'tag:yaml.org,2002:map', data, flow_style=True )

yaml.add_representer(flow_mapping, flow_mapping_representer)
    
###########
# Classes #
###########

class Silico_coords_ABC(Input_file):
    """
    ABC for classes that represents an input file in the silico input (.si) format.
    """

    @classmethod
    def from_xyz(self, geometry, **kwargs):
        """
        Create a Silico_coords object from a molecule in xyz format.
        
        :param geometry: The input geometry in xyz format.
        :param charge: The molecular charge.
        :param multiplicity: The molecular multiplicity (as an integer).
        :param name: Name of the system/molecule
        """
        raise NotImplementedError("Implement in subclass")
        
    @property
    def dict(self):
        """
        Get this input file as a dict.
        """
        return {
            'version': "2.0.0",
            'name': self.name,
            'charge': self.charge,
            'multiplicity': self.multiplicity,
            'atoms': self.atoms,
        }
        
    @property
    def yaml(self):
        """
        Get this input file in yaml format.
        """
        # Get in dict format.
        intermediate = self.dict
        
        # Wrap geometry in literal format (so it will appear line by line).
        intermediate['atoms'] = [flow_mapping(atom) for atom in intermediate['atoms']]
        
        # Convert.
        return yaml.dump(intermediate, sort_keys=False)
        
    @property
    def xyz(self):
        """
        Get the geometry of this input file in XYZ format.
        """
        return "{}\n\n{}".format(len(self.geometry.strip().split("\n")), self.geometry) 
        
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
        
        :param file: A file opened with open() (or similar).
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
            charge = self.charge if self.charge is not None else None
            multiplicity = self.multiplicity if self.multiplicity is not None else None
            return Openbabel_converter.get_cls("xyz")(input_file = self.xyz, input_file_path = self.implicit_name, input_file_type = "xyz").convert(file_type, file, charge = charge, multiplicity = multiplicity)

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
    
    def __eq__(self, other):
        """
        Check for equality with another coord object.
        """
        try:
            assert self.charge == other.charge
            assert self.multiplicity == other.multiplicity
            
            # TODO: This check should not care about the order of atoms.
            for index, atom in enumerate(self.atoms):
                assert atom == other.atoms[index]
            
            return True
        
        except AssertionError:
            return False
        


class Silico_coords_v2(Silico_coords_ABC):
    """
    Class that represents an input file in the silico input (.si) format (V2).
    
    The .si format is YAML based and stores atom positions in a dictionary format along with charge and multiplicity.
    """
    
    def __init__(self, atoms, *, charge = None, multiplicity = None, name = None, version = "2.0.0", file_name = None):
        """
        Constructor for .si files.
         
        :param atoms: A list of dictionaries specifying the geometry. Each dict should contain four keys, 'atom': The element, 'x': The x coord, 'y': The y coord, 'z': The z coord. All coordinates are in angstroms.
        :param charge: The molecular charge (as an integer).
        :param multiplicity: The molecular multiplicity (as an integer).
        :param name: Name of the system/molecule.
        :param version: The version of the .si file.
        :param file_name: The name of the file which was loaded. This can be used as a back-up file name.
        """
        super().__init__(charge, multiplicity, name, file_name)
        self.atoms = atoms
        self.version = version
        
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
        split_geom = geometry.strip().split("\n")
        
        # Remove the first line (the number of atoms) and the second (optional comment line)
        split_geom = split_geom[2:]
        
        atoms = []
        for geom_line in split_geom:
            atom, x_coord, y_coord, z_coord = geom_line.split()
            atoms.append({"atom": atom, "x": float(x_coord), "y": float(y_coord), "z": float(z_coord)})
        
        # Call our main constructor.
        return self(atoms, **kwargs)
        
    @property
    def geometry(self):
        """
        Get the geometry of this input file in XYZ format.
        """
        # Build the geometry section.
        # XYZ is a little bit loosely defined.
        # Obabel (defacto standard? uses ~10 chars for each column?).
        # That's 9 usable chars +1 obligate space.
        geometry = []
        for atom in self.atoms:
            geometry.append("{:<9} {:>14.8f} {:>14.8f} {:>14.8f}".format(atom['atom'], atom['x'], atom['y'], atom['z']).strip())
            
        return "\n".join(geometry)


class Silico_coords_v1(Silico_coords_ABC):
    """
    Class that represents an input file in the silico input (.si) format (V1).
     
    The .si format is yaml based, and stores the input geometry in a modified xyz format along with charge and multiplicity.
    """
    
    def __init__(self, geometry, *, charge = None, multiplicity = None, name = None, version = "1.0.0", file_name = None):
        """
        Constructor for .si files.
         
        :param geometry: The molecular geometry in .si format. Use from_xyz() instead if your format is in xyz.
        :param charge: The molecular charge (as an integer).
        :param multiplicity: The molecular multiplicity (as an integer).
        :param name: Name of the system/molecule.
        :param file_name: The name of the file which was loaded. This can be used as a back-up file name.
        """
        super().__init__(charge, multiplicity, name, file_name)
        self.geometry = geometry
        self.version = version

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
    
    @property
    def atoms(self):
        return Silico_coords_v2.from_xyz(self.xyz, charge = self.charge, multiplicity = self.multiplicity, name = self.name, file_name = self.file_name).atoms
        
#     @property
#     def dict(self):
#         """
#         Get this input file as a dict.
#         """
#         return {
#             'name': self.name,
#             'charge': self.charge,
#             'multiplicity': self.multiplicity,
#             'geometry': self.geometry,
#         }
#         
#     @property
#     def yaml(self):
#         """
#         Get this input file in yaml format.
#         """
#         # Get in dict format.
#         intermediate = self.dict
#         
#         # Wrap geometry in literal format (so it will appear line by line).
#         intermediate['geometry'] = literal_str(intermediate['geometry'])
#         
#         # Convert.
#         return yaml.dump(intermediate, sort_keys=False)


####################
# Helper utilities #
####################

Silico_coords = Silico_coords_v2

def si_from_yaml(yaml_dict, file_name = None, **kwargs):
    """
    """
    # In the first instance, look to see if there's a version string.
    if "version" in yaml_dict:
        version = packaging.version.parse(yaml_dict['version'])
        
        if version < packaging.version.parse("2"):
            # V1.
            cls = Silico_coords_v1
        
        elif version < packaging.version.parse("3"):
            # V2.
            cls = Silico_coords_v2
        
        else:
            raise Silico_exception("Unsupported .si file version '{}'".format(yaml_dict['version']))
    
    else:
        # No explicit version.
        # If we have the old 'geometry' section, assume v1.
        if "geometry" in yaml_dict:
            cls = Silico_coords_v1
            
        else:
            # Assume the latest format.
            cls = Silico_coords    
            
    return cls.from_yaml(yaml_dict, file_name, **kwargs)


def si_from_file(file_name, file_type = None, *, gen3D = None, **kwargs):
        """
        Create a Silico_coords object from a file in arbitrary format.
        
        :param file_name: Name/path of the input file to read from.
        :param file_type: The format of the file; a string recognised by openbabel. If not given, an attempt will be made to guess from the file name (see Openbabel_converter.type_from_file_name()).
        :param charge: The molecular charge.
        :param multiplicity: The molecular multiplicity (as an integer).
        :param name: Name of the system/molecule.
        """
        file_name = Path(file_name)
        silico.log.get_logger().info("Parsing coordinate file '{}'".format(file_name))
        
        # Get the file format.
        if file_type is None:
            file_type = Openbabel_converter.type_from_file_name(file_name)
                
        # Certain formats we support natively; others we convert to an intermediate format.
        if file_type in ["com", "gau", "gjc", "gjf"]:
            # Gaussian input format.
            with open(file_name, "rt") as com_file:
                return Silico_coords.from_com(com_file.read(), file_name = file_name, **kwargs)

        elif file_type == "si":
            # Silico input format.
            with open(file_name, "rt") as si_file:
                return si_from_yaml(yaml.safe_load(si_file.read()), file_name = file_name, **kwargs)
            
        else:
            # Generic input format.
            
            # We convert all formats to gaussian input formats (because this format contains charge and multiplicity, which we can extract).
            com_file = Openbabel_converter.from_file(file_name, file_type).convert("com", gen3D = gen3D)         
        
            # Continue with other constructors.
            return Silico_coords.from_com(com_file, file_name = file_name, **kwargs)
    