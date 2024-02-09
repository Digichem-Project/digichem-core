# General imports.
import subprocess
from subprocess import CalledProcessError
import re
import os
import copy
from pathlib import Path

# Silico imports.
from silico.exception.base import Silico_exception
import silico.log


# Try and load openbabel bindings.
HAVE_PYBEL = False

# Openbabel (and naturally also pybel) are GPL, and therefore incompatible with this library at present.
# try:
#     from openbabel import pybel
#     HAVE_PYBEL = True
# except ModuleNotFoundError:
#     # No bindings, carry on.
#     silico.log.get_logger().debug("Could not load python pybel bindings; falling back to obabel executable", exc_info = True)
# except Exception:
#     # Some other error occurred; print an error but continue.
#     silico.log.get_logger().error("Found but could not load python pybel bindings; falling back to obabel executable", exc_info = True)
    

class Openbabel_converter():
    """
    Top level class for openbabel wrappers.
    """
    
    def __init__(self, *, input_file = None, input_file_path = None, input_file_type = None):
        """
        Constructor for the OpenBabel converter.
        
        :param input_file: A string (unicode or byte) in the format given by input_file_type that should be converted.
        :param input_file_path: Alternatively, a Path to a file that should be converted.
        :param gen3D: If True and the loaded molecule does not have 3D coordinates, these will be generated (this will scramble atom coordinates).
        """
        # Do some arg checking.
        #if (input_file is None and input_file_path is None) or (input_file is not None and input_file_path is not None):
        #    raise TypeError(type(self).__name__ + "; exactly one of input_file or input_file_path must be given as argument")
        
        if input_file_path is None and input_file_type is None:
            raise TypeError("Missing argument: input_file_type")
        
        self.input_file = input_file
        self.input_file_path = input_file_path
        self.input_file_type = input_file_type
        # Currently, we always use add H because certain formats (xyz) cannot have H added.
        self.add_H = True
        
    @classmethod
    def type_from_file_name(self, input_file_name, allow_none = False):
        """
        Get the type of a file based on its file name.
        
        This method largely uses the file extension (.com, .tmol etc), with a few other simple rules.
        
        :param input_file_name: The name of the file to check.
        :param allow_none: If the type of the file cannot be determined and allow_none is False (the default), an exception is raised. Otherwise, None is returned.
        """
        input_file_name = Path(input_file_name)
        
        # Get file extension (removing the dot character).
        extension = input_file_name.suffix[1:]
        
        if extension != "":
            # All done.
            return extension.lower()
        elif input_file_name.name == "coord":
            # This is a turbomole file.
            return "tmol"
        else:
            if allow_none:
                return None
            
            else:
                # Don't recognise the file format.
                raise Silico_exception("Could not determine file format of file '{}'; the file does not have an extension and is not recognised".format(input_file_name))
        
        
    @property
    def input_name(self):
        """
        A descriptive name of the file we are converting. Works even if converting from memory.
        """
        if self.input_file_path is not None:
            return self.input_file_path
        else:
            return "(file loaded from memory)"
        
    @classmethod
    def from_file(self, input_file_path, input_file_type = None, **kwargs):
        """
        A more powerful constructor that automatically decides which concrete class to use.
        
        :param input_file_path: A Path to a file that should be converted.
        :param input_file_type: The format of the file; a string recognised by openbabel. If not given, an attempt will be made to guess from the file name (see type_from_file_name()).
        :param gen3D: If True and the loaded molecule does not have 3D coordinates, these will be generated (this will scramble atom coordinates).
        """
        # First, get out file format if it wasn't given to us.
        if input_file_type is None:
            input_file_type = self.type_from_file_name(input_file_path)
        
        # Next, decide which class
        cls = self.get_cls(input_file_type)
        # Normally we use input_file_path, not input_file.
        input_file = None
                
        # And return, if gen3D is not set, we turn it on for cdx files (which have to use the naive Obabel_converter and are always 2D).
        return cls(input_file_path = input_file_path, input_file = input_file, input_file_type = input_file_type, **kwargs)
        
        
    def convert(self, output_file_type):
        """
        Convert the input file wrapped by this class to the designated output_file_type.
        
        Inheriting classes should write their own implementation.
        
        :param output_file_type: The file type to convert to.
        """
        raise NotImplementedError("Abstract class Babel_converter does not have a convert() method defined (inheriting classes should write their own)")
        
    @classmethod
    def get_cls(self, input_file_type):
        """
        Automatically get a concrete Babel_converter class that can be used to convert a file.
        
        If the pybel bindings are available and loaded successfully; then a Pybel_converter object will be returned,
        otherwise, the Obabel_converter will be returned (this requires openbabel to be installed and obabel to be in the path).
        
        The only exception is for the cdx format for which Obabel_converter is always returned (because of bug https://github.com/openbabel/openbabel/issues/1690 which still seems to be plaguing us in mid-2020).
        """
        return Obabel_converter
#         
#         if not HAVE_PYBEL or input_file_type.lower() == "cdx":
#             return Obabel_converter
#         else:
#             return Pybel_converter
        
            

# Babel_convert doesn't inherit from File_converter because we are only interested in reading to/from stdin/out (which File_converter doesn't support)
#TODO: Add stdin/stdout support to File_converter
# if HAVE_PYBEL:
#     class Pybel_converter(Openbabel_converter):
#         """
#         Wrapper class for pybel
#         
#         """            
#         
#         def convert(self, output_file_type, output_file = None, *, gen3D = None, charge = None, multiplicity = None):
#             """
#             Convert the input file wrapped by this class to the designated output_file_type.
#             
#             :param output_file_type: The file type to convert to.
#             :param output_file: Optional file name to write to. If not given, the converted file will be returned as a string (or binary string depending on format).
#             :param gen3D: If True and the loaded molecule does not have 3D coordinates, these will be generated (this will scramble atom coordinates).
#             :param charge: Optional charge of the output format.
#             :param multiplicity: Optional multiplicity of the output format.
#             :return: The converted file, or None if output_file is not None.
#             """
#             output_file = str(output_file) if output_file is not None else None
#             if output_file is None and output_file_type == "png":
#                 raise ValueError("output_file must not be None if format is png")
#             
#             try:
#                 # Stop logging from pybel, we'll handle that ourselves.
#                 pybel.ob.obErrorLog.StopLogging()
#                 
#                 # For Pybel, gen3D defaults to True, because we'll only use gen3d if not already in 3D.
#                 gen3D = gen3D if gen3D is not None else True
#                 
#                 # Get upset if input_file_type is empty (because openbabel acts weird when it is).
#                 if self.input_file_type is None or self.input_file_type == "":
#                     raise TypeError("Cannot convert file; input_file_type '{}' is None or empty".format(self.input_file_type))
#                 
#                 # Read in the molecule(s) in the given file.
#                 try:
#                     # This is a generator.
#                     # Use a different func depending on whether we're reading from file or memory.
#                     if self.input_file is not None:
#                         # Reading from memory.
#                         molecule = pybel.readstring(self.input_file_type, str(self.input_file))
#                     else: 
#                         # Readfile gives us an iterator of molecules...
#                         molecules = pybel.readfile(self.input_file_type, str(self.input_file_path))
#                         
#                         # ...but we're only ever interested in one.
#                         # Try and get the first molecule.
#                         try:
#                             molecule = next(molecules)
#                         except StopIteration:
#                             raise Silico_exception("Cannot read file '{}'; file does not contain any molecules".format(self.input_name)) from None
#                         
#                         
#                 except Exception as e:
#                     raise Silico_exception("Failed to parse file '{}'".format(self.input_name)) from e
#                 
#                 if charge is not None:
#                     molecule.OBMol.SetTotalCharge(charge)
#                     
#                 if multiplicity is not None:
#                     molecule.OBMol.SetTotalSpinMultiplicity(multiplicity)
#                 
#                 # If we got a 2D (or 1D) format, convert to 3D (but warn that we are doing so.)
#                 if molecule.dim != 3 and gen3D:
#                     # We're missing 3D coords.
#                     silico.log.get_logger().warning("Generating 3D coordinates from {}D file '{}'; this will scramble atom coordinates".format(molecule.dim, self.input_name))
#                     molecule.localopt()
#                     
#                 if self.add_H:
#                     # Add hydrogens.
#                     #silico.log.get_logger().info("Adding any missing hydrogens to structure loaded from file '{}'".format(self.input_name))
#                     molecule.addh()
#                 
#                 # Now convert and return
#                 # If the format is png, use the draw() method instead because write() is bugged.
#                 if output_file_type == "png":
#                     molecule.draw(False, output_file)
#                 else:
#                     return molecule.write(output_file_type, output_file, overwrite = True)
#             
#             finally:
#                 # Start logging again.
#                 pybel.ob.obErrorLog.StartLogging()
#                 #pybel.ob.obErrorLog.SetOutputLevel(log_level)


class Obabel_converter(Openbabel_converter):
    """
    Wrapper class for openbabel.
    
    Obabel_convert calls the obabel command directly.
    """
    
    # The regex we'll use to check obabel converted successfully.
    obabel_fail = re.compile(r"\b0 molecules converted")
    
    # 'Path' to the obabel executable.
    obabel_execuable = "obabel"        
    
    def convert(self, output_file_type, output_file = None, *, gen3D = None, charge = None, multiplicity = None):
        """
        Convert the input file wrapped by this class to the designated output_file_type.
         
        :param output_file_type: The file type to convert to.
        :param output_file: Optional file name to write to. If not given, the converted file will be returned as a string (or binary string depending on format).
        :param gen3D: If True and the loaded molecule does not have 3D coordinates, these will be generated (this will scramble atom coordinates).        
        :param charge:  Optional charge of the output format.
        :param multiplicity: Optional multiplicity of the output format.
        :return: The converted file, or None if output_file is not None.
        """
        output_file = str(output_file) if output_file is not None else None
        # For Obabel, gen3D defaults to False, because we can't determine ahead of time whether we're in 3D or not (unless format is cdx, which is always 2D).
        if gen3D is None:
            if self.input_file_type.lower() == "cdx":
                gen3D = True
            else:
                gen3D = False
        
        if charge is not None:
            # We can't set charge with obabel sadly.
            silico.log.get_logger().warning("Unable to set charge '{}' of molecule loaded from file '{}' with obabel converter".format(charge, self.input_name))
            
        if multiplicity is not None:
            # We can't set multiplicity with obabel sadly.
            silico.log.get_logger().warning("Unable to set multiplicity '{}' of molecule loaded from file '{}' with obabel converter".format(multiplicity, self.input_name))
        
        # Run
        return self.run_obabel(output_file_type, output_file, gen3D = gen3D)
        
    def run_obabel(self, output_file_type, output_file, *, gen3D):
        """
        Run obabel, converting the input file wrapped by this class to the designated output_file_type.
        
        :param output_file_type: The file type to convert to.
        :param output_file: Optional file name to write to. If not given, the converted file will be returned as a string (or binary string depending on format).
        :param gen3D: If True and the loaded molecule does not have 3D coordinates, these will be generated (this will scramble atom coordinates).
        :return: The converted file.
        """
        # The signature we'll use to run obabel.
        sig = [self.obabel_execuable]
        
        # Add the input path if we're reading from file.
        if self.input_file is None:
            sig.append(str(self.input_file_path))
            
        # Now add the input and output switches.
        sig.extend([
             "-o", output_file_type,
             "-i", self.input_file_type
        ])
        
        # Add gen3D command if we've been asked to.
        if gen3D:
            silico.log.get_logger().warning("Generating 3D coordinates from file '{}'; this will scramble atom coordinates".format(self.input_name))
            sig.append("--gen3D")
        
        # Add H if we've been asked.    
        if self.add_H:
            #silico.log.get_logger().info("Adding any missing hydrogens to structure loaded from file '{}'".format(self.input_name))
            sig.append("-h")
            
        # If a file to write to has been given, set it.
        if output_file is not None:
            sig.extend(['-O', output_file])
        
        # There are several openbabel bugs re. the chem draw format; one of them occurs when we are frozen and have set the BABEL_LIBDIR env variable.
        # The workaround is to temp unset BABEL_LIBDIR.
        # Get our current environment.
        env = copy.copy(os.environ)
        # Now delete BABEL_LIBDIR if we are frozen.
        if silico.frozen:
            try:
                pass
                #del env['BABEL_LIBDIR']
            except KeyError:
                # The BABEL_LIBDIR isn't set.
                pass
        
        # Give our input_file as stdin if we're not reading from file.
        inputs = self.input_file        
        
        # GO.
        done_process = subprocess.run(
             sig,
             input = inputs,
             stdout = subprocess.PIPE,
             stderr = subprocess.PIPE,
             # TODO: Using universal newlines is probably not safe here; some formats are binary (.cdx etc...)
             universal_newlines = True,
             check = True,
             env = env
         )
        
        # Sadly, openbabel doesn't appear to make use of return codes all the time.
        # We'll do basic error checking on whether our output contains a certain string.
        #if not self.obabel_success.search(done_process.stderr):
        if self.obabel_fail.search(done_process.stderr):
            raise Silico_exception("obabel command did not appear to complete successfully") from CalledProcessError(done_process.returncode, " ".join(done_process.args), done_process.stdout, done_process.stderr)
        
        # Return our output.
        return done_process.stdout if output_file is None else None


class Obabel_formats():
    """
    Class for retrieving the supported file formats from obabel.
    """
    
    # Bit of a hack.
    obabel_execuable = Obabel_converter.obabel_execuable
    
    def __init__(self):
        """
        """
    
    def run(self, readwrite):
        """
        """
        if readwrite not in ["read", "write"]:
            raise ValueError("readwrite must be one of either 'read' or 'write'")
        
        # The signature we'll use to run obabel.
        sig = [
            self.obabel_execuable,
            "-L",
            "formats",
            readwrite
        ]
        
        # GO.
        done_process = subprocess.run(
             sig,
             stdout = subprocess.PIPE,
             stderr = subprocess.PIPE,
             # TODO: Using universal newlines is probably not safe here; some formats are binary (.cdx etc...)
             universal_newlines = True,
             check = True, # Obabel doesn't use return codes properly, but this doesn't hurt.
         )
        
        # The output should look like this:
        # abinit -- ABINIT Output Format
        # acesout -- ACES output format
        # acr -- ACR format
        # ...
        try:
            formats = {}
            for line in done_process.stdout.splitlines():
                code, desc = (part.strip() for part in line.split("--"))
                formats[code] = desc
                
        
        except Exception as e:
            # Output was bad.
            raise Silico_exception("Failed to parse stdout from obabel command") from e
        
        return formats

    
    def read(self, refresh = False):
        """
        Retrieve supported input (read) formats.
        """
        if refresh:
            try:
                del self._read
            
            except AttributeError:
                pass
                
        try:
            return self._read
        
        except AttributeError:
            # Cache miss.
            self._read = self.run("read")
            return self._read
        
    def write(self, refresh = False):
        """
        Retrieve supported input (read) formats.
        """
        if refresh:
            try:
                del self._write
            
            except AttributeError:
                pass
        
        try:
            return self._write
        
        except AttributeError:
            # Cache miss.
            self._write = self.run("write")
            return self._write
    
    