"""This module contains translation tables between various different, but equivalent, terms expected by different calculation programs."""

import basis_set_exchange as bse
import csv
from pathlib import Path
import pkg_resources

import silico.log
from silico.misc.base import is_int


# Load functional conversion table from file.
functional_db = {}
with open(Path(pkg_resources.resource_filename('silico', 'data/functionals.csv'))) as csv_file:
    _reader = csv.reader(csv_file)
    headers = None
    for line_num, line in enumerate(_reader):
        if line_num == 0:
            headers = line
        
        else:
            line = [line_part if line_part != "" else None for line_part in line]
            
            # The various parts of the definition.
            functional_name = line[0]
            aliases = line[1].split(",") if line[1] is not None else []
            functional = {"name": functional_name, "aliases": aliases, "turbomole": line[2], "gaussian": line[3]}
            
            # Get all the names we can refer to this functional by.
            all_names = [functional_name] + aliases
            for program in ("turbomole", "gaussian"):
                # Make sure we don't overwrite existing values with program specific ones.
                if functional[program] is not None and functional[program] not in functional_db:
                    all_names.append(functional[program])
            
            # Set into the DB, using all possible names.
            for ref_name in all_names:
                functional_db[ref_name] = functional


# Load solvent conversion table from file.
solvent_db = {}
with open(Path(pkg_resources.resource_filename('silico', 'data/solvents.csv'))) as csv_file:
    _reader = csv.reader(csv_file)
    headers = None
    for line_num, line in enumerate(_reader):
        if line_num == 0:
            headers = line
        
        else:
            line = [line_part if line_part != "" else None for line_part in line]
            
            # The various parts of the definition.
            solvent = {"name": line[0], "aliases": [line[1]] if line[1] is not None else [], "gaussian": line[2], "epsilon": float(line[3]) if line[3] is not None else None, "refractive": float(line[4]) if line[4] is not None else None}
            
            solvent_db[solvent['name']] = solvent


class Translate():
    """
    ABC for translate classes.
    """
    
    def __init__(self, value):
        """
        Constructor for Translate classes.
        """
        if isinstance(value, type(self)):
            self.value = value.value
        
        else:
            self.value = value
    
    def translate(self, program):
        """
        Translate into a name appropriate for a given program.
        """
        raise NotImplementedError("Implement in subclass")
    
    def to_turbomole(self):
        """
        Translate into a name appropriate for Turbomole.
        """
        return self.translate("turbomole")
    
    def to_gaussian(self):
        """
        Translate into a name appropriate for Turbomole.
        """
        return self.translate("gaussian")
    
    def to_orca(self):
        """
        Translate into a name appropriate for ORCA.
        """
        return self.translate("orca")


class Basis_set(Translate):
    """
    A class for converting between the names of basis sets.
    """
    
    @property
    def basis_set(self):
        return self.value
    
    @classmethod
    def build_db(self):
        """
        Build a table of basis set names.
        
        :returns: A dictionary, where each key is the name of a basis set and each value is a BSE metadata dict.
        """
        # Build a database of basis set metadata, adding the display name and alternative names as additional keys for each basis set.
        # BSE's database has some semi-duplicate entries (such as '6-31g(d,p)' and '6-31g_st__st_').
        # Although these basis sets are identical, they have reversed metadata (ie, for '6-31g(d,p)', 'display_name' is set to '6-31G(d,p)' and 'other_names'
        # contains '6-31G**', while for '6-31g_st__st_' 'display_name' is '6-31G**' and 'other_names' is '6-31G(d,p)'.
        # This is a problem for us because we need to predictable choose from these names (ie the star form '6-31G**' is suitable for both
        # Gaussian and Turbomole), but the keys referring to these same names are different.
        #
        # Firstly, build a new dict using basename as the key, which appears to be the same even for duplicates.
        db = {value["basename"]: value for key, value in bse.get_metadata().items()}
        
        new_db = {}
        
        for basis_key, basis_set_meta in db.items():
            # First, add the normal key.
            basis_set = {
                "key": basis_set_meta['basename'],
                "name": basis_set_meta['display_name'],
                "aliases": basis_set_meta['other_names'],
                "meta": basis_set_meta
            }
            
            # Decide on our program specific names.
            basis_set["gaussian"] = basis_set['name']
            basis_set["turbomole"] = basis_set['name']
            
            # Turbomole doesn't like pople sets with polarisation function spelled out,
            # it prefers the star format.
            if basis_set['name'][-5:] == "(d,p)":
                basis_set["turbomole"] = basis_set['name'][:-5] + "**"
            
            new_db[basis_key] = basis_set
                
        return new_db
    
    @classmethod
    def find_in_db(self, hint):
        """
        Try and find the entry for a basis set in the basis set exchange.
        
        :raises ValueError: If the given hint could not be found.
        :param hint: The name of the basis set to search for.
        :return: A dictionary of BSE metadata.
        """
        # Get the basis set DB.
        db = self.build_db()
        
        # First, try and see if we can just use hint as an exact match.
        try:
            return db[hint]
            
        except KeyError:
            pass
        
        # Try again, this time ignoring case and looking through all possible names.
        str_hint = "".join(str(hint).lower().split())
        for basis_set in db.values():
            if str_hint in [name.lower() for name in basis_set['aliases'] + [basis_set['key'], basis_set['name'], basis_set['gaussian'], basis_set['turbomole']] if name is not None]:
                return basis_set
        
        # No luck.
        raise ValueError("Could not find basis set definition for '{}'".format(hint))
    
    def translate(self, program = "name"):
        """
        Translate the name of this basis set into one appropriate for a calculation program.
        
        Basis sets are currently handled the same for all programs.
        """
        try:
            basis_def = self.find_in_db(self.basis_set)
            return basis_def[program]
        
        except ValueError:
            # Just return as is.
            if self.basis_set != "auto":
                silico.log.get_logger().debug("Could not find basis set with name '{}' in the basis set exchange; using name unmodified".format(self.basis_set))
            return self.basis_set
    
        except KeyError:
            return basis_def['name']
    
    def __str__(self):
        return self.translate()


class Functional(Translate):
    """
    A class for converting between the names of DFT functionals.
    """
    
    @property
    def functional(self):
        return self.value
            
    
    @classmethod
    def find_in_db(self, hint):
        """
        Try and find an entry for a functional in the internal library.
        
        :raises ValueError: If the functional could not be found.
        :param hint: The name of the functional to look for.
        :returns: The corresponding functional dict.
        """
        # First try using exact name.
        try:
            return functional_db[hint]
            
        except KeyError:
            pass
        
        # Try again, this time ignoring case.
        functional_names = [functional_name.upper() for functional_name in functional_db]
        
        try:
            return functional_db[list(functional_db.keys())[functional_names.index(hint.upper())]]
        
        except IndexError:
            pass
        
        # No luck.
        raise ValueError("Could not find functional definition for  '{}'".format(hint))
        
    
    def translate(self, program):
        """
        Translate this functional name to one recognised by a given program.
        
        :param program: The program to translate for.
        """
        # Try and get a definition for the functional from our db.
        try:
            func_def = self.find_in_db(self.functional)
                
            # If there's an explicit value for our program, use that.
            if program in func_def and func_def[program] is not None:
                return func_def[program]
            
            # Otherwise, use the main common name.
            return func_def['name']
            
        except ValueError:
            pass
        
        # No definition found, use as is.
        return self.functional
    
    def to_turbomole(self):
        """
        Translate into a name appropriate for Turbomole.
        """
        func_name = self.translate("turbomole")
        # All turbomole functionals are lower case, except for the word 'Gaussian'.
        return func_name.lower().replace("gaussian", "Gaussian")
    
    def __str__(self):
        return self.translate("name")


class Solvent(Translate):
    """A class for converting between different representations of solvents."""
    
    @property
    def solvent(self):
        return self.value
    
    @classmethod
    def find_in_db(self, hint):
        """
        Try and find an entry for a solvent in the internal library.
        
        :raises ValueError: If the multiplicity could not be found.
        :param hint: The name, number or symbol of the multiplicity to look for.
        :returns: The corresponding multiplicity dict.
        """
        # First, just try to lookup exactly.
        try:
            return solvent_db[hint]
            
        except KeyError:
            # No luck, look through for a match.
            str_hint = "".join(str(hint).lower().split())
            for solvent in solvent_db.values():
                if hint == solvent['epsilon'] or str_hint in [name.lower() for name in solvent['aliases'] + [solvent['name']] + [solvent['gaussian']] if name is not None]:
                    return solvent
        
        # No luck.
        raise ValueError("Could not find solvent definition for '{}'".format(hint))
    
    @property
    def epsilon(self):
        return self.find_in_db(self.solvent)['epsilon']
    
    @property
    def refractive_index(self):
        return self.find_in_db(self.solvent)['refractive']
    
    def translate(self, program):
        """
        Translate this functional name to one recognised by a given program.
        
        :param program: The program to translate for.
        """
        if program == "turbomole":
            program = "epsilon"
        
        try:
            solvent_def = self.find_in_db(self.solvent)
            
            if program in solvent_def and solvent_def[program] is not None:
                return solvent_def[program]
            
            else:
                return solvent_def["name"]
        
        except ValueError:
            return self.solvent
        
    def __str__(self):
        return self.translate("name")


class Multiplicity(Translate):
    """A class for converting between different representations of multiplicity."""
    
    table = [
        {"name": "no multiplicity", "number": 0, "symbol": "?"},
        {"name": "Singlet", "number": 1, "symbol": "S"},
        {"name": "Doublet", "number": 2, "symbol": "D"},
        {"name": "Triplet", "number": 3, "symbol": "T"},
        {"name": "Quartet", "number": 4, "symbol": "Q"}
    ]
    
    @property
    def multiplicity(self):
        return self.value
        
    @classmethod
    def find_in_db(self, hint):
        """
        Try and find an entry for a multiplicity in the internal library.
        
        :raises ValueError: If the multiplicity could not be found.
        :param hint: The name, number or symbol of the multiplicity to look for.
        :returns: The corresponding multiplicity dict.
        """
        for row in self.table:
            name, number, symbol = row.values()
            if hint == number or str(hint).upper() in (str(number), name.upper(), symbol.upper()):
                return row
            
        # No luck.
        raise ValueError("Could not find multiplicity definition for '{}'".format(hint))
    
    @property
    def symbol(self):
        """
        Get this multiplicity as a symbol.
        """
        # Get a shorthand symbol if we can.
        try:
            return self.find_in_db(self.multiplicity)['symbol']
            
        except ValueError:
            if self.multiplicity % 1 == 0:
                # Multiplicity is an integer, so return as a stringy whole number.
                return str(int(self.multiplicity))
            else:
                return str(self.multiplicity)
            
    @property
    def string(self):
        """
        Get this multiplicity as a string.
        """
        return self.translate("name")
    
    def __str__(self):
        return self.string
    
    @property
    def number(self):
        """
        Get this multiplicity as a number.
        """
        try:
            return self.find_in_db(self.multiplicity)['number']
        
        except ValueError:
            # No pre-defined number, see if it is an int.
            if is_int(self.multiplicity):
                return int(self.multiplicity)
            
            else:
                return float(self.multiplicity)
    
    def translate(self, to_type):
        """
        Translate into a name appropriate for a given program.
        """
        if to_type == "gaussian":
            to_type = "name"
        
        elif to_type == "turbomole":
            to_type = "number"
        
        
        try:
            return self.find_in_db(self.multiplicity)[to_type]
        
        except ValueError:
            # Couldn't find in conversion table.
            return self.multiplicity
