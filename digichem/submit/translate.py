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
            solvent = {"name": line[0], "aliases": [line[1]] if line[1] is not None else [], "gaussian": line[2], "epsilon": float(line[3])}
            
            solvent_db[solvent['name']] = solvent


class Translate():
    """
    ABC for translate classes.
    """
    
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


class Basis_set(Translate):
    """
    A class for converting between the names of basis sets.
    """
    
    def __init__(self, basis_set):
        """
        Constructor for Basis_set translator objects.
        
        :param basis_set: The basis set to convert.
        """
        self.basis_set = basis_set
    
    @classmethod
    def build_db(self):
        """
        Build a table of basis set names.
        
        :returns: A dictionary, where each key is the name of a basis set and each value is a BSE metadata dict.
        """
        # Build a database of basis set metadata, adding the display name and alternative names as additional keys for each basis set.
        bse_db = bse.get_metadata()
        new_db = {}
        
        for basis_key, basis_set in bse_db.items():
            # First, add the normal key.
            new_db[basis_key] = basis_set
            
            # Add display names.
            new_db[basis_set['display_name']] = basis_set
            
            # Add each alternative name.
            for alt_name in basis_set['other_names']:
                new_db[alt_name] = basis_set
                
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
        
        
        # Try again, this time ignoring case.
        basis_names = [basis_name.upper() for basis_name in db]
        
        try:
            return db[list(db.keys())[basis_names.index(hint.upper())]]
        
        except IndexError:
            pass
        
        # Couldn't find it.
        raise ValueError("Could not find basis set definition for '{}'".format(hint))
    
    def translate(self, program = None):
        """
        Translate the name of this basis set into one appropriate for a calculation program.
        
        Basis sets are currently handled the same for all programs.
        """
        try:
            return self.find_in_db(self.basis_set)['display_name']
        
        except ValueError:
            # Just return as is.
            if self.basis_set != "auto":
                silico.log.get_logger().debug("Could not find basis set with name '{}' in the basis set exchange; using name unmodified".format(self.basis_set))
            return self.basis_set
    
    def __str__(self):
        return self.translate()


class Functional(Translate):
    """
    A class for converting between the names of DFT functionals.
    """
    
    def __init__(self, functional):
        """
        """
        self.functional = functional
            
    
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
            if func_def[program] is not None:
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
    
    def __init__(self, solvent):
        """
        """
        self.solvent = solvent
    
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
    
    def translate(self, program):
        """
        Translate this functional name to one recognised by a given program.
        
        :param program: The program to translate for.
        """
        if program == "turbomole":
            program = "epsilon"
        
        try:
            solvent_def = self.find_in_db(self.solvent)
            
            if solvent_def[program] is not None:
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
    
    def __init__(self, multiplicity):
        self.multiplicity = multiplicity
        
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
            if hint == number or hint.upper() in (name.upper(), symbol.upper()):
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
            self.find_in_db(self.multiplicity)['number']
        
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
        if to_type == "Gaussian":
            to_type = "name"
        
        elif to_type == "Turbomole":
            to_type = "number"
        
        
        try:
            return self.find_in_db(self.multiplicity)[to_type]
        
        except ValueError:
            # Couldn't find in conversion table.
            return self.multiplicity
