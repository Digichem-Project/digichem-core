"""Tests for the native silico input format .si"""

import pytest
from pathlib import Path

from silico.test.util import pyridine_si_v2, pyridine_si_v1, pyridine_cml
from silico.file.input.coord import si_from_file

@pytest.mark.parametrize("file_path", [
        pyridine_si_v2,
        pyridine_si_v1
     ])
def test_si_reading(file_path):
    """
    Test whether we can correctly read from an existing si file.
    """
    si_file = si_from_file(file_path)
    
    # Check charge and mult.
    assert si_file.charge == 0
    assert si_file.multiplicity == 1
    
    # Check geometry.
    # All files should contain data for pyridine.
    assert len(si_file.atoms) == 11
    assert si_file.formula == "C5H5N"


def test_si_writing(tmp_path):
    """
    Test whether we can correctly write a new .si file.
    """
    # Read in an input file.
    si_file = si_from_file(pyridine_cml)
    
    # Write to our temp dir.
    new_si_path = Path(tmp_path, "Pyrdine.si")
    with open(new_si_path, "wt") as out_file:
        si_file.to_file(out_file)
        
    # Parse the newly written file.
    new_si_file = si_from_file(new_si_path)
    
    # Check they're the same.
    assert new_si_file == si_file