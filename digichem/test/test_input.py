"""Tests for parsing input files"""

import pytest
from pathlib import Path

from digichem.test.util import pyridine_si_v2, pyridine_si_v1, pyridine_cml,\
    result_files
from digichem.input.digichem_input import si_from_file

@pytest.mark.parametrize("file_path", [
        pyridine_si_v2,
        pyridine_si_v1
     ])
def test_input_reading(file_path):
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
    assert str(si_file.formula) == "NC5H5"


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


@pytest.mark.parametrize("file_path, sha", [
        [result_files['gaussian'][0], "ebd4f4e9f81cdb57ec2d2f2e1fba9cef0698f4c0"],
        [result_files['turbomole'][0], "60a8ebd9e701b849cfccd9cbb41684519a7fdf0b"],
        [result_files['orca'][0], "e48e7f653f4e67c1bd4c5c4bb76405fad2d441d0"],
     ])
def test_input_history(file_path, sha):
    """
    Test whether the history attribute is set properly.
    """
    si_file = si_from_file(file_path)
    
    assert si_file.history == sha
    
    dump = si_file.dict
    
    assert dump['history'] == sha
