"""Tests for result parsing and processing"""

from pathlib import Path
import itertools
import pytest

from silico.parser import parse_calculation
from silico.test import data_directory
from silico.result import Result_set
from silico.parser import parse_multiple_calculations, parse_calculations


gaussian_files = [Path(data_directory(), datum) for datum in [
        'Naphthalene/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p)',
        'Naphthalene/Gaussian 16 Excited States TDA 10 Singlets 10 Triplets PBE1PBE (GD3BJ) Toluene 6-31G(d,p)',
        'Naphthalene/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p)',  
    ]]

turbomole_files = [Path(data_directory(), datum) for datum in [
        'Naphthalene/Turbomole Optimisation ADC(2) cc-pVDZ',
        'Naphthalene/Turbomole Excited States ADC(2) S(1) and S(2) cc-pVDZ',
        'Naphthalene/Turbomole Excited States ADC(2) T(1) and T(2) cc-pVDZ',
    ]]


@pytest.mark.parametrize("result_data", list(itertools.chain(gaussian_files, turbomole_files)))
def test_parsing(result_data):
    """Test whether we can parse various calc results."""
    result = parse_calculation(result_data)
    assert isinstance(result, Result_set)

        
def test_multi_parsing():
    """Test whether we can parse in parallel."""
    
    results = parse_multiple_calculations(*list(itertools.chain(gaussian_files, turbomole_files)))
    
    # Check length.
    assert len(results) == len(gaussian_files) + len(turbomole_files)


@pytest.mark.parametrize("data_set", [gaussian_files, turbomole_files])
def test_merged_parsing(data_set):
    """Test whether we can parse and merge calc results."""
    result = parse_calculations(*data_set)
    
    assert isinstance(result, Result_set)
    assert len(result.results) == 3