"""Tests for checking the values of parsed results."""

import pytest
from pathlib import Path

from silico.parser import parse_calculation
from silico.test import data_directory


@pytest.fixture
def gaussian_opt_result():
    return parse_calculation(Path(data_directory(), "Naphthalene/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p)"))

@pytest.fixture
def gaussian_ES_result():
    return parse_calculation(Path(data_directory(), "Naphthalene/Gaussian 16 Excited States TDA 10 Singlets 10 Triplets PBE1PBE (GD3BJ) Toluene 6-31G(d,p)"))

@pytest.fixture
def gaussian_opt_ES_result():
    return parse_calculation(Path(data_directory(), "Naphthalene/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p)"))


def test_energy(gaussian_opt_result, gaussian_ES_result, gaussian_opt_ES_result):
    """Test the parsed energies are correct."""
    
    # These are DFT calcs, so only SCF energy is available.
    # Check length.
    assert len(gaussian_opt_result.SCF_energies) == 5
    assert len(gaussian_ES_result.SCF_energies) == 1
    assert len(gaussian_opt_ES_result.SCF_energies) == 4
    
    # Check values.
    assert gaussian_opt_result.SCF_energies.final == pytest.approx(-10488.995711)
    assert gaussian_ES_result.SCF_energies.final == pytest.approx(-10488.995711)
    assert gaussian_opt_ES_result.SCF_energies.final == pytest.approx(-10488.888883)
    
    # Check other energies are empty.
    for result in [gaussian_opt_result, gaussian_ES_result, gaussian_opt_ES_result]:
        assert len(result.MP_energies) == 0
        assert len(result.CC_energies) == 0