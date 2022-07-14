"""Common testing utilities and convenience functions."""

import pytest
from pathlib import Path
import pkg_resources

from silico.config.parser import Config_file_parser


def data_directory():
    """
    Get a path to the test data directory.
    """
    return Path(pkg_resources.resource_filename('silico', 'test/data'))

@pytest.fixture(scope="package")
def silico_options():
    return Config_file_parser.silico_options()

def check_float_list(test_list, expected_list, abs = 1e-4):
    """
    Helper function to compare to lists of floats.
    """
    # Pytest bug #9921 prevents this comparison from working for now...
    #assert test_list == pytest.approx(expected_list, abs=abs)
    for list_index in range(0,len(expected_list)):
        assert test_list[list_index] == pytest.approx(expected_list[list_index], abs=abs)


def check_dipole(dipole_moment, coords):
    """
    Helper function to check a dipole moment (PDM, TEDM or TMDM) matches some expected values.
    """
    check_float_list(dipole_moment.vector_coords, coords)
        
    # Check total
    assert dipole_moment.total == pytest.approx((coords[0] **2 + coords[1] **2 + coords[2] **2) **0.5, abs=1e-4)


def check_orbitals(orbitals, num_occ, num_unocc, homo, lumo):
    """Helper function for checking orbital energies"""
    # Check numbers.
    assert len(orbitals) == num_occ + num_unocc
    assert len(orbitals.occupied) == num_occ
    assert len(orbitals.virtual) == num_unocc
    
    # Check energies.
    assert orbitals.HOMO_energy == pytest.approx(homo)
    assert orbitals.LUMO_energy == pytest.approx(lumo)
    assert orbitals.HOMO_LUMO_energy == pytest.approx(lumo - homo)
    
    # Check ordering.
    assert [orbital.level for orbital in orbitals] == list(range(1, num_occ + num_unocc +1))