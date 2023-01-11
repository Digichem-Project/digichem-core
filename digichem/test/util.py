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


# Input files.
benzene_cdx = Path(data_directory(), "Input", "Benzene.cdx")
pyridine_cml = Path(data_directory(), "Input", "Pyridine.cml")
cyclopentane_com = Path(data_directory(), "Input", "Cyclopentane.com")
water_xyz = Path(data_directory(), "Input", "Water.xyz")
ethane_xyz = Path(data_directory(), "Input", "Ethane.xyz")
pyridine_si_v2 = Path(data_directory(), "Input/Pyridine.v2.si")
pyridine_si_v1 = Path(data_directory(), "Input/Pyridine.v1.si")
pyridine_resume_pickle = Path(data_directory(), "Input/silico.resume.pickle")


# Result files for testing parsing.
result_files = {
    "gaussian": [Path(data_directory(), datum) for datum in [
        'Naphthalene/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p)',
        'Naphthalene/Gaussian 16 Excited States TDA 10 Singlets 10 Triplets PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz',
        'Naphthalene/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz',  
    ]],
    "turbomole": [Path(data_directory(), datum) for datum in [
        'Naphthalene/Turbomole Optimisation ADC(2) cc-pVDZ.tar.gz',
        'Naphthalene/Turbomole Excited States ADC(2) S(1) and S(2) cc-pVDZ.tar.gz',
        'Naphthalene/Turbomole Excited States ADC(2) T(1) and T(2) cc-pVDZ.tar.gz',
    ]]
}

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


def check_dipole(dipole_moment, coords, abs = 1e-4):
    """
    Helper function to check a dipole moment (PDM, TEDM or TMDM) matches some expected values.
    """
    check_float_list(dipole_moment.vector_coords, coords, abs)
        
    # Check total
    assert dipole_moment.total == pytest.approx((coords[0] **2 + coords[1] **2 + coords[2] **2) **0.5, abs = abs)


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