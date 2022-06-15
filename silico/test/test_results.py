"""Tests for checking the values of parsed results."""

import pytest
from pathlib import Path

from silico.parser import parse_calculation
from silico.test import data_directory


@pytest.fixture(scope="module")
def gaussian_SP_result():
    return parse_calculation(Path(data_directory(), "Naphthalene/Gaussian 16 Single Point (Singlet) PBE1PBE (GD3BJ) Toluene 6-31G(d,p)"))

@pytest.fixture(scope="module")
def gaussian_opt_result():
    return parse_calculation(Path(data_directory(), "Naphthalene/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p)"))

@pytest.fixture(scope="module")
def gaussian_ES_result():
    return parse_calculation(Path(data_directory(), "Naphthalene/Gaussian 16 Excited States TDA 10 Singlets 10 Triplets PBE1PBE (GD3BJ) Toluene 6-31G(d,p)"))

@pytest.fixture(scope="module")
def gaussian_opt_ES_result():
    return parse_calculation(Path(data_directory(), "Naphthalene/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p)"))

@pytest.fixture(scope="module")
def turbomole_ADC2_opt_result():
    return parse_calculation(Path(data_directory(), "Naphthalene/Turbomole Optimisation ADC(2) cc-pVDZ"))

@pytest.fixture(scope="module")
def turbomole_ADC2_singlets_result():
    return parse_calculation(Path(data_directory(), "Naphthalene/Turbomole Excited States ADC(2) S(1) and S(2) cc-pVDZ"))

@pytest.fixture(scope="module")
def turbomole_ADC2_triplets_result():
    return parse_calculation(Path(data_directory(), "Naphthalene/Turbomole Excited States ADC(2) T(1) and T(2) cc-pVDZ"))


@pytest.mark.parametrize("result_set, num, final", [
        (pytest.lazy_fixture("gaussian_SP_result"), 1, -10488.995711),
        (pytest.lazy_fixture("gaussian_opt_result"), 5, -10488.995711),
        (pytest.lazy_fixture("gaussian_ES_result"), 1, -10488.995711),
        (pytest.lazy_fixture("gaussian_opt_ES_result"), 4, -10488.888883),
        (pytest.lazy_fixture("turbomole_ADC2_opt_result"), 7, -10432.315825879925),
        (pytest.lazy_fixture("turbomole_ADC2_singlets_result"), 1, -10432.315825879925),
        (pytest.lazy_fixture("turbomole_ADC2_triplets_result"), 1, -10432.315825879925)
    ])
def test_SCF_energy(result_set, num, final):
    """Test the parsed energy is correct"""
    assert result_set.SCF_energies.final == pytest.approx(final)
    
    # Check length, which will be 1 for SP, and >1 for the opts.
    assert len(result_set.SCF_energies) == num


@pytest.mark.parametrize("result_set, num, final", [
        (pytest.lazy_fixture("turbomole_ADC2_opt_result"), 7, -10467.162645663258),
        (pytest.lazy_fixture("turbomole_ADC2_singlets_result"), 1, -10467.162645663258),
        (pytest.lazy_fixture("turbomole_ADC2_triplets_result"), 1, -10467.162645663258)
    ])
def test_MP_energy(result_set, num, final):
    """Test the parsed energy is correct"""
    assert result_set.MP_energies.final == pytest.approx(final)
    
    # Check length, which will be 1 for SP, and >1 for the opts.
    assert len(result_set.SCF_energies) == num
    assert len(result_set.MP_energies) == num


@pytest.mark.parametrize("result_set, charge, mult, energy", [
        (pytest.lazy_fixture("gaussian_opt_result"), 0, 1, -10488.995711),
        (pytest.lazy_fixture("turbomole_ADC2_opt_result"), 0, 1, -10467.162645663258),
    ])
def test_ground_state(result_set, charge, mult, energy):
    """Test the ground state properties"""
    
    assert result_set.ground_state.charge == charge
    assert result_set.ground_state.multiplicity == mult
    # Ground states have zero excited state energy by definition.
    assert result_set.ground_state.energy == 0.0
    # This is the total energy.
    assert result_set.ground_state.absolute_energy == pytest.approx(energy)


@pytest.mark.parametrize("result_set", [
        pytest.lazy_fixture("gaussian_SP_result"),
        pytest.lazy_fixture("gaussian_opt_result"),
        pytest.lazy_fixture("gaussian_ES_result"),
        pytest.lazy_fixture("gaussian_opt_ES_result"),
        pytest.lazy_fixture("turbomole_ADC2_opt_result"),
        pytest.lazy_fixture("turbomole_ADC2_singlets_result"),
        pytest.lazy_fixture("turbomole_ADC2_triplets_result"),
    ])
def test_atoms(result_set):
    """Test the parsed (unaligned) atoms are correct."""
    # The molecule is Naphthalene; C10H8
    # Check num atoms.
    assert len(result_set.atoms) == 18
    
    # Check atom types.
    assert result_set.atoms.element_dict['C'] == 10
    assert result_set.atoms.element_dict['H'] == 8
    
    # Check mass.
    assert result_set.atoms.molar_mass == pytest.approx(128.1705)
    
    # We don't check positions here because these can vary from calc to calc due to reorientation...


@pytest.mark.parametrize("result_set, length, width, height", [
        (pytest.lazy_fixture("gaussian_SP_result"), 6.739394, 4.972394, 0.000004),
        (pytest.lazy_fixture("gaussian_opt_result"), 6.739394, 4.972394, 0.000004),
        (pytest.lazy_fixture("gaussian_ES_result"), 6.739394, 4.972394, 0.000004),
        (pytest.lazy_fixture("turbomole_ADC2_opt_result"), 6.8000171599, 5.0161488999, 0.0000475622),
        (pytest.lazy_fixture("turbomole_ADC2_singlets_result"), 6.8000171599, 5.0161488999, 0.0000475622),
        (pytest.lazy_fixture("turbomole_ADC2_triplets_result"), 6.8000171599, 5.0161488999, 0.0000475622),
    ])
def test_alignment(result_set, length, width, height):
    """Test the aligned atom positions are correct."""
    # Rounding errors plague parsing and conversion of atomic coordinates.
    # Much of this error comes from the accuracy printed by the CC programs themselves,
    # so there is not much we can do about it.
    tolerance = 1e-5
    
    # Dimensions in angstrom.
    try:
        linear_ratio = 1-(width/length)
    except (FloatingPointError, ZeroDivisionError):
        linear_ratio = 1.0
        
    try:
        planar_ratio = 1-(height/width)
    except (FloatingPointError, ZeroDivisionError):
        planar_ratio = 1.0
    
    # Check dimensions.
    assert result_set.alignment.X_length == pytest.approx(length, abs=tolerance)
    assert result_set.alignment.Y_length == pytest.approx(width, abs=tolerance)
    assert result_set.alignment.Z_length == pytest.approx(height, abs=tolerance)
    
    # Check params.
    assert result_set.alignment.get_linear_ratio() == pytest.approx(linear_ratio, abs=tolerance)
    assert result_set.alignment.get_planar_ratio() == pytest.approx(planar_ratio, abs=tolerance)


@pytest.mark.parametrize("result_set, num_occ, num_unocc, homo, lumo", [
        (pytest.lazy_fixture("gaussian_opt_result"), 34, 156, -6.13072481, -0.92437071),
        (pytest.lazy_fixture("turbomole_ADC2_opt_result"), 34, 146, -7.78345591, 2.37048311)
    ])
def test_orbitals(result_set, num_occ, num_unocc, homo, lumo):
    """Test the energies of the parsed MOs"""
    # Check numbers.
    assert len(result_set.molecular_orbitals) == num_occ + num_unocc
    assert len(result_set.molecular_orbitals.occupied) == num_occ
    assert len(result_set.molecular_orbitals.virtual) == num_unocc
    
    # Check energies.
    assert result_set.molecular_orbitals.HOMO_energy == pytest.approx(homo)
    assert result_set.molecular_orbitals.LUMO_energy == pytest.approx(lumo)
    assert result_set.molecular_orbitals.HOMO_LUMO_energy == pytest.approx(lumo - homo)
    
    # Check ordering.
    assert [orbital.level for orbital in result_set.molecular_orbitals] == list(range(1, num_occ + num_unocc +1))

