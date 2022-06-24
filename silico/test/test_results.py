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
def gaussian_radical_anion_result():
    return parse_calculation(Path(data_directory(), "Benzene Anion/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Gas Phase 6-31G(d,p)"))

@pytest.fixture(scope="module")
def gaussian_ES_result():
    return parse_calculation(Path(data_directory(), "Naphthalene/Gaussian 16 Excited States TDA 10 Singlets 10 Triplets PBE1PBE (GD3BJ) Toluene 6-31G(d,p)"))

@pytest.fixture(scope="module")
def gaussian_opt_ES_result():
    return parse_calculation(Path(data_directory(), "Naphthalene/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p)"))

@pytest.fixture(scope="module")
def gaussian_PDM_result():
    return parse_calculation(Path(data_directory(), "Pyridine/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p)"))

@pytest.fixture(scope="module")
def gaussian_TDM_result():
    return parse_calculation(Path(data_directory(), "Pyridine/Gaussian 16 Excited States TDA 10 Singlets 10 Triplets PBE1PBE (GD3BJ) Toluene 6-31G(d,p)"))

@pytest.fixture(scope="module")
def gaussian_PDM_ES_result():
    return parse_calculation(Path(data_directory(), "Pyridine/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p)"))

@pytest.fixture(scope="module")
def turbomole_radical_anion_result():
    return parse_calculation(Path(data_directory(), "Benzene Anion/Turbomole Optimisation Frequency PBE0 (GD3BJ) 6-31G**"))

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
    
    if num > 1:
        pytest.skip("A cclib bug currently parses the wrong number of MP energies")
    # Check length, which will be 1 for SP, and >1 for the opts.
    assert len(result_set.SCF_energies) == num
    assert len(result_set.MP_energies) == num


@pytest.mark.parametrize("result_set, charge, mult, energy", [
        (pytest.lazy_fixture("gaussian_opt_result"), 0, 1, -10488.995711),
        (pytest.lazy_fixture("turbomole_ADC2_opt_result"), 0, 1, -10467.162645663258),
        # Radical anions.
        (pytest.lazy_fixture("gaussian_radical_anion_result"), -1, 2, -6310.537803021298),
        (pytest.lazy_fixture("turbomole_radical_anion_result"), -1, 2, -6310.4569545463946),
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


@pytest.mark.parametrize("result_set, num_occ, num_unocc, homo, lumo", [
        (pytest.lazy_fixture("gaussian_opt_result"), 34, 156, -6.13072481, -0.92437071),
        (pytest.lazy_fixture("turbomole_ADC2_opt_result"), 34, 146, -7.78345591, 2.37048311)
    ])
def test_orbitals(result_set, num_occ, num_unocc, homo, lumo):
    """Test the energies of the parsed MOs"""
    check_orbitals(result_set.molecular_orbitals, num_occ, num_unocc, homo, lumo)
    
    # Check spin labels
    assert result_set.molecular_orbitals.spin_type == "none"


@pytest.mark.parametrize("result_set, num_occ, num_unocc, homo, lumo", [
        (pytest.lazy_fixture("gaussian_radical_anion_result"), (22, 21), (98, 99), (3.832451319, -0.0157826027), (6.268142298, 6.34270149)),
        (pytest.lazy_fixture("turbomole_radical_anion_result"), (22, 21), (92, 93), (3.825336056180260069, 0.0034865153882969327284), (6.275672672175910627, 6.387468734183322283))
    ])
def test_unrestricted_orbitals(result_set, num_occ, num_unocc, homo, lumo):
    """Test the energies of unrestricted orbitals"""
    # Check overall length.
    assert len(result_set.molecular_orbitals) == len(result_set.beta_orbitals)
    
    # Check each set of orbitals.
    check_orbitals(result_set.molecular_orbitals, num_occ[0], num_unocc[0], homo[0], lumo[0])
    check_orbitals(result_set.beta_orbitals, num_occ[1], num_unocc[1], homo[1], lumo[1])
    
    # Check spin labels
    assert result_set.molecular_orbitals.spin_type == "alpha"
    assert result_set.beta_orbitals.spin_type == "beta"


@pytest.mark.parametrize("result_set, num, index, frequency, intensity", [
        (pytest.lazy_fixture("gaussian_opt_result"), 48, 13, 806.2853, 121.0015),
        (pytest.lazy_fixture("gaussian_SP_result"), 0, None, None, None),
    ])
def test_frequencies(result_set, num, index, frequency, intensity):
    """Test the vibrational frequencies are correct"""
    
    # Check lengths.
    assert len(result_set.vibrations) == num
    assert len(result_set.vibrations.negative_frequencies) == 0

    if index != None:
        # Check one of the energies.
        assert result_set.vibrations[index].frequency == pytest.approx(frequency)
        assert result_set.vibrations[index].intensity == pytest.approx(intensity)
        
        # Check the frequencies are ordered.
        indices = [freq.level for freq in result_set.vibrations]
        assert indices == list(range(1, num+1))
        
    else:
        with pytest.raises(IndexError):
            result_set.vibrations[-1].frequency


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


@pytest.mark.parametrize("result_set, coords, axis_angle, plane_angle", [
        (pytest.lazy_fixture("gaussian_PDM_result"), (0.0, 2.5103, 0.0), 90.0, 0.0),
        (pytest.lazy_fixture("gaussian_PDM_ES_result"), (0.0001, -0.6147, 0.0001), 90.0, 0.0),
    ])
def test_pdm(result_set, coords, axis_angle, plane_angle):
    """Test the permanent dipole moment"""
     
    check_dipole(result_set.dipole_moment, coords)
    
    # Check angles
    assert float(result_set.dipole_moment.X_axis_angle) == pytest.approx(axis_angle, abs=1e-2)
    assert float(result_set.dipole_moment.XY_plane_angle) == pytest.approx(plane_angle, abs=1e-2)


@pytest.mark.parametrize("result_set, number, S1_TEDM, S1_TMDM", [
        (pytest.lazy_fixture("gaussian_TDM_result"), 20, (0.0, -0.0003, 0.5879), (0.6909, 0.0, 0.0))
    ])
def test_tdm(result_set, number, S1_TEDM, S1_TMDM):
    """Test transition dipole moments"""
    # Check number of dipoles.
    assert len([excited_state.transition_dipole_moment for excited_state in result_set.excited_states if excited_state.transition_dipole_moment is not None]) == number
    
    # Check the S1 moments are correct.
    S1_TDM = result_set.excited_states.get_state("S(1)").transition_dipole_moment
    
    # Check coords and magnitude.
    check_dipole(S1_TDM.electric, S1_TEDM)
    check_dipole(S1_TDM.magnetic, S1_TMDM)


@pytest.mark.parametrize("result_set, state, g_lum", [
        (pytest.lazy_fixture("gaussian_TDM_result"), "S(5)", -0.995)
    ])
def test_g_lum(result_set, state, g_lum):
    """Test calculation of G(lum) values for circularly polarized light."""
    
    TDM = result_set.excited_states.get_state(state).transition_dipole_moment
    
    TEDM = TDM.electric.gaussian_cgs_vector
    TMDM = TDM.magnetic.gaussian_cgs_vector
    
    # Test values in Gaussian-CGS units.))
    check_float_list(TEDM, tuple(254.35e-20 * i for i in TEDM))
    check_float_list(TMDM, tuple(-0.97401e-20 * i for i in TMDM))
    
    # Check angle.
    angle = (TEDM[0] * TMDM[0] + TEDM[1] * TMDM[1] + TEDM[2] * TMDM[2]) / ( ((TEDM[0]**2 + TEDM[1] **2 + TEDM[2] **2)**0.5) * ((TMDM[0]**2 + TMDM[1] **2 + TMDM[2] **2)**0.5) )
    assert TDM.cos_angle() == pytest.approx(angle)
    
    # Check g_lum
    assert TDM.g_value == pytest.approx(g_lum, abs=1e-3)