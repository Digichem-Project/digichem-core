"""Tests for checking the values of parsed results."""

import pytest
from pathlib import Path
import scipy.constants

from silico.parser import parse_calculation
from silico.test import data_directory
from silico.test.util import check_float_list
from silico.parser.util import parse_calculations


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
def turbomole_opt_result():
    return parse_calculation(Path(data_directory(), "Naphthalene/Turbomole Optimisation Frequency PBE0 (GD3BJ) 6-31G**"))

@pytest.fixture(scope="module")
def turbomole_ES_singlets_result():
    return parse_calculation(Path(data_directory(), "Naphthalene/Turbomole Excited States TDA 10 Singlets PBE0 (GD3BJ) 6-31G**"))

@pytest.fixture(scope="module")
def turbomole_ES_triplets_result():
    return parse_calculation(Path(data_directory(), "Naphthalene/Turbomole Excited States TDA 10 Triplets PBE0 (GD3BJ) 6-31G**"))

@pytest.fixture(scope="module")
def turbomole_ES_result():
    return parse_calculations(
        Path(data_directory(), "Naphthalene/Turbomole Excited States TDA 10 Singlets PBE0 (GD3BJ) 6-31G**"),
        Path(data_directory(), "Naphthalene/Turbomole Excited States TDA 10 Triplets PBE0 (GD3BJ) 6-31G**")
    )

@pytest.fixture(scope="module")
def turbomole_PDM_result():
    return parse_calculation(Path(data_directory(), "Pyridine/Turbomole Optimisation Frequency PBE0 (GD3BJ) 6-31G**"))

@pytest.fixture(scope="module")
def turbomole_TDM_result():
    return parse_calculation(Path(data_directory(), "Pyridine/Turbomole Excited States TDA 10 Singlets PBE0 (GD3BJ) 6-31G**"))

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
        (pytest.lazy_fixture("turbomole_opt_result"), 48, 13, 805.28, 103.77),
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
        (pytest.lazy_fixture("turbomole_PDM_result"), (-0.000743401, -2.20405438, -0.000494493), 90.0, 0.0),
    ])
def test_pdm(result_set, coords, axis_angle, plane_angle):
    """Test the permanent dipole moment"""
     
    check_dipole(result_set.dipole_moment, coords)
    
    # Check angles
    assert float(result_set.dipole_moment.X_axis_angle) == pytest.approx(axis_angle, abs=1e-1)
    assert float(result_set.dipole_moment.XY_plane_angle) == pytest.approx(plane_angle, abs=1e-1)


@pytest.mark.parametrize("result_set, number, TDM", [
        (pytest.lazy_fixture("gaussian_TDM_result"), 20, (0.0, -0.0003, 0.5879)),
        (pytest.lazy_fixture("turbomole_TDM_result"), 10, (-0.000628, -0.000132, 0.558096))
    ])
def test_tedm(result_set, number, TDM):
    """Test transition dipole moments"""
    # Check number of dipoles.
    assert len([excited_state.transition_dipole_moment for excited_state in result_set.excited_states if excited_state.transition_dipole_moment.electric is not None]) == number
    
    # Check the S1 moments are correct.
    S1_TDM = result_set.excited_states.get_state("S(1)").transition_dipole_moment
    
    # Check coords and magnitude.
    check_dipole(S1_TDM.electric, TDM)


@pytest.mark.parametrize("result_set, number, TDM", [
        (pytest.lazy_fixture("gaussian_TDM_result"), 20, (0.6909, 0.0, 0.0)),
        (pytest.lazy_fixture("turbomole_TDM_result"), 10, (0.002457, -0.000001, -0.000001))
    ])
def test_tmdm(result_set, number, TDM):
    """Test transition dipole moments"""
    if result_set.metadata.package == "Turbomole":
        pytest.skip("Turbomole does not yet support parsing TMDM")
    
    # Check number of dipoles.
    assert len([excited_state.transition_dipole_moment for excited_state in result_set.excited_states if excited_state.transition_dipole_moment.magnetic is not None]) == number
    
    # Check the S1 moments are correct.
    S1_TDM = result_set.excited_states.get_state("S(1)").transition_dipole_moment
    
    # Check coords and magnitude.
    check_dipole(S1_TDM.magnetic, TDM)


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


@pytest.mark.parametrize("result_set, num_singlets, num_triplets, dest, state_label, state_index, state_symmetry, state_energy, state_f", [
        (pytest.lazy_fixture("gaussian_ES_result"), 10, 10, 1.6231, "S(2)", 7, "Singlet-B1U", 4.7387, 0.1168),
        (pytest.lazy_fixture("turbomole_ES_result"), 10, 10, 1.6298, "S(2)", 7, "Singlet-A", 4.790868037954573, 0.0805)
    ])
def test_excited_states(result_set, num_singlets, num_triplets, dest, state_label, state_index, state_symmetry, state_energy, state_f):
    """Test excited states"""
    
    # Check number of states.
    assert result_set.excited_states.num_singlets == num_singlets
    assert result_set.excited_states.num_triplets == num_triplets
    
    # Check singlet-triplet splitting energy (dE(ST)).
    assert result_set.excited_states.singlet_triplet_energy == pytest.approx(dest, abs=1e-4)
    
    # Now check a specific state.
    state = result_set.excited_states.get_state(state_label)
    assert state.level == state_index
    assert state.symmetry == state_symmetry
    assert state.energy == pytest.approx(state_energy, abs=1e-4)
    assert state.oscillator_strength == pytest.approx(state_f, abs=1e-4)
    
    # Check wavelength.
    # e in J
    joule_energy = state_energy * scipy.constants.physical_constants['electron volt'][0]
    # wavelength in m
    meter_wavelength = (scipy.constants.Planck * scipy.constants.c) / joule_energy
    # wavelength in nm
    wavelength = meter_wavelength *1000000000
    
    assert state.wavelength == pytest.approx(wavelength, abs=1e-4)


@pytest.mark.parametrize("result_set, state_index, transition_index, start_orbital, end_orbital, coefficient", [
        (pytest.lazy_fixture("gaussian_ES_result"), 1, 1, 34, 35, 0.957238734),
        (pytest.lazy_fixture("gaussian_ES_result"), 1, 2, 32, 37, 0.186520627),
        (pytest.lazy_fixture("gaussian_ES_result"), 1, 3, 33, 36, 0.175928167),
    ])
def test_excited_state_transitions(result_set, state_index, transition_index, start_orbital, end_orbital, coefficient):
    """Test excited state transitions."""
    
    state = result_set.excited_states[state_index -1]
    transition = state.transitions[transition_index -1]
    
    assert transition.level == transition_index
    assert transition.starting_mo.level == start_orbital
    assert transition.ending_mo.level == end_orbital
    assert transition.coefficient == pytest.approx(coefficient)
    assert transition.probability == pytest.approx(coefficient **2)
    
