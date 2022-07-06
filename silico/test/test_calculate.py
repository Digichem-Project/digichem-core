"""Tests for calculated results"""

import pytest



def test_soc(result_set, s0_t1, s1_t1):
    """Test calculation of spin-orbit coupling values."""
    
    # Check we have SOC between each singlet state (including the ground) and each triplet state.
    assert len(result_set.spin_orbit_coupling) == (result_set.excited_states.num_singlets +1) * result_set.excited_states.num_triplets
    
    # Check some SOC values.
    assert result_set.spin_orbit_coupling.between("S(0)", "T(1)") == pytest.approx(s0_t1)
    assert result_set.spin_orbit_coupling.between("S(1)", "T(1)") == pytest.approx(s1_t1)