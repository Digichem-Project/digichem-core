"""Common testing utilities and convenience functions."""

import pytest


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