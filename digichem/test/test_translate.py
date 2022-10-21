"""Tests for translation classes."""
from silico.submit.translate import Solvent
import pytest


def test_solvent():
    """
    Test conversion for DCM.
    """
    # Get some different ways to represent DCM.
    solvents = (Solvent("dcm"), Solvent("DCM"), Solvent("dichloroMethane"), Solvent("Dichloro Methane"))
    
    # Check all have equivalent string representations.
    assert all([str(solvents[0]) == str(solvent) for solvent in solvents])
    
    # And the correct gaussian name.
    assert all(["Dichloromethane" == solvent.to_gaussian() for solvent in solvents])
    
    # And epsilon
    assert all([pytest.approx(8.93) == solvent.translate("epsilon") for solvent in solvents])