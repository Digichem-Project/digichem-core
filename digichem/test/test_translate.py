"""Tests for translation classes."""
from silico.submit.translate import Solvent, Basis_set
import pytest


def test_solvent():
    """
    Test conversions for DCM.
    """
    # Get some different ways to represent DCM.
    solvents = (Solvent("dcm"), Solvent("DCM"), Solvent("dichloroMethane"), Solvent("Dichloro Methane"))
    
    # Check all have equivalent string representations.
    assert all([str(solvents[0]) == str(solvent) for solvent in solvents])
    
    # And the correct gaussian name.
    assert all(["Dichloromethane" == solvent.to_gaussian() for solvent in solvents])
    
    # And epsilon
    assert all([pytest.approx(8.93) == solvent.translate("epsilon") for solvent in solvents])
    

def test_basis_set():
    """
    Test conversions for basis sets.
    """
    # Test both types of 6-31G**.
    basis_sets = (Basis_set("6-31g(d,p)"), Basis_set("6-31G**"))
    
    # Check both are equivalent.
    assert str(basis_sets[0]) == str(basis_sets[1])
    
    # Check Gaussian version.
    assert all(["6-31G**" == basis_set.to_gaussian() for basis_set in basis_sets])
    
    # Check Turbomole version.
    assert all(["6-31G**" == basis_set.to_turbomole() for basis_set in basis_sets])