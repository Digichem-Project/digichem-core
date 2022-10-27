"""Tests for translation classes."""
from silico.submit.translate import Solvent, Basis_set, Multiplicity, Functional
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
    assert all(["DiChloroMethane" == solvent.to_gaussian() for solvent in solvents])
    
    # And epsilon
    assert all([pytest.approx(8.93) == solvent.translate("epsilon") for solvent in solvents])
    

def test_631g_str_str():
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
    
    
def test_ccpvdz():
    assert Basis_set("cc-PvDz").to_gaussian() == "cc-pVDZ"
    assert str(Basis_set("cc-PvDz")) == "cc-pVDZ"
    

def test_b3lyp():
    """
    Test conversions for DFT functionals.
    """
    funcs = (Functional("b3-LYp"), Functional("B3lyp"))
    
    assert all(["b3-lyp" == func.to_turbomole() for func in funcs])
    
    assert all(["B3LYP" == func.to_gaussian() for func in funcs])


def test_pbe0gaussian():
    # Turbomole functionals are case sensitive. This is particularly annoying because functional
    # names are all lower-case, except for the word 'Gaussian'. Check this.
    assert Functional("PBE0 gaussian").to_turbomole() == "pbe0 Gaussian"


def test_multiplicity():
    """
    Test conversion for multiplicity.
    """
    mults = (Multiplicity(1), Multiplicity("1"), Multiplicity("singlet"), Multiplicity("sInGlEt"), Multiplicity("S"))
    
    # Check conversions.
    assert all(["Singlet" == mult.to_gaussian() for mult in mults])
    
    assert all([1 == mult.to_turbomole() for mult in mults])
    
    assert all(["S" == mult.symbol for mult in mults])
    
    assert all([1 == mult.number for mult in mults])
    
    assert all(["Singlet" == mult.string for mult in mults])
    