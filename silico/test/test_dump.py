"""Tests for resulting dumping (ie, output from `silico result`)."""

import pytest

from silico.test.util import silico_options
from silico.test.test_result import gaussian_opt_result, gaussian_ES_result
from silico.result.format.filter import Result_filter


# Tests for checking the filtering mechanism is working correctly.
@pytest.mark.parametrize("filter_string", [
        # Basic filtering tests.
        "atoms",
        "energies",
        "orbitals",
        # Nested filtering of single items (energies:scf:final, orbitals:ΔE(HOMO-LUMO) etc)
        "atoms:x-extension",
        "energies:scf:final",
        "orbitals:ΔE(HOMO-LUMO)",
        # Nested filtering of items from lists (energies:scf:0, energies:scf:values:0, orbitals:0, orbitals:values:0)
        "atoms:0",
        "energies:scf:0",
        "energies:scf:values:0",
        "orbitals:-1",
        # Nested filtering of single items from lists (energies:scf:0:units, orbitals:values:0:symmetry)
        "atoms:0:element",
        "energies:scf:0:units",
        "orbitals:values:0:symmetry",
    ])
def test_filter(gaussian_opt_result, filter_string, silico_options):
    """Basic filtering tests."""
    Result_filter(filter_string, silico_options).filter(gaussian_opt_result)