"""Tests for resulting dumping (ie, output from `silico result`)."""

import pytest
import yaml

from silico.test.util import silico_options
from silico.test.test_result import gaussian_opt_result, gaussian_ES_result
from silico.result.format.filter import Result_filter
from silico.result.format.yaml import Yaml_dumper


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


@pytest.mark.parametrize("result_set, filter_string, value", [
        (pytest.lazy_fixture("gaussian_opt_result"), "energies:scf:final:value", -10488.990333747),
    ])
def test_numeric_result(result_set, filter_string, value, silico_options):
    """Dumping of numeric values."""
    # Dump the result to text (YAML).
    data = Yaml_dumper.from_text(filter_string, silico_options = silico_options).dump(result_set)
    
    # Parse with yaml.
    parsed = yaml.safe_load(data)
    
    # Check the value
    assert list(parsed.items())[0][1] == pytest.approx(value, abs=1e-5)