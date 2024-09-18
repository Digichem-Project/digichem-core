"""Tests for checking the molecular alignment procedures"""
from pathlib import Path
import pytest

from digichem.parse import parse_calculation
from digichem.test.util import data_directory, result_files, digichem_options

@pytest.mark.parametrize(
    "method",
    ["MIN", "FAP", "AA", "AAA", "GRID", "NEST"]
)
def test_alignment_method(method, digichem_options):
    result = parse_calculation(Path(data_directory(), "Naphthalene/Gaussian 16 Single Point (Singlet) PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"), options = digichem_options, ornt = method)

    # Check axes are in the correct order.
    assert result.atoms.X_length >= result.atoms.Y_length
    assert result.atoms.Y_length >= result.atoms.Z_length