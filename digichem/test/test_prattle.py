"""Tests for functioning of oprattle"""

import pytest

from digichem.file.prattle import Oprattle_formats, Openprattle_converter

from digichem.test.util import ethane_xyz, benzene_cdx

@pytest.mark.parametrize("readwrite", ["read", "write"])
def test_formats(readwrite):
    if readwrite == "read":
        assert len(Oprattle_formats().read()) > 0
    
    else:
        assert len(Oprattle_formats().write()) > 0

@pytest.mark.parametrize("file_path", [
    ethane_xyz,
    benzene_cdx
])
def test_convert(file_path):
    Openprattle_converter(input_file_path = file_path).convert("com")


@pytest.mark.parametrize("file_path", [
    ethane_xyz,
    benzene_cdx
])
def test_gen3d(file_path):
    first = Openprattle_converter(input_file_path = file_path).convert("xyz", gen3D = False)

    same = Openprattle_converter(input_file_path = file_path).convert("xyz", gen3D = False)
    different = Openprattle_converter(input_file_path = file_path).convert("xyz", gen3D = True)

    assert first == same
    assert first != different