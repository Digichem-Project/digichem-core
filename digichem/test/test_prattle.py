"""Tests for functioning of oprattle"""
from pathlib import Path

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
def test_from_path(file_path):
    """
    Can we convert a file found on the filesystem?
    """
    Openprattle_converter(input_file_path = file_path).convert("com")


@pytest.mark.parametrize("file_path, mode", [
    (ethane_xyz, "r"),
    (benzene_cdx, "rb")
])
def test_from_file(file_path, mode):
    """
    Can we convert from an open file?
    """
    with open(file_path, mode) as input_file:
        Openprattle_converter(input_file = input_file, input_file_type = Path(file_path).suffix[1:]).convert("com")


@pytest.mark.parametrize("file_path, mode", [
    (ethane_xyz, "r"),
    (benzene_cdx, "rb")
])
def test_from_buffer(file_path, mode):
    """
    Can we convert from an open file?
    """
    with open(file_path, mode) as input_file:
        data = input_file.read()
        Openprattle_converter(input_file_buffer = data, input_file_type = Path(file_path).suffix[1:]).convert("com")


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