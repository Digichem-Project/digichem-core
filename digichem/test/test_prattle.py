"""Tests for functioning of oprattle"""

import pytest

from digichem.file.prattle import Oprattle_formats

@pytest.mark.parametrize("readwrite", ["read", "write"])
def test_formats(readwrite):
    if readwrite == "read":
        assert len(Oprattle_formats().read()) > 0
    
    else:
        assert len(Oprattle_formats().write()) > 0