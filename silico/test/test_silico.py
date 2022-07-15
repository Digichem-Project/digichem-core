#!/usr/bin/env python3

import pytest
from pytest_lazyfixture import lazy_fixture

if __name__ == '__main__':
    pytest.main([], plugins = [lazy_fixture])