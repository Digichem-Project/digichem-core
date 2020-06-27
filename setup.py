#!/usr/bin/env python3

from setuptools import setup
from setuptools import find_packages
from pathlib import Path

setup(
	name = "silico",
	version = "0.6",
	description = "Silico Computational Chemistry Package",
	author = "Oliver Lee",

	# These 'paths' are relative to this script.
	packages = find_packages('src'),
	package_dir = {'': 'src'},
	package_data = {'silico': [str(path.relative_to('src/silico')) for path in Path('src/silico/data').glob('**/*') if path.is_file()]},
	#data_files = [("", [str(path) for path in Path('data').glob('**/*') if path.is_file()])],
	scripts = ['src/silico/program/creport.py', 'src/silico/program/cresult.py'],
	
)
