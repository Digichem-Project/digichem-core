"""Tests for result parsing and processing"""

from pathlib import Path
import itertools
import numpy
import pytest

from silico.parse import parse_calculation
from silico.test.util import data_directory, result_files, silico_options
from silico.result import Result_set
from silico.parse import parse_multiple_calculations, parse_and_merge_calculations
from silico.result.format.yaml import Yaml_dumper

@pytest.mark.parametrize("result_data", list(itertools.chain(result_files['gaussian'], result_files['turbomole'])))
def test_parsing(result_data):
    """Test whether we can parse various calc results."""
    result = parse_calculation(result_data)
    assert isinstance(result, Result_set)

        
def test_multi_parsing():
    """Test whether we can parse in parallel."""
    
    results = parse_multiple_calculations(*list(itertools.chain(result_files['gaussian'], result_files['turbomole'])))
    
    # Check length.
    assert len(results) == len(result_files['gaussian']) + len(result_files['turbomole'])


@pytest.mark.parametrize("data_set", [result_files['gaussian'], result_files['turbomole']])
def test_merged_parsing(data_set):
    """Test whether we can parse and merge calc results."""
    result = parse_and_merge_calculations(*data_set)
    
    assert isinstance(result, Result_set)
    assert len(result.results) == 3


@pytest.mark.parametrize("result_file", [result_file for program_files in result_files.values() for result_file in program_files])
def test_dump_and_parse(result_file, tmp_path, silico_options):
    """Test whether we can dump and reload a result to get the same data back again."""
    # Set numpy errors (not sure why this isn't the default...)
    numpy.seterr(invalid = 'raise', divide = 'raise')
    
    # Parse the raw data.
    raw = parse_calculation(result_file)
    
    # Instead of just dumping to memory we'll dump to a file to test the full mechanism.
    with open(Path(tmp_path, "dump.sir"), "wt") as dump_file:
        dump_file.write(Yaml_dumper(silico_options = silico_options).dump(raw))
    
    # Now load the result back again.
    parsed = parse_calculation(Path(tmp_path, "dump.sir"))
    
    # Dump both the raw and parsed result sets and compare to make sure they're the same.
    assert raw.dump(silico_options) == parsed.dump(silico_options)
    