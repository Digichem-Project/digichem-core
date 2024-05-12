"""Tests for result parsing and processing"""

from pathlib import Path
import itertools
import numpy
import pytest
import yaml

from digichem.parse import parse_calculation, parse_multiple_calculations, parse_and_merge_calculations
from digichem.test.util import data_directory, result_files, digichem_options
from digichem.result import Result_set

@pytest.mark.parametrize("result_data", list(itertools.chain(*list(result_files.values()))))
def test_parsing(result_data, digichem_options):
    """Test whether we can parse various calc results."""
    result = parse_calculation(result_data, options = digichem_options)
    assert isinstance(result, Result_set)

        
def test_multi_parsing(digichem_options):
    """Test whether we can parse in parallel."""
    
    results = parse_multiple_calculations(*list(itertools.chain(*list(result_files.values()))), options = digichem_options)
    
    # Check length.
    assert len(results) == len(list(itertools.chain(*list(result_files.values()))))


@pytest.mark.parametrize("data_set", list(result_files.values()))
def test_merged_parsing(data_set, digichem_options):
    """Test whether we can parse and merge calc results."""
    result = parse_and_merge_calculations(*data_set, options = digichem_options)
    
    assert isinstance(result, Result_set)
    assert len(result.results) == len(data_set)


@pytest.mark.parametrize(
    "result_files",
    itertools.chain(
        [(result_file,) for program_files in result_files.values() for result_file in program_files],
        [(
            Path(data_directory(), "Naphthalene/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p)"),
            Path(data_directory(), "Naphthalene/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz")
        )]
    )
)
def test_dump_and_parse(result_files, tmp_path, digichem_options):
    """Test whether we can dump and reload a result to get the same data back again."""
    # Set numpy errors (not sure why this isn't the default...)
    numpy.seterr(invalid = 'raise', divide = 'raise')
    
    # Parse the raw data.
    if len(result_files) == 1:
        raw = parse_calculation(result_files[0], options = digichem_options)
    
    else:
        raw = parse_and_merge_calculations(*result_files, options = digichem_options)
    
    # Instead of just dumping to memory we'll dump to a file to test the full mechanism.
    with open(Path(tmp_path, "dump.sir"), "wt") as dump_file:
        dump_file.write(yaml.safe_dump(raw.dump(digichem_options = digichem_options)))
    
    # Now load the result back again.
    parsed = parse_calculation(Path(tmp_path, "dump.sir"), options = digichem_options)
    
    # Dump both the raw and parsed result sets and compare to make sure they're the same.
    raw_dump = raw.dump(digichem_options)
    parsed_dump = parsed.dump(digichem_options)
    
    # Aux files are not preserved, so don't compare these.
    raw_dump['metadata'].pop("auxiliary_files")
    raw_dump['metadata'].pop("log_files")
    parsed_dump['metadata'].pop("auxiliary_files")
    parsed_dump['metadata'].pop("log_files")
    
    
    assert raw_dump == parsed_dump
        