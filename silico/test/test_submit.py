"""Tests for calculation submission."""

import pytest

from silico.test.util import benzene_cdx, pyridine_cml, cyclopentane_com, silico_options
from silico.input.silico import si_from_file
from silico.submit.calculation.base import Calculation_target

# The destination to use.
destination = "Series"

@pytest.mark.slow
@pytest.mark.parametrize("coordinate_files, method_codes", [
    # Multiple Submission.
    ([benzene_cdx, pyridine_cml, cyclopentane_com], [f"{destination}/Turbomole/[Turbomole Single Point, Turbomole PBE0 (GD3BJ), Turbomole 6-31G**]"]),
    # Gaussian organics.
    ([benzene_cdx], [f"{destination}/Gaussian 16/Gaussian Auto Organic TDA Emission"]),
    # Gaussian organometallics.
    ([benzene_cdx], [f"{destination}/Gaussian 16/Gaussian Auto Organometallic TD-DFT Unrestricted Triplet"]),
    # Turbomole.
    ([benzene_cdx], [f"{destination}/Turbomole/Turbomole Auto ADC(2)"]),
])
def test_submit(coordinate_files, method_codes, tmp_path, silico_options):
    
    # Resolve methods.
    methods = [silico_options.methods.resolve_method_string(method_id) for method_id in method_codes]
    
    # Load files.
    coordinates = [si_from_file(coordinate_file) for coordinate_file in coordinate_files]
    
    # The first method.
    first = Calculation_target.link(methods, global_silico_options = silico_options)
    
    # Submit.
    for coordinate in coordinates:
        # Prepare.
        first.prepare(tmp_path, coordinate)
        
        # Go.
        first.submit()
    