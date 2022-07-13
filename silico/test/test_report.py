"""Test report writing."""

import pytest
from pathlib import Path

from silico.report.pdf import PDF_report
from silico.test.util import silico_options
from silico.test.test_result import gaussian_opt_result, gaussian_radical_anion_result, gaussian_ES_result, gaussian_emission_result, \
                                    turbomole_opt_result, turbomole_ES_result, turbomole_ADC2_opt_result


@pytest.mark.parametrize("result_set", [pytest.lazy_fixture(fixture_name) for fixture_name in [
        "gaussian_opt_result", "gaussian_radical_anion_result", "gaussian_ES_result", "gaussian_emission_result", "turbomole_opt_result", "turbomole_ES_result", "turbomole_ADC2_opt_result"
    ]])
def test_report_writing(result_set, tmp_path, silico_options):
    # Create the report.
    report = PDF_report(result_set, options = silico_options)
    
    # And write.
    report.write(Path(tmp_path, "report.pdf"))
    
    # Check the pdf is there
    assert Path(tmp_path, "report.pdf").exists()
    