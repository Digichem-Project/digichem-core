"""Test report writing."""

import pytest
from pathlib import Path
import shutil
import time

from silico.report.pdf import PDF_report
from silico.test.util import silico_options
from silico.test.test_result import gaussian_opt_result, gaussian_radical_anion_result, gaussian_ES_result, gaussian_emission_result, \
                                    turbomole_opt_result, turbomole_ES_result, turbomole_ADC2_opt_result, gaussian_PDM_ES_result, \
                                    turbomole_PDM_result, turbomole_TDM_result

@pytest.mark.skip("Report writing from archives is currently broken")
@pytest.mark.slow
@pytest.mark.parametrize("result_set", [pytest.lazy_fixture(fixture_name) for fixture_name in [
        "gaussian_opt_result",
        "gaussian_radical_anion_result",
        "gaussian_ES_result",
        "gaussian_emission_result",
        "turbomole_opt_result",
        "turbomole_ES_result",
        "turbomole_ADC2_opt_result"
    ]])
def test_report_complete(result_set, tmp_path, silico_options):
    """Test creating a report from scratch"""
    # Create the report.
    report = PDF_report(result_set, options = silico_options)
    
    # And write.
    report.write(Path(tmp_path, "report.pdf"))
    
    # Check the pdf is there
    assert Path(tmp_path, "report.pdf").exists()
    
    # Check certain images are there.
    for file_path in report.images['structure'].file_path.values():
        assert file_path.exists()


@pytest.mark.skip("Report writing from archives is currently broken")
@pytest.mark.parametrize("result_set", [pytest.lazy_fixture(fixture_name) for fixture_name in [
        "turbomole_PDM_result",
#        "gaussian_PDM_ES_result",
#         "turbomole_TDM_result",
    ]])
def test_pdf_writing(result_set, tmp_path, silico_options):
    """Test writing a report from existing image files."""
    
    
    time.sleep(20)
    # Copy the existing report to the temp dir.
    print(result_set.metadata.log_files[0])
    print(result_set.metadata.log_files[0].exists())
    shutil.copytree(Path(result_set.metadata.log_files[0].parents[1], "Report"), Path(tmp_path, "Report"),  copy_function = shutil.copy)
    
    
    