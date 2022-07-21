"""Test report writing."""

import pytest
from pathlib import Path
import shutil

from silico.report.pdf import PDF_report
from silico.test.util import silico_options, data_directory
from silico.parse.util import open_for_parsing, parse_calculation


def make_report(result_path, tmp_path, silico_options):
    """Helper function which unpacks a (number of) data archives and writes a report from them."""
    with open_for_parsing(Path(data_directory(), result_path)) as log_files:
        result_set = parse_calculation(*log_files)
        report = PDF_report(result_set, options = silico_options)
     
        # And write.
        report.write(Path(tmp_path, "report.pdf"))
        return report


@pytest.mark.slow
@pytest.mark.parametrize("result_path", [
        "Benzene Anion/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Gas Phase 6-31G(d,p).tar.gz"
     ])
def test_report_complete(result_path, tmp_path, silico_options):
    """Test creating a report from scratch"""
    report = make_report(result_path, tmp_path, silico_options)
         
    # Check the pdf is there
    assert Path(tmp_path, "report.pdf").exists()
     
    # Check certain images are there.
    for file_path in report.images['structure'].file_path.values():
        assert file_path.exists()
            

@pytest.mark.slow
@pytest.mark.parametrize("result_path", [
        [
            "Pyridine/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz",
            "Pyridine/Gaussian 16 Excited States TDA 10 Singlets 10 Triplets PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz",
            "Pyridine/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz",
        ],
        [
            "Pyridine/Turbomole Optimisation Frequency PBE0 (GD3BJ) 6-31G**.tar.gz"
            "Pyridine/Turbomole Excited States TDA 10 Triplets PBE0 (GD3BJ) 6-31G**.tar.gz",
            "Pyridine/Turbomole Excited States TDA 10 Singlets PBE0 (GD3BJ) 6-31G**.tar.gz",
        ]
     ])
def test_report_merged(result_path, tmp_path, silico_options):
    """Test creating a report from scratch"""
    report = make_report(result_path, tmp_path, silico_options)
         
    # Check the pdf is there
    assert Path(tmp_path, "report.pdf").exists()
     
    # Check certain images are there.
    for file_path in report.images['structure'].file_path.values():
        assert file_path.exists()


@pytest.mark.parametrize("result_path", [
        "Pyridine/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz",
        "Benzene Anion/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Gas Phase 6-31G(d,p).tar.gz",
        "Pyridine/Turbomole Optimisation Frequency PBE0 (GD3BJ) 6-31G**.tar.gz",
        "Pyridine/Turbomole Excited States TDA 10 Triplets PBE0 (GD3BJ) 6-31G**.tar.gz",
        "Pyridine/Turbomole Excited States TDA 10 Singlets PBE0 (GD3BJ) 6-31G**.tar.gz"

    ])
def test_pdf_writing(result_path, tmp_path, silico_options):
    """Test writing a report from existing image files."""
    # Unpack the archive.
    with open_for_parsing(Path(data_directory(), result_path)) as log_files:
        result_set = parse_calculation(*log_files)
        report = PDF_report(result_set, options = silico_options)
    
        # The location to write the pdf to.
        output = Path(tmp_path, "Report/Report.pdf")
    
        # Copy the existing report to the temp dir.
        shutil.copytree(Path(result_set.metadata.log_files[0].parents[1], "Report"), output.parent,  copy_function = shutil.copy)
    
        # Write the report, using the newly copied directory as output.
        report.write(output)
    
    