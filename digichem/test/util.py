"""Common testing utilities and convenience functions."""

import pytest
from pathlib import Path
import pkg_resources

from silico.config.parse import Config_file_parser


def data_directory():
    """
    Get a path to the test data directory.
    """
    return Path(pkg_resources.resource_filename('silico', 'test/data'))


# Input files.
benzene_cdx = Path(data_directory(), "Input", "Benzene.cdx")
pyridine_cml = Path(data_directory(), "Input", "Pyridine.cml")
cyclopentane_com = Path(data_directory(), "Input", "Cyclopentane.com")
water_xyz = Path(data_directory(), "Input", "Water.xyz")
ethane_xyz = Path(data_directory(), "Input", "Ethane.xyz")
pyridine_si_v2 = Path(data_directory(), "Input/Pyridine.v2.si")
pyridine_si_v1 = Path(data_directory(), "Input/Pyridine.v1.si")
pyridine_resume_pickle = Path(data_directory(), "Input/silico.resume.pickle")


# Result files for testing parsing.
result_files = {
    "gaussian": [Path(data_directory(), datum) for datum in [
        'Naphthalene/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p)',
        'Naphthalene/Gaussian 16 Excited States TDA 10 Singlets 10 Triplets PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz',
        'Naphthalene/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz',  
    ]],
    "turbomole": [Path(data_directory(), datum) for datum in [
        'Naphthalene/Turbomole Optimisation ADC(2) cc-pVDZ.tar.gz',
        'Naphthalene/Turbomole Excited States ADC(2) S(1) and S(2) cc-pVDZ.tar.gz',
        'Naphthalene/Turbomole Excited States ADC(2) T(1) and T(2) cc-pVDZ.tar.gz',
    ]]
}

@pytest.fixture(scope="package")
def silico_options(tmpdir_factory):
    config = Config_file_parser.silico_options()
    # Change the default DB location to prevent SPAM.
    #config['database']['db_type']
    config['database']['path'] = Path(tmpdir_factory.mktemp("silico_database"), "tmp.db")
    return config


def check_float_list(test_list, expected_list, abs = 1e-4):
    """
    Helper function to compare to lists of floats.
    """
    # Pytest bug #9921 prevents this comparison from working for now...
    #assert test_list == pytest.approx(expected_list, abs=abs)
    for list_index in range(0,len(expected_list)):
        assert test_list[list_index] == pytest.approx(expected_list[list_index], abs=abs)


def check_dipole(dipole_moment, coords, abs = 1e-4):
    """
    Helper function to check a dipole moment (PDM, TEDM or TMDM) matches some expected values.
    """
    check_float_list(dipole_moment.vector_coords, coords, abs)
        
    # Check total
    assert dipole_moment.total == pytest.approx((coords[0] **2 + coords[1] **2 + coords[2] **2) **0.5, abs = abs)


def check_orbitals(orbitals, num_occ, num_unocc, homo, lumo):
    """Helper function for checking orbital energies"""
    # Check numbers.
    assert len(orbitals) == num_occ + num_unocc
    assert len(orbitals.occupied) == num_occ
    assert len(orbitals.virtual) == num_unocc
    
    # Check energies.
    assert orbitals.HOMO_energy == pytest.approx(homo)
    assert orbitals.LUMO_energy == pytest.approx(lumo)
    assert orbitals.HOMO_LUMO_energy == pytest.approx(lumo - homo)
    
    # Check ordering.
    assert [orbital.level for orbital in orbitals] == list(range(1, num_occ + num_unocc +1))


def check_rendered_image(base_name):
    """Check whether the eight separate images of a given rendered image have been created."""
    for extension in ["jpg", "png"]:
        for orientation in ["x0y0z0", "x0y90z0", "x45y45z45", "x90y0z0"]:
            assert Path(base_name.parent, f"{base_name.name}.{orientation}.{extension}").exists()


def check_report(report_folder, base_name, *,
    optimisation = False,
    vibrations = False,
    excited_states = False,
    absorption = False,
    s1_tdm = False,
    s1_diff = False,
    t1_diff = False,
    s1_nto = False,
    t1_nto = False,
    vertical_emission = False,
    adiabatic_emission = False):
    """
    Check that certain report and image files have been correctly rendered.
    """
    # First, check the report file itself is there
    assert (Path(report_folder, f"{base_name}.journal.pdf").exists() and Path(report_folder, f"{base_name}.traditional.pdf").exists()) or Path(report_folder, f"{base_name}.pdf").exists()
    
    # Now check for images.
    image_folder = Path(report_folder, "image")
    # HOMO and LUMO
    for image in ["HOMO", "LUMO"]:
        check_rendered_image(Path(image_folder, f"{image}", f"{base_name}.{image}"))
    
    # HOMO-LUMO overlap
    check_rendered_image(Path(image_folder, "HOMO LUMO", f"{base_name}.HOMO_LUMO"))
    
    # Structure and density.
    check_rendered_image(Path(image_folder, "Structure", f"{base_name}.structure"))
    check_rendered_image(Path(image_folder, "Density", f"{base_name}.SCF"))
    
    # Dipole.
    check_rendered_image(Path(image_folder, "Dipole Moment", f"{base_name}.dipole"))
    
    if optimisation:
        assert Path(image_folder, f"{base_name}.SCF_graph.png").exists()
    
    # HOMO-LUMO diagrams.
    assert Path(image_folder, "Orbital Diagram", f"{base_name}.HOMO_LUMO.png").exists()
    assert Path(image_folder, "Orbital Diagram", f"{base_name}.orbitals.png").exists()
    
    # Excited state stuff.
    if s1_tdm:
        check_rendered_image(Path(image_folder, "S(1) Transition Dipole Moment", f"{base_name}.S(1)_dipole"))
    
    if s1_diff:
        check_rendered_image(Path(image_folder, "S(1)", f"{base_name}.S(1)_difference_density"))
    
    if t1_diff:
        check_rendered_image(Path(image_folder, "T(1)", f"{base_name}.T(1)_difference_density"))
        
    if s1_nto:
        check_rendered_image(Path(image_folder, "S(1)", f"{base_name}.S(1)_NTO"))
    
    if t1_nto:
        check_rendered_image(Path(image_folder, "T(1)", f"{base_name}.T(1)_NTO"))
        
    if excited_states:
        assert Path(image_folder, f"{base_name}.excited_states.png").exists()
    
    if absorption:
        assert Path(image_folder, f"{base_name}.simulated_absorption_spectrum.png").exists()
        
    if vertical_emission:
        assert Path(image_folder, f"{base_name}.vertical_S(1)_emission_states.png").exists()
        assert Path(image_folder, f"{base_name}.simulated_vertical_S(1)_emission_spectrum.png").exists()
    
    if adiabatic_emission:
        assert Path(image_folder, f"{base_name}.adiabatic_S(1)_emission_states.png").exists()
        assert Path(image_folder, f"{base_name}.simulated_adiabatic_S(1)_emission_spectrum.png").exists()
    
    # Vibrations.
    if vibrations:
        assert Path(image_folder, f"{base_name}.simulated_frequencies.png").exists()