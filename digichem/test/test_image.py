import pytest
from pathlib import Path
import shutil

from digichem.image.excited_states import Excited_states_diagram_maker
from digichem.image.graph import Convergence_graph_maker
from digichem.image.orbitals import Orbital_diagram_maker
from digichem.image.spectroscopy import Absorption_graph_maker, Emission_graph_maker, \
    Frequency_graph_maker, NMR_graph_maker, NMR_graph_zoom_maker
from digichem.image.structure import Skeletal_image_maker
from digichem.image.vmd import Structure_image_maker, Orbital_image_maker, Density_image_maker, \
    Combined_orbital_image_maker, Permanent_dipole_image_maker, Transition_dipole_image_maker
from digichem.parse.util import open_for_parsing, parse_calculation

from digichem.test.util import digichem_options, data_directory
from digichem.test.test_result import gaussian_ES_result, turbomole_ES_result, orca_ES_result, \
    gaussian_opt_result, turbomole_opt_result, orca_opt_result, orca_opt_freq_result, \
    orca_nmr_result

HAS_VMD = shutil.which("vmd")

@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("gaussian_ES_result"),
    pytest.lazy_fixture("turbomole_ES_result"),
    pytest.lazy_fixture("orca_ES_result")
])
def test_es_diagram(result_set, tmp_path, digichem_options):
    """Can we make an excited states diagram?"""
    maker = Excited_states_diagram_maker.from_options(
        tmp_path / "tmp.png",
        excited_states = result_set.excited_states,
        ground_state = result_set.ground_state,
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("gaussian_opt_result"),
    pytest.lazy_fixture("turbomole_opt_result"),
    pytest.lazy_fixture("orca_opt_result")
])
def test_convergence_diagram(result_set, tmp_path, digichem_options):
    """Can we make an optimisation graph?"""
    maker = Convergence_graph_maker.from_options(
        tmp_path / "tmp.png",
        energies = result_set.energies.scf,
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("gaussian_opt_result"),
    pytest.lazy_fixture("turbomole_opt_result"),
    pytest.lazy_fixture("orca_opt_result")
])
def test_orbital_diagram(result_set, tmp_path, digichem_options):
    """Can we make an orbital diagram?"""
    maker = Orbital_diagram_maker.from_options(
        tmp_path / "tmp.png",
        orbitals = result_set.orbitals,
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("gaussian_ES_result"),
    pytest.lazy_fixture("turbomole_ES_result"),
    pytest.lazy_fixture("orca_ES_result")
])
def test_abs_diagram(result_set, tmp_path, digichem_options):
    """Can we make an absorption spectrum?"""
    maker = Absorption_graph_maker.from_options(
        tmp_path / "tmp.png",
        excited_states = result_set.excited_states,
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("gaussian_ES_result"),
    pytest.lazy_fixture("turbomole_ES_result"),
    pytest.lazy_fixture("orca_ES_result")
])
def test_pl_diagram(result_set, tmp_path, digichem_options):
    """Can we make an emission spectrum?"""
    maker = Emission_graph_maker.from_options(
        tmp_path / "tmp.png",
        excited_states = result_set.excited_states,
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("gaussian_opt_result"),
    pytest.lazy_fixture("turbomole_opt_result"),
    pytest.lazy_fixture("orca_opt_freq_result")
])
def test_freq_diagram(result_set, tmp_path, digichem_options):
    """Can we make an IR spectrum?"""
    maker = Frequency_graph_maker.from_options(
        tmp_path / "tmp.png",
        vibrations = result_set.vibrations,
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("orca_nmr_result"),
])
def test_nmr_diagram(result_set, tmp_path, digichem_options):
    """Can we make an NMR spectrum?"""
    # Get a suitable spectrometer object.
    spectrometer = result_set.nmr.spectrometer

    maker = NMR_graph_maker.from_options(
        tmp_path / "tmp.png",
        graph = spectrometer.get_graph(1, 1),
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("orca_nmr_result"),
])
def test_focus_nmr_diagram(result_set, tmp_path, digichem_options):
    """Can we make an NMR spectrum?"""
    # Get a suitable spectrometer object.
    spectrometer = result_set.nmr.spectrometer

    maker = NMR_graph_zoom_maker.from_options(
        tmp_path / "tmp.png",
        graph = spectrometer.get_graph(1, 1),
        focus = result_set.atoms.groups[(2, "H")],
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("gaussian_opt_result"),
    pytest.lazy_fixture("turbomole_opt_result"),
    pytest.lazy_fixture("orca_opt_freq_result")
])
def test_2d_diagram(result_set, tmp_path, digichem_options):
    """Can we make an IR spectrum?"""
    maker = Skeletal_image_maker.from_options(
        tmp_path / "tmp.png",
        atoms = result_set.atoms,
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.skipif(not HAS_VMD,
                    reason="No VMD available")
@pytest.mark.parametrize("cube_file", [
    Path(data_directory(), "Cubes/Pyridine.HOMO.cube")
], ids = ["Gaussian"])
@pytest.mark.parametrize("maker_cls", [
    Structure_image_maker,
    Orbital_image_maker,
])
def test_3d_image(cube_file, maker_cls, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    maker = maker_cls.from_options(
        tmp_path / "tmp.png",
        cube_file = cube_file,
        options = digichem_options
    )

    maker.get_file('x0y0z0')
    assert Path(tmp_path, "tmp.x0y0z0.png").exists()


# @pytest.mark.parametrize("result_path", [
#     Path(data_directory(), "Pyridine/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz")
# ], ids = ["Gaussian"])
# @pytest.mark.parametrize("density", [
#     "SCF"
# ])

@pytest.mark.skipif(not HAS_VMD,
                    reason="No VMD available")
@pytest.mark.parametrize("cube_file", [
    Path(data_directory(), "Cubes/Pyridine.SCF.cube")
], ids = ["Gaussian"])
def test_density_image(cube_file, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    maker = Density_image_maker.from_options(
        tmp_path / "tmp.png",
        cube_file = cube_file,
        options = digichem_options
    )

    maker.get_file('x0y0z0')
    assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.skipif(not HAS_VMD,
                    reason="No VMD available")
@pytest.mark.parametrize("cube_file", [
    Path(data_directory(), "Cubes/Benzene.anion.spin.cube")
], ids = ["Gaussian"])
def test_spin_density_image(cube_file, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    maker = Density_image_maker.from_options(
        tmp_path / "tmp.png",
        cube_file = cube_file,
        options = digichem_options
    )

    maker.get_file('x0y0z0')
    assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.skipif(not HAS_VMD,
                    reason="No VMD available")
@pytest.mark.parametrize("homo_cube, lumo_cube", [
    [Path(data_directory(), "Cubes/Pyridine.HOMO.cube"), Path(data_directory(), "Cubes/Pyridine.LUMO.cube")]
], ids = ["Gaussian"])
def test_combined_orbital_image(homo_cube, lumo_cube, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    maker = Combined_orbital_image_maker.from_options(
        tmp_path / "tmp.png",
        HOMO_cube_file = homo_cube,
        LUMO_cube_file = lumo_cube,
        options = digichem_options
    )

    maker.get_file('x0y0z0')
    assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.skipif(not HAS_VMD,
                    reason="No VMD available")
@pytest.mark.parametrize("cube_file", [
    Path(data_directory(), "Cubes/Benzene.anion.AHOMO.cube")
], ids = ["Gaussian"])
def test_unrestricted_orbital_image(cube_file, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    maker = Orbital_image_maker.from_options(
        tmp_path / "tmp.png",
        cube_file = cube_file,
        options = digichem_options
    )

    maker.get_file('x0y0z0')
    assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.skipif(not HAS_VMD,
                    reason="No VMD available")
@pytest.mark.parametrize("result_path, cube_file", [
    [  
        Path(data_directory(), "Pyridine/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"),
        Path(data_directory(), "Cubes/Pyridine.HOMO.cube"),
    ]
])
def test_pdm_image(result_path, cube_file, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    with open_for_parsing(result_path) as open_files:
        result_set = parse_calculation(*open_files, options = digichem_options)
        maker = Permanent_dipole_image_maker.from_options(
            tmp_path / "tmp.png",
            cube_file = cube_file,
            dipole_moment = result_set.dipole_moment,
            options = digichem_options
        )

        maker.get_file('x0y0z0')
        assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.skipif(not HAS_VMD,
                    reason="No VMD available")
@pytest.mark.parametrize("result_path, cube_file", [
    [
        Path(data_directory(), "Pyridine/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"),
        Path(data_directory(), "Cubes/Pyridine.HOMO.cube"),
    ]
])
def test_tdm_image(result_path, cube_file, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    with open_for_parsing(result_path) as open_files:
        result_set = parse_calculation(*open_files, options = digichem_options)
        maker = Transition_dipole_image_maker.from_options(
            tmp_path / "tmp.png",
            cube_file = cube_file,
            dipole_moment = result_set.excited_states.find("S(1)").transition_dipole_moment.electric,
            magnetic_dipole_moment = result_set.excited_states.find("S(1)").transition_dipole_moment.magnetic,
            options = digichem_options
        )

        maker.get_file('x0y0z0')
        assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.skipif(not HAS_VMD,
                    reason="No VMD available")
@pytest.mark.parametrize("cube_file", [
    Path(data_directory(), "Cubes/Pyridine.HOMO.cube"),
])
def test_nto_image(cube_file, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    maker = Orbital_image_maker.from_options(
        tmp_path / "tmp.png",
        cube_file = cube_file,
        options = digichem_options
    )

    maker.get_file('x0y0z0')
    assert Path(tmp_path, "tmp.x0y0z0.png").exists()