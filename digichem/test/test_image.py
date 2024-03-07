import pytest
from pathlib import Path

from digichem.image.excited_states import Excited_states_diagram_maker
from digichem.image.graph import Convergence_graph_maker
from digichem.image.orbitals import Orbital_diagram_maker
from digichem.image.spectroscopy import Absorption_graph_maker, Emission_graph_maker

from digichem.test.util import digichem_options
from digichem.test.test_result import gaussian_ES_result, turbomole_ES_result, orca_ES_result, \
    gaussian_opt_result, turbomole_opt_result, orca_opt_result


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