"""
Tests for calculation submission.

These tests run actual calculations. This requires a significant number of computational resources and time.
They should probably only be run from inside a calculation server environment.
"""

import pytest
import itertools
import os

from silico.test.util import benzene_cdx, pyridine_cml, cyclopentane_com, water_xyz, ethane_xyz, silico_options
from silico.input.silico import si_from_file
from silico.submit.calculation.base import Calculation_target

# The destination to use.
destination = "Series"

# Basis sets to test.
turbomole_basis_sets = [
    "6-31G*",
    "6-31G**",
    "6-31+G*",
    "6-31+G**",
    "6-31++G*",
    "6-31++G**",
    "6-311G*",
    "6-311G**",
    "6-311+G*",
    "6-311+G**",
    "6-311++G*",
    "6-311++G**",
    "cc-pVDZ",
    "cc-pVTZ",
    "cc-pVQZ",
    "cc-pV5Z",
    "cc-pV6Z",
]
gaussian_basis_sets = list(itertools.chain(turbomole_basis_sets, [
    "6-31+G** SBKJC-VDZ (ECP)",
    "6-31+G** LANL2DZ (ECP)",
]))

# Solvents to test.
gaussian_solvents = [
    "Gas Phase",
    "1,4-Dioxane",
    "2-Propanol",
    "Acetone",
    "Acetonitrile",
    "Benzene",
    "CycloHexane",
    "Dichloromethane",
    "DiethylEther",
    "n,n-DiMethylFormamide",
    "DiMethylSulfoxide",
    "Ethanol",
    "EthylEthanoate",
    "Methanol",
    "Heptane",
    "n-Hexane",
    "TetraHydroFuran",
    "Toluene",
    "Water"
]
turbomole_solvents = list(gaussian_solvents)

# Methods to test.
gaussian_methods = [
    "PBE0 (GD3BJ)",
    "B3LYP (GD3BJ)",
    "CAM-B3LYP (GD3BJ)",
    "M062X (GD3)",
    "wB97XD",
    "MP2",
    "MP3",
    "MP4",
    "MP5",
    "CCD",
    "CCSD"
]
turbomole_methods = [
    "Standard DFT: PBE0 (GD3BJ)",
    "Standard DFT: B3LYP (GD3BJ)",
    "Standard DFT: M062X (GD3)",
    "RI-DFT: PBE0 (GD3BJ)",
    "RI-DFT: B3LYP (GD3BJ)",
    "RI-DFT: M062X (GD3)",
#     "RIJK-DFT: PBE0 (GD3BJ)",
#     "RIJK-DFT: B3LYP (GD3BJ)",
#     "RIJK-DFT: M062X (GD3)",
]
turbomole_ri_methods = [
    "RI-MP2",
    "RI-MP3",
    "RI-MP4",
    "SCS-ADC(2)",
    "SCS-CC2",
    "RI-CCSD",
    "RI-CCSD(T)"
]

# Calculations/properties to test.
gaussian_properties = [
    "Single Point Singlet",
    "Gradient",
    "Optimisation",
    "Optimisation Frequencies",
    "Optimisation Frequencies Unrestricted Triplet",
    "Frequencies",
    "Excited States: TDA: 10 Singlets 10 Triplets",
    "Excited States: TDA: 5 Triplets",
    "Excited States: TDA: 100 Singlets",
    "Excited States: TDA: Optimised S(1)",
    "Excited States: TDA: Optimised T(1)",
    "Excited States: TD-DFT: 10 Singlets 10 Triplets",
    "Excited States: TD-DFT: 5 Triplets",
    "Excited States: TD-DFT: 100 Singlets",
    "Excited States: TD-DFT: Optimised S(1)",
    "Excited States: TD-DFT: Optimised T(1)",
]
turbomole_properties = [
    "Single Point",
    "Gradient",
    "Optimisation",
    "Optimisation Frequencies",
    "Frequencies",
    "Excited States: TDA: 10 Singlets",
    "Excited States: TDA: 10 Triplets",
    "Excited States: TDA: 100 Singlets",
    "Excited States: TDA: Optimised S(1)",
    "Excited States: TDA: Optimised T(1)",
    "Excited States: TD-DFT: 10 Singlets",
    "Excited States: TD-DFT: 10 Triplets",
    "Excited States: TD-DFT: 100 Singlets",
    "Excited States: TD-DFT: Optimised S(1)",
    "Excited States: TD-DFT: Optimised T(1)",
]
turbomole_posthf_properties = [
    "Excited States: Post-HF Excited States: Post-HF ES Singlet",
    "Excited States: Post-HF Excited States: Post-HF ES Triplet",
    "Excited States: Post-HF Excited States: Post-HF Optimised S(1)",
    "Excited States: Post-HF Excited States: Post-HF Optimised T(1)",
]



def run_submission_test(coordinate_files, method_codes, tmp_path, silico_options):
    """Run a submission test."""
    
    # Resolve methods.
    methods = [silico_options.methods.resolve_method_string(method_code) for method_code in method_codes]
    
    # We expect to be run in a SLURM context.
    # Get our available resources from the context.
    try:
        num_cpus = int(os.environ['SLURM_CPUS_ON_NODE'])
        # Weirdly, this (and other memory related) variable is not set by srun..
        #memory = int(os.environ['SLURM_MEM_PER_CPU']) * num_cpus
        # 2GB per CPU
        memory = 2 *1000 *1000 *1000 * num_cpus
        
    except KeyError:
        # Not running in SLURM.
        pytest.skip("Not running in a SLURM context")
    
    # Change necessary params.
    for method in methods:
        calc = method[2]
        # Set performance options.
        calc.performance['num_cpu'] = num_cpus
        # in MB
        calc.performance['memory'] = memory
        
        # Disable reports
        calc.post_process['write_report'] = False
        
        # Validate to set types.
        calc.validate()
    
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


@pytest.mark.slow
@pytest.mark.parametrize("coordinate_files, method_codes", [
    # Multiple Submission.
    ([benzene_cdx, pyridine_cml, cyclopentane_com], [f"{destination}/Turbomole/Turbomole:: Single Point: Standard DFT: PBE0 (GD3BJ): Gas Phase: 6-31G**"]),
    # Gaussian organics.
    ([pyridine_cml], [f"{destination}/Gaussian 16/Gaussian Auto Organic TDA Emission"]),
    # Gaussian organometallics.
    ([ethane_xyz], [f"{destination}/Gaussian 16/Gaussian Auto Organometallic TD-DFT Unrestricted Triplet"]),
    # Turbomole.
    ([ethane_xyz], [f"{destination}/Turbomole/Turbomole Auto SCS-ADC(2)"]),
])
def test_submit(coordinate_files, method_codes, tmp_path, silico_options):
    """Test some specific calcs"""
    # Run the test.
    run_submission_test(coordinate_files, method_codes, tmp_path, silico_options)


@pytest.mark.slow
@pytest.mark.parametrize("basis_set", gaussian_basis_sets)
def test_gaussian_basis(basis_set, tmp_path, silico_options):
    """Test different gaussian basis sets."""
    # Get our method code.
    method_code = f"{destination}/Gaussian 16/Gaussian:: Single Point Singlet: PBE0 (GD3BJ): Gas Phase: {basis_set}"
    
    # Run the test.
    run_submission_test([ethane_xyz], [method_code], tmp_path, silico_options)


@pytest.mark.slow
@pytest.mark.parametrize("basis_set", turbomole_basis_sets)
def test_turbomole_basis(basis_set, tmp_path, silico_options):
    """Test different gaussian basis sets."""
    # Get our method code.
    method_code = f"{destination}/Turbomole/Turbomole:: Single Point: Standard DFT: PBE0 (GD3BJ): Gas Phase: {basis_set}"
    
    # Run the test.
    run_submission_test([water_xyz], [method_code], tmp_path, silico_options)


@pytest.mark.slow
@pytest.mark.parametrize("solvent", gaussian_solvents)
def test_gaussian_solvents(solvent, tmp_path, silico_options):
    """Test different gaussian solvents."""
    # Get our method code.
    method_code = f"{destination}/Gaussian 16/Gaussian:: Single Point Singlet: PBE0 (GD3BJ): {solvent}: 6-31G*"
    
    # Run the test.
    run_submission_test([ethane_xyz], [method_code], tmp_path, silico_options)


@pytest.mark.slow
@pytest.mark.parametrize("solvent", turbomole_solvents)
def test_turbomole_solvents(solvent, tmp_path, silico_options):
    """Test different gaussian solvents."""
    # Get our method code.
    method_code = f"{destination}/Turbomole/Turbomole:: Single Point: Standard DFT: PBE0 (GD3BJ): {solvent}: 6-31G*"
    
    # Run the test.
    run_submission_test([water_xyz], [method_code], tmp_path, silico_options)


@pytest.mark.slow
@pytest.mark.parametrize("method", gaussian_methods)
def test_gaussian_methods(method, tmp_path, silico_options):
    """Test different gaussian methods."""
    # Get our method code.
    method_code = f"{destination}/Gaussian 16/Gaussian:: Single Point Singlet: {method}: Gas Phase: 6-31G*"
    
    # Run the test.
    run_submission_test([ethane_xyz], [method_code], tmp_path, silico_options)


@pytest.mark.slow
@pytest.mark.parametrize("method", turbomole_methods)
def test_turbomole_methods(method, tmp_path, silico_options):
    """Test different turbomole methods."""
    # Get our method code.
    method_code = f"{destination}/Turbomole/Turbomole:: Single Point: {method}: Gas Phase: cc-pVDZ"
    
    # Run the test.
    run_submission_test([water_xyz], [method_code], tmp_path, silico_options)


@pytest.mark.slow
@pytest.mark.parametrize("method", turbomole_ri_methods)
def test_turbomole_ri_methods(method, tmp_path, silico_options):
    """Test different turbomole methods."""
    # Get our method code.
    method_code = f"{destination}/Turbomole/Turbomole:: Single Point: {method}: Gas Phase RI: cc-pVDZ"
    
    # Run the test.
    run_submission_test([water_xyz], [method_code], tmp_path, silico_options)


# For properties, we use the larger pyridine test molecule.
# This is because lots of properties that we test for are not properly expressed by the smaller test molecules.
# For example, water doesn't have 100 electronic excitations, and neither does it have a stable S(1) or T(1) opt structure.

@pytest.mark.slow
@pytest.mark.parametrize("prop", gaussian_properties)
def test_gaussian_properties(prop, tmp_path, silico_options):
    """Test different gaussian properties."""
    # Get our method code.
    method_code = f"{destination}/Gaussian 16/Gaussian:: {prop}: PBE0 (GD3BJ): Gas Phase: 6-31G*"
    
    # Run the test.
    run_submission_test([pyridine_cml], [method_code], tmp_path, silico_options)


@pytest.mark.slow
@pytest.mark.parametrize("prop", turbomole_properties)
def test_turbomole_properties(prop, tmp_path, silico_options):
    """Test different turbomole properties."""
    # Get our method code.
    method_code = f"{destination}/Turbomole/Turbomole:: {prop}: Standard DFT: PBE0 (GD3BJ): Gas Phase: 6-31G*"
    
    # Run the test.
    run_submission_test([pyridine_cml], [method_code], tmp_path, silico_options)

@pytest.mark.slow
@pytest.mark.parametrize("prop", turbomole_posthf_properties)
def test_turbomole_posthf_properties(prop, tmp_path, silico_options):
    """Test different turbomole properties (with post HF methods)."""
    # Get our method code.
    method_code = f"{destination}/Turbomole/Turbomole:: {prop}: SCS-ADC(2): Gas Phase RI: cc-pVDZ"
    
    # Run the test.
    run_submission_test([pyridine_cml], [method_code], tmp_path, silico_options)


