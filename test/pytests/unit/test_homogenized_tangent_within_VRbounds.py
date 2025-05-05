import os
import numpy as np
import json
import pytest
from fans_dashboard.core.utils import identify_hierarchy, extract_and_organize_data
from fans_dashboard.core.tensortools import Ciso, compute_VoigtReuss_bounds, is_spd, compute_volume_fractions


@pytest.fixture(
    params=[
        (os.path.join("input_files", "test_J2Plasticity.json"), "test_J2Plasticity.h5"),
        (
            os.path.join("input_files", "test_LinearElastic.json"),
            "test_LinearElastic.h5",
        ),
        (
            os.path.join("input_files", "test_LinearThermal.json"),
            "test_LinearThermal.h5",
        ),
        (
            os.path.join("input_files", "test_PseudoPlastic.json"),
            "test_PseudoPlastic.h5",
        ),
    ]
)
def test_files(request):
    base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../..")

    json_path = os.path.join(base_dir, request.param[0])
    h5_path = os.path.join(base_dir, request.param[1])

    if os.path.exists(json_path) and os.path.exists(h5_path):
        return json_path, h5_path
    pytest.skip(f"Required test files not found: {json_path} or {h5_path}")


def test_homogenized_tangent_within_VRbounds(test_files):
    """
    This test verifies that the homogenized tangent is within Voigt and Reuss bounds
    for all microstructures and load cases for linear elastic and thermal problems.

    Parameters
    ----------
    test_files : tuple
        A tuple containing (input_json_file, results_h5_file) paths.
        - input_json_file: Path to the JSON file containing configuration data
        - results_h5_file: Path to the HDF5 file containing simulation results

    Returns
    -------
    None
        The test passes if the homogenized tangent is within bounds.
    """
    input_json_file, results_h5_file = test_files
    print(f"Running test on: {input_json_file} and {results_h5_file}")

    # Load the input json file to check material model and properties
    with open(input_json_file, "r") as f:
        input_data = json.load(f)
    
    # Try to access results or top-level fields
    if "results" not in input_data:
        print(f"ERROR: 'results' field not found in input_data")
        for key in input_data:
            print(f"- {key}: {type(input_data[key])}")
    results = input_data.get("results", [])
    
    # Get material model
    mat_model = input_data.get("matmodel")
    
    # Skip if not a linear model we're interested in
    if mat_model not in ["LinearElasticIsotropic", "LinearThermalIsotropic"]:
        pytest.skip(f"Skipping test: Material model {mat_model} is not supported for bounds check")
        return

    # Check if homogenized_tangent field is available
    if "homogenized_tangent" not in results:
        pytest.skip(f"Skipping test: No homogenized_tangent field found in {input_json_file}")
        return
    
    # Check if microstructure field is available
    if "microstructure" not in results:
        pytest.skip(f"Skipping test: No microstructure field found in {input_json_file}")
        return

    # Extract hierarchy information from the h5 file
    hierarchy = identify_hierarchy(results_h5_file)

    # Load the data from the HDF5 file
    microstructures_to_load = list(hierarchy.keys())

    quantities_to_load = ["homogenized_tangent", "microstructure"]

    time_steps_to_load = []
    load_cases_to_load = []

    # Get all unique load cases across all microstructures
    for microstructure in microstructures_to_load:
        for load_case in hierarchy[microstructure].keys():
            if load_case not in load_cases_to_load:
                load_cases_to_load.append(load_case)

    # Extract the specified data, organized and sorted by time steps
    data = extract_and_organize_data(
        results_h5_file,
        hierarchy,
        quantities_to_load,
        microstructures_to_load,
        load_cases_to_load,
        time_steps_to_load,
    )

    print(f"\nVerifying homogenized tangent is within Voigt and Reuss bounds...")

    # Get material properties based on material model
    material_properties = input_data.get("material_properties", {})
    
    for microstructure_name in microstructures_to_load:
        # Get the microstructure data and compute volume fractions
        if "microstructure" not in data[microstructure_name][list(data[microstructure_name].keys())[0]]:
            print(f"Skipping {microstructure_name}: Missing microstructure data")
            continue
            
        microstructure_data = data[microstructure_name][list(data[microstructure_name].keys())[0]]["microstructure"]
        volume_fractions = compute_volume_fractions(microstructure_data)
        
        # Create phase tensors based on material model
        phase_tensors = []
        if mat_model == "LinearThermalIsotropic":
            conductivities = material_properties.get("conductivity")
            for k in conductivities:
                # Create conductivity tensor (diagonal matrix with conductivity values)
                phase_tensors.append(k * np.eye(3))
        elif mat_model == "LinearElasticIsotropic":
            bulk_moduli = material_properties.get("bulk_modulus")
            shear_moduli = material_properties.get("shear_modulus")
            for k, g in zip(bulk_moduli, shear_moduli):
                phase_tensors.append(Ciso(k, g))
                
        # Compute Voigt and Reuss bounds
        voigt, reuss = compute_VoigtReuss_bounds(phase_tensors, volume_fractions)
        
        for load_case in load_cases_to_load:
            if load_case in hierarchy[microstructure_name]:
                if "homogenized_tangent" not in data[microstructure_name][load_case]:
                    print(f"Skipping {microstructure_name}/{load_case}: Missing homogenized_tangent field")
                    continue

                tangent_data = data[microstructure_name][load_case]["homogenized_tangent"]

                # Check each time step
                for time_idx, homogenized_tangent in enumerate(tangent_data):
                    # Check if homogenized_tangent is within bounds
                    # Voigt - homogenized_tangent should be SPD (positive eigenvalues)
                    voigt_minus_hom, voigt_eigs = is_spd(voigt - homogenized_tangent)
                    
                    # homogenized_tangent - Reuss should be SPD (positive eigenvalues)
                    hom_minus_reuss, reuss_eigs = is_spd(homogenized_tangent - reuss)
                    
                    assert voigt_minus_hom, (
                        f"For microstructure {microstructure_name}, load case {load_case}, time step {time_idx}: "
                        f"Homogenized tangent exceeds Voigt bound."
                        f"\nMin eigenvalue of (Voigt - Homogenized): {min(voigt_eigs)}"
                    )

                    assert hom_minus_reuss, (
                        f"For microstructure {microstructure_name}, load case {load_case}, time step {time_idx}: "
                        f"Homogenized tangent is below Reuss bound."
                        f"\nMin eigenvalue of (Homogenized - Reuss): {min(reuss_eigs)}"
                    )

                print(
                    f"Verified: {microstructure_name}, load case {load_case} - homogenized tangent is within bounds "
                    f"({tangent_data.shape[0]} time steps)"
                )


if __name__ == "__main__":
    pytest.main(["-v", "-s", __file__])