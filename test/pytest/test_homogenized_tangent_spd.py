import os
import numpy as np
import json
import pytest
from fans_dashboard.core.utils import identify_hierarchy, extract_and_organize_data
from scipy.linalg import eigvalsh


@pytest.fixture(
    params=[
        "test_J2Plasticity",
        "test_LinearElastic",
        "test_LinearThermal",
        "test_PseudoPlastic",
    ]
)
def test_files(request):
    json_base_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "../input_files/"
    )
    h5_base_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "../../build/test/"
    )

    json_path = os.path.join(json_base_dir, f"{request.param}.json")
    h5_path = os.path.join(h5_base_dir, f"{request.param}.h5")

    if os.path.exists(json_path) and os.path.exists(h5_path):
        return json_path, h5_path
    pytest.skip(f"Required test files not found: {json_path} or {h5_path}")


def test_homogenized_tangent_spd(test_files):
    """
    This test verifies that the homogenized tangent is strictly Symmetric Positive Definite (SPD)
    for all microstructures and load cases.

    Parameters
    ----------
    test_files : tuple
        A tuple containing (input_json_file, results_h5_file) paths.
        - input_json_file: Path to the JSON file containing configuration data
        - results_h5_file: Path to the HDF5 file containing simulation results
    """
    input_json_file, results_h5_file = test_files

    # Load the input json file to check which fields are requested
    with open(input_json_file, "r") as f:
        input_data = json.load(f)

    # Check which fields are available in the results
    results = input_data.get("results", [])

    # Check if homogenized_tangent field is available
    if "homogenized_tangent" not in results:
        pytest.skip(
            f"Skipping test: No homogenized_tangent field found in {input_json_file}"
        )
        return

    # Extract hierarchy information from the h5 file
    hierarchy = identify_hierarchy(results_h5_file)

    # Load the data from the HDF5 file
    microstructures_to_load = list(hierarchy.keys())

    quantities_to_load = ["homogenized_tangent"]

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

    print(f"\nVerifying homogenized tangent is strictly SPD...")

    for microstructure in microstructures_to_load:
        for load_case in load_cases_to_load:
            if load_case in hierarchy[microstructure]:
                if "homogenized_tangent" not in data[microstructure][load_case]:
                    print(
                        f"Skipping {microstructure}/{load_case}: Missing homogenized_tangent field"
                    )
                    continue

                tangent_data = data[microstructure][load_case]["homogenized_tangent"]

                # Check each time step
                for time_idx, tangent in enumerate(tangent_data):
                    # Check symmetry
                    is_symmetric = np.allclose(tangent, tangent.T, rtol=1e-5, atol=1e-8)

                    # Check positive definiteness by computing eigenvalues
                    eigenvalues = eigvalsh(tangent)
                    is_positive_definite = np.all(eigenvalues > 0)

                    assert is_symmetric, (
                        f"For microstructure {microstructure}, load case {load_case}, time step {time_idx}: "
                        f"Homogenized tangent is not symmetric."
                        f"\nTangent shape: {tangent.shape}, Max asymmetry: {np.max(np.abs(tangent - tangent.T))}"
                    )

                    assert is_positive_definite, (
                        f"For microstructure {microstructure}, load case {load_case}, time step {time_idx}: "
                        f"Homogenized tangent is not positive definite."
                        f"\nEigenvalues: {eigenvalues}, Min eigenvalue: {np.min(eigenvalues)}"
                    )

                print(
                    f"Verified: {microstructure}, load case {load_case} - homogenized tangent is strictly SPD "
                    f"({tangent_data.shape[0]} time steps)"
                )


if __name__ == "__main__":
    pytest.main(["-v", "-s", __file__])
