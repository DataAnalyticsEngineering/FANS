import os
import numpy as np
import json
import pytest
from fans_dashboard.core.utils import identify_hierarchy, extract_and_organize_data


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


def test_homogenization_consistency(test_files):
    """
    This test verifies that the relationship stress_average = homogenized_tangent * strain_average
    holds for all microstructures and load cases.

    Parameters
    ----------
    test_files : tuple
        A tuple containing (input_json_file, results_h5_file) paths.
        - input_json_file: Path to the JSON file containing configuration data
        - results_h5_file: Path to the HDF5 file containing simulation results

    Returns
    -------
    None
        The test passes if stress_average equals homogenized_tangent * strain_average.
    """
    input_json_file, results_h5_file = test_files

    # Load the input json file to check which fields are requested
    with open(input_json_file, "r") as f:
        input_data = json.load(f)

    # Check which fields are available in the results
    results = input_data.get("results", [])

    # Check if required fields are available
    required_fields = ["strain_average", "stress_average", "homogenized_tangent"]
    missing_fields = [field for field in required_fields if field not in results]

    if missing_fields:
        pytest.skip(
            f"Skipping test: Missing required fields {', '.join(missing_fields)} in {input_json_file}"
        )

    # Extract hierarchy information from the h5 file
    hierarchy = identify_hierarchy(results_h5_file)

    # Load the data from the HDF5 file
    microstructures_to_load = list(hierarchy.keys())

    quantities_to_load = ["strain_average", "stress_average", "homogenized_tangent"]

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

    print(f"\nVerifying stress_average = homogenized_tangent * strain_average...")

    for microstructure in microstructures_to_load:
        for load_case in load_cases_to_load:
            if load_case in hierarchy[microstructure]:
                # Check if all required fields exist for this microstructure and load case
                required_data = {
                    field: field in data[microstructure][load_case]
                    for field in quantities_to_load
                }

                if not all(required_data.values()):
                    missing = [
                        field for field, exists in required_data.items() if not exists
                    ]
                    print(
                        f"Skipping {microstructure}/{load_case}: Missing fields: {', '.join(missing)}"
                    )
                    continue

                strain_avg = data[microstructure][load_case]["strain_average"]
                stress_avg = data[microstructure][load_case]["stress_average"]
                h_tangent = data[microstructure][load_case]["homogenized_tangent"]

                # Loop through time steps
                for t in range(len(strain_avg)):
                    # Compute stress using homogenized tangent and strain
                    computed_stress = np.matmul(h_tangent[t], strain_avg[t])

                    # Check if computed stress matches the reported average stress
                    assert np.allclose(
                        computed_stress, stress_avg[t], rtol=1e-4, atol=1e-1
                    ), (
                        f"For microstructure {microstructure}, load case {load_case}, time step {t}: "
                        f"stress_average != homogenized_tangent * strain_average"
                        f"\nComputed stress: {computed_stress}"
                        f"\nReported stress: {stress_avg[t]}"
                    )

                print(
                    f"Verified: {microstructure}, load case {load_case} - "
                    f"stress_average = homogenized_tangent * strain_average for all time steps"
                )


if __name__ == "__main__":
    pytest.main(["-v", "-s", __file__])
