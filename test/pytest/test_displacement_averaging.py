import numpy as np
import json
import pytest
from MSUtils.fans_dashboard.utils import identify_hierarchy, extract_and_organize_data


def test_displacement_averaging(test_files):
    """
    This test verifies that the average of displacement fluctuations is zero for all
    microstructures and load cases.

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

    # Check if displacement_fluctuation field is available
    if "displacement_fluctuation" not in results:
        pytest.skip(
            f"Skipping test: No displacement_fluctuation field found in {input_json_file}"
        )
        return

    # Extract hierarchy information from the h5 file
    hierarchy = identify_hierarchy(results_h5_file)

    # Load the data from the HDF5 file
    microstructures_to_load = list(hierarchy.keys())

    quantities_to_load = ["displacement_fluctuation"]

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

    print(f"\nVerifying displacement fluctuation averages are zero...")

    for microstructure in microstructures_to_load:
        for load_case in load_cases_to_load:
            if load_case in hierarchy[microstructure]:
                if "displacement_fluctuation" not in data[microstructure][load_case]:
                    print(
                        f"Skipping {microstructure}/{load_case}: Missing displacement_fluctuation field"
                    )
                    continue

                displacement_data = data[microstructure][load_case][
                    "displacement_fluctuation"
                ]

                # Compute average manually by averaging over spatial dimensions (1, 2, 3)
                # displacement_data shape: time_steps x Nx x Ny x Nz x components
                computed_average = np.mean(displacement_data, axis=(1, 2, 3))

                # Check if the computed average is approximately zero
                zero_array = np.zeros_like(computed_average)

                assert np.allclose(
                    computed_average, zero_array, rtol=1e-5, atol=1e-8
                ), (
                    f"For microstructure {microstructure}, load case {load_case}: "
                    f"Average of displacement fluctuations is not zero."
                    f"\nComputed average shape: {computed_average.shape}, Value: {computed_average}"
                )

                print(
                    f"Verified: {microstructure}, load case {load_case} - displacement fluctuation average is zero"
                )


if __name__ == "__main__":
    pytest.main(["-v", "-s", __file__])
