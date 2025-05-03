import os
import numpy as np
import json
import pytest
from fans_dashboard.core.utils import identify_hierarchy, extract_and_organize_data

@pytest.fixture(params=[
    (os.path.join("input_files", "test_J2Plasticity.json"), "test_J2Plasticity.h5"),
    (os.path.join("input_files", "test_LinearElastic.json"), "test_LinearElastic.h5"),
    (os.path.join("input_files", "test_LinearThermal.json"), "test_LinearThermal.h5"),
    (os.path.join("input_files", "test_PseudoPlastic.json"), "test_PseudoPlastic.h5")
])
def test_files(request):
    base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../..")

    json_path = os.path.join(base_dir, request.param[0])
    h5_path = os.path.join(base_dir, request.param[1])
    
    if os.path.exists(json_path) and os.path.exists(h5_path):
        return json_path, h5_path
    pytest.skip(f"Required test files not found: {json_path} or {h5_path}")


def test_loading_to_strain_average(test_files):
    """
    This test verifies that the strain_average field in the results matches the macroscale_loading
    specified in the input JSON file for all microstructures and load cases.
    Parameters
    ----------
    test_files : tuple
        A tuple containing (input_json_file, results_h5_file) paths.
        - input_json_file: Path to the JSON file containing macroscale_loading data
        - results_h5_file: Path to the HDF5 file containing simulation results
    Returns
    -------
    None
        The test passes if strain_average equals macroscale_loading for all cases.
        The test is skipped if 'strain_average' is not found in the results field.
    """
    input_json_file, results_h5_file = test_files
    
    # Load the macroscale_loading field from the input json file
    with open(input_json_file, 'r') as f:
        input_data = json.load(f)    
    
    # Check if 'strain_average' exists in the results field
    results = input_data.get('results', [])
    if 'strain_average' not in results:
        pytest.skip(f"Skipping test: 'strain_average' not requested in {input_json_file}")
        return
        
    macroscale_loading = np.array(input_data.get('macroscale_loading', {}))
    
    # Extract hierarchy information from the h5 file
    hierarchy = identify_hierarchy(results_h5_file)
    
    # Load the data from the HDF5 file
    microstructures_to_load = list(hierarchy.keys())
    quantities_to_load = ['strain_average']
    time_steps_to_load = []

    load_cases_to_load = []
    # Get all unique load cases across all microstructures
    for microstructure in microstructures_to_load:
        for load_case in hierarchy[microstructure].keys():
            if load_case not in load_cases_to_load:
                load_cases_to_load.append(load_case)

    # Extract the specified data, organized and sorted by time steps
    data = extract_and_organize_data(results_h5_file, hierarchy, 
                                    quantities_to_load, 
                                    microstructures_to_load, 
                                    load_cases_to_load, 
                                    time_steps_to_load)
    
    # Comprehensive check for all microstructures and load cases
    for i, microstructure in enumerate(microstructures_to_load):
        for j, load_case in enumerate(load_cases_to_load):
            if load_case in hierarchy[microstructure]:
                strain_average = data[microstructure][load_case]['strain_average']
                
                # Get corresponding macroscale loading
                # Assuming macroscale_loading is organized to match load cases
                current_loading = np.array(macroscale_loading)[j] if j < len(macroscale_loading) else None
                
                if current_loading is not None:
                    assert np.allclose(strain_average, current_loading), \
                        f"For microstructure {microstructure}, load case {load_case}: " \
                        f"strain_average and macroscale_loading are not close."
                    print(f"Verified: microstructure {microstructure}, load case {load_case}")
                    
                    
if __name__ == "__main__":
    
    pytest.main(["-v", "-s", __file__])