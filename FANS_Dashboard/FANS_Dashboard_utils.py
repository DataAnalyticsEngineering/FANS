import h5py
import numpy as np
from collections import defaultdict
from postprocessing import compute_rank2tensor_measures

def recursively_find_structure(group, current_path=""):
    """
    Recursively find the structure of the HDF5 file and return a dictionary
    containing the datasets' names and their shapes.
    """
    structure = {}
    if isinstance(group, h5py.Group):
        for key in group.keys():
            item = group[key]
            new_path = f"{current_path}/{key}"
            if isinstance(item, h5py.Group):
                nested_structure = recursively_find_structure(item, new_path)
                if nested_structure:
                    structure[key] = nested_structure
            elif isinstance(item, h5py.Dataset):
                structure[key] = item.shape  # Store the shape for summarization
    return structure

def identify_hierarchy(file_path):
    """
    Identify the hierarchy of microstructures, load cases, and time steps.
    """
    def find_microstructures(group, path=''):
        microstructures = {}
        for key in group.keys():
            item = group[key]
            new_path = f"{path}/{key}"
            if isinstance(item, h5py.Group):
                if any(subkey.startswith('load') for subkey in item.keys()):
                    microstructures[new_path] = recursively_find_structure(item, new_path)
                else:
                    deeper_microstructures = find_microstructures(item, new_path)
                    microstructures.update(deeper_microstructures)
        return microstructures

    with h5py.File(file_path, 'r') as h5file:
        hierarchy = find_microstructures(h5file)
    return hierarchy

def summarize_hierarchy(hierarchy):
    """
    Summarize the identified hierarchy by printing out the number of microstructures,
    load cases, the number of time steps, and the quantities available within each load case.
    """
    microstructure_count = len(hierarchy)
    print(f"Found {microstructure_count} microstructure(s):")
    
    for microstructure, loads in hierarchy.items():
        load_cases = {key: value for key, value in loads.items() if key.startswith('load')}
        load_case_count = len(load_cases)
        print(f"  Microstructure '{microstructure}' has {load_case_count} load case(s):")
        
        for load_case, time_steps in load_cases.items():
            time_step_count = len(time_steps)
            print(f"    Load case '{load_case}' has {time_step_count} time step(s) and the following quantities:")
            
            # Print available quantities in the first time step (assumed same for all time steps)
            first_time_step = next(iter(time_steps.values()))
            for quantity, shape in first_time_step.items():
                print(f"      - {quantity} with shape {shape}")

def extract_and_organize_data(file_path, hierarchy, quantities_to_load, microstructures_to_load=None, load_cases_to_load=None, time_steps_to_load=None):
    """
    Extract and organize the specified quantities into a structured format,
    keeping the original quantity names and aggregating them as batches along the 0th index.
    Users can also specify which microstructures, load cases, and time steps to load.
    """
    organized_data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    
    with h5py.File(file_path, 'r') as h5file:
        for microstructure, loads in hierarchy.items():
            if microstructures_to_load and microstructure not in microstructures_to_load:
                continue  # Skip microstructures that are not in the user-provided list
            
            for load_case, time_steps in loads.items():
                if load_cases_to_load and load_case not in load_cases_to_load:
                    continue  # Skip load cases that are not in the user-provided list
                
                extracted_data = {quantity: [] for quantity in quantities_to_load}
                time_step_values = []
                
                sorted_time_steps = sorted(time_steps.items(), key=lambda x: int(x[0].replace('time_step', '')))
                
                for time_step, _ in sorted_time_steps:
                    time_step_number = int(time_step.replace('time_step', ''))
                    if time_steps_to_load and time_step_number not in time_steps_to_load:
                        continue  # Skip time steps that are not in the user-provided list
                    
                    for quantity in quantities_to_load:
                        if quantity in time_steps[time_step]:
                            data_array = np.array(h5file[f"{microstructure}/{load_case}/{time_step}/{quantity}"])
                            extracted_data[quantity].append(data_array)
                    
                    time_step_values.append(time_step_number)
                
                # Aggregate the extracted quantities as batches along the 0th index
                for quantity, data_list in extracted_data.items():
                    if data_list:
                        organized_data[microstructure][load_case][quantity] = np.stack(data_list, axis=0)
                
                organized_data[microstructure][load_case]['time_steps'] = np.array(time_step_values)
    
    return organized_data






def write_processed_data_to_h5(file_path, processed_data, overwrite=True):
    """
    Writes processed data back into the HDF5 file at the correct locations.
    
    Parameters:
    - file_path: str, path to the HDF5 file.
    - processed_data: dict, structured data with computed measures to write back.
    - overwrite: bool, whether to overwrite existing datasets.
    
    Returns:
    - None
    """
    with h5py.File(file_path, 'a') as h5file:
        for microstructure, loads in processed_data.items():
            for load_case, data in loads.items():
                time_steps = data.pop('time_steps', [])
                for measure_name, data_array in data.items():
                    for time_step_idx, time_step in enumerate(time_steps):
                        time_step_group = f"{microstructure}/{load_case}/time_step{time_step}"
                        dataset_path = f"{time_step_group}/{measure_name}"
                        
                        # Ensure the group exists
                        if time_step_group not in h5file:
                            raise ValueError(f"Group '{time_step_group}' does not exist in the HDF5 file.")
                        
                        # Check if the dataset already exists
                        if dataset_path in h5file:
                            if overwrite:
                                del h5file[dataset_path]
                                h5file.create_dataset(dataset_path, data=data_array[time_step_idx])
                            else:
                                print(f"Dataset exists and overwrite=False: {dataset_path}")
                        else:
                            h5file.create_dataset(dataset_path, data=data_array[time_step_idx])


def postprocess_and_write_back(file_path, hierarchy, quantities_to_process, measures, 
                               microstructures_to_load=None, load_cases_to_load=None):
    """
    A higher-level function that extracts specific data from an HDF5 file, processes it to compute various tensor measures,
    and writes the results back into the HDF5 file using the naming convention 'quantity_measure'.

    Parameters:
    - file_path: str
        Path to the HDF5 file containing the data.
    - hierarchy: dict
        The structure of the HDF5 file as identified by the `identify_hierarchy` function.
    - quantities_to_process: list of str
        List of quantities (e.g., 'stress_average', 'strain_average') to extract from the HDF5 file and process.
    - measures: list of str
        List of tensor measures to compute for the extracted quantities. 
        Available options can be found in the `compute_rank2tensor_measures` function in the `postprocessing` module.
    - microstructures_to_load: list of str, optional
        List of microstructures to process. If None, all microstructures found in the hierarchy will be processed.
    - load_cases_to_load: list of str, optional
        List of load cases to process. If None, all load cases found in the hierarchy will be processed.

    Returns:
    - processed_data: dict
        A dictionary containing the processed data, organized by microstructure and load case. 
        The computed measures are stored under keys following the 'quantity_measure' naming convention.
        Additionally, the processed data is also written back into the HDF5 file at the corresponding locations.
    """
    # Extract the data (loads all time steps if time_steps_to_load is None or empty)
    extracted_data = extract_and_organize_data(file_path, hierarchy, quantities_to_process, 
                                               microstructures_to_load, load_cases_to_load)
    
    # Process the data and prepare for writing
    processed_data = defaultdict(lambda: defaultdict(dict))
    for microstructure, loads in extracted_data.items():
        for load_case, quantities in loads.items():
            time_steps = quantities.pop('time_steps', [])
            for quantity_name, tensor_data in quantities.items():
                measures_results = compute_rank2tensor_measures(tensor_data, measures)
                for measure_name, measure_data in measures_results.items():
                    key = f"{quantity_name}_{measure_name}"
                    processed_data[microstructure][load_case][key] = measure_data
            processed_data[microstructure][load_case]['time_steps'] = time_steps
    
    # Write the processed data back to the HDF5 file
    write_processed_data_to_h5(file_path, processed_data)

    return processed_data