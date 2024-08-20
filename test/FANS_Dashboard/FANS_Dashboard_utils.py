import h5py
import numpy as np
from collections import defaultdict

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












import plotly.graph_objects as go
from plotly.subplots import make_subplots

def plot_subplots(data1, data2, labels, title="Subplot Grid", nrows=None, ncols=None):
    """
    Plot a grid of subplots using Plotly, handling both single-component (scalar vs scalar) and multi-component data.
    
    Parameters:
    - data1: numpy array, first set of data to plot (e.g., strain, time)
    - data2: numpy array, second set of data to plot (e.g., stress)
    - labels: tuple of strings, labels for the x and y axes (e.g., ("Strain", "Stress"))
    - title: string, title of the overall plot
    - nrows: int, number of rows in the subplot grid (optional)
    - ncols: int, number of columns in the subplot grid (optional)
    """
    # Ensure data1 and data2 are 2D arrays
    if data1.ndim == 1:
        data1 = data1[:, np.newaxis]
    if data2.ndim == 1:
        data2 = data2[:, np.newaxis]
    
    # Set the number of components based on data shape
    n_components = data1.shape[1]
    
    # If nrows or ncols is not specified, determine an optimal grid layout
    if nrows is None or ncols is None:
        nrows = int(np.ceil(np.sqrt(n_components)))
        ncols = int(np.ceil(n_components / nrows))
    
    # Create the subplot figure
    fig = make_subplots(rows=nrows, cols=ncols, subplot_titles=[f'Component {i+1}' for i in range(n_components)])

    # Add traces for each component
    for i in range(n_components):
        row = i // ncols + 1
        col = i % ncols + 1
        fig.add_trace(go.Scatter(
            x=data1[:, i],
            y=data2[:, i],
            mode='lines+markers',
            marker=dict(symbol='x', size=4),
            line=dict(width=1),
            name=f'Component {i+1}'
        ), row=row, col=col)

    # Update layout with text labels and styling
    fig.update_layout(
        height=600,
        width=900,
        title_text=title,
        showlegend=False,
        template="plotly_white",
    )

    # Update axes with text labels
    for i in range(n_components):
        row = i // ncols + 1
        col = i % ncols + 1
        fig.update_xaxes(title_text=labels[0], row=row, col=col, showgrid=True)
        fig.update_yaxes(title_text=labels[1], row=row, col=col, showgrid=True)

    # Adjust layout for tight spacing
    fig.update_layout(margin=dict(l=20, r=20, t=50, b=20), title_x=0.5)

    # Show the plot
    fig.show()