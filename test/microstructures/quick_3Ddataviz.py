# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     custom_cell_magics: kql
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: pytorch
#     language: python
#     name: python3
# ---

# %%
import h5py
import numpy as np

# %%
def visualize_3d_data(file_path, group_name, dset_name, image_name, reduction_factor=1, output_filename=None): 
    with h5py.File(file_path, 'r') as file:
        group = file[group_name]
        dset = group[dset_name]
        image = dset[image_name]

        data = np.array(image, dtype=np.uint8)

        if reduction_factor > 1:
            data = data[::reduction_factor, ::reduction_factor, ::reduction_factor]
        
        # write to h5 file output_filename
        if output_filename:
            with h5py.File(output_filename, 'w') as f:
                dset = f.create_dataset(image_name, data=data)


# %%
# Use the function
visualize_3d_data(file_path='ulm_dataset.h5',
                  group_name='structures_test',
                  dset_name='dset_3456',
                  image_name='image',
                  reduction_factor=1,
                  output_filename='test.h5')
