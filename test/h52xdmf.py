#!/home/keshav/.pyenv/versions/3.12.4/bin/python3

import argparse
import h5py
import os
from lxml import etree as ET

VERBOSE = False

def set_verbose(value):
    global VERBOSE
    VERBOSE = value

def print_verbose(*args, **kwargs):
    if VERBOSE:
        print(*args, **kwargs)

def write_xdmf(h5_filepath, xdmf_filepath=None, cube_length=[1, 1, 1], ensemble=False):
    """
    Function to convert HDF5 files to XDMF format for visualization.
    Args:
        h5_filepath (str): The path to the input HDF5 file.
        xdmf_filepath (str): The path to the output XDMF file. If None, the name of the HDF5 file is used.
        cube_length (list): The total length of the cube in each dimension.
        ensemble (bool): Whether to visualize the ensemble side by side.
    """
    
    if xdmf_filepath is None:
        xdmf_filepath = os.path.splitext(h5_filepath)[0] + '.xdmf'

    def dataset_to_xdmf(dset, grid):
        dimensions = len(dset.shape)  # Number of dimensions in the dataset.

        if dimensions == 3:
            Nt = 1
            Nx, Ny, Nz = dset.shape
        elif dimensions == 4:
            Nx, Ny, Nz, Nt = dset.shape
            if Nt not in [1, 3, 6, 9]:
                print_verbose(f'Omitting dataset {dset.name} due to unsupported dimension: {Nt}')
                return
        else:
            print_verbose(f'Omitting dataset {dset.name} due to unexpected number of dimensions: {dimensions}')
            return

        element = {
            1: "Scalar",
            3: "Vector",
            6: "Tensor6",
            9: "Tensor",
        }.get(Nt)

        attr = ET.SubElement(
            grid,
            "Attribute",
            Name=dset.name,
            AttributeType=element,
            Center="Cell",
        )

        data_item = ET.SubElement(
            attr,
            "DataItem",
            Dimensions=" ".join(str(i) for i in dset.shape),
            NumberType="Float",
            Precision="4",
            Format="HDF",
        )

        data_item.text = h5_filepath + ":" + dset.name
        print_verbose(f'Converted dataset {dset.name} to XDMF and added to grid {grid.get("Name")}.')

    def crawl_h5_group(group, grid_dict):
        print_verbose(f'Crawling through group {group.name}')

        for dset in group.values():
            if isinstance(dset, h5py.Dataset):
                if len(dset.shape) in [3, 4]:
                    grid_size = dset.shape[:3]
                    if grid_size not in grid_dict:
                        grid_dict[grid_size] = []
                    grid_dict[grid_size].append(dset)
                else:
                    print_verbose(f'Omitting dataset {dset.name} due to unexpected number of dimensions: {len(dset.shape)}')

        for subgroup in group.values():
            if isinstance(subgroup, h5py.Group):
                crawl_h5_group(subgroup, grid_dict)

    with h5py.File(h5_filepath, "r") as h5_file:
        print_verbose(f'Opened HDF5 file: {h5_filepath}')

        xdmf_root = ET.Element("Xdmf", Version="2.0")
        domain = ET.SubElement(xdmf_root, "Domain")

        grid_dict = {}
        crawl_h5_group(h5_file, grid_dict)

        print_verbose(f'Found {len(grid_dict)} unique grid sizes.')

        grid_counter = 0
        for grid_size, dsets in grid_dict.items():
            Nx, Ny, Nz = grid_size
            DX, DY, DZ = cube_length[0] / Nx, cube_length[1] / Ny, cube_length[2] / Nz
            grid = ET.SubElement(domain, "Grid", Name=f"{Nx}x{Ny}x{Nz}", GridType="Uniform")

            print_verbose(f'Creating grid with size {Nx}x{Ny}x{Nz}.')

            topology = ET.SubElement(grid, "Topology", TopologyType="3DCoRectMesh", NumberOfElements=f"{Nx} {Ny} {Nz}")
            geometry = ET.SubElement(grid, "Geometry", GeometryType="ORIGIN_DXDYDZ")

            origin = f"{grid_counter * cube_length[0]} 0 0" if ensemble else "0 0 0"
            grid_counter += 1

            ET.SubElement(geometry, "DataItem", Dimensions="3", NumberType="Float", Precision="4",
                          Format="XML").text = origin
            ET.SubElement(geometry, "DataItem", Dimensions="3", NumberType="Float", Precision="4",
                          Format="XML").text = f"{DX} {DY} {DZ}"

            for dset in dsets:
                dataset_to_xdmf(dset, grid)

        xdmf_tree = ET.ElementTree(xdmf_root)
        xdmf_tree.write(xdmf_filepath, pretty_print=True)

        print_verbose(f'Wrote XDMF file: {xdmf_filepath}.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert HDF5 to XDMF format.")
    parser.add_argument('h5_filepath', type=str, help='Input HDF5 file path.')
    parser.add_argument('-x', '--xdmf_filepath', type=str, default=None, 
                        help='Output XDMF file path. Optional. If not given, uses the HDF5 file name with .xdmf extension.')
    parser.add_argument('-e', '--ensemble', action='store_true', help='Visualize ensemble side by side.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output.')
    parser.add_argument('-c', '--cube_length', type=float, nargs=3, default=[1.0, 1.0, 1.0], metavar=('Lx', 'Ly', 'Lz'),
                        help='Cube length in x, y, z dimensions. Provide three floats. Default is [1.0, 1.0, 1.0].')

    args = parser.parse_args()

    set_verbose(args.verbose)

    write_xdmf(args.h5_filepath, args.xdmf_filepath, args.cube_length, args.ensemble)