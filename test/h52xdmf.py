#!/home/keshav/.pyenv/versions/3.12.4/bin/python3
import argparse
import h5py
import os
from lxml import etree as ET
import re

VERBOSE = False

def set_verbose(value):
    global VERBOSE
    VERBOSE = value

def print_verbose(*args, **kwargs):
    if VERBOSE:
        print(*args, **kwargs)

def write_xdmf(h5_filepath, xdmf_filepath=None, cube_length=[1, 1, 1], time_series=False, time_keyword="load"):
    """
    Function to convert HDF5 files to XDMF format for visualization.
    Args:
        h5_filepath (str): The path to the input HDF5 file.
        xdmf_filepath (str): The path to the output XDMF file. If None, the name of the HDF5 file is used.
        cube_length (list): The total length of the cube in each dimension.
        time_series (bool): Whether to treat the data as a time series.
        time_keyword (str): Keyword used to identify temporal datasets.
    """
    
    if xdmf_filepath is None:
        xdmf_filepath = os.path.splitext(h5_filepath)[0] + '.xdmf'

    def dataset_to_xdmf(dset, grid, time=None):
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

        if time_series:
            # Remove the time_keyword{time_step}/ part from the attribute name
            attr_name_parts = dset.name.split('/')
            time_index = [i for i, part in enumerate(attr_name_parts) if part.startswith(time_keyword)]
            if time_index:
                time_index = time_index[0]
                attr_name = '/'.join(attr_name_parts[:time_index] + attr_name_parts[time_index + 1:])
            else:
                attr_name = dset.name
        else:
            attr_name = dset.name
        
        attr = ET.SubElement(
            grid,
            "Attribute",
            Name=attr_name,
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

    def crawl_h5_group(group, grid_dict, time=None):
        print_verbose(f'Crawling through group {group.name}')

        for dset in group.values():
            if isinstance(dset, h5py.Dataset):
                if len(dset.shape) in [3, 4]:
                    grid_size = dset.shape[:3]
                    if time_series and time is not None:
                        grid_key = (grid_size, time)
                    else:
                        if time_series and time is None:
                            # Skip datasets without a time step in time series mode
                            continue
                        grid_key = (grid_size, None)
                    
                    if grid_key not in grid_dict:
                        grid_dict[grid_key] = []
                    grid_dict[grid_key].append(dset)
                else:
                    print_verbose(f'Omitting dataset {dset.name} due to unexpected number of dimensions: {len(dset.shape)}')

        for subgroup in group.values():
            if isinstance(subgroup, h5py.Group):
                if time_series and re.search(rf'{time_keyword}\d*\.?\d+', subgroup.name):
                    time_step = float(re.search(rf'{time_keyword}(\d*\.?\d+)', subgroup.name).group(1))
                    crawl_h5_group(subgroup, grid_dict, time=time_step)
                else:
                    crawl_h5_group(subgroup, grid_dict, time=time)

    with h5py.File(h5_filepath, "r") as h5_file:
        print_verbose(f'Opened HDF5 file: {h5_filepath}')

        xdmf_root = ET.Element("Xdmf", Version="2.0")
        domain = ET.SubElement(xdmf_root, "Domain")

        grid_dict = {}
        crawl_h5_group(h5_file, grid_dict)

        print_verbose(f'Found {len(grid_dict)} unique grid sizes.')

        temporal_grids = {}

        for (grid_size, time), dsets in list(grid_dict.items()):
            Nx, Ny, Nz = grid_size
            DX, DY, DZ = cube_length[0] / Nx, cube_length[1] / Ny, cube_length[2] / Nz
            grid_name = f"{Nx}x{Ny}x{Nz}"

            if time_series:
                if time is None:
                    continue
                if grid_size not in temporal_grids:
                    grid = ET.SubElement(domain, "Grid", Name=grid_name, GridType="Collection", CollectionType="Temporal")
                    temporal_grids[grid_size] = grid
                else:
                    grid = temporal_grids[grid_size]
                
                subgrid = ET.SubElement(grid, "Grid", Name=grid_name, GridType="Uniform")
                ET.SubElement(subgrid, "Time", Type="Single", Value=str(time))
            else:
                grid = ET.SubElement(domain, "Grid", Name=grid_name, GridType="Uniform")
                subgrid = grid

            topology = ET.SubElement(subgrid, "Topology", TopologyType="3DCoRectMesh", NumberOfElements=f"{Nx} {Ny} {Nz}")
            geometry = ET.SubElement(subgrid, "Geometry", GeometryType="ORIGIN_DXDYDZ")

            origin = "0 0 0"

            ET.SubElement(geometry, "DataItem", Dimensions="3", NumberType="Float", Precision="4",
                          Format="XML").text = origin
            ET.SubElement(geometry, "DataItem", Dimensions="3", NumberType="Float", Precision="4",
                          Format="XML").text = f"{DX} {DY} {DZ}"

            for dset in dsets:
                dataset_to_xdmf(dset, subgrid, time)

        xdmf_tree = ET.ElementTree(xdmf_root)
        xdmf_tree.write(xdmf_filepath, pretty_print=True)

        print_verbose(f'Wrote XDMF file: {xdmf_filepath}.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Convert HDF5 files with 3D spaciotemporal microstructure data to XDMF representation for visualization.\n\n"
            "This script reads an input HDF5 file containing 3D microstructural field data and converts it into an XDMF file. "
            "The resulting XDMF file can be used for visualization in software like ParaView."
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'h5_filepath', 
        type=str, 
        help='Path to the input HDF5 file.'
    )
    parser.add_argument(
        '-x', '--xdmf_filepath', 
        type=str, 
        default=None, 
        help='Path to the output XDMF file. Optional. If not provided, the HDF5 filename with .xdmf extension is used.'
    )
    parser.add_argument(
        '-v', '--verbose', 
        action='store_true', 
        help='Enable verbose output.'
    )
    parser.add_argument(
        '-c', '--cube_length', 
        type=float, 
        nargs=3, 
        default=[1.0, 1.0, 1.0], 
        metavar=('Lx', 'Ly', 'Lz'),
        help=(
            "Cube length in x, y, z dimensions.\n"
            "Provide three floats. Default is [1.0, 1.0, 1.0]."
        )
    )
    parser.add_argument(
        '-t', '--time-series', 
        action='store_true', 
        help=(
            "Treat datasets as a time series based on load groups.\n"
            "If this flag is set, the script will search for datasets within groups named '<keyword>{number}' and creates a temporal collection in the XDMF file."
        )
    )
    parser.add_argument(
        '-k', '--time-keyword', 
        type=str, 
        default="load",
        help=(
            "Keyword used to identify temporal datasets. Default is 'load'.\n"
            "This keyword should be followed by a number to denote the time in the group names."
        )
    )

    args = parser.parse_args()
    set_verbose(args.verbose)

    write_xdmf(args.h5_filepath, args.xdmf_filepath, args.cube_length, args.time_series, args.time_keyword)