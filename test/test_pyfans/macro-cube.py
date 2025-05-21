"""
Run FANS as a Python callable library.
"""
import numpy as np
import precice


def main():
    np_axis = 2

    # preCICE setup
    participant = precice.Participant("macro-cube", "precice-config.xml", 0, 1)
    mesh_name = "cube"

    # Coupling mesh - unit cube
    x_coords, y_coords, z_coords = np.meshgrid(
        np.linspace(0, 1, np_axis),
        np.linspace(0, 1, np_axis),
        np.linspace(0, 1, np_axis),
    )

    nv = np_axis ** participant.get_mesh_dimensions(mesh_name)
    coords = np.zeros((nv, participant.get_mesh_dimensions(mesh_name)))

    # Define unit cube coordinates
    for z in range(np_axis):
        for y in range(np_axis):
            for x in range(np_axis):
                n = x + y * np_axis + z * np_axis * np_axis
                coords[n, 0] = x_coords[x, y, z]
                coords[n, 1] = y_coords[x, y, z]
                coords[n, 2] = z_coords[x, y, z]

    vertex_ids = participant.set_mesh_vertices(mesh_name, coords)

    participant.initialize()

    dt = participant.get_max_time_step_size()

    strains1to3 = np.full((nv, 3), 0.005)
    strains4to6 = np.full((nv, 3), 0.0005)

    # time loop
    while participant.is_coupling_ongoing():
        stress1to3 = participant.read_data(mesh_name, "stresses1to3", vertex_ids, dt)
        stress4to6 = participant.read_data(mesh_name, "stresses4to6", vertex_ids, dt)
        cmat1 = participant.read_data(mesh_name, "cmat1", vertex_ids, dt)
        cmat2 = participant.read_data(mesh_name, "cmat2", vertex_ids, dt)
        cmat3 = participant.read_data(mesh_name, "cmat3", vertex_ids, dt)
        cmat4 = participant.read_data(mesh_name, "cmat4", vertex_ids, dt)
        cmat5 = participant.read_data(mesh_name, "cmat5", vertex_ids, dt)
        cmat6 = participant.read_data(mesh_name, "cmat6", vertex_ids, dt)
        cmat7 = participant.read_data(mesh_name, "cmat7", vertex_ids, dt)

        participant.write_data(mesh_name, "strains1to3", vertex_ids, strains1to3)
        participant.write_data(mesh_name, "strains4to6", vertex_ids, strains4to6)

        participant.advance(dt)
        dt = participant.get_max_time_step_size()

    participant.finalize()


if __name__ == "__main__":
    main()
