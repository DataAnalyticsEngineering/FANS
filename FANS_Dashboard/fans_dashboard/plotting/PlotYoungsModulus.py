import numpy as np
import plotly.graph_objs as go
import meshio


def compute_YoungsModulus3D(C_batch):
    """
    Compute Young's modulus for all directions in 3D for a batch of stiffness tensors.

    Args:
        C_batch (ndarray): Batch of stiffness tensors in Mandel notation, shape (n, 6, 6).

    Returns:
        tuple: A tuple containing:
            - X_batch (ndarray): X-coordinates for plotting the modulus surface, shape (n, n_theta, n_phi).
            - Y_batch (ndarray): Y-coordinates for plotting the modulus surface, shape (n, n_theta, n_phi).
            - Z_batch (ndarray): Z-coordinates for plotting the modulus surface, shape (n, n_theta, n_phi).
            - E_batch (ndarray): Young's modulus in all directions, shape (n, n_theta, n_phi).
    """
    n = C_batch.shape[0]
    n_theta = 180
    n_phi = 360

    theta = np.linspace(0, np.pi, n_theta)
    phi = np.linspace(0, 2 * np.pi, n_phi)
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="ij")

    d_x = np.sin(theta_grid) * np.cos(phi_grid)  # Shape (n_theta, n_phi)
    d_y = np.sin(theta_grid) * np.sin(phi_grid)
    d_z = np.cos(theta_grid)

    N = np.stack(
        (
            d_x**2,
            d_y**2,
            d_z**2,
            np.sqrt(2) * d_x * d_y,
            np.sqrt(2) * d_x * d_z,
            np.sqrt(2) * d_y * d_z,
        ),
        axis=-1,
    )  # Shape (n_theta, n_phi, 6)

    N_flat = N.reshape(-1, 6)  # Shape (n_points, 6)

    # Invert stiffness tensors to get compliance tensors
    S_batch = np.linalg.inv(C_batch)  # Shape (n, 6, 6)

    # Compute E for each tensor in the batch
    NSN = np.einsum("pi,nij,pj->np", N_flat, S_batch, N_flat)  # Shape (n, n_points)
    E_batch = 1.0 / NSN  # Shape (n, n_points)

    # Reshape E_batch back to (n, n_theta, n_phi)
    E_batch = E_batch.reshape(n, *d_x.shape)

    X_batch = E_batch * d_x  # Shape (n, n_theta, n_phi)
    Y_batch = E_batch * d_y
    Z_batch = E_batch * d_z

    return X_batch, Y_batch, Z_batch, E_batch


def plot_YoungsModulus3D(C, title="Young's Modulus Surface"):
    """
    Plot a 3D surface of Young's modulus.

    Args:
        C (ndarray): Stiffness tensor in Mandel notation. Can be a single tensor of shape (6,6) or a batch of tensors of shape (n,6,6).
        title (str): Title of the plot.

    Raises:
        ValueError: If C is not of shape (6,6) or (1,6,6).
    """
    if C.shape == (6, 6):
        C_batch = C[np.newaxis, :, :]
    elif C.shape == (1, 6, 6):
        C_batch = C
    else:
        raise ValueError(
            "C must be either a (6,6) tensor or a batch with one tensor of shape (1,6,6)."
        )

    X_batch, Y_batch, Z_batch, E_batch = compute_YoungsModulus3D(C_batch)
    X, Y, Z, E = X_batch[0], Y_batch[0], Z_batch[0], E_batch[0]

    surface = go.Surface(x=X, y=Y, z=Z, surfacecolor=E, colorscale="Viridis")
    layout = go.Layout(
        title=title,
        scene=dict(
            xaxis=dict(title="X"),
            yaxis=dict(title="Y"),
            zaxis=dict(title="Z"),
            aspectmode="auto",
        ),
    )

    fig = go.Figure(data=[surface], layout=layout)
    fig.show()


def export_YoungsModulus3D_to_vtk(C, prefix="youngs_modulus_surface"):
    """
    Export the computed Young's modulus surfaces to VTK files for Paraview visualization.

    Args:
        C (ndarray): Stiffness tensor in Mandel notation. Can be a single tensor of shape (6,6) or a batch of tensors of shape (n,6,6).
        prefix (str): Prefix for the output files.

    Returns:
        None
    """
    X_batch, Y_batch, Z_batch, E_batch = compute_YoungsModulus3D(C)
    n, n_theta, n_phi = X_batch.shape

    for k in range(n):
        points = np.vstack(
            (X_batch[k].ravel(), Y_batch[k].ravel(), Z_batch[k].ravel())
        ).T
        cells = [
            (
                "quad",
                np.array(
                    [
                        [
                            i * n_phi + j,
                            (i + 1) * n_phi + j,
                            (i + 1) * n_phi + (j + 1),
                            i * n_phi + (j + 1),
                        ]
                        for i in range(n_theta - 1)
                        for j in range(n_phi - 1)
                    ],
                    dtype=np.int32,
                ),
            )
        ]
        mesh = meshio.Mesh(
            points=points,
            cells=cells,
            point_data={"Youngs_Modulus": E_batch[k].ravel()},
        )
        filename = f"{prefix}_{k}.vtk"
        meshio.write(filename, mesh)
        print(f"Exported {filename}")


def demoCubic():
    """
    Demonstrates the Young's modulus surface plotting routine for a cubic material (Copper).

    This function generates the stiffness tensor for a cubic material, specifically copper,
    and then plots the 3D Young's modulus surface using the generated tensor.

    Args:
        None

    Returns:
        None
    """
    P1 = np.zeros((6, 6))
    P1[:3, :3] = 1.0 / 3.0
    D = np.diag([1, 1, 1, 0, 0, 0])
    P2 = D - P1
    P3 = np.eye(6) - D

    # generate stiffness for a cubic material: copper
    l1, l2, l3 = 136.67, 46, 150
    C = 3 * l1 * P1 + l2 * P2 + l3 * P3

    # show the 3D Young's modulus plot for copper
    plot_YoungsModulus3D(C, title="Young's Modulus Surface for Copper")
