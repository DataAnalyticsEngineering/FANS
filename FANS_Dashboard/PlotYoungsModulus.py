import numpy as np
import plotly.graph_objs as go


def compute_3d_youngs_modulus(C):
    """
    Compute Young's modulus for all directions in 3D.

    Parameters:
    C : ndarray
        Stiffness tensor in Mandel notation.

    Returns:
    E: ndarray
        Young's modulus in all directions.
    X, Y, Z: ndarrays
        Coordinates for plotting the modulus surface.
    """

    n_theta = 180
    n_phi = 360
    theta = np.linspace(0, np.pi, n_theta)  # Polar angle
    phi = np.linspace(0, 2 * np.pi, n_phi)  # Azimuthal angle

    S = np.linalg.inv(C)

    E = np.zeros((n_theta, n_phi))

    for i in range(n_theta):
        for j in range(n_phi):
            d = np.array(
                [
                    np.sin(theta[i]) * np.cos(phi[j]),
                    np.sin(theta[i]) * np.sin(phi[j]),
                    np.cos(theta[i]),
                ]
            )

            N = np.array(
                [
                    d[0] ** 2,
                    d[1] ** 2,
                    d[2] ** 2,
                    np.sqrt(2.0) * d[0] * d[1],
                    np.sqrt(2.0) * d[0] * d[2],
                    np.sqrt(2.0) * d[2] * d[1],
                ]
            )

            E[i, j] = 1.0 / (N.T @ S @ N)

    X = E * np.sin(theta)[:, np.newaxis] * np.cos(phi)[np.newaxis, :]
    Y = E * np.sin(theta)[:, np.newaxis] * np.sin(phi)[np.newaxis, :]
    Z = E * np.cos(theta)[:, np.newaxis]

    return X, Y, Z, E


def plot_3d_youngs_modulus_surface(C, title="Young's Modulus Surface"):
    """
    Plot a 3D surface of Young's modulus.

    Parameters:
    C : ndarray
        Stiffness tensor in Mandel notation.
    title : str
        Title of the plot.

    """
    X, Y, Z, E = compute_3d_youngs_modulus(C)

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


def demoCubic():
    """
    Demonstrates the Young's modulus surface plotting routine for a cubic material (Copper)

    Returns
    -------
    None.

    """
    P1 = np.zeros((6, 6))
    P1[:3, :3] = 1.0 / 3.0
    D = np.diag([1, 1, 1, 0, 0, 0])
    P2 = D - P1
    P3 = np.eye(6) - D

    # generate stiffness for a cubic material: copper
    l1, l2, l3 = 136.67, 46, 150
    C = 3 * l1 * P1 + l2 * P2 + l3 * P3

    print(C)

    # show the 3D Young's modulus plot for copper
    plot_3d_youngs_modulus_surface(C, title="Young's Modulus Surface for Copper")
