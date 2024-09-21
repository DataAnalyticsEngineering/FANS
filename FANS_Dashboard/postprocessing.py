import numpy as np


def compute_rank2tensor_measures(tensor_matrix, measures_to_compute=None):
    """
    Computes various tensor measures from a given stress or strain tensor in Mandel notation.
    The user can specify which measures to compute. This function supports input tensors with arbitrary leading dimensions,
    as long as the last dimension is 6 (Mandel notation).

    Based on : https://doc.comsol.com/5.5/doc/com.comsol.help.sme/sme_ug_theory.06.16.html

    Parameters:
    - tensor_matrix: numpy array, tensor (stress or strain) in Mandel notation with shape (..., 6).
                     The tensor should be organized as follows:
                     [s11, s22, s33, s12, s13, s23] for stress or
                     [e11, e22, e33, e12, e13, e23] for strain.
    - measures_to_compute: list of strings, optional, specifying which measures to compute.
                           If not provided, default measures ['von_mises', 'hydrostatic', 'deviatoric'] will be computed.
                           Available options include:
                           - 'von_mises': Computes the von Mises stress/strain.
                           - 'hydrostatic': Computes the hydrostatic stress/strain.
                           - 'deviatoric': Computes the deviatoric stress/strain.
                           - 'principal': Computes the principal stresses/strains (eigenvalues).
                           - 'max_shear': Computes the maximum shear stress/strain.
                           - 'I_invariants': Computes the I1, I2, I3 invariants.
                           - 'J_invariants': Computes the J1, J2, J3 invariants of the deviatoric tensor.
                           - 'eigenvalues': Computes the eigenvalues of the stress/strain tensor.
                           - 'eigenvectors': Computes the eigenvectors of the stress/strain tensor.
                           - 'lode_angle': Computes the Lode angle, useful in advanced plasticity models.

    Returns:
    - result: dictionary, keys are the requested measure names and values are the computed measures.
              Each returned measure will have the same leading dimensions as the input tensor_matrix,
              with the last dimension adjusted based on the measure (e.g., eigenvalues will have an extra dimension for the 3 components).
    """
    if measures_to_compute is None:
        measures_to_compute = ["von_mises", "hydrostatic", "deviatoric"]

    original_shape = tensor_matrix.shape[:-1]  # All dimensions except the last one
    tensor_matrix = tensor_matrix.reshape(-1, 6)  # Flatten to (N, 6) for processing

    result = {}

    # Hydrostatic stress/strain (mean of the diagonal components)
    hydrostatic = np.mean(tensor_matrix[:, :3], axis=1)
    if "hydrostatic" in measures_to_compute:
        result["hydrostatic"] = hydrostatic.reshape(original_shape)

    # Deviatoric stress/strain and von Mises stress/strain
    deviatoric = tensor_matrix[:, :3] - hydrostatic[:, np.newaxis]
    deviatoric_shear = tensor_matrix[:, 3:6]
    deviatoric_tensor = np.hstack([deviatoric, deviatoric_shear])
    if "deviatoric" in measures_to_compute:
        result["deviatoric"] = deviatoric_tensor.reshape(original_shape + (6,))

    if "von_mises" in measures_to_compute:
        von_mises = np.sqrt(
            0.5
            * (
                (deviatoric[:, 0] - deviatoric[:, 1]) ** 2
                + (deviatoric[:, 1] - deviatoric[:, 2]) ** 2
                + (deviatoric[:, 2] - deviatoric[:, 0]) ** 2
                + 6
                * (
                    deviatoric_shear[:, 0] ** 2
                    + deviatoric_shear[:, 1] ** 2
                    + deviatoric_shear[:, 2] ** 2
                )
            )
        )
        result["von_mises"] = von_mises.reshape(original_shape)

    # Compute I1, I2, I3 invariants if requested
    if "I_invariants" in measures_to_compute:
        I1 = np.sum(tensor_matrix[:, :3], axis=1)
        I2 = (
            tensor_matrix[:, 0] * tensor_matrix[:, 1]
            + tensor_matrix[:, 1] * tensor_matrix[:, 2]
            + tensor_matrix[:, 2] * tensor_matrix[:, 0]
            - tensor_matrix[:, 3] ** 2
            - tensor_matrix[:, 4] ** 2
            - tensor_matrix[:, 5] ** 2
        )
        if "full_tensor" not in locals():
            full_tensor = mandel_to_matrix(tensor_matrix)
        I3 = np.linalg.det(full_tensor)
        result["I_invariants"] = np.stack([I1, I2, I3], axis=-1).reshape(
            original_shape + (3,)
        )

    # Compute J1, J2, J3 invariants if requested
    if "J_invariants" in measures_to_compute or "lode_angle" in measures_to_compute:
        J1 = np.sum(deviatoric_tensor[:, :3], axis=1)
        J2 = 0.5 * np.sum(deviatoric**2 + 2 * deviatoric_shear**2, axis=1)
        full_deviatoric_tensor = mandel_to_matrix(deviatoric_tensor)
        J3 = np.linalg.det(full_deviatoric_tensor)
        result["J_invariants"] = np.stack([J1, J2, J3], axis=-1).reshape(
            original_shape + (3,)
        )

    # Principal stresses/strains, maximum shear, eigenvalues, and eigenvectors
    if any(
        measure in measures_to_compute
        for measure in ["principal", "max_shear", "eigenvalues", "eigenvectors"]
    ):
        full_tensor = mandel_to_matrix(tensor_matrix)
        eigenvalues, eigenvectors = np.linalg.eigh(full_tensor)
        if "principal" in measures_to_compute:
            result["principal"] = eigenvalues.reshape(original_shape + (3,))
        if "max_shear" in measures_to_compute:
            max_shear = 0.5 * (eigenvalues[:, 2] - eigenvalues[:, 0])
            result["max_shear"] = max_shear.reshape(original_shape)
        if "eigenvalues" in measures_to_compute:
            result["eigenvalues"] = eigenvalues.reshape(original_shape + (3,))
        if "eigenvectors" in measures_to_compute:
            result["eigenvectors"] = eigenvectors.reshape(original_shape + (3, 3))

    # Lode angle calculation
    if "lode_angle" in measures_to_compute:
        if "J2" not in locals():  # Compute J2 if not already computed
            J2 = 0.5 * np.sum(deviatoric**2 + 2 * deviatoric_shear**2, axis=1)
        if "J3" not in locals():  # Compute J3 if not already computed
            full_deviatoric_tensor = mandel_to_matrix(deviatoric_tensor)
            J3 = np.linalg.det(full_deviatoric_tensor)
        # Handle very small J2 values to prevent division by zero
        safe_J2 = np.where(J2 > 1e-12, J2, 1e-12)
        sqrt_3_3 = (3 * np.sqrt(3)) / 2
        cos_3theta = np.clip(sqrt_3_3 * (J3 / safe_J2 ** (3 / 2)), -1, 1)
        lode_angle = (1.0 / 3.0) * np.arccos(cos_3theta)
        result["lode_angle"] = lode_angle.reshape(original_shape)

    return result


def mandel_to_matrix(mandel_tensor):
    """
    Convert a tensor from Mandel notation to a full 3x3 matrix.

    Parameters:
    - mandel_tensor: numpy array, tensor in Mandel notation with shape (n_steps, 6).
                     The tensor should be organized as follows:
                     [s11, s22, s33, s12, s13, s23] for stress or
                     [e11, e22, e33, e12, e13, e23] for strain.

    Returns:
    - full_tensor: numpy array, tensor in full 3x3 matrix form with shape (n_steps, 3, 3).
    """
    full_tensor = np.zeros((mandel_tensor.shape[0], 3, 3))
    full_tensor[:, 0, 0] = mandel_tensor[:, 0]  # s11 or e11
    full_tensor[:, 1, 1] = mandel_tensor[:, 1]  # s22 or e22
    full_tensor[:, 2, 2] = mandel_tensor[:, 2]  # s33 or e33
    full_tensor[:, 0, 1] = full_tensor[:, 1, 0] = mandel_tensor[:, 3] / np.sqrt(
        2
    )  # s12 or e12
    full_tensor[:, 0, 2] = full_tensor[:, 2, 0] = mandel_tensor[:, 4] / np.sqrt(
        2
    )  # s13 or e13
    full_tensor[:, 1, 2] = full_tensor[:, 2, 1] = mandel_tensor[:, 5] / np.sqrt(
        2
    )  # s23 or e23
    return full_tensor


def matrix_to_mandel(full_tensor, tolerance=1e-8):
    """
    Convert a full 3x3 symmetric tensor to Mandel notation in a vectorized and efficient way.

    Parameters:
    - full_tensor: numpy array, tensor in full 3x3 matrix form with shape (n_steps, 3, 3).
    - tolerance: float, optional, tolerance for checking symmetry. Default is 1e-8.

    Returns:
    - mandel_tensor: numpy array, tensor in Mandel notation with shape (n_steps, 6).
                     The tensor will be organized as follows:
                     [s11, s22, s33, s12, s13, s23] for stress or
                     [e11, e22, e33, e12, e13, e23] for strain.

    Raises:
    - ValueError: if any of the tensors in the batch are not symmetric within the specified tolerance.
    """
    # Check if the tensors are symmetric within the given tolerance
    if not np.allclose(full_tensor, full_tensor.transpose(0, 2, 1), atol=tolerance):
        raise ValueError(
            "One or more tensors are not symmetric within the specified tolerance."
        )

    # Efficiently extract and scale the relevant components
    mandel_tensor = np.zeros((full_tensor.shape[0], 6))
    mandel_tensor[:, 0] = full_tensor[:, 0, 0]  # s11 or e11
    mandel_tensor[:, 1] = full_tensor[:, 1, 1]  # s22 or e22
    mandel_tensor[:, 2] = full_tensor[:, 2, 2]  # s33 or e33
    mandel_tensor[:, 3] = np.sqrt(2) * full_tensor[:, 0, 1]  # s12 or e12
    mandel_tensor[:, 4] = np.sqrt(2) * full_tensor[:, 0, 2]  # s13 or e13
    mandel_tensor[:, 5] = np.sqrt(2) * full_tensor[:, 1, 2]  # s23 or e23

    return mandel_tensor
