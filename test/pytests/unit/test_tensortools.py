import numpy as np
import pytest
from fans_dashboard.core.tensortools import Full2Mandel, Mandel2Full
from fans_dashboard.core.tensortools import VoigtStrain2Mandel, VoigtStress2Mandel
from fans_dashboard.core.tensortools import Mandel2VoigtStrain, Mandel2VoigtStress
from fans_dashboard.core.tensortools import Ciso, IsoProjectionC, IsoProjectionKappa


def test_Conversion():
    """Test various conversions, e.g., Full tensor <-> Mandel, Mandel <-> Voigt (both orderings)"""
    n_test = 10
    for i in range(n_test):
        # generate random tensor
        A = np.random.uniform(-1, 1, size=(3, 3))
        A = A + A.T
        A_m = Full2Mandel(A)
        assert (
            np.linalg.norm(A - Mandel2Full(A_m)) < 1.0e-10
        ), "Full->Mandel->Full failed for " + str(A)
        A_m = np.random.uniform(-1, 1, size=(6,))
        assert (
            np.linalg.norm(A_m - Full2Mandel(Mandel2Full(A_m))) < 1.0e-10
        ), "Mandel->Full->Mandel failed for A_m=" + str(A_m)

        assert (
            np.linalg.norm(A_m - VoigtStrain2Mandel(Mandel2VoigtStrain(A_m))) < 1.0e-10
        ), "Mandel->Voigt_eps->Mandel failed for A_m=" + str(A_m)
        assert (
            np.linalg.norm(A_m - VoigtStress2Mandel(Mandel2VoigtStress(A_m))) < 1.0e-10
        ), "Mandel->Voigt_sig->Mandel failed for A_m=" + str(A_m)

        assert (
            np.linalg.norm(
                A_m
                - VoigtStrain2Mandel(Mandel2VoigtStrain(A_m, order="abq"), order="abq")
            )
            < 1.0e-10
        ), "Mandel->Voigt_eps,ABQ->Mandel failed for A_m=" + str(A_m)
        assert (
            np.linalg.norm(
                A_m
                - VoigtStress2Mandel(Mandel2VoigtStress(A_m, order="abq"), order="abq")
            )
            < 1.0e-10
        ), "Mandel->Voigt_sig,ABQ->Mandel failed for A_m=" + str(A_m)


def test_Ciso():
    """test the isotropic projection for 4-tensors using randomized data"""
    n_test = 10
    G = np.random.uniform(0, 10, size=n_test)
    K = np.random.uniform(0, 10, size=n_test)
    for k, g in zip(K, G):
        C = Ciso(k, g)
        k_fit, g_fit = IsoProjectionC(C)
        assert (
            np.abs(k - k_fit) + np.abs(g - g_fit) < 1.0e-10
        ), f"Error in isotropic projection (4-tensor): K, G true = [{k}, {g}]; K, G projected = [{k_fit}, {g_fit}]"


def test_kappaiso():
    """test the isotropic projection for 2-tensors using randomized data"""
    n_test = 10
    Id_m = np.array((1.0, 1.0, 1.0, 0.0, 0.0, 0.0))
    for kappa in np.random.uniform(0, 10, size=n_test):
        kappa_fit = IsoProjectionKappa(kappa * Id_m)
        assert (
            np.abs(kappa - kappa_fit) < 1.0e-10
        ), f"Error in isotropic projection (2-tensor): kappa = {kappa}; kappa projected = {kappa_fit}"


if __name__ == "__main__":
    pytest.main(["-v", "-s", __file__])
