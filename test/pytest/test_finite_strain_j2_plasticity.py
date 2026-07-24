import json
import os

import numpy as np
import pytest
from MSUtils.fans_dashboard.utils import identify_hierarchy, extract_and_organize_data


@pytest.fixture
def finite_strain_files(request):
    test_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(
        test_dir,
        "../input_files/test_FiniteStrainJ2Plasticity.json",
    )
    if request.config.getoption("--from-pixi"):
        h5_path = os.path.join(
            test_dir,
            "../output/test_FiniteStrainJ2Plasticity.h5",
        )
    else:
        h5_path = os.path.join(
            test_dir,
            "../../build/test/test_FiniteStrainJ2Plasticity.h5",
        )

    assert os.path.exists(json_path), f"Required test input not found: {json_path}"
    assert os.path.exists(h5_path), f"Required test output not found: {h5_path}"
    return json_path, h5_path


def _symmetric_tensor_function(tensor, function):
    tensor = 0.5 * (tensor + tensor.T)
    eigenvalues, eigenvectors = np.linalg.eigh(tensor)
    return (eigenvectors * function(eigenvalues)) @ eigenvectors.T


def _constitutive_oracle(load_path, bulk, shear, yield_stress, hardening):
    """Independent transcription of GooseFFT's finite-strain return map."""
    identity = np.eye(3)
    previous_F = identity.copy()
    previous_be = identity.copy()
    previous_ep = 0.0

    first_piola = []
    equivalent_plastic_strain = []
    elastic_finger = []

    for flattened_F in load_path:
        F = np.asarray(flattened_F, dtype=float).reshape(3, 3)
        F_increment = F @ np.linalg.inv(previous_F)
        be_predictor = F_increment @ previous_be @ F_increment.T
        log_be_predictor = _symmetric_tensor_function(be_predictor, np.log)

        trace_log_be = np.trace(log_be_predictor)
        dev_log_be = log_be_predictor - trace_log_be / 3.0 * identity
        tau_predictor = 0.5 * bulk * trace_log_be * identity + shear * dev_log_be
        dev_tau = tau_predictor - np.trace(tau_predictor) / 3.0 * identity
        equivalent_tau = np.sqrt(1.5 * np.sum(dev_tau * dev_tau))
        yield_function = equivalent_tau - (yield_stress + hardening * previous_ep)

        tau = tau_predictor
        log_be = log_be_predictor
        ep = previous_ep
        tolerance = 1.0e-12 * max(1.0, equivalent_tau, yield_stress + hardening * ep)
        if yield_function > tolerance:
            flow_direction = 1.5 * dev_tau / equivalent_tau
            plastic_increment = yield_function / (hardening + 3.0 * shear)
            ep += plastic_increment
            tau = tau_predictor - 2.0 * shear * plastic_increment * flow_direction
            log_be = log_be_predictor - 2.0 * plastic_increment * flow_direction

        be = _symmetric_tensor_function(log_be, np.exp)
        P = tau @ np.linalg.inv(F).T

        first_piola.append(P.reshape(9))
        equivalent_plastic_strain.append(ep)
        elastic_finger.append(be.reshape(9))
        previous_F = F
        previous_be = be
        previous_ep = ep

    return (
        np.asarray(first_piola),
        np.asarray(equivalent_plastic_strain),
        np.asarray(elastic_finger),
    )


def test_finite_strain_response_matches_oracle(
    finite_strain_files,
):
    json_path, h5_path = finite_strain_files
    with open(json_path, encoding="utf-8") as input_file:
        input_data = json.load(input_file)

    properties = input_data["materials"][0]["material_properties"]
    material_parameters = {
        "bulk": properties["bulk_modulus"][0],
        "shear": properties["shear_modulus"][0],
        "yield_stress": properties["yield_stress"][0],
        "hardening": properties["isotropic_hardening_parameter"][0],
    }

    load_paths = input_data["macroscale_loading"]
    hierarchy = identify_hierarchy(h5_path)
    microstructures = list(hierarchy)
    load_cases = [f"load{i}" for i in range(len(load_paths))]
    results = extract_and_organize_data(
        h5_path,
        hierarchy,
        [
            "strain_average",
            "stress_average",
            "equivalent_plastic_strain",
            "elastic_finger",
        ],
        microstructures,
        load_cases,
        [],
    )

    for load_index, load_path in enumerate(load_paths):
        expected_stress, expected_ep, expected_be = _constitutive_oracle(
            load_path, **material_parameters
        )
        load_case = f"load{load_index}"

        for microstructure in microstructures:
            actual = results[microstructure][load_case]
            actual_ep = actual["equivalent_plastic_strain"][..., 0]
            actual_be = actual["elastic_finger"]

            assert np.allclose(
                actual["strain_average"],
                load_path,
                rtol=1.0e-10,
                atol=1.0e-12,
            )
            assert np.allclose(
                actual["stress_average"],
                expected_stress,
                rtol=1.0e-10,
                atol=1.0e-11,
            )
            assert np.allclose(
                actual_ep,
                expected_ep[:, None, None, None],
                rtol=1.0e-10,
                atol=1.0e-12,
            )
            assert np.allclose(
                actual_be,
                expected_be[:, None, None, None, :],
                rtol=1.0e-10,
                atol=1.0e-12,
            )
