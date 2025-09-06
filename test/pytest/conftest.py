import pytest
import os

def pytest_addoption(parser):
    parser.addoption(
        "--from-pixi", action="store_true", default=False,
        help="Run tests from pixi task."
    )

@pytest.fixture(
    params=[
        "test_J2Plasticity",
        "test_LinearElastic",
        "test_LinearThermal",
        "test_PseudoPlastic",
    ]
)
def test_files(request):
    json_base_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "../input_files/"
    )
    # Determine the HDF5 base directory based on the --from-pixi option
    # if run from pixi, use output directory; otherwise, use build/test directory, as it is the one used by ctest
    if request.config.getoption("--from-pixi"):
        h5_base_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "../output"
        )
    else:
        h5_base_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "../../build/test/"
        )

    json_path = os.path.join(json_base_dir, f"{request.param}.json")
    h5_path = os.path.join(h5_base_dir, f"{request.param}.h5")

    if os.path.exists(json_path) and os.path.exists(h5_path):
        return json_path, h5_path
    pytest.skip(f"Required test files not found: {json_path} or {h5_path}")