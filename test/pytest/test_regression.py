import json
import subprocess
import os
import pytest
from pathlib import Path

# tmp_path leads to all file being cleaned up
def test_datasetname_without_leading_slash(tmp_path):
    """
    Regression test:
    Removing the leading '/' from microstructure.datasetname
    must NOT cause FANS to crash.
    """

    # Locate original working JSON
    here = Path(__file__).resolve().parent
    input_json = (here / "../input_files/test_LinearElastic.json").resolve()
    assert input_json.exists(), "Reference LinearElastic JSON not found"

    # Load JSON
    with open(input_json) as f:
        data = json.load(f)

    # 1) Remove leading slash from datasetname (the bug trigger)
    ds = data["microstructure"]["datasetname"]
    assert ds.startswith("/"), "Test assumes original datasetname starts with '/'"
    data["microstructure"]["datasetname"] = ds.lstrip("/")

    # 2) Make microstructure filepath absolute
    ms_path = Path(data["microstructure"]["filepath"])
    if not ms_path.is_absolute():
        ms_path = ( here.parent / ms_path).resolve()
    assert ms_path.exists(), "Referenced microstructure HDF5 file not found"

    data["microstructure"]["filepath"] = str(ms_path)

    # Write modified JSON to tmp_path
    modified_json = tmp_path / "LinearElastic_no_slash.json"
    with open(modified_json, "w") as f:
        json.dump(data, f, indent=2)

    # Output file (also in tmp_path)
    output_h5 = tmp_path / "LinearElastic_no_slash.h5"

    # Locate executable
    fans_exec = (
        Path(os.environ["CONDA_PREFIX"]) / "bin/FANS"
        if "CONDA_PREFIX" in os.environ
        else Path("./FANS")
    )
    assert fans_exec.exists(), "FANS executable not found"

    # Run CLI
    cmd = [
        "mpiexec",
        "-n",
        "1",
        str(fans_exec),
        str(modified_json),
        str(output_h5),
    ]

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    # After the fix, this must succeed
    assert result.returncode == 0, (
        "FANS failed when datasetname lacked leading '/'\n"
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )

    assert output_h5.exists(), "Output HDF5 file was not created"
