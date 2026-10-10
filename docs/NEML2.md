# NEML2 materials in FANS

Any [NEML2](https://github.com/applied-material-modeling/neml2) material model can be the material of a phase: `"matmodel": "NEML2"`.

## Capabilities

- Thermal, small-strain and large-strain problems: the model maps the temperature gradient to the heat flux, the strain to the stress, or the deformation gradient to the first Piola stress.
- Any model `neml2-compile` can compile, e.g. from NEML2's own building blocks (J2 or crystal plasticity with an implicit update).
- Your own model is a few lines of [PyTorch](https://pytorch.org/), automatic differentiation included (e.g. a hyperelastic law from just its strain energy): once compiled, it runs in FANS as is.
- CPU or GPU. On a GPU all MPI ranks share the card.
- History variables are found, stored and written by FANS (`"internal_variables"` in `results`).
- Further model inputs ("fields", e.g. a grain orientation) come from the microstructure file, per voxel or per phase.
- Time-dependent models get the time step of the input file (`time_step`).
- Native and NEML2 materials can be mixed in one simulation.
- NEML2 and libtorch never reach FANS itself: the material is a plugin library, `libfans_neml2.so`, loaded at run time.

On the CPU, expect a NEML2 material to run about 2 to 6 times slower than the same model written natively in FANS.

## Build

```bash
pixi run -e dev-neml2 build-fans_neml2
```

This gives `test/FANS_neml2`. Without pixi, configure with `-DFANS_NEML2=ON` in a Python environment that has `neml2` (3.1 or newer) and PyTorch installed.

## Compile a model

```bash
pixi run -e dev-neml2 neml2-compile model.i --model model --device cpu cuda --dtype float64 --output-dir compiled_models -d state/S:forces/E
```

This compiles the block `[model]` of `model.i` for both devices into `compiled_models/model`. `--load a.py` (repeatable) imports a Python file the model needs. `-d flux:gradient` also compiles the model's tangent, which the `homogenized_tangent` result needs.

## Input file

```json
{
    "phases": [0, 1],
    "matmodel": "NEML2",
    "material_properties": {
        "artifact": "compiled_models/model",
        "gradient": "forces/E",
        "flux": "state/S",
        "device": "cuda",
        "fields": {"orientation": "rotation_matrices"}
    }
}
```

| Property | Meaning | Default |
| --- | --- | --- |
| `artifact` | Folder of the compiled model | required |
| `gradient`, `flux` | Names of the model's gradient input and flux output | required |
| `device` | `"cpu"`, `"cuda"`, ... | `"cpu"` |
| `batch_size` | Material points per evaluation | 1024 on the CPU, 65536 on a GPU |
| `fields` | For every further model input, the dataset to read it from: next to the microstructure, or an absolute path; shaped `[Z][Y][X][...]` (per voxel) or `[n_phase][...]` (per phase) | required if the model has such inputs |
| `linear` | `true` if the flux is linear in the gradient: FANS then solves with its linear CG. | `false` |

The input file also needs a top-level `reference_material`, the reference stiffness of the fundamental solution.
