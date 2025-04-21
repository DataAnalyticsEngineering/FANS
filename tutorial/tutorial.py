

import marimo

__generated_with = "0.13.0"
app = marimo.App(app_title="AiiDA-FANS Tutorial")


@app.cell(hide_code=True)
def _():
    import marimo as mo
    from pathlib import Path
    return Path, mo


@app.cell(hide_code=True)
def _(mo):
    nav_menu = mo.nav_menu(
        {
            "#aiida-setup": "AiiDA Setup",
            "#fans-rundown": "FANS Rundown",
            "#submitting-jobs": "Submitting Jobs",
            "#analysing-the-results": "Analysing the Results",
            "Links": {
                "https://github.com/ethan-shanahan/aiida-fans": "aiida-fans",
                "https://github.com/DataAnalyticsEngineering/FANS": "FANS",
                "https://www.aiida.net/": "AiiDA",
            },
        }
    )

    _tip = mo.md("""
    **Requirements:**

    The rest of this tutorial assumes you have read the attached README and have installed the requirements described therein. Although not foolproof, you can run the following commands to check if AiiDA, the plugin, and FANS are installed correctly.

    ```
    verdi plugin list aiida.calculations fans
    FANS
    ```

    Notice that we assume FANS is located on your PATH (or at least it is in your active  environment). While this is not necessary in general practice, the tutorial will continue under this assumption.
    """).callout("warn")

    mo.md(rf"""{nav_menu}

    ---

    # AiiDA-FANS Tutorial

    The goal of this tutorial is to give you an idea of how to utilise the `aiida-fans` plugin as well as an introduction to `AiiDA` and `FANS`. By the end of this tutorial, you should know how to:

    - Setup your AiiDA profile, computer, and code.
    - Define FANS options and prepare a parameter space study. 
    - Write a `submit.py` script to run your jobs.
    - Query and read the results.

    {_tip}
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    _note = mo.md(r"""
    **Note:** _not your first profile..._

    This section assumes you have not already set up an appropriate profile, computer, and code for using AiiDA and FANS. If you have already done this, you may wish to skip to the next section.

    However, this tutorial is designed to work with a blank profile specifically.
    """).callout("info")

    mo.md(rf"""
    ## AiiDA Setup

    Before we can truly begin, we must set up AiiDA on your machine. This means three things.

    1. Create a Profile
    2. Specify a Computer
    3. Define a Code

    AiiDA has multiple user interfaces but their CLI, `verdi`, is particularly well suited to these three steps since they need to be performed only rarely. Therefore, you will need access to the terminal to proceed.

    {_note}
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### 1. Create a Profile

        By default, AiiDA stores app data at the user level. Even when AiiDA is installed in a virtual environment, it will still read and write to `.aiida` in your home directory. However, AiiDA provides users a way to seperate their data into "profiles". Let's create a profile for this tutorial.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    profile_settings = mo.hstack(
        [
            mo.vstack([
                "Profile Name:",
                "First Name:",
                "Last Name:",
                "Email:",
                "Institution:",
            ], align="start", heights="equal", gap=0.8),
            mo.vstack([
                "{profile_name}",
                "{first_name}",
                "{last_name}",
                "{email}",
                "{institution}",
            ], align="start", heights="equal", gap=0.5)
        ],
        justify="center", align="stretch", gap=2.0,
    ).batch(
        profile_name=mo.ui.text("aiida-fans-tutorial"),
        first_name=mo.ui.text("Max"),
        last_name=mo.ui.text("Mustermann"),
        email=mo.ui.text("example@nomail.com"),
        institution=mo.ui.text("MIB"),
    ).form(
        show_clear_button=True, clear_button_label="Reset", bordered=True
    )

    mo.vstack([
        mo.md("**Fill in the details below to generate your custom profile configuration.**"),
        profile_settings
    ], align="center")
    return (profile_settings,)


@app.cell(hide_code=True)
def _(mo, profile_settings):
    mo.stop(
        profile_settings.value is None,
        mo.status.spinner(title="Awaiting input above ...", remove_on_exit=False)
    )

    profile_config = \
    rf"""profile: {profile_settings.value["profile_name"]}
    first_name: {profile_settings.value["first_name"]}
    last_name: {profile_settings.value["last_name"]}
    email: {profile_settings.value["email"]}
    institution: {profile_settings.value["institution"]}
    use_rabbitmq: false
    set_as_default: true
    non_interactive: true
    """

    with open("configure_profile.yaml", "w") as _f:
        _f.write(profile_config)

    mo.md(f"""
    With the values you input above, a `configure_profile.yaml` has been automatically written to the working directory. It contains the following data:

    ```yaml
    {profile_config}
    ```
    """).callout(kind="success")
    return


@app.cell(hide_code=True)
def _(mo, profile_settings):
    _note = mo.md(rf"""
    **Note:** _on default profiles..._

    We have made the new profile our default profile. This means that any further calls to `verdi` will implicitly use the {"<profile_name>" if profile_settings.value is None else profile_settings.value["profile_name"]} profile. You can change the profile on a per-call basis with the `-p/--profile` option. To change the default profile use:

    ```
    verdi profile set-default <profile_name>
    ```
    """).callout(kind="info")

    mo.md(f"""
    To create your new profile from this file run:

    ```
    verdi profile setup core.sqlite_dos --config configure_profile.yaml
    ```

    Hopefully, that completed successfully. Using these commands, you should see your new profile listed (alone if this is your first profile) and a report on it also:

    ```
    verdi profile list
    verdi profile show {"<profile_name>" if profile_settings.value is None else profile_settings.value["profile_name"]}
    ```

    {_note}
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### 2. Specify a Computer

        Before you proceed, ensure that your local computer satisfies the following requirements:

        - it runs a Unix-like operating system (Linux distros and MacOS should work fine)
        - it has `bash` installed

        AiiDA does not assume what computer you wish to run jobs on, so even if you are only using your local machine, you must tell it as much. That is what we will do here; specify the localhost computer.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    computer_settings = mo.hstack(
        [
            mo.vstack([
                "Computer Label:",
                "MPI processes:",
                "Description:",
            ], align="start", heights="equal", gap=0.8),
            mo.vstack([
                "{label}",
                "{mpiprocs}",
                "{description}",
            ], align="start", heights="equal", gap=0.5)
        ],
        justify="center", align="stretch", gap=2.0,
    ).batch(
            label=mo.ui.text("localhost"),
            mpiprocs=mo.ui.text("2"),
            description=mo.ui.text_area("This is my local machine."),
    ).form(
        show_clear_button=True, clear_button_label="Reset", bordered=True
    )

    mo.vstack([
        mo.md("**Fill in the details below to generate your custom computer configuration.**"),
        computer_settings
    ], align="center")
    return (computer_settings,)


@app.cell(hide_code=True)
def _(Path, computer_settings, mo):
    mo.stop(
        computer_settings.value is None,
        mo.status.spinner(title="Awaiting input above ...", remove_on_exit=False)
    )

    computer_config = \
    rf"""label: {computer_settings.value["label"]}
    description: {computer_settings.value["description"]}
    hostname: localhost
    transport: core.local
    scheduler: core.direct
    shebang: #!/bin/bash
    work_dir: {Path.cwd()}/.aiida_run""" + r"""
    mpirun_command: mpiexec -n {tot_num_mpiprocs}""" + rf"""
    mpiprocs_per_machine: {computer_settings.value["mpiprocs"]}
    default_memory_per_machine: null
    use_double_quotes: false
    prepend_text: ' '
    append_text: ' '
    non_interactive: true
    """

    with open("configure_computer.yaml", "w") as _f:
        _f.write(computer_config)

    mo.md(f"""
    With the values you input above, a `configure_computer.yaml` has been automatically written to the working directory. It contains the following data:

    ```yaml
    {computer_config}
    ```
    """).callout(kind="success")
    return


@app.cell(hide_code=True)
def _(computer_settings, mo):
    mo.md(rf"""
    To specify your new computer from this file run:

    ```
    verdi computer setup --config configure_computer.yaml
    ```

    Then you must configure the computer with the following command:

    ```
    verdi computer configure core.local {"<computer_label>" if computer_settings.value is None else computer_settings.value["label"]}
    ```

    The default options should be suitable.

    Hopefully, that completed successfully. Using this command, you should test that AiiDA can connect to the machine:

    ```
    verdi computer test {"<computer_label>" if computer_settings.value is None else computer_settings.value["label"]}
    ```
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### 3. Define a Code

        The final step to setup AiiDA is to define the "code" you wish to utilise. Here, the "code" refers to FANS. This step is important as it tells AiiDA how to execute FANS and which plugin should handle its jobs. AiiDA provides many ways of handling the "code" of your project. Since we installed FANS in the environment, we can simply make use of it there.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    code_settings = mo.hstack(
        [
            mo.vstack([
                mo.vstack([
                    "Code Label:",
                    "Code Executable:",
                    "Description:",
                ], align="start", heights="equal", gap=0.8),
                    "Environment Activation Script:",
            ], align="start", heights="equal", gap=5.55),    
            mo.vstack([
                "{label}",
                "{executable}",
                "{description}",
                "{environment}",
            ], align="start", heights="equal", gap=0.5)
        ],
        justify="center", align="stretch", gap=2.0,
    ).batch(
        label=mo.ui.text("FANS"),
        executable=mo.ui.text("FANS"),
        description=mo.ui.text_area("The FANS executable."),
        environment=mo.ui.text_area("eval \"$(conda shell.bash hook)\"\nconda activate aiida-fans-tutorial"),
    ).form(
        show_clear_button=True, clear_button_label="Reset", bordered=True
    )

    mo.vstack([
        mo.md("**Fill in the details below to generate your custom code configuration.**"),
        code_settings
    ], align="center")
    return (code_settings,)


@app.cell(hide_code=True)
def _(code_settings, computer_settings, mo):
    mo.stop(
        code_settings.value is None or computer_settings.value is None,
        mo.status.spinner(title="Awaiting input above ...", remove_on_exit=False)
    )

    code_config = \
    rf"""label: {code_settings.value["label"]}
    description: {code_settings.value["description"]}
    default_calc_job_plugin: fans
    use_double_quotes: false
    with_mpi: true
    computer: {computer_settings.value["label"]}
    filepath_executable: {code_settings.value["executable"]}
    prepend_text: |
    {"\n".join([f"    {ln}" for ln in code_settings.value["environment"].split("\n")])}
    append_text: ' '
    non_interactive: true
    """

    with open("configure_code.yaml", "w") as _f:
        _f.write(code_config)

    mo.md(f"""
    With the values you input above, a `configure_code.yaml` has been automatically written to the working directory. It contains the following data:

    ```yaml
    {code_config}
    ```
    """).callout(kind="success")
    return


@app.cell(hide_code=True)
def _(code_settings, mo):
    _note = mo.md(r"""
    **Note:** _your first node..._

    You should also note that the code is saved by AiiDA as a node, and thus we have created our first node. Any calculation jobs we perform will be connected to this code node in the provenance graph.

    To list all the nodes stored in your profile, run:

    ```
    verdi node list
    ```
    """).callout(kind="info")

    mo.md(rf"""
    To define your new code from this file run:

    ```
    verdi code create core.code.installed --config configure_code.yaml
    ```

    Hopefully, that completed successfully. Using these commands, you can show the details of your new code and verify that AiiDA can connect to it:

    ```
    verdi code show {"<code_label>" if code_settings.value is None else code_settings.value["label"]}
    verdi code test {"<code_label>" if code_settings.value is None else code_settings.value["label"]}
    ```

    {_note}
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## FANS Rundown

        FANS requires a JSON input file. The input file can be thought of in 5 sections, each specifying the various problem parameters as well as runtime settings. Each setting also notes the appropriate AiiDA datatype. This is the type of node that you must give AiiDA when running jobs, as we will see later.

        ### Microstructure Definition

        ```json
        "ms_filename": "microstructures/sphere32.h5",
        "ms_datasetname": "/sphere/32x32x32/ms",
        "ms_L": [1.0, 1.0, 1.0]
        ```

        - `ms_filename`: This specifies the path to the HDF5 file that contains the microstructure data. (AiiDA type: `SinglefileData`)
        - `ms_datasetname`: This is the path within the HDF5 file to the specific dataset that represents the microstructure. (AiiDA type: `Str`)
        - `ms_L`: Microstructure length defines the physical dimensions of the microstructure in the x, y, and z directions. (AiiDA type: `List`)

        ### Problem Type and Material Model

        ```json
        "problem_type": "mechanical",
        "matmodel": "LinearElasticIsotropic",
        "material_properties": {
            "bulk_modulus": [62.5000, 222.222],
            "shear_modulus": [28.8462, 166.6667]
        }
        ```

        - `problem_type`: This defines the type of physical problem you are solving. Common options include "thermal" problems and "mechanical" problems. (AiiDA type: `Str`)
        - `matmodel`: This specifies the material model to be used in the simulation. Examples include `LinearThermalIsotropic` for isotropic linear thermal problems, `LinearElasticIsotropic` for isotropic linear elastic mechanical problems, `PseudoPlasticLinearHardening`/`PseudoPlasticNonLinearHardening` for plasticity mimicking model with linear/nonlinear hardening, and `J2ViscoPlastic_LinearIsotropicHardening`/ `J2ViscoPlastic_NonLinearIsotropicHardening` for rate dependent J2 plasticity model with linear/nonlinear isotropic hardening. (AiiDA type: `Str`)
        - `material_properties`: This provides the necessary material parameters for the chosen material model. For thermal problems, you might specify `conductivity`, while mechanical problems might require `bulk_modulus`, `shear_modulus`, and more properties for advanced material models. These properties can be defined as arrays to represent multiple phases within the microstructure. (AiiDA type: `Dict`)

        ### Solver Settings

        ```json
        "method": "cg",
        "error_parameters":{
            "measure": "Linfinity",
            "type": "absolute",
            "tolerance": 1e-10
        },
        "n_it": 100
        ```

        - `method`: This indicates the numerical method to be used for solving the system of equations. `cg` stands for the Conjugate Gradient method, and `fp` stands for the Fixed Point method. (AiiDA type: `Str`)
        - `error_parameters`: This section defines the error parameters for the solver. Error control is applied on the finite element nodal residual of the problem.
            - `measure`: Specifies the norm used to measure the error. Options include `Linfinity`, `L1`, or `L2`. (AiiDA type: `Str`)
            - `type`: Defines the type of error measurement. Options are `absolute` or `relative`. (AiiDA type: `Str`)
            - `tolerance`: Sets the tolerance level for the solver, defining the convergence criterion based on the chosen error measure. The solver iterates until the solution meets this tolerance. (AiiDA type: `Float`)
        - `n_it`: Specifies the maximum number of iterations allowed for the FANS solver. (AiiDA type: `Int`)


        ### Macroscale Loading Conditions

        ```json
        "macroscale_loading":   [
                                    [
                                        [0.004, -0.002, -0.002, 0, 0, 0],
                                        [0.008, -0.004, -0.004, 0, 0, 0],
                                        [0.012, -0.006, -0.006, 0, 0, 0],
                                        [0.016, -0.008, -0.008, 0, 0, 0],
                                    ],
                                    [
                                        [0, 0, 0, 0.002, 0, 0],
                                        [0, 0, 0, 0.004, 0, 0],
                                        [0, 0, 0, 0.006, 0, 0],
                                        [0, 0, 0, 0.008, 0, 0],
                                    ]
                                ]
        ```

        - `macroscale_loading`: This defines the external loading applied to the microstructure. It is an array of arrays, where each sub-array represents a loading condition applied to the system. The format of the loading array depends on the problem type (AiiDA type: `ArrayData`):
            - For `thermal` problems, the array typically has 3 components, representing the temperature gradients in the x, y, and z directions.
            - For `mechanical` problems, the array must have 6 components, corresponding to the components of the strain tensor in Mandel notation (e.g., $[[ε_{11}, ε_{22}, ε_{33}, \sqrt{2} ε_{12}, \sqrt{2} ε_{13}, \sqrt{2} ε_{23}]]$).

        In the case of path/time-dependent loading as shown, for example as in plasticity problems, the `macroscale_loading` array can include multiple steps with corresponding loading conditions.

        ### Results Specification

        ```json
        "results": [
            "stress", "strain",
            "stress_average", "strain_average",
            "phase_stress_average", "phase_strain_average",
            "microstructure",
            "displacement",
            "absolute_error",
        ]
        ```

        - `results`: This array lists the quantities that should be stored into the results HDF5 file during the simulation. Each string in the array corresponds to a specific result (AiiDA type: `List`):
            - `stress` and `strain`: The stress and strain fields at each voxel in the microstructure.
            - `stress_average` and `strain_average`: Volume averaged- homogenized stress and strain over the entire microstructure.
            - `phase_stress_average` and `phase_strain_average`: Volume averaged- homogenized stress and strain for each phase within the microstructure.
            - `microstructure`: The original microstructure data.
            - `displacement`: The displacement fluctuation field (for mechanical problems) and temperature fluctuation field (for thermal problems).
            - `absolute_error`: The L-infinity error of finite element nodal residual at each iteration.

        Additional material model specific results can be included depending on the problem type and material model.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Submitting Jobs

        Now that AiiDA is suitably prepared and we're familiar with the FANS parameter specifications, its time to get to work. We will conduct a mock experiment to demonstrate the simplicity and flexibility that using the plugin offers. Breaking down the submission of jobs into two steps makes for a clean workflow.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    import_button = mo.ui.run_button(label="RUN")

    mo.md(rf"""
    ### Creating Input Parameters

    We will create all the input parameters we wish to study today at once. To begin, we import everything that will be needed and call the `load_profile()` function to activate the default profile within this script.

    Press this button only after you have created a default profile as described above in [AiiDA Setup](#aiida-setup).

    {import_button}

    ```py
    from aiida.engine import run                  # run jobs
    from aiida.plugins import CalculationFactory  # generates fans calculation
    from aiida.orm import (
        Group,                                    # node organisation tool
        SinglefileData,                           # \
        Str,                                      # |
        Float,                                    # |
        Int,                                      # |- AiiDA datatypes
        List,                                     # |
        Dict,                                     # |
        ArrayData,                                # /
        CalcJobNode,                              # node type for calculation jobs
        QueryBuilder,                             # advanced query tool
        load_node,                                # basic query tool for nodes
        load_code,                                # basic query tool for codes
    )
    from numpy import array                       # numpy array
    from itertools import product                 # for parameter space generation
    from random import uniform                    # for parameter space generation

    from aiida import load_profile                # injects profile context into script
    load_profile()
    ```
    """)
    return (import_button,)


@app.cell(hide_code=True)
def imports(import_button, mo):
    mo.stop(not import_button.value)  # run on click

    try:
        from aiida.common.exceptions import ProfileConfigurationError, ConfigurationError

        from aiida.engine import run
        from aiida.plugins import CalculationFactory
        from aiida.orm import (
            Group,
            SinglefileData,
            Str,
            Float,
            Int,
            List,
            Dict,
            ArrayData,
            CalcJobNode,
            QueryBuilder,
            load_node,
            load_code,
        )
        from numpy import array
        from itertools import product
        from random import uniform

        from aiida import load_profile
        load_profile()

    except ImportError:
        mo.stop(True, output=mo.md("**Imports failed to load properly!**").style(text_align="center").callout(kind="danger"))

    except ProfileConfigurationError:
        mo.stop(True, output=mo.md("**Your profile failed to load properly!**").style(text_align="center").callout(kind="danger"))

    mo.md("**Success!**").style(text_align="center").callout(kind="success")
    return (
        ArrayData,
        CalcJobNode,
        CalculationFactory,
        ConfigurationError,
        Dict,
        Float,
        Group,
        Int,
        List,
        QueryBuilder,
        SinglefileData,
        Str,
        array,
        load_code,
        product,
        run,
        uniform,
    )


@app.cell(hide_code=True)
def _(mo):
    _code = r"""
    groups = QueryBuilder(                # we will use the advanced query method
    ).append(                             # it consists of a series of `.append` methods
        Group, filters={
            Group.fields.label: "inputs"  # here we filter by groups labeled "inputs"
        }
    ).all(                                # conclude with a method to fetch the results
        flat=True
    )

    if len(groups) == 0:                  # if the "inputs" group does not exist...
        inputs = Group(                   # make it
            label="inputs",
            description="Herein are all the manually defined inputs for FANS."
        ).store()                           # the `.store` method saves it in the database

    elif len(groups) == 1:                # otherwise don't make it
        inputs = groups.pop()
    """

    mo.md(rf"""
    Next, we will create a "group". This is purely an organisational tool that AiiDA provides. It may come in handy later to see what nodes belong to the inputs we are creating today.

    We use the `QueryBuilder` to find all groups with label "inputs". If none exist, we create one and provide it a short description.

    ```py
    {_code}
    ```
    """)
    return


@app.cell(hide_code=True)
def group(Group, QueryBuilder, import_button, mo):
    mo.stop(not import_button.value)  # run on click

    groups = QueryBuilder(                # we will use the advanced query method
    ).append(                             # it consists of a series of `.append` methods
        Group, filters={
            Group.fields.label: "inputs"  # here we filter by groups labeled "inputs"
        }
    ).all(                                # conclude with a method to fetch the results
        flat=True
    )

    if len(groups) == 0:                  # if the "inputs" group does not exist...
        inputs = Group(                   # make it
            label="inputs",
            description="Herein are all the manually defined inputs for FANS."
        ).store()                           # the `.store` method saves it in the database

    elif len(groups) == 1:                # otherwise don't make it
        inputs = groups.pop()

    else:
        raise
    return (inputs,)


@app.cell
def _(Path, mo):
    try:
        dataset_path = Path("tutorial_dataset.h5").absolute()
    except:
        mo.stop(True)
    return (dataset_path,)


@app.cell
def _(dataset_path, mo):
    _code = r"""
    microstructurequery = QueryBuilder(
    ).append(
        SinglefileData, filters={
            SinglefileData.fields.label: "microstructure"
        }
    ).all(
        flat=True
    )

    if len(microstructurequery) == 0:
        microstructurefile = SinglefileData(
            Path('""" + str(dataset_path) + r"""'),
            label="microstructure"
        ).store()
    elif len(microstructurequery) == 1:
        microstructurefile = microstructurequery.pop()
    else:
        raise

    inputs.add_nodes(microstructurefile)           # add the node to the "inputs" group
    """

    mo.md(rf"""
    Next, we store the microstructure file in the database. Using a similar strategy as with the group definition, the `QueryBuilder` first searches for existing microstructures. If none are found, we define a new one in the form of a `SinglefileData` node. This built-in AiiDA datatype points to a file via a path. Finally, the microstructure node is included in our "inputs" group.

    ```py
    {_code}
    ```
    """)
    return


@app.cell
def microstructure(
    QueryBuilder,
    SinglefileData,
    dataset_path,
    import_button,
    inputs,
    mo,
):
    mo.stop(not import_button.value)  # run on click

    microstructurequery = QueryBuilder(
    ).append(
        SinglefileData, filters={
            SinglefileData.fields.label: "microstructure"
        }
    ).all(flat=True)

    if len(microstructurequery) == 0:
        microstructurefile = SinglefileData(
            dataset_path,
            label="microstructure"
        ).store()
    elif len(microstructurequery) == 1:
        microstructurefile = microstructurequery.pop()
    else:
        raise

    inputs.add_nodes(microstructurefile)
    return


@app.cell
def _(mo):
    mo.md(r"""
    **Note:** _more nodes..._

    The microstructure file node we just created is saved by AiiDA as a node. Just as before, we can list all the nodes we've created thus far; and it may be helpful to do so every once in a while to ensure everything is proceeding as expected.

    To list all the nodes stored in your profile, run:

    ```
    verdi node list
    ```
    """).callout(kind="info")
    return


@app.cell
def _(mo):
    def_nodes_button = mo.ui.run_button(label="RUN", kind="warn")
    def_nodes_code_switch = mo.ui.switch(label="*show full code...*")

    _code = r"""
    # Microstructure Definition
    Str("/dset_0/image", label="ms_datasetname"),
    Str("/dset_1/image", label="ms_datasetname"),
    Str("/dset_2/image", label="ms_datasetname"),
    List([1.0, 1.0, 1.0], label="ms_L"),

    ...

    # Problem Type and Material Model: Moduli
    Dict({"bulk_modulus": bulk, "shear_modulus": shear}, label="material_properties")
     for bulk, shear in product(
        [[uniform(50, 75), uniform(200, 250)] for _ in range(2)],
        [[uniform(25, 50), uniform(150, 200)] for _ in range(2)]
    )

    ...
    """

    mo.md(rf"""
    Now, we will define the rest of our parameters. This is mostly straightforward, but we treat `ms_datasetname` and `material_properties` a little differently.

    - `ms_datasetname`: Three different datasets are chosen from the sample microstructure file provided.
    - `material_properties`: A mock parameter space study is realised by randomly picking bulk and shear moduli from within a range.

    When it comes time to run our calculations, we will run the "product" of all these parameters.

    ```py
    {_code}
    ```
    {def_nodes_code_switch}
    """)
    return def_nodes_button, def_nodes_code_switch


@app.cell(hide_code=True)
def _(def_nodes_button, def_nodes_code_switch, mo):
    def gatekeep1():
        mo.stop(not def_nodes_button.value and def_nodes_code_switch.value, output=mo.show_code())
        mo.stop(not def_nodes_button.value)
    return (gatekeep1,)


@app.cell
def node_definition(
    ArrayData,
    Dict,
    Float,
    Int,
    List,
    Str,
    array,
    gatekeep1,
    product,
    uniform,
):
    gatekeep1() # Ignore this line.

    nodes = [

    # Microstructure Definition
    Str("/dset_0/image", label="ms_datasetname"),
    Str("/dset_1/image", label="ms_datasetname"),
    Str("/dset_2/image", label="ms_datasetname"),
    List([1.0, 1.0, 1.0], label="ms_L"),

    # Problem Type and Material Model
    Str("mechanical", label="problem_type"),
    Str("LinearElasticIsotropic", label="matmodel")
    ] + [
    Dict({"bulk_modulus": bulk, "shear_modulus": shear}, label="material_properties")
     for bulk, shear in product(
        [[uniform(50, 75), uniform(200, 250)] for _ in range(2)],
        [[uniform(25, 50), uniform(150, 200)] for _ in range(2)]
    )] + [

    # Solver Settings
    Str("cg", label="method"),
    Str("Linfinity", label="error_parameters.measure"),
    Str("absolute", label="error_parameters.type"),
    Float(1e-10, label="error_parameters.tolerance"),
    Int(100, label="n_it"),

    # Macroscale Loading Conditions
    ArrayData({
        "0": array([[0,0,0,0,0,0]])
    }, label="macroscale_loading"),

    # Results Specification
    List(["stress", "strain", "stress_average", "strain_average", 
          "absolute_error", "phase_stress_average", "phase_strain_average", 
          "microstructure", "displacement"], label="results")

    ]
    return (nodes,)


@app.cell
def _(def_nodes_button, mo):
    mo.md(rf"""
    While the cell above defined all the parameters, they still need to be stored in the database. Otherwise, they will be lost when the session ends. AiiDA automatically stores nodes when submitting them to a job, but it is good practice to handle this yourself. Moreover, you get to see your database grow step by step. After clicking the button below, try running `verdi node list` in your terminal to see all the new additions we've made so far, and `verdi node show <id>` for more information about specific nodes.

    It is important to note that this time we did not make any checks through the QueryBuilder to ensure that indentical nodes don't already exist. This means that if you click the button below repeatedly, you *may* cause duplicate nodes to be created. Since these are some the first nodes we're making, it is not so critical, but in practice you would want to first fetch existing nodes you want to reuse before creating the remainder of the nodes you wish to study.

    {def_nodes_button}

    ```py
    for node in nodes:             # iterate over the list of new node
        node.store()               # store each one in the database
        inputs.add_nodes(node)     # assign each one to the "inputs" group
    ```
    """)
    return


@app.cell
def node_storage(def_nodes_button, inputs, mo, nodes):
    mo.stop(not def_nodes_button.value) # Ignore this line.

    for node in nodes:             # iterate over the list of new node
        node.store()               # store each one in the database
        inputs.add_nodes(node)     # assign each one to the "inputs" group
    return


@app.cell
def _(mo):
    mk_params_code_switch = mo.ui.switch(label="*show full code...*")

    _code = r"""
    some_params = [{
        "problem_type": fetch("problem_type", "mechanical"),
        "matmodel": fetch("matmodel", "LinearElasticIsotropic"),
        ...
    }]

    ...

    ms_datasetname_params = [
        {"microstructure":{"file": ms_file,"L": ms_L,
            "datasetname": fetch("ms_datasetname", "/dset_0/image"),}
        },
        {"microstructure":{"file": ms_file,"L": ms_L,
            "datasetname": fetch("ms_datasetname", "/dset_1/image"),}
        },
        {"microstructure":{"file": ms_file,"L": ms_L,
            "datasetname": fetch("ms_datasetname", "/dset_2/image"),}
        }
    ]

    material_properties_params = [
        {"material_properties": mp.pop()}
        for mp in QueryBuilder().append(
            Dict, filters={
                Dict.fields.label: "material_properties"
            }
        ).iterall()
    ]
    """

    mo.md(rf"""
    ### Executing Calculations

    Now that all the input parameters have been specified, it it time to run some calculations. We create lists of dictionaries for each set of paramaters we wish to vary. In our case, `microsctructure` needs a list, as does `material_properties`. Everything else falls into a list of length one. The keys of the dictionaries here are important and are specified by the plugin. More information is available in the documentation, but efforts are being made to synchronise these with the FANS parameter specification.

    Below, some nodes are fetched using a helper function (see [Appendix A](#appendix)) which essentially queries the database for a single node with a particular label and value. You could also use the nodes we created above instead, passing them forward as variables, but here we demonstrate how you might run calculations using a either new or old nodes at once.

    Click the button bellow when you are sure that all the nodes above have been successfully created and stored. Try `verdi node list` to see them all.


    ```py
    {_code}
    ```

    {mk_params_code_switch}
    """)
    return (mk_params_code_switch,)


@app.cell
def _(def_nodes_button, mk_params_code_switch, mo):
    def gatekeep2():
        mo.stop(not def_nodes_button.value and mk_params_code_switch.value, output=mo.show_code())
        mo.stop(not def_nodes_button.value)

    def gatekeep3():
        return mo.show_code() if mk_params_code_switch.value else None
    return


@app.cell
def parameter_definition(
    ArrayData,
    Dict,
    QueryBuilder,
    SinglefileData,
    fetch,
    mk_params_code_switch,
    mo,
):
    some_params = [{
        "problem_type": fetch("problem_type", "mechanical"),
        "matmodel": fetch("matmodel", "LinearElasticIsotropic"),
        "method": fetch("method", "cg"),
        "error_parameters": {
            "measure": fetch("error_parameters.measure", "Linfinity"),
            "type": fetch("error_parameters.type", "absolute"),
            "tolerance": fetch("error_parameters.tolerance", 1e-10)
        },
        "n_it": fetch("n_it", 100),
        "macroscale_loading": QueryBuilder().append(
            ArrayData, filters={
                ArrayData.fields.label: "macroscale_loading"
            }
        ).first(flat=True),
        "results": fetch("results", ["stress", "strain", "stress_average", "strain_average", "absolute_error", "phase_stress_average", "phase_strain_average", "microstructure", "displacement"])
    }]

    ms_file = QueryBuilder().append(
        SinglefileData, filters={
            SinglefileData.fields.label: "microstructure"
        }
    ).first(flat=True)

    ms_L = fetch("ms_L", [1.0, 1.0, 1.0])

    ms_datasetname_params = [
        {"microstructure":{"file": ms_file,"L": ms_L,
            "datasetname": fetch("ms_datasetname", "/dset_0/image"),}
        },
        {"microstructure":{"file": ms_file,"L": ms_L,
            "datasetname": fetch("ms_datasetname", "/dset_1/image"),}
        },
        {"microstructure":{"file": ms_file,"L": ms_L,
            "datasetname": fetch("ms_datasetname", "/dset_2/image"),}
        }
    ]

    material_properties_params = [
        {"material_properties": mp.pop()}
        for mp in QueryBuilder().append(
            Dict, filters={
                Dict.fields.label: "material_properties"
            }
        ).iterall()
    ]

    mo.show_code() if mk_params_code_switch.value else None # Ignore this line.
    return material_properties_params, ms_datasetname_params, some_params


@app.cell
def _(code_settings, mo):
    calculate_button = mo.ui.run_button(label="RUN", kind="warn")

    get_calc_state, set_calc_state = mo.state(False)

    _code = r"""
    FANSCalculation = CalculationFactory("fans")      # get the plugin's process class
    code = {"code": load_code('""" + f"{"<code_label>')}" if code_settings.value is None else code_settings.value["label"] + "')}" : <22}" + """ # get the existing code node

    for sp, dsp, mpp in product(some_params, ms_datasetname_params, material_properties_params):
        all_params = sp | dsp | mpp                   # merge this permutation of params
        run(FANSCalculation, all_params | code)       # finally run the job
    """

    mo.md(rf"""
    Once these lists are defined, we use the `product` function to explore every permutation of their contents. Each permutation is coupled with the code node, defined earlier, and given to the `run` function with the plugin specific `FANSCalculation` process class.

    Much like last time, we aren't checking if these calculations have already been run, so clicking the button below repeatedly will request duplicate calulations to be run and duplicate results will be generated.

    {calculate_button}

    ```py
    {_code}
    ```
    """)
    return calculate_button, get_calc_state, set_calc_state


@app.cell
def calculations(
    CalculationFactory,
    ConfigurationError,
    calculate_button,
    code_settings,
    load_code,
    material_properties_params,
    mo,
    ms_datasetname_params,
    product,
    run,
    set_calc_state,
    some_params,
):
    mo.stop(not calculate_button.value)

    FANSCalculation = CalculationFactory("fans")      # get the plugin's process class
    try:                                              # get the existing code node
        code = {"code": load_code(code_settings.value["label"])}
    except ConfigurationError:
        mo.stop(True, output=mo.md("**Your code failed to load properly!**\n\nPlease submit the 'Define a Code' form in the [AiiDA Setup](aiida-setup) section.").style(text_align="center").callout(kind="danger"))

    for sp, dsp, mpp in mo.status.progress_bar(
        list(product(some_params, ms_datasetname_params, material_properties_params)),
        title="Calculating Jobs...", completion_title="Finished!"
    ):
        all_params = sp | dsp | mpp                   # merge this permutation of params
        run(FANSCalculation, all_params | code)       # finally run the job
    else:
        set_calc_state(True)
    return


@app.cell(hide_code=True)
def _(mo):
    query_button = mo.ui.run_button(label="RUN")

    get_query_state, set_query_state = mo.state(False)

    mo.md(rf"""
    ## Analysing the Results

    Once our calculations are complete, we can make use of the QueryBuilder again to find and analyse the results.

    {query_button}
    """)
    return get_query_state, query_button, set_query_state


@app.cell(hide_code=True)
def _(mo, query_button, set_query_state):
    #! DO NOT DELETE
    # This cell saves the query_button state to allow for user confirmation!
    mo.stop(not query_button.value)  # run on click
    set_query_state(True)
    return


@app.cell(hide_code=True)
def _(get_calc_state, set_calc_state):
    #! DO NOT DELETE
    # This cell triggers the following cell upon confirmation! 
    set_calc_state(get_calc_state())
    return


@app.cell
def _(
    CalcJobNode,
    Int,
    QueryBuilder,
    Str,
    get_calc_state,
    get_query_state,
    mo,
    set_calc_state,
    set_query_state,
):
    confirm = mo.ui.button(label="Are you sure?", on_click=lambda _: set_calc_state(True))
    are_you_sure = mo.md(rf"""
    It seems the jobs were not calculated in this session. If you are sure that they have been completed, you may proceed.

    {confirm}
    """).callout(kind="danger")
    mo.stop(not get_query_state())
    mo.stop(not get_calc_state(), output=are_you_sure)
    set_query_state(False)


    # QUERY:

    calc = QueryBuilder().append(CalcJobNode).first(flat=True)
    # Inputs
    ins = list(calc.inputs._get_keys())
    ins = "<br>".join(ins)
    # Microstructure Dataset Name
    ms_datasetname = calc.inputs.microstructure.datasetname.value
    # Material Properties
    mat_props = {
        "b": (
            calc.inputs.material_properties["bulk_modulus"][0],
            calc.inputs.material_properties["bulk_modulus"][1]
        ),
        "s": (
            calc.inputs.material_properties["shear_modulus"][0],
            calc.inputs.material_properties["shear_modulus"][1]
        )
    }
    # Outputs
    outs = list(calc.outputs._get_keys())
    outs = ", ".join(outs)
    # Stresses and Strains
    log = calc.outputs.retrieved.get_object_content("input.json.log").split("\n")
    stresses = []
    strains = []
    for ln in log:
        if "Effective Stress" in ln:
            stresses.append(list(map(
                lambda n: round(float(n), ndigits=3),
                ln.lstrip("# Effective Stress .. ")
                .replace("(", "").replace(")", "").strip(" ")
                .split(" ")
            )))
        if "Effective Strain" in ln:
            strains.append(list(map(
                lambda n: round(float(n), ndigits=3),
                ln.lstrip("# Effective Strain .. ")
                .replace("(", "").replace(")", "").strip(" ")
                .split(" ")
            )))
    stress_strains = [{"stress": stress, "strain": strain} for stress, strain in zip(stresses, strains)]
    # Filtered Query
    filtered_calcs = \
    QueryBuilder(
    ).append(                              # In the first `.append` we look for nodes
        Str,                               # of the `Str` AiiDA datatype,
        filters={                          # then apply the filters for:
            Int.fields.label: "ms_datasetname", # 
            Int.fields.value: {"==": "/dset_0/image"} #
        },
        tag="ms_datasetname"               # The `tag` is an internal reference.

    ).append(                              # In the second `.append` we look for nodes
        CalcJobNode,                       # of the `CalcJobNode` AiiDA datatype,
        with_incoming="ms_datasetname"     # and specify required incoming nodes with
                                           # the `tag` we defined above.
    ).all(flat=True)


    # DISPLAY:

    _code = r"""
    QueryBuilder(
    ).append(                              # In the first `.append` we look for nodes
        Str,                               # of the `Str` AiiDA datatype,
        filters={                          # then apply the filters for:
            Int.fields.label: "ms_datasetname", # 
            Int.fields.value: {"==": "dset_0"}  #
        },
        tag="ms_datasetname"               # The `tag` is an internal reference.

    ).append(                              # In the second `.append` we look for nodes
        CalcJobNode,                       # of the `CalcJobNode` AiiDA datatype,
        with_incoming="ms_datasetname"     # and specify required incoming nodes with
                                           # the `tag` we defined above.
    ).all(flat=True)
    """

    mo.md(rf"""
    ### Fetch a single calculation...

    We will begin by querying the database for the first `CalcJobNode` present. This is the AiiDA datatype given to nodes that represent the exectution of an individual job.

    ```py
    calc = QueryBuilder().append(CalcJobNode).first(flat=True)
    ```

    From this calculation job node we can gleam some identifying information, such as the type of calculation job (i.e. the process label) or its primary key in the database. Additionally, we can list the available inputs and outputs provided by this kind of job.

    |                    |                                       |
    |--------------------|---------------------------------------|
    | **Process Label:** | {calc.process_label}                  |
    | **Primary Key:**   | {calc.pk}                             |
    | **Inputs:**        | {ins}                             |
    | **Outputs:**       | {outs} |

    ### Identify some input parameters...

    Of course, it would be helpful to know exactly what inputs were used in the calculation of this particular job. The inputs can be accessed via dot notation which provides the respective values as AiiDA datatypes.

    When it comes to the microstructure dataset name, the inputs's value is accessed through the `value` attribute.

    ```py
    calc.inputs.microstructure.datasetname.value
    ```

    | | |
    |-|-|
    | **Microstructure Dataset Name:** | {ms_datasetname} |

    In the case of the material properties, this attribute takes the form an AiiDA `Dict` which has methods just like an ordinary `dict`.

    ```py
    calc.inputs.material_properties.items()
    ```

    | | | |
    |-|-|-|
    |Bulk Modulus: | {mat_props["b"][0]} | {mat_props["b"][1]} |
    |Shear Modulus: | {mat_props["s"][0]} | {mat_props["s"][1]} |

    ### Effective stress and strain...

    To extract the effective stress and strain per loading condition from the output of FANS, we can use the `std_out` it produces. This text is stored in the `retrieved` folder output. We can get its contents and parse it to determine our results.

    ```py
    log = calc.outputs.retrieved.get_object_content("input.json.log")
    for ln in log:
        ...
    ```

    | Loading <br> Condition: | Stress: | Strain:                       |
    |---|-------------------------------|-------------------------------|
    | **1** | {stress_strains[0]["stress"]} | {stress_strains[0]["strain"]} |

    ### Perform a filtered query...

    Aside from manually examining the inputs and outputs of individual calculation jobs, the `QueryBuilder` offers the ability to filter your query based on a variety of criteria. In this instance, we query for all jobs that used the "dset_0" microstructure dataset. This time, we are given back a list of calculation job nodes to do with as we please.

    ```py
    {_code}
    ```

    | | | | | |
    |-|-|-|-|-|
    | **Primary Keys:** | {filtered_calcs[0].pk} | {filtered_calcs[1].pk} | {filtered_calcs[2].pk} | {filtered_calcs[3].pk} |

    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Appendix""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## A. `fetch()`

        This is a helper function to simplify the querying of individual nodes when the label and value are known.
        """
    )
    return


@app.cell
def _(Dict, Float, Int, List, QueryBuilder, Str, mo):
    def fetch(label : str, value):
        """Helper function to return a node whose label and value are known.

        Returns an error if more or less than 1 suitable node is found.
        """
        match value:
            case str():
                datatype = Str
            case int():
                datatype = Int
            case float():
                datatype = Float
            case list():
                datatype = List
            case dict():
                datatype = Dict
            case _:
                raise NotImplementedError

        bone = QueryBuilder().append(
            datatype,
            filters={
                datatype.fields.label: label,
                "attributes.value": value
            } if datatype is not List else {
                datatype.fields.label: label,
                "attributes.list": value
            },
        ).all(flat=True)

        if len(bone) != 1:
            raise RuntimeError

        return bone.pop()

    mo.show_code()
    return (fetch,)


if __name__ == "__main__":
    app.run()
