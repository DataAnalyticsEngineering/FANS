# aiida-fans-tutorial
Learn how to use aiida-fans in this marimo powered tutorial.

## Usage

Assuming you have FANS, python 3.13, venv, and pip installed on a linux system, the recommended way to use this tutorial is by creating a virtual environment in this directory with the following command:

```
python -m venv .venv
```
Then activate this environment like so:

```
source .venv/bin/activate
```

You can ensure the the environment was succesfully activated with `which python` and ultimately deactivate the environment with `deactivate` when you're finished.

You may need to install/upgrade pip now with your virtual environment activated. Run the following command:

```
python -m pip install --upgrade pip
```

Once pip is up to date, run the following command to install the tutorial's dependencies:

```
python -m pip install -r requirements.txt
```

Now you are ready to launch the notebook and begin the tutorial. Run the following command and access the marimo notebook at the port provided:

```
marimo run tutorial.py
```

## Alternative Usage

### 1. Conda

> [!WARNING]  
> This method is a work-in-progress!

### 2. Pixi

> [!WARNING]  
> This method is a work-in-progress!

You can use pixi to install everything you need as defined by the pyproject.toml file. It should bundle python, FANS, AiiDA, aiida-fans, and marimo all into a virtual environment located in a .pixi directory. You can proceed to directly begin the tutorial with:

```
marimo run tutorial.py
```

Activating the environment may look something like this:

```
pixi shell --manifest-path ~/FANS/tutorial/pyproject.toml
```
