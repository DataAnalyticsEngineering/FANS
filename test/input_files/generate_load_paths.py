import numpy as np


def generate_strain_path(steps):
    """
    Generate a strain path with smooth transitions based on user-defined steps,
    without repeating strain values at the transition points.

    Parameters:
        steps (list of dict): Each dictionary contains:
            - 'final_strain': A list of length 6 specifying the final strain position to reach
            - 'n_steps': The number of time steps for this loading segment

    Returns:
        strain_path (np.ndarray): An array of shape (total_timesteps, 6) representing the strain path.
    """

    strain_path = []
    current_strain = np.zeros(6)  # Start from zero strain

    for step in steps:
        final_strain = np.array(step["final_strain"])
        n_steps = step["n_steps"]

        # Generate incremental steps excluding the start point
        step_strain = np.linspace(
            current_strain, final_strain, n_steps + 1, endpoint=True
        )[1:]

        # Append the step strain path
        strain_path.append(step_strain)

        # Update current strain to the last value of this step
        current_strain = final_strain

    # Concatenate all steps to form the complete strain path
    strain_path = np.vstack(strain_path)

    return strain_path


def generate_random_walk_path(n_steps, step_config, initial_strain=None):
    """
    Generate a strain path based on a random walk.

    Parameters:
        n_steps (int): The total number of time steps for the random walk.
        step_config (list of dict): A list of dictionaries for each strain component (length 6).
            Each dictionary can contain:
            - 'range': A tuple (min_step, max_step) indicating the range of random step sizes for this component.
            - OR 'fixed_step': A fixed step size for this component.
        initial_strain (list or np.ndarray): The starting strain position. If None, starts from zero.

    Returns:
        random_walk_path (np.ndarray): An array of shape (n_steps, 6) representing the random walk strain path.
    """

    if initial_strain is None:
        current_strain = np.zeros(6)
    else:
        current_strain = np.array(initial_strain)

    random_walk_path = [current_strain]

    for _ in range(n_steps):
        step = np.zeros(6)
        for i, config in enumerate(step_config):
            if "range" in config:
                # Generate a random step within the specified range
                step[i] = np.random.uniform(config["range"][0], config["range"][1])
            elif "fixed_step" in config:
                # Use the fixed step size
                step[i] = config["fixed_step"]
            # Randomly decide whether to add or subtract the step
            step[i] *= np.random.choice([-1, 1])

        # Update the current strain
        current_strain = current_strain + step
        random_walk_path.append(current_strain)

    random_walk_path = np.array(random_walk_path)

    return random_walk_path


def print_formatted_strain_path(strain_path):
    formatted_strain = ",\n".join(
        [f"[{', '.join(f'{v:.4g}' for v in strain)}]" for strain in strain_path]
    )
    print(formatted_strain)


# Example usage for regular strain path
steps = [
    {"final_strain": [0.005, 0, 0, 0, 0, 0], "n_steps": 50},
    {"final_strain": [-0.005, 0, 0, 0, 0, 0], "n_steps": 50},
    {"final_strain": [0.005, 0, 0, 0, 0, 0], "n_steps": 50},
]
strain_path = generate_strain_path(steps)
print("Generated strain path:")
print_formatted_strain_path(strain_path)


# # Example usage for random walk strain path
# step_config = [
#     {'range': (-0.001, 0.001)},  # Component 1 with a range
#     {'range': (-0.001, 0.001)},      # Component 2 with a fixed step size
#     {'range': (-0.001, 0.001)},# Component 3 with a different range
#     {'range': (-0.0001, 0.0001)},# Component 4 with another range
#     {'range': (-0.0001, 0.0001)},      # Component 5 with a fixed step size
#     {'range': (-0.0001, 0.0001)} # Component 6 with another range
# ]

# random_walk_path = generate_random_walk_path(n_steps=100, step_config=step_config)
# print("\nGenerated random walk strain path:")
# print_formatted_strain_path(random_walk_path)
