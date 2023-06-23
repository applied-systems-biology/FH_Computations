# Copyright by Christoph Saffer
# Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
# https://www.leibniz-hki.de/en/applied-systems-biology.html
# HKI-Center for Systems Biology of Infection
# Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
# Adolf-Reichwein-Straße 23, 07745 Jena, Germany
#
# This code is licensed under BSD 2-Clause
# See the LICENSE file provided with this code for the full license.


import numpy as np
import plotly.graph_objs as go
import pandas as pd
import itertools


def angle_between_vecs(vector1, vector2):
    """
    Calculates angle between two vectors v1 and v2 in degrees.
    Args:
        vector1 (numpy.ndarray): first input vector
        vector2 (numpy.ndarray): second input vector
    Returns:
        theta (float): angle between two vectors in degrees
    """

    # Normalize the vectors
    normalized_vector1 = vector1 / np.linalg.norm(vector1)
    normalized_vector2 = vector2 / np.linalg.norm(vector2)

    # Compute the angle in radians
    angle_radians = np.arccos(np.clip(np.dot(normalized_vector1, normalized_vector2), -1.0, 1.0))

    # Convert the angle from radians to degrees
    angle_degrees = np.degrees(angle_radians)

    return angle_degrees


def generate_spheres(number_of_spheres, radius, angles, distribution=True):
    """
    Generates a specified number of spheres in 3D space, where each sphere touches the previous one.
    Args:
        number_of_spheres (int): number of spheres to generate
        radius (float): radius of each sphere
        angles (list of float): list of angles defining the direction of each new sphere
        distribution (bool): if True, the angle of the new sphere will be randomly chosen from a uniform distribution
                              between -angle and angle. If False, the angle of the new sphere will be exactly the value
                              provided in the angles list
    Returns:
        centers (list of numpy.ndarray): list of 3D coordinates of the center of each sphere
    """
    assert len(angles) == number_of_spheres - 1, "Angles array should have N-1 elements"

    # Initialize the validity of the chain
    valid_chain_found = False

    while not valid_chain_found:
        # Start with one sphere at origin plus radius in z-direction
        centers = [np.array([0.0, 0.0, radius])]
        # Initial direction vector for the first sphere
        direction = np.array([radius * 2, 0.0, 0.0])

        # Generate the chain of spheres
        for i in range(number_of_spheres - 1):
            valid_sphere_found = False
            attempts = 0

            # Keep trying until a valid sphere is found or max attempts are reached
            while not valid_sphere_found and attempts < 100:
                attempts += 1

                # Generate a random vector and make it perpendicular to the current direction
                random_vector = np.random.randn(3)
                perpendicular_vector = random_vector - direction * np.dot(random_vector, direction) / np.linalg.norm(
                    direction) ** 2
                # Make the vector of length equal to the diameter of a sphere
                perpendicular_vector = perpendicular_vector / np.linalg.norm(perpendicular_vector) * radius * 2

                # Select the angle based on the distribution flag
                angle = np.random.uniform(-angles[i], angles[i]) if distribution else angles[i]
                factor_b = np.tan(np.deg2rad(angle))
                factor_a = 1.0 if abs(angle) < 90 else -1.0

                # Compute the new direction vector
                new_direction = factor_a * direction + factor_b * perpendicular_vector
                new_direction = new_direction / np.linalg.norm(new_direction) * radius * 2

                # Compute the new center
                new_center = centers[-1] + new_direction

                # Check if the new sphere is valid (does not overlap and is not below z = radius)
                if new_center[2] >= radius and not any(
                        np.linalg.norm(new_center - center) < 2 * radius for center in centers):
                    valid_sphere_found = True
                    # Calculate distance to previous sphere and angle between the direction vectors
                    distance_to_previous_sphere = np.linalg.norm(new_center - centers[-1])
                    angle_between = angle_between_vecs(direction, new_direction)

                    direction = new_direction
                    centers.append(new_center)

                    # Commented out print statement for debugging purpose
                    # print(f"Sphere {i + 2}: Distance to previous sphere = {distance_to_previous_sphere}, Angle with previous direction = {angle_between} vs actual angle = {angle}")

        # A valid chain is found when the number of centers equals the number of spheres
        if len(centers) == number_of_spheres:
            valid_chain_found = True

    return centers


def draw_sphere(radius, center, color='blue', dens=10, opacity=1.0, hemisphere=False):
    """
        Approximates sphere from dens * dens dots to ensure size of sphere
    """
    # Define the range of theta and phi
    fac = 0.5 if hemisphere else 1.0
    theta = np.linspace(0, 2. * np.pi, dens * 2)
    phi = np.linspace(0, fac * np.pi, dens)
    theta, phi = np.meshgrid(theta, phi)

    # Convert to Cartesian coordinates and add the center point
    x = radius * np.sin(phi) * np.cos(theta) + center[0]
    y = radius * np.sin(phi) * np.sin(theta) + center[1]
    z = radius * np.cos(phi) + center[2]

    # Create a surface trace
    surface = go.Surface(
        x=x, y=y, z=z,
        showscale=False,
        surfacecolor=np.ones_like(z) * radius,
        colorbar=dict(thickness=20, tickvals=[radius], ticktext=[color]),
        colorscale=[[0, color], [1, color]],
        opacity=opacity
    )

    # Package the trace into a figure and return it
    return surface


def draw_rotated_half_cylinder(radius, height, center, color='blue', opacity=1.0):
    # Define the range of theta and z
    theta = np.linspace(0, np.pi, 100)
    z = np.linspace(-height/2, height/2, 50)
    theta, z = np.meshgrid(theta, z)

    # Convert to Cartesian coordinates
    x = z
    y = radius*np.cos(theta)
    z = radius*np.sin(theta)

    # Create a surface trace
    surface = go.Surface(
        x=x+center[0], y=y+center[1], z=z+center[2], # Add the center here
        showscale=False,
        surfacecolor=np.ones_like(z) * radius,
        colorscale=[[0, color], [1, color]],
        opacity=opacity
    )

    angles = np.linspace(0, np.pi, 100)
    radii = np.linspace(0, radius, 50)

    # Create a grid of angles and radii
    grid_angles, grid_radii = np.meshgrid(angles, radii)

    # Convert to Cartesian coordinates
    y = grid_radii*np.cos(grid_angles)
    z = grid_radii*np.sin(grid_angles)
    x = np.zeros_like(y) - height/2

    # Create a surface trace for the first end
    end1 = go.Surface(
        x=x+center[0], y=y+center[1], z=z+center[2], # Add the center here
        showscale=False,
        surfacecolor=np.ones_like(z) * radius,
        colorscale=[[0, color], [1, color]],
        opacity=opacity
    )

    # For the second end, change x to height/2
    x = np.zeros_like(y) + height/2

    end2 = go.Surface(
        x=x+center[0], y=y+center[1], z=z+center[2], # Add the center here
        showscale=False,
        surfacecolor=np.ones_like(z) * radius,
        colorscale=[[0, color], [1, color]],
        opacity=opacity
    )

    return [surface, end1, end2]


def generate_fhr1(radius, angleX, angleZ, distr=False):
    """
    Generates a set of points in a 3D space that follows a certain pattern with rotations applied.
    Args:
        radius (float): The base distance unit for positioning the points.
        angleX (float): The angle of rotation around the x-axis (in degrees).
        angleZ (float): The angle of rotation around the y-axis (in degrees).
        distr (bool): If True, the angles of rotation are randomly chosen from a uniform distribution between -angle and 0.
                      If False, the angles of rotation are exactly the values provided.
    Returns:
        positions (numpy.ndarray): Array of the 3D coordinates of each point.
        real_angle (float): The actual angle between vectors formed by first point and seventh point before and after the rotation.
    """

    def rotate_around_x_axis(point, delta):
        """
        Rotates a given point around the x-axis by a specified angle.
        Args:
            point (numpy.ndarray): The 3D coordinates of the point to rotate.
            delta (float): The angle of rotation (in degrees).
        Returns:
            numpy.ndarray: The 3D coordinates of the rotated point.
        """
        rad_delta = np.deg2rad(delta)
        y = point[1] * np.cos(rad_delta) - point[2] * np.sin(rad_delta)
        z = point[1] * np.sin(rad_delta) + point[2] * np.cos(rad_delta)
        return np.array([point[0], y, z])

    def rotate_around_y_axis(point, delta):
        """
        Rotates a given point around the y-axis by a specified angle.
        Args:
            point (numpy.ndarray): The 3D coordinates of the point to rotate.
            delta (float): The angle of rotation (in degrees).
        Returns:
            numpy.ndarray: The 3D coordinates of the rotated point.
        """
        rad_delta = np.deg2rad(delta)
        x = point[0] * np.cos(rad_delta) - point[2] * np.sin(rad_delta)
        z = point[0] * np.sin(rad_delta) + point[2] * np.cos(rad_delta)
        return np.array([x, point[1], z])

    # Initialize array for storing the positions
    positions = np.zeros((10, 3))

    # Manually set some initial positions
    positions[0, 0] = radius * 2
    positions[1, 0] = radius * 4
    positions[2, 0] = radius * 6
    positions[3, 0] = -radius * 2
    positions[4, 0] = -radius * 4
    positions[5, 0] = -radius * 6
    positions[6, 0] = radius
    positions[7, 0] = radius
    positions[8, 0] = -radius
    positions[9, 0] = -radius
    # sin(pi/3) = 60 degree turned binding units
    positions[6, 2] = radius * 2 * np.sin(np.pi/3)
    positions[7, 2] = radius * 4 * np.sin(np.pi/3)
    positions[8, 2] = radius * 2 * np.sin(np.pi/3)
    positions[9, 2] = radius * 4 * np.sin(np.pi/3)

    # Define initial vector between position 6 and position 0
    vec1 = positions[6] - positions[0]

    # Determine the angles of rotation based on the 'distr' flag
    delta_x = np.random.uniform(-angleX, angleX) if distr else angleX
    delta_z = np.random.uniform(-angleZ, 0) if distr else -angleZ

    # Rotate positions around the y-axis
    for i in range(len(positions)):
        if i > 2:
            positions[i][0] -= 2 * radius
            positions[i] = rotate_around_y_axis(positions[i], delta_z)
            positions[i][0] += 2 * radius

    # Rotate positions around the x-axis
    for i in range(len(positions)):
        if i > 2:
            positions[i] = rotate_around_x_axis(positions[i], delta_x)

    # Define the vector between position 6 and position 0 after rotation
    vec2 = positions[6] - positions[0]

    # Compute the actual angle between the vectors before and after rotation
    real_angle = angle_between_vecs(vec1, vec2)

    # Shift all points in the z-direction by 'radius' units
    for point in positions:
        point[2] += radius

    return positions, real_angle


def create_grid(radius, n_spheres, dens=2):
    """
    Creates a grid of points in 3D space.
    Args:
        radius (float): The distance between the spheres.
        n_spheres (int): The number of spheres.
        dens (float): The density of the grid, i.e., the spacing between the grid points is radius/dens.
    Returns:
        grid (numpy.ndarray): A numpy array of the 3D coordinates of each grid point.
    """

    # Initialize an empty list to store the grid points
    grid = []

    # Calculate the range for the grid points
    grid_range = np.arange(-n_spheres * radius * 2, n_spheres * radius * 2, dens * radius)

    # Generate grid points only in the upper half (z >= 0) of the 3D space
    for x, y, z in itertools.product(grid_range, grid_range, grid_range):
        if z >= 0:
            grid.append(np.array([x, y, z]))

    return np.array(grid)


def find_nearest_grid_point(position, grid):
    """
    Finds the nearest grid point to a given position in the grid.
    Args:
        position (numpy.ndarray): The 3D coordinates of the position.
        grid (numpy.ndarray): A numpy array of the 3D coordinates of each grid point.
    Returns:
        int: The index of the grid point that is closest to the given position.
    """

    # Calculate the Euclidean distances from the given position to each grid point
    distances = np.linalg.norm(grid - position, axis=1)

    # Return the index of the grid point with the smallest distance to the given position
    return np.argmin(distances)


def monte_carlo_simulation(n_spheres, radius, delta, n_runs):
    """
    Performs a Monte Carlo simulation to compute the frequency of spheres ending up at different locations in a 3D grid.
    Args:
        n_spheres (int): The number of spheres.
        radius (float): The radius of the spheres.
        delta (numpy.ndarray): The array of angle deviations.
        n_runs (int): The number of Monte Carlo runs.
    Returns:
        grid (numpy.ndarray): A numpy array of the 3D coordinates of each grid point.
        grid_counter (numpy.ndarray): A numpy array storing the count of spheres at each grid point.
    """

    # Create a 3D grid of points
    grid = create_grid(radius, n_spheres)

    # Initialize a counter for each grid point
    grid_counter = np.zeros(len(grid))

    for _ in range(n_runs):
        # Generate sphere positions
        positions = generate_spheres(n_spheres, radius, delta)

        # Find the nearest grid point for the last 4 sphere positions and increment the counter
        for pos in positions[-4:]:
            nearest_idx = find_nearest_grid_point(pos, grid)
            grid_counter[nearest_idx] += 1

    return grid, grid_counter


def run_factorh_simulations(n_runs, n_spheres, radius, deltas):
    """
    Runs a Monte Carlo simulation for different angle deviations and collates the results into a DataFrame.
    Args:
        n_runs (int): The number of Monte Carlo runs.
        n_spheres (int): The number of spheres.
        radius (float): The radius of the spheres.
        deltas (list): The list of angle deviations.
    Returns:
        pandas.DataFrame: A DataFrame containing the x, y, and z coordinates of each grid point,
                          the count of spheres at each point, and the corresponding angle deviation.
    """

    # Initialize an empty list to store the results
    result_data = []

    for d in deltas:
        print("delta: ", d)

        # Initialize an array of angles with d and set the first two elements to 0
        angles = np.ones(n_spheres - 1) * d
        angles[0], angles[1] = 0.0, 0.0

        # Run the Monte Carlo simulation
        grid, grid_counter = monte_carlo_simulation(n_spheres, radius, angles, n_runs)

        # Append the grid point coordinates, the count, and the delta to the results
        for grid_point, count in zip(grid, grid_counter):
            result_data.append([grid_point[0], grid_point[1], grid_point[2], count, d])

    # Convert the results to a DataFrame and remove all gridpoints with counter zero
    return pd.DataFrame(result_data, columns=["x", "y", "z", "counter", "delta"]).query("counter != 0.0")


def monte_carlo_simulation_fhr1(n_spheres, radius, deltai, deltaj, n_runs):
    """
    Performs a Monte Carlo simulation for 'Fhr1' function to compute the frequency of spheres ending up at different
    locations in a 3D grid.
    Args:
        n_spheres (int): The number of spheres.
        radius (float): The radius of the spheres.
        deltai, deltaj (float): The range of angle deviations.
        n_runs (int): The number of Monte Carlo runs.
    Returns:
        grid (numpy.ndarray): A numpy array of the 3D coordinates of each grid point.
        grid_counter (numpy.ndarray): A numpy array storing the count of spheres at each grid point.
    """
    # Create a 3D grid of points with a specified density
    grid = create_grid(radius, n_spheres, dens=0.33333)

    # Initialize a counter for each grid point
    grid_counter = np.zeros(len(grid))

    # Run the simulation for the specified number of times
    for _ in range(n_runs):
        # Generate positions for the spheres with a specified angle deviation
        positions, a = generate_fhr1(radius, deltai, deltaj, distr=True)

        # For the last four sphere positions, find the nearest grid point and increment its count
        for pos in positions[-4:]:
            nearest_idx = find_nearest_grid_point(pos, grid)
            grid_counter[nearest_idx] += 1

    return grid, grid_counter


def run_fhr1_simulations(n_runs, n_spheres, radius, deltas, bend_both_dimensions=False):
    """
    Runs a Monte Carlo simulation for different angle deviations and collates the results into a DataFrame for 'Fhr1' function.
    Args:
        n_runs (int): The number of Monte Carlo runs.
        n_spheres (int): The number of spheres.
        radius (float): The radius of the spheres.
        deltas (list): The list of angle deviations.
        bend_both_dimensions (bool): Whether to bend both dimensions.
    Returns:
        pandas.DataFrame: A DataFrame containing the x, y, and z coordinates of each grid point,
                          the count of spheres at each point, and the corresponding angle deviation.
    """
    # Initialize an empty list to store the results
    result_data = []

    for i in deltas:
        print("delta: ", i)
        j = i if bend_both_dimensions else 0

        # Generate the positions of the spheres and calculate the actual angle deviation
        vecc, _ = generate_fhr1(radius, i, j, distr=False)
        vecc1 = vecc[7] - vecc[6]
        vecc2 = np.array([0., 0., 2*radius*np.sin(np.pi/3)])
        actual_delta = int(np.round(angle_between_vecs(vecc1, vecc2), 0))

        # Run the Monte Carlo simulation
        grid, grid_counter = monte_carlo_simulation_fhr1(n_spheres, radius, i, j, n_runs)

        # Append the grid point coordinates, the count, and the delta to the results
        for grid_point, count in zip(grid, grid_counter):
            result_data.append([grid_point[0], grid_point[1], grid_point[2], count, actual_delta])

    # Convert the results to a DataFrame and remove all gridpoints with counter zero
    return pd.DataFrame(result_data, columns=["x", "y", "z", "counter", "delta"]).query("counter != 0.0")
