# Copyright by Christoph Saffer
# Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
# https://www.leibniz-hki.de/en/applied-systems-biology.html
# HKI-Center for Systems Biology of Infection
# Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
# Adolf-Reichwein-Straße 23, 07745 Jena, Germany
#
# This code is licensed under BSD 2-Clause
# See the LICENSE file provided with this code for the full license.


import plotly.graph_objs as go
import numpy as np
import pandas as pd
from utils import generate_spheres, generate_fhr1


def draw_sphere(radius, center, color='blue', dens=10, opacity=1.0, hemisphere=False):
    """
    Create a plotly surface object that represents a sphere or a hemisphere.

    Args:
        radius (float): The radius of the sphere.
        center (list or np.array): The 3D coordinates of the center of the sphere.
        color (str): The color of the sphere.
        dens (int): The density of the surface grid, higher values provide a smoother appearance.
        opacity (float): The opacity of the sphere, 1 is fully opaque, 0 is fully transparent.
        hemisphere (bool): If True, only the top hemisphere is drawn.

    Returns:
        go.Surface: A plotly surface object that represents the sphere.
    """
    # Factor to determine if the entire sphere or just the top hemisphere should be drawn
    fac = 0.5 if hemisphere else 1.0

    # Define the range of theta (angle from the x-axis) and phi (angle from the z-axis)
    theta = np.linspace(0, 2. * np.pi, 2 * dens)
    phi = np.linspace(0, fac * np.pi, dens)

    # Create a grid of points on the surface of the sphere
    theta, phi = np.meshgrid(theta, phi)

    # Convert the spherical coordinates to Cartesian coordinates
    x = radius * np.sin(phi) * np.cos(theta) + center[0]
    y = radius * np.sin(phi) * np.sin(theta) + center[1]
    z = radius * np.cos(phi) + center[2]

    # Create a surface trace for the sphere
    surface = go.Surface(
        x=x, y=y, z=z,
        showscale=False,
        surfacecolor=np.ones_like(z) * radius,
        colorbar=dict(thickness=20, tickvals=[radius], ticktext=[color]),
        colorscale=[[0, color], [1, color]],
        opacity=opacity
    )

    return surface


def draw_rotated_half_cylinder(radius, height, center, color='blue', opacity=1.0):
    """
    Create a list of plotly surface objects that represent a half-cylinder.

    Args:
        radius (float): The radius of the half-cylinder.
        height (float): The height of the half-cylinder.
        center (list or np.array): The 3D coordinates of the center of the half-cylinder.
        color (str): The color of the half-cylinder.
        opacity (float): The opacity of the half-cylinder, 1 is fully opaque, 0 is fully transparent.

    Returns:
        list of go.Surface: A list of plotly surface objects that represent the half-cylinder.
    """

    # Define the range of theta (angle from the x-axis) and z (height)
    theta = np.linspace(0, np.pi, 100)
    z = np.linspace(-height / 2, height / 2, 50)
    # Create a grid of points on the surface of the half-cylinder
    theta, z = np.meshgrid(theta, z)

    # Convert the cylindrical coordinates to Cartesian coordinates
    x = z
    y = radius * np.cos(theta)
    z = radius * np.sin(theta)

    # Create a surface trace for the half-cylinder
    surface = go.Surface(
        x=x + center[0], y=y + center[1], z=z + center[2],
        showscale=False,
        surfacecolor=np.ones_like(z) * radius,
        colorscale=[[0, color], [1, color]],
        opacity=opacity
    )

    # Define the range of radii and angles for the ends of the half-cylinder
    angles = np.linspace(0, np.pi, 100)
    radii = np.linspace(0, radius, 50)

    # Create a grid of angles and radii for the ends
    grid_angles, grid_radii = np.meshgrid(angles, radii)

    # Convert to Cartesian coordinates for the ends
    y = grid_radii * np.cos(grid_angles)
    z = grid_radii * np.sin(grid_angles)
    x = np.zeros_like(y) - height / 2

    # Create a surface trace for the first end of the half-cylinder
    end1 = go.Surface(
        x=x + center[0], y=y + center[1], z=z + center[2],
        showscale=False,
        surfacecolor=np.ones_like(z) * radius,
        colorscale=[[0, color], [1, color]],
        opacity=opacity
    )

    # For the second end, change x to height/2
    x = np.zeros_like(y) + height / 2

    # Create a surface trace for the second end of the half-cylinder
    end2 = go.Surface(
        x=x + center[0], y=y + center[1], z=z + center[2],
        showscale=False,
        surfacecolor=np.ones_like(z) * radius,
        colorscale=[[0, color], [1, color]],
        opacity=opacity
    )

    return [surface, end1, end2]


def generate_plot_layout(delta="", radius=1.5, n_spheres=20):
    """
    Generate the layout for a 3D plot with specified parameters.

    Args:
        delta (str, optional): Maximum bending angle in degrees, default is an empty string.
        radius (float, optional): Radius of the spheres in the plot, default is 1.5.
        n_spheres (int, optional): Number of spheres in the plot, default is 20.

    Returns:
        go.Layout: A plotly layout object for a 3D plot.
    """
    title = "<i>Maximum bending angle: 𝛅 = " + str(delta) + "°</i>" if delta != "" else ""
    layout = go.Layout(font=dict(size=15),
        scene=dict(
            xaxis=dict(nticks=8, range=[-n_spheres * radius * 2 , n_spheres * radius * 2], title="x-axis (nm)"),
            yaxis=dict(nticks=8, range=[-n_spheres * radius * 2 , n_spheres * radius * 2],  title="y-axis (nm)"),
            zaxis=dict(nticks=4, tickvals=[20, 40, 60], range=[0, n_spheres * radius * 2], title="z-axis (nm)"),
            aspectratio=dict(x=1, y=1, z=0.5),
            camera=dict(
                up=dict(x=0, y=0, z=1),  # this is the 'up' direction for the camera
                center=dict(x=0, y=0, z=-0.1),  # this will move the center of the plot
                eye=dict(x=0.82, y=1.27, z=0.47)  # this is the position of the camera 'eye'
            )
        ),
        width=1300,
        height=1000,
        title=dict(
            text=title,
            y=0.83,
            font=dict(
                size=25,  # font size
            )
        )
    )

    return layout


def draw_factorh(centers, radius=1.5):
    """
    Creates a list of spheres for Factor H visualization. We distinguish between surface units and binding units.

    Args:
        centers (array-like): An array of center coordinates for each sphere.
        radius (float, optional): The radius of the spheres. Defaults to 1.5.

    Returns:
        list: A list of plotly Surface objects representing the spheres.
    """
    traces = []
    for i, pos in enumerate(centers):
        if i < 3:
            sphere = draw_sphere(radius, pos, 'blue', dens=10)
        elif i < 16:
            sphere = draw_sphere(radius, pos, 'grey', dens=10)
        else:
            sphere = draw_sphere(radius, pos, 'yellow', dens=10)

        traces.append(sphere)
    return traces


def draw_fhr1(centers, radius=1.5):
    """
    Creates a list of spheres for FHR1 visualization. We distinguish between surface units and binding units.

    Args:
        centers (array-like): An array of center coordinates for each sphere.
        radius (float, optional): The radius of the spheres. Defaults to 1.5.

    Returns:
        list: A list of plotly Surface objects representing the spheres.
    """
    traces = []
    for i, pos in enumerate(centers):
        if i < 6:
            sphere = draw_sphere(radius, pos, 'blue', dens=20)
        else:
            sphere = draw_sphere(radius, pos, 'yellow', dens=20)

        traces.append(sphere)
    return traces


def draw_probability_cloud(Dfff, scaling_f=2.3, radius=1.5):
    # Function to draw a 3D probability cloud of binding units' positions
    # Initialize an empty list to store 3D scatter plot objects
    traces = []

    # Loop over the unique counts in the 'counter' column of the dataframe
    for count in np.unique(Dfff["counter"]):
        # Filter the dataframe to keep only rows with the current count
        Dfn = Dfff.query('counter == @count').copy()

        # Create a 3D scatter plot object for the filtered dataframe
        # Marker size is proportional to radius and scaling factor
        # Opacity is proportional to the square root of the ratio of the count to the total sum of counts
        # This effectively highlights the regions where the count is high (higher opacity)
        traces.append(go.Scatter3d(
            x=Dfn['x'], y=Dfn['y'], z=Dfn['z'],
            mode='markers',
            marker=dict(
                size=3 * scaling_f * radius,
                color='yellow',
                opacity=np.sqrt(count / Dfff["counter"].sum()),
                line=dict(color='black', width=0)  # Set width to 0 to remove the border
            ),
            showlegend=False  # Don't show a legend for these traces
        ))

    # Generate a dummy surface plot for the colorbar
    dummy_z = np.linspace(0, 1, 10).reshape(10, 1)
    colorscale = [[0, 'rgba(255, 255, 0, 0)'], [1, 'rgba(255, 255, 0, 1)']]

    # The dummy surface plot is fully transparent (opacity=0)
    dummy_trace = go.Surface(
        z=dummy_z,
        colorscale=colorscale,
        showscale=True,
        opacity=0,
        colorbar=dict(
            title='Probability of<br>binding units\'<br>position:',  # Colorbar title
            thickness=25,  # Thickness of the colorbar
            len=0.65,  # Length of the colorbar as a fraction of plot height
        )
    )

    # Add the dummy trace to the list of traces
    traces.append(dummy_trace)

    # Return the list of traces
    return traces


def create_figure_factorh(delta, radius, n_spheres, df=[], show_proability_cloud=False):
    # Generate layout for plot
    layout = generate_plot_layout(delta, radius, n_spheres)
    fig = go.Figure(layout=layout)

    # Generate random factorh configuration
    angles = np.ones(19)
    angles[0], angles[1] = 0.0, 0.0
    centers = generate_spheres(20, radius, angles * delta)

    # Generate spheres
    spheres = draw_factorh(centers, radius)
    for sphere in spheres:
        fig.add_trace(sphere)

    # Generate analytically computed action area and volume
    fig.add_trace(draw_sphere(35 * radius, [4.0 * radius, 0, radius], color='orange', dens=50, opacity=0.1, hemisphere=True))

    if show_proability_cloud:
        df_delta = df.query('counter != 0.0').query('delta == @delta').copy()
        scaling_f = 5
        traces = draw_probability_cloud(df_delta, scaling_f, radius)
        for trace in traces:
            fig.add_trace(trace)

    return fig


def create_figure_fhr(delta, radius, df=[], show_proability_cloud=False, two_bending=False):
    # Generate layout for plot
    layout = generate_plot_layout(delta, radius, n_spheres=20)
    fig = go.Figure(layout=layout)

    # Generate random factorh configuration
    deltay = delta if two_bending else 0
    centers, _ = generate_fhr1(radius, delta, deltay, distr=True)

    # Generate spheres
    spheres = draw_fhr1(centers, radius)
    for sphere in spheres:
        fig.add_trace(sphere)

    # Generate analytically computed action area and volume
    if deltay == 0:
        areas = draw_rotated_half_cylinder(5 * radius, 4 * radius, [0, 0, radius], color='orange', opacity=0.2)
        for f in areas:
            fig.add_trace(f)
    else:
        areas = draw_rotated_half_cylinder(5 * radius, 8 * radius, [2 * radius, 0, radius], color='orange',
                                              opacity=0.2)
        for f in areas:
            fig.add_trace(f)

    if show_proability_cloud:
        df_delta = df.query('counter != 0.0').query('delta == @delta').copy()
        scaling_f = 2.3
        traces = draw_probability_cloud(df_delta, scaling_f, radius)
        for trace in traces:
            fig.add_trace(trace)

    return fig