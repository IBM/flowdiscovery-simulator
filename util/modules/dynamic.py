import json
import os

import cv2
import h5py
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.mplot3d import Axes3D


def read_network(out_folder, voxel_size, file_name='centerlines.json'):
    """
    Reads the 'centerlines.json' file and returns the network geometry.

    Parameters
    ----------
    out_folder : str
        folder containing '/snapshots/'.
    voxel_size : float
        size of the voxel in [m].
    filename : str
        name of the input JSON file.

    Returns
    -------
    origins : np.array
        (N, 3) array containing source node coordinates for all 'N' capillaries.
    lengths : np.array
        (N, 3) array containing the source-to-target distance for all 'N' capillaries.
    link_diameters : np.array
        (N,) array containing the diameters for all 'N' capillaries.
    node_diameters : np.array
        (N,) array containing the diameters for all 'N' capillaries.
    x : np.array
        (N,) array containing the node coordinates for all 'N' capillaries.
    y : np.array
        (N,) array containing the node coordinates for all 'N' capillaries.
    z : np.array
        (N,) array containing the node coordinates for all 'N' capillaries.
    """

    # Load centerlines input file
    with open(os.path.join(out_folder, file_name), mode='r') as file:
        data = json.load(file)

    # Extract node geometry arrays from JSON
    nodes = sorted(data['graph']['nodes'], key=lambda node: int(node['id']))
    x = np.array([node['metadata']['node_coordinates']['x'] for node in nodes], dtype=np.double)
    y = np.array([node['metadata']['node_coordinates']['y'] for node in nodes], dtype=np.double)
    z = np.array([node['metadata']['node_coordinates']['z'] for node in nodes], dtype=np.double)
    node_squared_radius = np.array([node['metadata']['node_squared_radius'] for node in nodes])

    # Extract link geometry arrays from JSON
    edges = sorted(data['graph']['edges'], key=lambda edge: int(edge['id']))
    sources = np.array([int(edge['source']) for edge in edges])
    targets = np.array([int(edge['target']) for edge in edges])
    link_squared_radius = np.array([edge['metadata']['link_squared_radius'] for edge in edges])
    link_length = np.array([edge['metadata']['link_length'] for edge in edges])

    # Extract capillary coordinates
    origins = np.column_stack([x[sources], y[sources], z[sources]])
    lengths = np.column_stack([x[targets] - x[sources],
                               y[targets] - y[sources],
                               z[targets] - z[sources]])

    # Convert return variable to appropriate units
    origins *= (voxel_size / 1.0e-6)
    lengths *= (voxel_size / 1.0e-6)
    link_diameters = 2.0 * np.sqrt(link_squared_radius) * (voxel_size / 1.0e-6)
    node_diameters = 2.0 * np.sqrt(node_squared_radius) * (voxel_size / 1.0e-6)
    link_lengths = link_length * (voxel_size / 1.0e-6)

    return sources, targets, origins, lengths, link_diameters, node_diameters, link_lengths, x, y, z


def list_output_files(out_folder):
    """
    Lists the names of all output files.

    Parameters
    ----------
    out_folder : str
        folder containing '/snapshots/'.

    Returns
    -------
    output_files : list(str)
        list of output file names.
    """

    # Read all file names
    folder = out_folder + '/snapshots/'
    output_files = sorted([folder + f for f in os.listdir(folder) if f.endswith('.h5')])

    return output_files


def get_dataset(hdf5_file, dataset_name):
    """
    Retrieves a scalar or an array from an HDF5 file.

    Parameters
    ----------
    hdf5_file : h5py.File
        an h5py.File object previously loaded into memory.
    dataset_name : str
        name of the dataset to be retrieved from the h5py.File

    Returns
    -------
    value : scalar or np.array
        dataset value as a np.array or scalar (if value.size == 1).
    """

    # Get dataset value as np.array
    value = hdf5_file[dataset_name][(0)]
    return value


def get_time_step(output_files):
    """
    Retrieves the time step between output files.

    Parameters
    ----------
    output_files : list(str)
        list of output file names.

    Returns
    -------
    time_step : double
        time step between output files.
    """

    # Read first output file as HDF5
    first_output_h5 = h5py.File(output_files[0], 'r')
    time_step = get_dataset(first_output_h5, 'current_time')[0]
    return time_step


def draw_interfaces(parameters):
    """
    Draws the entire network from a snapshot.

    Parameters
    ----------
    parameters : tuple
        various parameters to be unpacked inside the function.
    """

    # Unpack arguments tuple
    filename, origins, lengths, diameters, out_folder, plane, scaling_factor = parameters

    # Read HDF5 files
    hdf5_file = h5py.File(filename, 'r')
    snapshot_number = filename[-9:-3:]

    # Get datasets from h5py.File
    time = get_dataset(hdf5_file, 'current_time')[0]
    fluid_saturation = get_dataset(hdf5_file, 'fluid_saturation')
    fluid_at_source = get_dataset(hdf5_file, 'fluid_at_source')
    interface_offsets = get_dataset(hdf5_file, 'interface_offsets')
    interface_positions = get_dataset(hdf5_file, 'interface_positions')

    # Initialise figure
    fig = plt.figure(dpi=200)

    if plane is None:
        ax = fig.add_subplot(111, projection=Axes3D.name)
        ax.text2D(0.04, 0.95, 't = {0:.6f} s'.format(time), transform=ax.transAxes)
        ax.set_xlabel('X [$\\mu$m]')
        ax.set_ylabel('Y [$\\mu$m]')
        ax.set_zlabel('Z [$\\mu$m]')
        ax.view_init(30, int(snapshot_number))
    else:
        ax = fig.add_subplot(111)
        ax.text(0.04, 0.96, 't = {0:.6f} s'.format(time), transform=ax.transAxes)
        ax.set_aspect(aspect='equal', adjustable='box')
        if plane == 'xy':
            ax.set_xlabel('X [$\\mu$m]')
            ax.set_ylabel('Y [$\\mu$m]')
        elif plane == 'yz':
            ax.set_xlabel('Y [$\\mu$m]')
            ax.set_ylabel('Z [$\\mu$m]')
        else:
            ax.set_xlabel('Z [$\\mu$m]')
            ax.set_ylabel('X [$\\mu$m]')

    # Iterate over capillaries
    for cap_i in np.arange(fluid_at_source.size):
        # Determine line colours based on which fluid is at the source of the capillary
        colours = ['blue', 'darkorange'] if fluid_at_source[cap_i] else ['darkorange', 'blue']

        # Build an array with the interface positions for a given capillary
        index_range = np.arange(interface_offsets[cap_i], interface_offsets[cap_i + 1], dtype=int)
        position = np.array([0.0, 1.0])
        if index_range.size != 0:
            position = np.insert(position, 1, interface_positions[index_range])

        # Get geometrical information for a given capillary
        origin = origins[cap_i]
        length = lengths[cap_i]
        diameter = diameters[cap_i] / scaling_factor

        # Iterate over interfaces of a given capillary and draw them
        for int_i in np.arange(position.size - 1):
            colour = colours[0] if (int_i % 2 == 0) else colours[1]
            interface_span = np.array([position[int_i], position[int_i + 1]])
            P = np.column_stack([origin, origin]) + np.outer(length, interface_span)

            if plane is None:
                ax.plot(P[0], P[1], P[2], color=colour, linewidth=diameter, linestyle='-',
                        solid_capstyle='round', solid_joinstyle='round')
            elif plane == 'xy':
                ax.plot(P[0], P[1], color=colour, linewidth=diameter, linestyle='-',
                        solid_capstyle='round', solid_joinstyle='round')
            elif plane == 'yz':
                ax.plot(P[1], P[2], color=colour, linewidth=diameter, linestyle='-',
                        solid_capstyle='round', solid_joinstyle='round')
            else:
                ax.plot(P[2], P[0], color=colour, linewidth=diameter, linestyle='-',
                        solid_capstyle='round', solid_joinstyle='round')

    inset_ax = inset_axes(ax, width='100%', height='100%', bbox_to_anchor=[1100, 750, 125, 125])
    inset_ax.set_aspect(aspect='equal', adjustable='box')
    patches, _ = inset_ax.pie(fluid_saturation,
                              radius=1.2,
                              startangle=90,
                              colors=['darkorange', 'blue'])
    inset_ax.legend(handles=patches,
                    loc='center',
                    frameon=False,
                    fontsize='x-small',
                    bbox_to_anchor=(0.5, -0.25),
                    labels=['Resident', 'Injected'])

    plt.savefig(out_folder + '/dyn_interfaces_{0}.png'.format(snapshot_number), bbox_inches='tight')
    plt.close()


def create_video(folder, prefix, time_step, speed, clean):
    """
    Creates a movie out of individual images.

    Parameters
    ----------
    folder : str
        folder containing the image files.

    prefix : str
        prefix common to all file names.

    time_step : double
        time step between output files.

    speed : double
        speed multiplier with respect to real time.

    clean : boolean
        remove individual images after creating movie.
    """

    # List individual images
    images = sorted([os.path.join(folder, f) for f in os.listdir(folder)
                     if f.startswith(prefix) and f.endswith('.png')])

    # Define video properties
    fps = speed / time_step
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    height, width, channels = cv2.imread(images[0]).shape
    video = cv2.VideoWriter(os.path.join(folder, prefix + '.mp4'), fourcc, fps, (width, height))

    # Build movie out of individual frames
    [video.write(cv2.imread(image)) for image in images]
    video.release()

    # Clean all image files after saving movie
    if clean:
        [os.remove(image) for image in images]
