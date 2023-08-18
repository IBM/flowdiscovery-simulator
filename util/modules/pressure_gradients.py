import os

import h5py
import numpy as np


def pressure_gradients_simple(axis,
                              centerlines_data,
                              out_folder,
                              static_results_filename='static_results.h5'):
    """
    Calculates pressure gradients, using pressures at source and target capilaries as inputs.

    Parameters
    ----------
    out_folder : str
        directory containing output files.
    axis : either int (0, 1, 2) or str (x, y, z)
        flow axis.
    voxel : float
        size of the voxel in [m].
    filename : int
        name of the input JSON file.
    Returns
    -------
    delta_P : np.array
        dataset value as a np.array.
    zero_delta_P : np.array
        dataset value as a np.array.
    P : np.array
        dataset value as a np.array.
    """

    # Import fields from static_results.h5 at first time step
    with h5py.File(os.path.join(out_folder, static_results_filename), 'r') as file:
        # Load pressure input file
        P = np.squeeze(np.array(file[f'pressures_{axis}'], dtype=np.double))

    # Extract link geometry arrays from centerlines file
    edges = centerlines_data['graph']['edges']
    sources = np.array([int(edge['source']) for edge in edges])
    targets = np.array([int(edge['target']) for edge in edges])

    # Obtain number of links
    number_of_links = sources.size

    # Calculate pressure gradients
    pressure_gradients = np.array(
        [P[sources[link]] - P[targets[link]] for link in range(number_of_links)])
    pressure_gradients[(pressure_gradients < np.finfo(float).eps)] = np.finfo(float).eps

    # Identify capillaries with delta_P == np.finfo(float).eps for fluid-solid interaction analysis
    zero_pressure_gradients = np.array(pressure_gradients == np.finfo(float).eps, dtype=float)

    return (pressure_gradients, zero_pressure_gradients, P)
