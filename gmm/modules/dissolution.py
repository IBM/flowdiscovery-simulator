import os

import h5py
import numpy as np


def dissolution_geometry(case_data, link_radius, out_folder, time_step_index):
    """
    Calculates geometry evolution of a capillary network due to dissolution using
    phenomenological correlations.

    Parameters
    ----------
    case_data : .json file
        set of parameters of pore-scale processes.
    link_radius : np.array
        array containing link radius of 'n' capillaries [m].
    out_folder : str
        directory containing output files.
    time_step_index : int
        current time step of the simulation.

    Returns
    -------
    dissolution_results_n.h5 : h5 files
        dataset containing distribution of dissolution parameters.
    link_radius_variation_dissolution : np.array
        dataset value as a numpy array.
    """

    # Import parameters from case_data
    case_sample = case_data['simulation']['sample']
    voxel_size = case_sample['parameters']['voxel_size']

    if case_data['setup']['processes']['dissolution']:
        # Dissolution data
        dissolution_data = case_data['simulation']['processes']['dissolution']

        if dissolution_data['model']['name'] == 'molins2021':
            # Define dict name
            dissolution_parameters = dissolution_data['parameters']

            # Import time parameters from case_m.json
            dissolution_reaction_time = dissolution_parameters['dissolution_reaction_time']

            # Import solid phase parameters
            solid_properties = case_data['simulation']['phases']['solid']['properties']
            solid_density = solid_properties['solid_density']
            molar_weight_s = solid_properties['solid_molar_weight']

            # Calculate dissolution rate
            # METHOD: INPUT
            if dissolution_data['model']['dissolution_rate_method'] == 'input':
                # Import dissolution rate from case_data
                dissolution_rate = dissolution_parameters['dissolution_rate']

            # METHOD: CALCULATE
            elif dissolution_data['model']['dissolution_rate_method'] == 'calculation':
                # Import parameters from case_data
                activity_coefficient = dissolution_parameters['activity_coefficient']
                inlet_concentration = dissolution_parameters['inlet_concentration']
                dissolution_rate_constant = dissolution_parameters['dissolution_rate_constant']

                # Calculate dissolution rate
                dissolution_rate = dissolution_rate_constant * activity_coefficient * \
                    inlet_concentration

            # Calculate modified link squared radius due to dissolution
            modified_link_squared_radius_dissolution = link_radius**2 \
                + (2 * link_radius * dissolution_rate * molar_weight_s *
                   dissolution_reaction_time) / (solid_density * voxel_size)

            # Calculate modified link radius due to dissolution
            modified_link_radius_dissolution = np.sqrt(modified_link_squared_radius_dissolution)

        # Define dissolution_results folder
        dissolution_results_folder = os.path.join(out_folder, 'dissolution_results')

        # Create dissolution_results folder
        try:
            os.mkdir(dissolution_results_folder)
        except FileExistsError:
            pass

        # Define dissolution_results_path
        dissolution_results_path = os.path.join(dissolution_results_folder,
                                                f'dissolution_results_{time_step_index:06}.h5')

        # Create HDF5 files for dissolution results and initialize them
        dissolution_results = h5py.File(dissolution_results_path, 'w')

        # Save array containing dissolution rate to dissolution_results_00000n.h5
        dissolution_results.create_dataset('dissolution_rate', data=dissolution_rate)

        # Close dissolution_results_00000n.h5
        dissolution_results.close()

    else:
        # Modified link radius due to dissolution
        modified_link_radius_dissolution = link_radius

    # Variation of link radius
    link_radius_variation_dissolution = modified_link_radius_dissolution - link_radius

    return link_radius_variation_dissolution
