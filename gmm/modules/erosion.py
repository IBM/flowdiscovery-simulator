import os

import h5py
import numpy as np


def erosion_geometry(case_data,
                     link_length,
                     link_radius,
                     out_folder,
                     pressure_gradients_prev,
                     print_progress,
                     time_step_index):
    """
    Calculates geometry evolution of a capillary network due to erosion using
    phenomenological correlations.

    Parameters
    ----------
    case_data : .json file
        set of parameters of pore-scale processes.
    link_length : np.array
        array containing link lengths of 'n' capillaries [m].
    link_radius : np.array
        array containing link radius of 'n' capillaries [m].
    out_folder : str
        directory containing output files.
    pressure_gradients_prev : np.array
        array containing previous iteration pressure gradients of 'n' capillaries [Pa].
    time_step_index : int
        current time step of the simulation.
    voxel_size : float
        size of the voxel in [m].

    Returns
    -------
    erosion_results_n.h5 : h5 files
        dataset containing distribution of erosion parameters.
    link_radius_variation_erosion : np.array
        dataset value as a numpy array.
    """

    # Import parameters from case_data
    case_sample = case_data['simulation']['sample']
    voxel_size = case_sample['parameters']['voxel_size']

    # If erosion is considered
    if case_data['setup']['processes']['erosion']:
        # Import erosion data
        erosion_data = case_data['simulation']['processes']['erosion']

        # Data from centerlines_00000n.json
        number_of_links = link_length.size

        # Calculate wall shear stress
        wall_shear_stress = abs(pressure_gradients_prev) * link_radius / (2 * link_length)

        # Erosion correlation of Jager et al. (2017)
        if erosion_data['model']['name'] == 'jager2017':
            # Define erosion_parameters dictionary
            erosion_parameters = erosion_data['parameters']

            # Import time parameters from case_m.json
            erosion_reaction_time = erosion_parameters['erosion_reaction_time']

            # Import solid phase parameters
            solid_properties = case_data['simulation']['phases']['solid']['properties']
            solid_density = solid_properties['solid_density']

            # Import erosion parameters
            erosion_coefficient = erosion_parameters['erosion_coefficient']

            # Calculate erosion time scale
            erosion_time_scale = 2 * solid_density * link_length * voxel_size \
                / (erosion_coefficient * abs(pressure_gradients_prev))

            # Calculate erosion shear stress
            if erosion_data['model']['erosion_method'] == 'geometry_flow_dependent':
                modified_link_radius_erosion = \
                    link_radius * (1 + (wall_shear_stress / abs(pressure_gradients_prev))
                                   * (np.exp(erosion_reaction_time / erosion_time_scale) - 1))

            # Calculate erosion shear stress
            erosion_shear_stress = abs(pressure_gradients_prev) * modified_link_radius_erosion \
                / link_radius

            if all(wall_shear_stress > erosion_shear_stress):
                # Print erosion information
                if print_progress:
                    print('Uniform erosion within the capillary network')

                # Condition for erosion
                erosion_onset = np.full(number_of_links, 1)

                # Calculate erosion rate
                erosion_rate = -erosion_coefficient * (wall_shear_stress - erosion_shear_stress)

                # Calculate modified link radius due to erosion for no_theshold method
                if erosion_data['model']['erosion_method'] == 'geometry_flow_dependent':
                    # Calculate modified link radius due to erosion
                    modified_link_radius_erosion = (link_radius) * \
                        np.exp(erosion_reaction_time / erosion_time_scale)

            elif any(wall_shear_stress < erosion_shear_stress):
                # Locate links with wall_shear_stress > erosion_shear_stress
                erosion_link_location = \
                    np.where(wall_shear_stress > erosion_shear_stress)[0]

                # Obtain number of links under erosion
                number_erosion_link_location = erosion_link_location.size

                # Print erosion information
                if print_progress:
                    print(f'Erosion within {number_erosion_link_location} of \
{number_of_links} capillaries')

                # Calculate erosion onset
                erosion_onset = np.array(wall_shear_stress > erosion_shear_stress, dtype=int)

                # Calculate erosion rate
                erosion_rate = -erosion_coefficient * (wall_shear_stress - erosion_shear_stress)
                erosion_rate[np.logical_not(erosion_onset)] = 0

                # Calculate modified link radius due to erosion for no_theshold method
                if erosion_data['model']['erosion_method'] == 'geometry_flow_dependent':
                    modified_link_radius_erosion = \
                        link_radius * (1 + (wall_shear_stress / abs(pressure_gradients_prev))
                                       * (np.exp(erosion_reaction_time / erosion_time_scale) - 1))
                    modified_link_radius_erosion[np.logical_not(erosion_onset)] = \
                        link_radius[np.logical_not(erosion_onset)]

        # Define erosion_results folder
        erosion_results_folder = os.path.join(out_folder, 'erosion_results')

        # Create erosion_results folder
        try:
            os.mkdir(erosion_results_folder)
        except FileExistsError:
            pass

        # Create HDF5 files for erosion results and initialize them
        erosion_results = h5py.File(os.path.join(erosion_results_folder,
                                    f'erosion_results_{time_step_index:06}.h5'), 'w')

        # Save array containing erosion rate to erosion_results_00000n.h5
        erosion_results.create_dataset('erosion_rate', data=erosion_rate)

        # Save array containing onset of erosion to erosion_results_00000n.h5
        erosion_results.create_dataset('erosion_onset', data=erosion_onset)

        # Save array containing erosion shear stress to erosion_results_00000n.h5
        erosion_results.create_dataset('erosion_shear_stress', data=erosion_shear_stress)

        # Save array containing erosion time scale to erosion_results_00000n.h5
        erosion_results.create_dataset('erosion_time_scale', data=erosion_time_scale)

        # Close erosion_results_0000n.h5
        erosion_results.close()

    else:
        # Modified link radius due to erosion
        modified_link_radius_erosion = link_radius

    # Link radius variation
    link_radius_variation_erosion = modified_link_radius_erosion - link_radius

    return link_radius_variation_erosion
