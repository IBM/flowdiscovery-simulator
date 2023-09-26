import os

import h5py
import numpy as np
from modules.clogging import squared_link_radius_clogging
from modules.links_set_modification import links_set_modification


def deposition_geometry(case_data,
                        inlet_link_location,
                        link_length,
                        link_radius,
                        number_of_inlet_face_adjacent_links,
                        out_folder,
                        pressure_gradients_prev,
                        print_progress,
                        time_step_index):
    """
    Calculates geometry evolution of a capillary network due to deposition using
    phenomenological correlations.

    Parameters
    ----------
    case_data : .json file
        set of parameters of pore-scale processes.
    link_length : np.array
        array containing link lengths of 'n' capillaries [voxels].
    link_radius : np.array
        array containing link radius of 'n' capillaries [voxels].
    out_folder : str
        directory containing output files.
    pressure_gradients_prev : np.array
        array containing previous iterations pressure gradients of 'n' capillaries [Pa].
    time_step_index : int
        current time step of the simulation.

    Returns
    -------
    deposition_results_n.h5 : h5 files
        dataset containing distribution of deposition parameters.
    link_radius_variation_deposition : np.array
        dataset value as a numpy array.
    """

    # Import parameters from case_data
    case_sample = case_data['simulation']['sample']
    voxel_size = case_sample['parameters']['voxel_size']

    # Import sample data
    sample_data = case_data['simulation']['sample']['parameters']
    minimum_link_radius = sample_data['minimum_link_radius']

    if case_data['setup']['processes']['deposition']:
        # Import deposition data
        deposition_data = case_data['simulation']['processes']['deposition']

        # Obtain number of links
        number_of_links = link_length.size

        # Deposition based on Erosion and Deposition Law
        if deposition_data['model']['name'] == 'jager2017':
            # Import time parameters from case_m.json
            deposition_reaction_time = deposition_data['parameters']['deposition_reaction_time']
            deposition_coefficient = deposition_data['parameters']['deposition_coefficient']
            relative_concentration = deposition_data['parameters']['relative_concentration']
            experiment_condition_data = case_data['simulation']['flow']['experiment']
            driving_force = experiment_condition_data['boundary_condition']['value']

            # Wall shear stress
            wall_shear_stress = abs(pressure_gradients_prev) * link_radius / (2 * link_length)

            # Calculate modified link radius due to deposition
            modified_link_radius_deposition = \
                (((link_radius * voxel_size - 2 * wall_shear_stress / driving_force)
                    * np.exp(driving_force * deposition_coefficient * relative_concentration
                             * deposition_reaction_time / 2))
                    - (2 * wall_shear_stress / driving_force)) / voxel_size

            # Calculate erosion shear stress
            deposition_shear_stress = abs(pressure_gradients_prev) \
                * modified_link_radius_deposition / link_radius

            # Threshold deposition
            if deposition_data['model']['deposition_method'] == 'geometry_flow_dependent':
                # Define number of clogged links
                deposition_number_clogged_links = 0
                deposition_number_clogged_inlet_links = 0

                if all(wall_shear_stress < deposition_shear_stress):
                    # Print deposition information
                    if print_progress:
                        print("Uniform deposition within the capillary network")

                    # Calculate condition for deposition
                    deposition_onset = np.full(number_of_links, 1)

                    # Calculate deposition rate
                    deposition_rate = relative_concentration * deposition_coefficient \
                        * (deposition_shear_stress - wall_shear_stress)

                    # Calculate modified link radius due to deposition
                    modified_link_radius_deposition = \
                        (((link_radius * voxel_size - 2 * wall_shear_stress / driving_force)
                         * np.exp(driving_force * deposition_coefficient * relative_concentration
                                  * deposition_reaction_time / 2))
                         - (2 * wall_shear_stress / driving_force)) / voxel_size

                elif any(wall_shear_stress > deposition_shear_stress):
                    # Locate links with wall_shear_stress < deposition_shear_stress
                    deposition_link_location = \
                        np.where(wall_shear_stress < deposition_shear_stress)[0]

                    # Obtain number of links under deposition
                    number_deposition_link_location = deposition_link_location.size

                    # Print deposition information
                    if print_progress:
                        print(f'Deposition within {number_deposition_link_location} of \
{number_of_links} capillaries')

                    # Calculate deposition onset
                    deposition_onset = np.array(wall_shear_stress < deposition_shear_stress,
                                                dtype=int)
                    # Calculate deposition rate
                    deposition_rate = (deposition_shear_stress - wall_shear_stress) \
                        * relative_concentration * deposition_coefficient
                    deposition_rate[np.logical_not(deposition_onset)] = 0

                    # Calculate modified link radius due to deposition
                    modified_link_radius_deposition = \
                        (((link_radius * voxel_size - 2 * wall_shear_stress / driving_force)
                         * np.exp(driving_force * deposition_coefficient * relative_concentration
                                  * deposition_reaction_time / 2))
                         - (2 * wall_shear_stress / driving_force)) / voxel_size

                    modified_link_radius_deposition[np.logical_not(deposition_onset)] = \
                        link_radius[np.logical_not(deposition_onset)]

                # MODIFIED LINK RADIUS COMPUTATION based on minimum radius
                if any(modified_link_radius_deposition < minimum_link_radius):
                    # LOCATE LINKS WITH X < minimum_link_radius**2
                    clogged_link_location = \
                        np.where(modified_link_radius_deposition < minimum_link_radius)[0]

                    # Calulate deposition_number_clogged_links
                    deposition_number_clogged_links = clogged_link_location.size

                    # Calculate modified_link_radius_deposition
                    modified_link_squared_radius_deposition = \
                        squared_link_radius_clogging(minimum_link_radius,
                                                     modified_link_radius_deposition**2,
                                                     clogged_link_location)

                    # Calculate modified modified_link_radius_deposition
                    modified_link_radius_deposition = \
                        np.sqrt(modified_link_squared_radius_deposition)

        # Calculate updated modified_inlet_link_radius due to combined processes
        modified_inlet_link_radius_deposition = \
            links_set_modification(modified_link_radius_deposition,
                                   inlet_link_location)

        # INLET LINKS CLOGGING
        # Calculate number_clogged_inlet_links
        deposition_number_clogged_inlet_links = \
            np.count_nonzero(modified_inlet_link_radius_deposition <= minimum_link_radius)

        # Print clogging evolution
        if print_progress:
            # Print number of clogged links
            print(f'{deposition_number_clogged_links} of {number_of_links} TOTAL clogged links')

            # Print number_clogged_inlet_links
            print(f'{deposition_number_clogged_inlet_links} of \
{number_of_inlet_face_adjacent_links} clogged INLET links')

        # Define deposition_results folder
        deposition_results_folder = os.path.join(out_folder, 'deposition_results')

        # Create deposition_results folder
        try:
            os.mkdir(deposition_results_folder)
        except FileExistsError:
            pass

        # Define deposition_results_path
        deposition_results_path = os.path.join(deposition_results_folder,
                                               f'deposition_results_{time_step_index:06}.h5')

        # Create HDF5 files for deposition results and initialize them
        deposition_results = h5py.File(deposition_results_path, 'w')

        # SAVE RESULTS
        # Save array containing number of clogged links
        deposition_results.create_dataset('deposition_number_clogged_links',
                                          data=deposition_number_clogged_links)

        # Save array containing number of clogged inlet links
        deposition_results.create_dataset('deposition_number_clogged_inlet_links',
                                          data=deposition_number_clogged_inlet_links)

        # Save array containing deposition rate to deposition_results_00000n.h5
        deposition_results.create_dataset('deposition_rate', data=deposition_rate)

        # Save array containing condition for deposition to deposition_results_00000n.h5
        deposition_results.create_dataset('deposition_onset', data=deposition_onset)

        # Save array containing deposition shear stress to deposition_results_00000n.h5
        deposition_results.create_dataset('deposition_shear_stress', data=deposition_shear_stress)

        # Close deposition_results_00000n.h5
        deposition_results.close()

    else:
        # Modified link radius due to deposition
        modified_link_radius_deposition = link_radius

    # Link radius variation
    link_radius_variation_deposition = modified_link_radius_deposition - link_radius

    return link_radius_variation_deposition
