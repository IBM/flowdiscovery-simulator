import os
import time

import h5py
import numpy as np
from modules.clogging import squared_link_radius_clogging
from modules.links_set_modification import links_set_modification


def precipitation_geometry(case_data,
                           default_extract_axis,
                           initial_Q,
                           inlet_link_location,
                           link_length,
                           link_radius,
                           number_of_inlet_face_adjacent_links,
                           out_folder,
                           prev_static_results_filename,
                           print_progress,
                           time_step_index):
    """
    Calculates geometry evolution of a capillary network due to precipitation using
    phenomenological correlations.

    Parameters
    ----------
    case_data : .json file
        set of parameters of pore-scale processes.
    link_radius : np.array
        array containing link radius of 'n' links [m].
    out_folder : str
        directory containing output files.
    time_step_index : int
        current time step of the simulation.

    Returns
    -------
    precipitation_results_n.h5 : .h5 files
        dataset containing distribution of precipitation parameters.
    link_radius_variation_precipitation : np.array
        dataset value as a numpy array [voxels].
    """

    if case_data['setup']['processes']['precipitation']:
        # Start process time
        start_precipitation_time = time.process_time()

        # Define precipitation_results folder
        precipitation_results_folder = os.path.join(out_folder, 'precipitation_results')

        # Create precipitation_results folder
        try:
            os.mkdir(precipitation_results_folder)
        except FileExistsError:
            pass

        # Create HDF5 files for precipitation results and initialize them
        precipitation_results = h5py.File(os.path.join(precipitation_results_folder,
                                                       f'precipitation_results_\
{time_step_index:06}.h5'), 'w')

        # Define geometry_results and static_results folder names
        static_results_folder = os.path.join(out_folder, 'static_results')

        # Define flow_rate path
        flow_rate_path = os.path.join(static_results_folder, prev_static_results_filename)

        # Import flow experiment parameters from case.json
        case_experiment = case_data['simulation']['flow']['experiment']
        flow_axis = case_experiment['flow_axis']

        # Define case_sample
        case_sample = case_data['simulation']['sample']
        voxel_size = case_sample['parameters']['voxel_size']

        # Import precipitation data
        precipitation_data = case_data['simulation']['processes']['precipitation']

        # Define number of links
        number_of_links = link_radius.size

        if precipitation_data['model']['name'] == 'noiriel2012':
            # Define precipitation_parameters dictionary
            precipitation_parameters = precipitation_data['parameters']

            # Define solid_properties dictionary
            solid_properties = case_data['simulation']['phases']['solid']['properties']

            # Import parameters of the solid phase
            solid_density = solid_properties['solid_density']
            solid_molar_weight = solid_properties['solid_molar_weight']

            # Import fields from static_results_00000n.h5 at previous time step
            with h5py.File(flow_rate_path, 'r') as static_file_data:
                # Load flow speed input file
                Q = np.squeeze(np.array(static_file_data[f'flow_rate_\
{default_extract_axis[flow_axis]}'], dtype=np.double))

            # PRECIPITATION RATE CALCULATION
            # Define precipitation rate method
            precipitation_rate_method = precipitation_data['model']['precipitation_rate_method']

            # Import experiment factor
            experiment_factor = precipitation_parameters['experiment_factor']

            # Set precipitation rate using parameters from the literature
            if precipitation_rate_method == 'constant_precipitation_rate_input':
                # Import precipitation rate
                precipitation_rate_input = precipitation_parameters['precipitation_rate']

                # Calculate precipitation rate for experiments
                precipitation_rate = np.full(number_of_links,
                                             experiment_factor * precipitation_rate_input)

            # Calculate precipitation rate using parameters from the literature
            elif precipitation_rate_method == 'constant_precipitation_rate_parameters':
                # Import precipitation parameters
                m_coefficient = precipitation_parameters['m_coefficient']
                n_coefficient = precipitation_parameters['n_coefficient']
                precipitation_rate_constant = \
                    precipitation_parameters['precipitation_rate_constant']
                saturation_index = precipitation_parameters['saturation_index']

                # Calculate precipitation rate
                precipitation_rate_calculated = precipitation_rate_constant \
                    * (saturation_index**m_coefficient - 1)**n_coefficient

                # Calculate precipitation rate
                precipitation_rate = np.full(number_of_links, experiment_factor
                                             * precipitation_rate_calculated)

            elif any([precipitation_rate_method == 'constant_flow_rate',
                      precipitation_rate_method == 'flow_dependent_precipitation_rate']):
                # Import precipitation parameters
                calcium_concentration_variation = \
                    precipitation_parameters['calcium_concentration_variation']

                # FLOW DEPENDENT (PRECIPITATION RATE AS A FUNCTION OF FLOW RATE) Constant flow rate
                if precipitation_rate_method == 'constant_flow_rate':
                    # Calculate precipitation rate
                    precipitation_rate_calculated = - abs(initial_Q) \
                        * calcium_concentration_variation \
                        / (2 * np.pi * link_radius * link_length * voxel_size**2)

                # FLOW DEPENDENT (PRECIPITATION RATE AS A FUNCTION OF FLOW RATE)
                elif precipitation_rate_method == 'flow_dependent_precipitation_rate':
                    # Calculate precipitation rate
                    precipitation_rate_calculated = - abs(Q) * calcium_concentration_variation \
                        / (2 * np.pi * link_radius * link_length * voxel_size**2)

                # Calculate precipitation rate for experiments
                precipitation_rate = experiment_factor * precipitation_rate_calculated

            elif precipitation_rate_method == 'flow_dependent_precipitation_rate_input':
                # Import precipitation parameters
                calcium_concentration_variation = \
                    precipitation_parameters['calcium_concentration_variation']

                # Import precipitation rate
                precipitation_rate = precipitation_parameters['precipitation_rate']

                # Calculate of scaled_link_radius
                modified_link_radius_precipitation_calculated = - abs(Q) \
                    * calcium_concentration_variation \
                    / (2 * np.pi * (experiment_factor * precipitation_rate)
                        * link_length * voxel_size)

                if (print_progress):
                    print(f'precipitation_rate: {precipitation_rate}')
                    print(f'modified_link_radius_precipitation_calculated: \
{modified_link_radius_precipitation_calculated}')

                modified_link_radius_precipitation = \
                    modified_link_radius_precipitation_calculated / voxel_size

                # Calculate of scaled_link_squared_radius
                modified_link_squared_radius_precipitation = \
                    modified_link_radius_precipitation**2

                # PRINT PROGRESS (Remove after testing one-phase flow precipitation)
                if (print_progress):
                    # Print MAX modified_link_radius due to precipitation
                    print(f'PREV MAX modified_link_radius_precipitation (in voxels) = \
{np.max(modified_link_radius_precipitation)}')
                    print(f'PREV MAX link_radius (in voxels) = \
{np.max(link_radius)}')

            # # REGRESSION
            if any([precipitation_rate_method == 'constant_precipitation_rate_input',
                    precipitation_rate_method == 'constant_precipitation_rate_parameters',
                    precipitation_rate_method == 'constant_flow_rate',
                    precipitation_rate_method == 'flow_dependent_precipitation_rate']):

                # Import time parameters from case_m.json
                precipitation_reaction_time = \
                    precipitation_parameters['precipitation_reaction_time']

                # Calculate precipitation constant
                precipitation_coefficient = 2 * precipitation_rate * solid_molar_weight * \
                    precipitation_reaction_time / (solid_density * voxel_size)

                if precipitation_data['model']['regression_method'] == 'no_regression':
                    # Calculate of scaled_link_squared_radius
                    modified_link_squared_radius_precipitation = \
                        precipitation_coefficient * link_radius \
                        + link_radius**2

                # Evaluatiuon of scaled_link_squared_radius
                elif precipitation_data['model']['regression_method'] == 'quadratic':
                    # Quadratic regression
                    fit_model = quadratic_model(precipitation_coefficient,
                                                link_radius,
                                                number_of_links)

                    # Calculate of scaled_link_squared_radius
                    modified_link_squared_radius_precipitation = \
                        fit_model(link_radius)

                elif precipitation_data['model']['regression_method'] == 'exponential':
                    # Get exponential equation terms
                    fit_model = exponential_model_polyfit(precipitation_coefficient,
                                                          link_radius)

                    # Calculate of scaled_link_squared_radius
                    modified_link_squared_radius_precipitation = \
                        np.exp(fit_model[1]) * np.exp(fit_model[0] * link_radius)

            # Define number of clogged links
            precipitation_number_clogged_links = 0
            precipitation_number_clogged_inlet_links = 0

            # Import sample data
            sample_data = case_data['simulation']['sample']['parameters']
            minimum_link_radius = sample_data['minimum_link_radius']

            # MODIFIED LINK RADIUS COMPUTATION based on minimum radius
            if all(x > minimum_link_radius**2 for x in
                    modified_link_squared_radius_precipitation):

                # Modified link radius due to precipitation
                modified_link_radius_precipitation = \
                        np.sqrt(modified_link_squared_radius_precipitation)

            elif any(x < minimum_link_radius**2 for x in
                     modified_link_squared_radius_precipitation):

                # Locate links with X < minimum_link_radius**2
                clogged_link_location = \
                    np.where(modified_link_squared_radius_precipitation <
                             minimum_link_radius**2)[0]

                # Calulate precipitation_number_clogged_links
                precipitation_number_clogged_links = clogged_link_location.size

                # Calculate modified_link_radius_precipitation
                modified_link_squared_radius_precipitation = \
                    squared_link_radius_clogging(minimum_link_radius,
                                                 modified_link_squared_radius_precipitation,
                                                 clogged_link_location)

                # Modified link radius due to precipitation
                modified_link_radius_precipitation = \
                    np.sqrt(modified_link_squared_radius_precipitation)

        # Calculate updated modified_inlet_link_radius due to combined processes
        modified_inlet_link_radius_precipitation = \
            links_set_modification(modified_link_radius_precipitation,
                                   inlet_link_location)

        # INLET LINKS CLOGGING
        # Calculate number_clogged_inlet_links
        precipitation_number_clogged_inlet_links = \
            np.count_nonzero(modified_inlet_link_radius_precipitation <= minimum_link_radius)

        # Print clogging evolution
        if (print_progress):
            # Print number of clogged links
            print(f'{precipitation_number_clogged_links} of {number_of_links} \
TOTAL clogged links')

            # Print number_clogged_inlet_links
            print(f'{precipitation_number_clogged_inlet_links} of \
{number_of_inlet_face_adjacent_links} clogged INLET links')

        # SAVE RESULTS
        # Save array containing number of clogged links
        precipitation_results.create_dataset('precipitation_number_clogged_links',
                                             data=precipitation_number_clogged_links)

        # Save array containing number of clogged inlet links
        precipitation_results.create_dataset('precipitation_number_clogged_inlet_links',
                                             data=precipitation_number_clogged_inlet_links)

        # Save array containing distribution of clogged links
        precipitation_results.create_dataset('precipitation_rate',
                                             data=precipitation_rate)

        # Save array containing distribution of clogged links
        precipitation_results.create_dataset('maximum_precipitation_rate',
                                             data=np.max(abs(precipitation_rate)))

        # Close precipitation_results_0000n.h5
        precipitation_results.close()

        # PRINT PROGRESS (Remove after testing one-phase flow precipitation)
        if (print_progress):
            # Print MAX modified_link_radius due to precipitation
            print(f'MAX modified_link_radius_precipitation (in voxels) = \
{np.max(modified_link_radius_precipitation)}')

            # Print MIN precipitation rate
            print(f'MAX precipitation_rate = {abs(np.min(precipitation_rate))}')

        # Print process time
        if (print_progress):
            print(f'Geometry modification (precipitation) time: \
{time.process_time() - start_precipitation_time} s')

    else:
        # Modified link radius due to precipitation
        modified_link_radius_precipitation = link_radius

    # Variation of link radius
    link_radius_variation_precipitation = modified_link_radius_precipitation - link_radius

    return link_radius_variation_precipitation


def exponential_model_polyfit(precipitation_coefficient, link_radius):
    """
    Creates a exponential function for clogging and clogging due to precipitation.

    Parameters
    ----------
    precipitation_coefficient : float
        coefficient that combines phasic properties and precipitation parameters.
    link_radius : np.array
        array containing link radius of 'n' links [voxel].

    Returns
    -------
    fit: np.array
        a and b terms of the exponenitial model y = np.exp(a_term) * np.exp(b_term)**x
    """

    # Minimum and maximum values for quadratic regression
    maximum_value = np.max(link_radius)

    # Dataset for regression
    x = np.linspace(precipitation_coefficient + 1, maximum_value, num=100)
    y = x**2 + precipitation_coefficient*x

    # Polyfit the arrays
    fit = np.polyfit(x, np.log(y), 1, w=np.sqrt(y))

    return fit


def quadratic_model(precipitation_coefficient, link_radius, number_of_links):
    """
    Creates a quadratic function for clogging and clogging due to precipitation.

    Parameters
    ----------
    precipitation_coefficient : float
        coefficient that combines phasic properties and precipitation parameters.
    link_radius : np.array
        array containing link radius of 'n' links [voxel].

    Returns
    -------
    fit: np.poly1d
        a and b terms of the quadratic model y = ax^2 + bx + c
    """

    # Minimum and maximum values for quadratic regression
    min_value = np.min(link_radius)
    maximum_value = np.max(link_radius)

    # Dataset for regression
    x = np.linspace(min_value, maximum_value, num=number_of_links)
    y = precipitation_coefficient * x + x**2

    fit_model = np.poly1d(np.polyfit(x, y, 2))

    return fit_model
