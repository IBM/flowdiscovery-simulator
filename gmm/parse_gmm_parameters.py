#!/usr/bin/env python3

"""
This program parses the results of GMM simulations.
"""

import argparse
import json
import os

import h5py
import numpy as np

if __name__ == '__main__':

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Parse instance parameters.',
        epilog="""\
                Returns for simulation instance:
                ----------------------------------------
                parse_ti_tf.h5 : h5 files
                    gmm simulation parameters from 'ti' initial time step to 'tf' final
                    time step.""")
    # Required arguments
    parser.add_argument('--case_filename',
                        action='store',
                        metavar='CASE.json',
                        type=str,
                        required=True,
                        help='Name of the input CASE.json file, containing flow parameters.')
    parser.add_argument('--case_folder',
                        action='store',
                        metavar='CASE_FOLDER',
                        type=str,
                        required=True,
                        help='Directory containing CASE.json files.')
    parser.add_argument('in_folder',
                        action='store',
                        metavar='IN_FOLDER',
                        type=str,
                        help='Directory containing CONFIG.json and CENTERLINES.json files.')

    # Optional arguments
    parser.add_argument('--flow_results_foldername',
                        action='store',
                        metavar='FLOW_RESULTS',
                        type=str,
                        required=False,
                        default='flow_results',
                        help='Name of the output FLOW_RESULTS folder.')
    parser.add_argument('--geometry_results_foldername',
                        action='store',
                        metavar='GEOMETRY_RESULTS',
                        type=str,
                        required=False,
                        default='geometry_results',
                        help='Name of the output GEOMETRY_RESULTS folder.')
    parser.add_argument('--parse_parameters',
                        action='store',
                        metavar='parse_parameters',
                        type=str,
                        required=False,
                        default=['maximum_flow_rate',
                                 'maximum_precipitation_rate',
                                 'maximum_reynolds_number',
                                 'permeability',
                                 'pore_volume',
                                 'porosity',
                                 'porosity_time_step',
                                 'precipitation_number_clogged_links',
                                 'precipitation_number_clogged_inlet_links',
                                 'total_accumulated_volume_simulation',
                                 'total_accumulated_volume_time_step',
                                 'void_space_volume_ratio'],
                        help='Default parse parameters. This list can be modified.')
    parser.add_argument('--parse_results_filename',
                        action='store',
                        metavar='PARSE_RESULTS.h5',
                        type=str,
                        required=False,
                        default='parse_results',
                        help='Name of the output PARSE_RESULTS.h5 file.')
    parser.add_argument('--print_progress',
                        action='store_true',
                        default=False,
                        help='Print information during each iteration.')
    parser.add_argument('--remove_local_parameter_files',
                        action='store_true',
                        default=False,
                        help='After parsing, remove local parameter results files.')
    parser.add_argument('--static_results_foldername',
                        action='store',
                        metavar='STATIC_RESULTS',
                        type=str,
                        required=False,
                        default='static_results',
                        help='Name of the output STATIC_RESULTS folder.')
    parser.add_argument('--parsing_iterations',
                        action='store',
                        choices=['count_files',
                                 'time_steps_input'],
                        metavar='parsing_iterations',
                        type=str,
                        required=False,
                        default='count_files',
                        help='Method for calculating number of parsing iterations')
    parser.add_argument('--precipitation_results_foldername',
                        action='store',
                        metavar='PRECIPITATION_RESULTS',
                        type=str,
                        required=False,
                        default='precipitation_results',
                        help='Name of the output PRECIPITATION_RESULTS folder.')

    arg = parser.parse_args()

    # Default parameters
    default_flow_axis = ['x', 'y', 'z']

    # Start Geometry Modification Module (GMM)
    if (arg.print_progress):
        print('FLOWSIMULATOR::GMM SAYS:')

    # Load case parameters from case_m.json
    with open(os.path.join(arg.case_folder, f'{arg.case_filename}.json'),
              mode='r') as case_file_data:
        case_data = json.load(case_file_data)

    # Create out_folder for results
    out_folder = os.path.join(arg.in_folder, arg.case_filename)

    # COUNT FILES
    # Define static results folder name
    static_results_folder = os.path.join(out_folder, arg.static_results_foldername)

    print(f'static_results_folder: {static_results_folder}')

    # DEFINE TIME STEPS FOR PARSING
    # Import time_steps data from case_data
    if arg.parsing_iterations == 'count_files':
        # Count files
        _, _, files = next(os.walk(os.path.join(static_results_folder, '')))

        total_iterations = len(files)
        time_steps = total_iterations - 1

    elif arg.parsing_iterations == 'time_steps_input':
        time_steps = case_data['setup']['time_steps']

        # Define total iterations
        total_iterations = time_steps + 1

    # Define h5 file
    parse_results = h5py.File(os.path.join(out_folder, f'{arg.parse_results_filename}.h5'), 'w')

    # Start iterations
    for parameter in arg.parse_parameters:
        # Print progress
        if (arg.print_progress):
            print(f'Parsing: {parameter}')

        # Create array containing parameter value
        parameter_array = np.zeros(total_iterations, dtype=float)

        # Start iterations
        time_step_index = 0

        # Create dataset for 2d plot horizontal axis
        while time_step_index <= time_steps:
            # Print progress
            if (arg.print_progress) and (time_step_index == time_steps):
                print(f'Iteration {time_step_index} of {time_steps}')

            # Parameters in geometry_results dataset
            if any([parameter == 'maximum_flow_rate',
                    parameter == 'maximum_flow_speed',
                    parameter == 'maximum_reynolds_number']):
                # Define dataset to extract arrays
                dataset = 'flow'

                # Import h5 file with results
                results_filename = os.path.join(out_folder, f'{dataset}_results',
                                                f'{dataset}_results_{time_step_index:06}.h5')

                # Import parameter
                with h5py.File(results_filename, 'r') as flow_file_data:
                    # Load void_space_volume_ratio input file
                    array = np.array(flow_file_data[f'{parameter}'], dtype=np.double)
                    parameter_array[time_step_index] = np.squeeze(array)

            # Parameters in geometry_results dataset
            elif any([parameter == 'pore_volume',
                      parameter == 'porosity',
                      parameter == 'porosity_time_step',
                      parameter == 'total_accumulated_volume_simulation',
                      parameter == 'total_accumulated_volume_time_step',
                      parameter == 'void_space_volume_ratio']):
                # Define dataset to extract arrays
                dataset = 'geometry'

                # Import h5 file with results
                results_filename = os.path.join(out_folder, f'{dataset}_results',
                                                f'{dataset}_results_{time_step_index:06}.h5')

                # Import parameter
                with h5py.File(results_filename, 'r') as geometry_file_data:
                    # Define geometry data array
                    geometry_array = np.array(geometry_file_data[f'{parameter}'], dtype=np.double)

                    # Load void_space_volume_ratio input file
                    parameter_array[time_step_index] = np.squeeze(geometry_array)

            # Parameters in static_results dataset
            elif parameter == 'permeability':
                # Define dataset to extract arrays
                dataset = 'static'

                # Import flow experiment parameters from case.json
                case_experiment = case_data['simulation']['flow']['experiment']
                flow_axis = case_experiment['flow_axis']

                # Import h5 file with results
                results_filename = os.path.join(out_folder, f'{dataset}_results',
                                                f'{dataset}_results_{time_step_index:06}.h5')

                # Import parameter
                with h5py.File(results_filename, 'r') as static_file_data:
                    # Define permeability filename
                    permeability_filename = f'{parameter}_{default_flow_axis[flow_axis]}'

                    # Define array
                    permeability_array = np.array(static_file_data[permeability_filename],
                                                  dtype=np.double)

                    # Load void_space_volume_ratio input file
                    parameter_array[time_step_index] = np.squeeze(permeability_array)

            # Parameters in precipitation_results dataset
            elif any([parameter == 'maximum_precipitation_rate',
                      parameter == 'precipitation_number_clogged_links',
                      parameter == 'precipitation_number_clogged_inlet_links']):

                # If precipitation_results dataset does not exist:
                if os.path.exists(os.path.join(out_folder,
                                               arg.precipitation_results_foldername)) is True:
                    # Define dataset to extract arrays
                    dataset = 'precipitation'

                    # Import h5 file with results
                    results_filename = os.path.join(out_folder, f'{dataset}_results',
                                                    f'{dataset}_results_{time_step_index:06}.h5')

                    # Define parameter value at the initial conditions
                    if time_step_index == 0:
                        parameter_array[time_step_index] = 0

                    # Define parameter value during the iterations
                    else:
                        # Import precipitation_number_clogged_links
                        if parameter == 'precipitation_number_clogged_links':
                            with h5py.File(results_filename, 'r') as file:
                                # Load precipitation_number_clogged_links input file
                                array = np.array(file[f'{parameter}'], dtype=np.double)
                                parameter_array[time_step_index] = np.squeeze(array)

                        elif parameter == 'precipitation_number_clogged_inlet_links':
                            with h5py.File(results_filename, 'r') as file:
                                # Load precipitation_number_clogged_inlet_links input file
                                array = np.array(file[f'{parameter}'], dtype=np.double)
                                parameter_array[time_step_index] = np.squeeze(array)

                        # Import precipitation_maximum_rate
                        elif parameter == 'maximum_precipitation_rate':
                            with h5py.File(results_filename, 'r') as file:
                                # Load precipitation_maximum_rate input file
                                array = np.array(file['maximum_precipitation_rate'],
                                                 dtype=np.double)
                                parameter_array[time_step_index] = np.squeeze(array)

            # Next iteration
            time_step_index += 1

        # Create parameter datasets for flow_results_00000n.h5
        parse_results.create_dataset(f'{parameter}', data=parameter_array)

    # Close h5 file
    parse_results.close()

    # After parsing, remove files corresponding to local parameter results
    if (arg.remove_local_parameter_files):
        # Remove flow_results folder
        os.remove(os.path.join(out_folder, arg.flow_results_foldername))

        # Remove geometry_results folder
        os.remove(os.path.join(out_folder, arg.geometry_results_foldername))

        # Remove static_results folder
        os.remove(os.path.join(out_folder, arg.static_results_foldername))

        # Remove precipitation_results folder
        if os.path.exists(os.path.join(out_folder, arg.precipitation_results_foldername)) is True:
            os.remove(os.path.join(out_folder, arg.static_results_foldername))
