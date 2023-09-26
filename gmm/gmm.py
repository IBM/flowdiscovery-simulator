#!/usr/bin/env python3

"""
This program calculates the geometry evolution of porous media under
simultaneous mineral pore-scale processes (i.e., erosion, deposition,
dissolution, and precipitation).
"""

import argparse
import json
import os
import subprocess
import sys
import time

import h5py
import numpy as np
from modules.create_config import create_config_file
from modules.geometry_modification import (calculate_flow_geometry_parameters,
                                           geometry_evolution_single_process,
                                           save_flow_geometry_parameters)
from modules.links_set_modification import links_set_modification

if __name__ == '__main__':

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Geometry Modification Module (GMM) for geometry evolution prediction.',
        epilog="""\
                Returns for 'n' simulation time steps:
                ----------------------------------------
                centerlines_n.json : JSON file
                    updated centerlines.
                static_results_n.h5 : h5 files
                    fluid flow parameters.
                erosion_results_n.h5 : h5 files
                    erosion parameters.
                deposition_results_n.h5 : h5 files
                    deposition parameters.
                dissolution_results_n.h5 : h5 files
                    dissolution parameters.
                precipitation_results_n.h5 : h5 files
                    precipitation parameters.
                flow_results_n.h5 : h5 files
                    flow parameters of solid-fluid interactions.
                geometry_results_n.h5 : h5 files
                    geometry parameters of solid-fluid interactions.""")

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
    parser.add_argument('--binary_image_filename',
                        action='store',
                        metavar='BINARY_IMAGE.raw',
                        type=str,
                        required=False,
                        default='binary_image',
                        help='Name of the input BINARY_IMAGE.raw file.')
    parser.add_argument('--centerlines_filename',
                        action='store',
                        metavar='CENTERLINES.json',
                        type=str,
                        required=False,
                        default='centerlines',
                        help='Name of the input CENTERLINES.json file.')
    parser.add_argument('--config_filename',
                        action='store',
                        metavar='CONFIG.json',
                        type=str,
                        required=False,
                        default='config',
                        help='Name of the input CONFIG.json file, containing flow parameters.')
    parser.add_argument('--config_template_filename',
                        action='store',
                        metavar='CONFIG_TEMPLATE.json',
                        type=str,
                        required=False,
                        default='config_template',
                        help='Name of the input CONFIG_TEMPLATE.json fiLE.')
    parser.add_argument('--geometry_results_filename',
                        action='store',
                        metavar='GEOMETRY_RESULTS.h5',
                        type=str,
                        required=False,
                        default='geometry_results',
                        help='Name of the input GEOMETRY_RESULTS.h5 file.')
    parser.add_argument('--print_flow_info',
                        action='store_true',
                        required=False,
                        default=False,
                        help='Print flow additional information.')
    parser.add_argument('--print_progress',
                        action='store_true',
                        required=False,
                        default=False,
                        help='Print information during each iteration.')
    parser.add_argument('--save_flow_parameters',
                        action='store_true',
                        required=False,
                        default=True,
                        help='Save aditional flow parameters in FLOW_RESULTS_0000n.h5.')
    parser.add_argument('--save_geometry_parameters',
                        action='store_true',
                        required=False,
                        default=True,
                        help='Save aditional geometry parameters in GEOMETRY_RESULTS_0000n.h5.')
    parser.add_argument('--save_modified_centerlines',
                        action='store_true',
                        required=False,
                        default=False,
                        help='Save modified CENTERLINES_00000n.json.')
    parser.add_argument('--static_results_filename',
                        action='store',
                        metavar='STATIC_RESULTS.h5',
                        type=str,
                        required=False,
                        default='static_results',
                        help='Name of the input STATIC_RESULTS.h5 file.')

    arg = parser.parse_args()

    # Default parameters
    default_flow_axis = ['x', 'y', 'z']
    default_processes = ['deposition', 'dissolution', 'erosion', 'precipitation']

    # Start Geometry Modification Module (GMM)
    if (arg.print_progress):
        print('FLOWSIMULATOR::GMM SAYS:')
        print('Loading case parameters and centerlines data:')

    # Load case parameters from case_m.json
    with open(os.path.join(arg.case_folder, f'{arg.case_filename}.json'),
              mode='r') as case_file_data:
        case_data = json.load(case_file_data)

    # Define simulated pore-scale processes
    processes = np.array(list(case_data['setup']['processes'].values()))

    # Warning for missing input
    if not np.any(processes):
        print('Missing pore-scale process as input.')
        sys.exit(11)

    # IMPORT CENTERLINES.JSON
    # Open centerlines file
    centerlines_filename_input = os.path.join(arg.in_folder,
                                              f'{arg.centerlines_filename}.json')
    with open(centerlines_filename_input, mode='r') as centerlines_file_data:
        centerlines_data = json.load(centerlines_file_data)

    # Extract node geometry arrays from JSON of rock sample centerlines
    nodes = sorted(centerlines_data['graph']['nodes'], key=lambda node: int(node['id']))
    number_of_nodes = centerlines_data['graph']['metadata']['number_of_nodes']

    x = np.array([node['metadata']['node_coordinates']['x'] for node in nodes], dtype=float)
    y = np.array([node['metadata']['node_coordinates']['y'] for node in nodes], dtype=float)
    z = np.array([node['metadata']['node_coordinates']['z'] for node in nodes], dtype=float)
    coordinates = [x, y, z]

    # IDENTIFY INLET LINKS
    # Import flow experiment parameters from case.json
    case_experiment = case_data['simulation']['flow']['experiment']
    flow_axis = case_experiment['flow_axis']

    # Extract link geometry arrays from centerlines file
    edges = centerlines_data['graph']['edges']
    sources = np.array([int(edge['source']) for edge in edges])
    targets = np.array([int(edge['target']) for edge in edges])

    # Obtain number of links
    number_of_links = sources.size

    # Print number of links
    if (arg.print_progress):
        print(f'Number of links: {number_of_links}')

    # CREATE CONFIG.JSON for Flow Simulator
    create_config_file(arg.in_folder,
                       arg.config_template_filename,
                       arg.config_filename,
                       arg.centerlines_filename,
                       case_data,
                       centerlines_data)

    # INITIAL CONDITIONS
    # Print Initial Conditions
    if (arg.print_progress):
        print('Imposing initial conditions and setting up phasic properties:')

    # Start process time
    start_initial_conditions = time.process_time()

    # CALCULATE STATIC_RESULTS.H5
    # Run Flow Simulator at initial conditions in in_folder if static_results.h5 does not exist
    if os.path.exists(os.path.join(arg.in_folder,
                                   f'{arg.static_results_filename}.h5')) is False:
        # Run Flow Simulator
        if (arg.print_progress):
            print('Running Flow Simulator at the initial conditions:')

        # Calculate flow fields at the initial conditions
        subprocess.run([os.path.join('bin', f'flow-simulator.x --run_simulation {arg.in_folder}',
                                     f'{arg.config_filename}.json')], shell=True)

    # COPY CENTERLINES.JSON AND CONFIG.JSON FILES INTO OUT_FOLDER
    # Define out_folder for results
    out_folder = os.path.join(arg.in_folder, arg.case_filename)

    # Create out_folder for results
    try:
        os.mkdir(out_folder)
    except FileExistsError:
        pass

    # Define centerlines_filename_output
    centerlines_filename_output = os.path.join(out_folder, f'{arg.centerlines_filename}.json')

    # Import time_steps from case_data
    time_steps = case_data['setup']['time_steps']

    # Start geometry evolution computation
    time_step_index = 0

    # Print current time step
    if (arg.print_progress):
        print(f'Time step {time_step_index} of {time_steps}')

    # Copy centerlines.json and config.json to out_folder
    subprocess.run(os.path.join(f'cp {arg.in_folder}',
                                f'{arg.centerlines_filename}.json {arg.in_folder}',
                                f'{arg.config_filename}.json {out_folder}'),
                   shell=True)

    if arg.save_modified_centerlines is True:
        # Define geometry_results folder
        modified_centerlines_folder = os.path.join(out_folder, 'modified_centerlines')

        # Create geometry_results folder
        try:
            os.mkdir(modified_centerlines_folder)
        except FileExistsError:
            pass

        # Copy centerlines.json with 00000n format to out_folder
        subprocess.run(os.path.join(
            f'cp {arg.in_folder}',
            f'{arg.centerlines_filename}.json {modified_centerlines_folder}',
            f'{arg.centerlines_filename}_{time_step_index:06}.json'),
                        shell=True)

    # Define static_results_folder
    static_results_folder = os.path.join(out_folder, arg.static_results_filename)

    # Define geometry_results_folder
    geometry_results_folder = os.path.join(out_folder, arg.geometry_results_filename)

    # create static_results_folder
    try:
        os.mkdir(static_results_folder)
    except FileExistsError:
        pass

    # Copy static_results to static_results_folder
    subprocess.run(os.path.join(f'cp {arg.in_folder}',
                                f'{arg.static_results_filename}.h5 {static_results_folder}',
                                f'{arg.static_results_filename}_{time_step_index:06}.h5'),
                   shell=True)

    # Load config.json in out_folder
    config_filename_output = os.path.join(out_folder, f'{arg.config_filename}.json')
    with open(config_filename_output, mode='r') as config_file_data:
        config_data = json.load(config_file_data)

    # Update dict folder of config.json file
    config_data['setup']['folder'] = out_folder

    # Update config.json file with out_folder dictionary
    with open(config_filename_output, mode='w') as config_file_data:
        json.dump(config_data, config_file_data, sort_keys=False, indent=2)

    # Extract link geometry arrays from centerlines file
    edges = centerlines_data['graph']['edges']

    # Obtain link radius
    link_squared_radius = np.array([edge['metadata']['link_squared_radius'] for edge in edges])
    initial_link_radius = np.sqrt(link_squared_radius)

    # Obtain link length and voxel size
    link_length = np.array([edge['metadata']['link_length'] for edge in edges])
    voxel_size = case_data['simulation']['sample']['parameters']['voxel_size']

    # Calculate initial void space volume
    initial_void_space_volume = np.pi * initial_link_radius**2 * link_length * voxel_size**3

    # Define initial void space volume
    previous_void_space_volume = initial_void_space_volume

    # Current static_results_filename
    current_static_results = f'{arg.static_results_filename}_{time_step_index:06}.h5'

    # Define static results path
    results_path = os.path.join(out_folder,
                                arg.static_results_filename,
                                current_static_results)

    # Define initial pressures and flow rate array names
    initial_P_array = f'pressures_{default_flow_axis[flow_axis]}'
    initial_Q_array = f'flow_rate_{default_flow_axis[flow_axis]}'

    # Extract data
    with h5py.File(results_path, 'r') as initial_static_file:
        initial_P = np.squeeze(np.array(initial_static_file[initial_P_array],
                                        dtype=np.double))

        # Load flow speed input file
        initial_Q = np.squeeze(np.array(initial_static_file[initial_Q_array],
                                        dtype=np.double))

    # Create array for source_pressures
    source_pressures = initial_P[sources]

    # Define inlet links identification
    inlet_links_identification = case_data['simulation']['sample']['inlet_links_identification']

    # Absolute pressure method
    if inlet_links_identification['pressure_threshold_method'] == 'absolute_pressure':

        # Calculate pressure grandient
        pressure_difference = max(initial_P) - min(initial_P)

        # Calculate minimum inlet nodes pressure
        minimum_inlet_nodes_pressure = case_experiment['absolute_pressure'] + pressure_difference \
            * voxel_size / 2

        # Compute total inlet links
        total_inlet_links = (source_pressures >= minimum_inlet_nodes_pressure)

    # Maximum pressure method
    elif inlet_links_identification['pressure_threshold_method'] == 'maximum_pressure':
        # Calculate minimum inlet nodes pressure
        minimum_inlet_nodes_pressure = min(initial_P)

    # Inlet nodes method
    elif inlet_links_identification['pressure_threshold_method'] == 'inlet_nodes':
        # Obtain inlet nodes location
        inlet_nodes_location = (np.where(coordinates[flow_axis] == 0)[0])

        # Calculate number of inlet nodes
        number_of_inlet_nodes = inlet_nodes_location.size

        # Sort source pressures (descending)
        sorted_source_pressures = sorted(source_pressures, reverse=True)

        # Calculate minimum inlet nodes pressure
        minimum_inlet_nodes_pressure = sorted_source_pressures[number_of_inlet_nodes]

    # Compute total inlet links
    total_inlet_links = (source_pressures > minimum_inlet_nodes_pressure)

    # Define inlet links location
    inlet_link_location = np.where(total_inlet_links != 0)[0]

    # Define current total_inlet_links
    number_of_inlet_face_adjacent_links = np.count_nonzero(total_inlet_links)

    # Obtain initial_inlet_link_radius
    initial_inlet_link_radius = links_set_modification(initial_link_radius,
                                                       inlet_link_location)

    # Print inlet links number
    if (arg.print_progress):
        print(f'{number_of_inlet_face_adjacent_links} inlet links of {number_of_links} total links \
for capillary network clogging.')

    # Import sample data
    sample_data = case_data['simulation']['sample']['parameters']
    minimum_link_radius = sample_data['minimum_link_radius']

    # Calculate flow and geometry parameters
    (accumulated_volume_simulation,
     accumulated_volume_time_step,
     aspect_ratio,
     maximum_flow_rate,
     pore_volume,
     porosity,
     porosity_time_step,
     pressure_gradients,
     reactive_area,
     reynolds_number,
     total_accumulated_volume_simulation,
     total_accumulated_volume_time_step,
     void_space_volume,
     void_space_volume_ratio,
     wall_shear_stress) = calculate_flow_geometry_parameters(arg.binary_image_filename,
                                                             case_data,
                                                             centerlines_data,
                                                             default_flow_axis,
                                                             arg.geometry_results_filename,
                                                             geometry_results_folder,
                                                             arg.in_folder,
                                                             initial_link_radius,
                                                             initial_void_space_volume,
                                                             link_length,
                                                             initial_link_radius,
                                                             out_folder,
                                                             previous_void_space_volume,
                                                             arg.print_progress,
                                                             arg.static_results_filename,
                                                             time_step_index)

    # SAVE INITIAL CONDITIONS
    # Create FLOW_RESULTS_000000.H5 and GEOMETRY_RESULTS_000000.H5
    save_flow_geometry_parameters(accumulated_volume_simulation,
                                  accumulated_volume_time_step,
                                  aspect_ratio,
                                  initial_inlet_link_radius,
                                  link_length,
                                  initial_link_radius,
                                  maximum_flow_rate,
                                  out_folder,
                                  pore_volume,
                                  porosity,
                                  porosity_time_step,
                                  pressure_gradients,
                                  arg.print_progress,
                                  reactive_area,
                                  reynolds_number,
                                  arg.save_flow_parameters,
                                  arg.save_geometry_parameters,
                                  time_step_index,
                                  total_accumulated_volume_simulation,
                                  total_accumulated_volume_time_step,
                                  void_space_volume,
                                  void_space_volume_ratio,
                                  wall_shear_stress)

    # Print process time
    if (arg.print_progress):
        print(f'Time initial_conditions: {time.process_time() - start_initial_conditions} s')

    # Start process time
    start_transport_reaction = time.process_time()

    # Print Geometry Evolution
    if (arg.print_progress):
        print('Calculating GEOMETRY EVOLUTION')

    # CALCULATE GEOMETRY EVOLUTION DUE TO SINGLE PROCESSES
    if np.count_nonzero(processes) != 1:
        # Print pore-scale process information
        print('Simultaneous pore-scale processes')
        sys.exit(11)

    # Compute geometry evolution and save centerlines and h5 files
    geometry_evolution_single_process(arg.binary_image_filename,
                                      case_data,
                                      arg.centerlines_filename,
                                      centerlines_filename_output,
                                      minimum_link_radius,
                                      arg.config_filename,
                                      config_filename_output,
                                      default_flow_axis,
                                      arg.geometry_results_filename,
                                      geometry_results_folder,
                                      arg.in_folder,
                                      initial_link_radius,
                                      initial_Q,
                                      initial_void_space_volume,
                                      inlet_link_location,
                                      link_length,
                                      number_of_inlet_face_adjacent_links,
                                      out_folder,
                                      previous_void_space_volume,
                                      arg.print_flow_info,
                                      arg.print_progress,
                                      arg.save_flow_parameters,
                                      arg.save_geometry_parameters,
                                      arg.save_modified_centerlines,
                                      arg.static_results_filename,
                                      start_transport_reaction,
                                      time_step_index,
                                      total_inlet_links)
