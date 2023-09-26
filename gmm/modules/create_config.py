import json
import os
import subprocess

import numpy as np


def create_config_file(in_folder,
                       config_template_filename,
                       config_filename,
                       centerlines_filename,
                       case_data,
                       centerlines_data):
    """
    Creates a config file for Flow Simulator.

    Parameters
    ----------
    in_folder : str
        folder containing centerlines.json file.
    config_template_filename : str
        name of config_template.json file.
    config_filename : str
        name of the output config.json file.
    centerlines_filename : str
        name of centerlines.json file.
    case_data : .json file
        set of parameters of pore-scale processes.
    centerlines_data : str
        set of centerlines parameters.

    Returns
    -------
    config.json : JSON file
    """

    # Define modules_folder
    modules_folder = os.path.join('gmm', 'modules')

    # Copy centerlines.json and config.json to arg.in_folder/gmm_case_m
    subprocess.run(os.path.join(f'cp {modules_folder}',
                                f'{config_template_filename}.json {in_folder}',
                                f'{config_filename}.json'),
                   shell=True)

    # Define config_filename_input path
    config_filename_input = os.path.join(in_folder, f'{config_filename}.json')

    # Load config.json
    with open(config_filename_input, mode='r') as config_file_data:
        config_data = json.load(config_file_data)

    # Extract node geometry arrays from JSON of rock sample centerlines
    nodes = sorted(centerlines_data['graph']['nodes'], key=lambda node: int(node['id']))
    x = np.array([node['metadata']['node_coordinates']['x'] for node in nodes], dtype=float)
    y = np.array([node['metadata']['node_coordinates']['y'] for node in nodes], dtype=float)
    z = np.array([node['metadata']['node_coordinates']['z'] for node in nodes], dtype=float)

    # UPDATE CONFIG.JSON SETUP PARAMETERS
    # Update folder
    config_data['setup']['folder'] = in_folder

    # Update input_file
    config_data['setup']['input_file'] = f'{centerlines_filename}.json'

    # Update shape
    (config_data['setup']['shape']['x'],
     config_data['setup']['shape']['y'],
     config_data['setup']['shape']['z']) = check_coordinates(x, y, z)

    # Update voxel_size
    case_sample = case_data['simulation']['sample']
    config_data['setup']['voxel_size'] = case_sample['parameters']['voxel_size']

    # UPDATE CONFIG.JSON SIMULATION PARAMETERS
    # Define fluid dictionaries in config.json
    config_fluid = config_data['simulation']['fluid']
    config_properties = config_fluid[0]['properties']

    # Define liquid dictionaries in case.json
    case_liquid = case_data['simulation']['phases']['liquid']

    # Update fluid in config.json
    config_fluid[0]['name'] = case_liquid['liquid_name']
    config_fluid[0]['viscosity_behaviour'] = case_liquid['viscosity_behaviour']

    # Update fluid properties in config.json
    config_properties['dynamic_viscosity'] = case_liquid['properties']['liquid_dynamic_viscosity']

    # Define algorithm dictionaries in config.json
    config_algorithm = config_data['simulation']['algorithm']

    # Define algorithm dictionaries in case.json
    case_algorithm = case_data['simulation']['flow']['algorithm']

    # Update algorithm dictionaries in config.json
    config_algorithm['name'] = case_algorithm['name']
    config_algorithm['model'] = case_algorithm['model']

    # Define experiment dictionaries in config.json
    config_experiment = config_data['simulation']['experiment']
    config_boundary_condition = config_experiment['boundary_condition']

    # Define experiment dictionaries in case.json
    case_experiment = case_data['simulation']['flow']['experiment']
    case_boundary_condition = case_data['simulation']['flow']['experiment']['boundary_condition']

    # Update experiment dictionaries in config.json
    config_experiment['absolute_pressure'] = case_experiment['absolute_pressure']
    config_experiment['boundary_thickness'] = case_experiment['boundary_thickness']
    config_experiment['flow_axis'] = case_experiment['flow_axis']
    config_experiment['temperature'] = case_experiment['temperature']

    # Update boundary_condition dictionaries in config.json
    config_boundary_condition['driving_force'] = case_boundary_condition['driving_force']
    config_boundary_condition['value'] = case_boundary_condition['value']

    # Save modified config.json file
    with open(config_filename_input, mode='w') as updated_config_data:
        json.dump(config_data, updated_config_data, sort_keys=False, indent=4)

    return


def check_coordinates(x, y, z):
    """
    Calculates centerlines shape.

    Parameters
    ----------
    x : np.array
        (N,) array containing the node coordinates for all 'N' capillaries.
    y : np.array
        (N,) array containing the node coordinates for all 'N' capillaries.
    z : np.array
        (N,) array containing the node coordinates for all 'N' capillaries.

    Returns
    ----------
    coord_x : float
        x axis shape.
    coord_y : float
        y axis shape.
    coord_z : float
        z axis shape.
    """

    # Calculate max coordinates
    xmax = np.amax(x)
    ymax = np.amax(y)
    zmax = np.amax(z)

    # Calculate shape
    coord_x = int(xmax)+1
    coord_y = int(ymax)+1
    coord_z = int(zmax)+1

    return (coord_x, coord_y, coord_z)
