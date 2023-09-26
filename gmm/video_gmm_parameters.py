#!/usr/bin/env python3


"""
This program creates videos based on 3d plots of geometry, fluid flow, and pore-scale
process parameters during the geometry evolution of porous media.
"""

import argparse
import glob
import json
import os

import cv2
import numpy as np
from modules.dynamic import read_network
from modules.visualization import (create_datasets, plot_3d_parameter, plot_histogram,
                                   scale_parameters_video)

if __name__ == '__main__':

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Geometry evolution parametric videos',
        epilog="""\
            Returns
            -------
            parameter_video : .mp4
                3d-plot-based video containing the parameter spatiotemporal evolution.
            parameter_video_snapshots (optional) :
                folder containing the snapshots of the parameter spatiotemporal evolution.""")

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
    parser.add_argument('--dataset',
                        action='store',
                        choices=['centerlines', 'dissolution', 'erosion', 'flow',
                                 'geometry', 'precipitation', 'static'],
                        metavar='dataset',
                        type=str,
                        required=True,
                        help='dataset of centerlines, flow, geometry, and pore-scale parameters.')
    parser.add_argument('--final_time_step_index',
                        action='store',
                        metavar='final_time_step_index',
                        type=int,
                        required=True,
                        help='Final time step index. Used for parameter evolution analysis.')
    parser.add_argument('in_folder',
                        action='store',
                        metavar='IN_FOLDER',
                        type=str,
                        help='Directory containing input files.')
    parser.add_argument('--video_type',
                        action='store',
                        choices=['3d_plot', 'histogram'],
                        default='3d_plot',
                        metavar='video_type',
                        type=str,
                        required=True,
                        help='Video type.')

    # Optional arguments
    parser.add_argument('--centerlines_parameter',
                        action='store',
                        choices=['node_diameters'],
                        metavar='centerlines_parameter',
                        type=str,
                        required=False,
                        help='Centerlines parameter.')
    parser.add_argument('--deposition_parameter',
                        action='store',
                        choices=['deposition_onset',
                                 'deposition_rate',
                                 'deposition_shear_stress'],
                        metavar='deposition_parameter',
                        type=str,
                        required=False,
                        help='Deposition parameter.')
    parser.add_argument('--erosion_parameter',
                        action='store',
                        choices=['erosion_onset',
                                 'erosion_rate',
                                 'erosion_time_scale',
                                 'erosion_shear_stress'],
                        metavar='erosion_parameter',
                        type=str,
                        required=False,
                        help='Erosion parameter.')
    parser.add_argument('--flow_parameter',
                        action='store',
                        choices=['maximum_flow_rate',
                                 'pressure_gradients',
                                 'reynolds_number',
                                 'wall_shear_stress'],
                        metavar='flow_parameter',
                        type=str,
                        required=False,
                        help='Flow parameter.')
    parser.add_argument('--geometry_parameter',
                        action='store',
                        choices=['accumulated_volume_simulation',
                                 'accumulated_volume_time_step',
                                 'aspect_ratio',
                                 'link_diameters',
                                 'link_length',
                                 'pore_volume',
                                 'porosity',
                                 'porosity_time_step',
                                 'reactive_area',
                                 'total_accumulated_volume_simulation',
                                 'total_accumulated_volume_time_step',
                                 'void_space_volume',
                                 'void_space_volume_ratio'],
                        metavar='geometry_parameter',
                        type=str,
                        required=False,
                        help='Geometry parameter.')
    parser.add_argument('--initial_time_step_index',
                        action='store',
                        metavar='initial_time_step_index',
                        type=int,
                        required=False,
                        help='Initial time step. Default at initial conditions.',
                        default=0)
    parser.add_argument('--mayavi3d_type',
                        action='store',
                        choices=['quiver3d', 'points3d'],
                        metavar='mayavi3d_type',
                        type=str,
                        required=False,
                        help='Mayavi 3d plot type for link parameters',
                        default='quiver3d')
    parser.add_argument('--transparency',
                        action='store_true',
                        default=False,
                        help='Enable transparency.')
    parser.add_argument('--precipitation_parameter',
                        action='store',
                        choices=['maximum_precipitation_rate',
                                 'precipitation_number_clogged_links',
                                 'precipitation_number_clogged_inlet_links',
                                 'precipitation_clogging_evolution',
                                 'precipitation_rate'],
                        metavar='precipitation_parameter',
                        type=str,
                        required=False,
                        help='Precipation parameter.')
    parser.add_argument('--save_snapshots',
                        action='store_true',
                        default=False,
                        help='Save plots of each time step in parameter_video_snapshots folder.')
    parser.add_argument('--scale_plot',
                        action='store',
                        choices=['set', 'video'],
                        metavar='scale_plot',
                        type=str,
                        required=False,
                        default='video',
                        help='Scale type of 3d plots.')
    parser.add_argument('--static_parameter',
                        action='store',
                        choices=['pressures', 'flow_rate', 'flow_speed', 'permeability'],
                        metavar='static_parameter',
                        type=str,
                        required=False,
                        help='Static parameter.')
    parser.add_argument('--video_duration',
                        action='store',
                        metavar='video_duration',
                        type=int,
                        required=False,
                        default=10,
                        help='Video duration in seconds.')
    parser.add_argument('--video_size',
                        action='store',
                        metavar='video_size',
                        type=int,
                        required=False,
                        default=(800, 700),
                        help='Video size (width, height).')
    parser.add_argument('--vmax',
                        action='store',
                        metavar='vmax',
                        type=float,
                        required=False,
                        help='Maximum scale value.',
                        default=False)
    parser.add_argument('--vmin',
                        action='store',
                        metavar='vmin',
                        type=float,
                        required=False,
                        help='Minimum scale value.',
                        default=False)

    arg = parser.parse_args()

    # Default parameters
    default_flow_axis = ['x', 'y', 'z']

    # print information on geometry evolution due to pore-scale processes
    print('FLOWSIMULATOR::VIDEO GMM PARAMETERS SAYS:')

    # Define h5 files folder
    out_folder = os.path.join(arg.in_folder, arg.case_filename, f'{arg.dataset}_results')

    # Time evolution condition
    if (arg.final_time_step_index):
        if arg.final_time_step_index < arg.initial_time_step_index:
            print('Final time step must be greater than initial time step.')
            exit()

    # Load case parameters from case_m.json
    with open(os.path.join(arg.case_folder, f'{arg.case_filename}.json'),
              mode='r') as case_file_data:
        case_data = json.load(case_file_data)

    # Define case_sample
    case_sample = case_data['simulation']['sample']
    voxel_size = case_sample['parameters']['voxel_size']

    # Use read_network(out_folder, voxel_size):
    x = read_network(arg.in_folder, voxel_size)[7]
    y = read_network(arg.in_folder, voxel_size)[8]
    z = read_network(arg.in_folder, voxel_size)[9]

    # Use read_network(out_folder, voxel_size)
    sources = read_network(arg.in_folder, voxel_size)[0]
    targets = read_network(arg.in_folder, voxel_size)[1]

    # Number of links
    number_of_links = sources.size

    # Calculate node origins (X, Y, Z) and node-to-node vectors (U[0], U[1], U[2])
    X = x[sources]
    Y = y[sources]
    Z = z[sources]
    U = np.array([x[targets] - x[sources],
                  y[targets] - y[sources],
                  z[targets] - z[sources]])

    # Invert origin and components if (q < 0)
    u = U / np.linalg.norm(U, axis=0)

    # Calculate scale for video
    print('Calculating scale for video.')

    # Define axis
    flow_axis = case_data['simulation']['flow']['experiment']['flow_axis']

    # Dfine minimum and maximum scale values during the simulation time steps
    scale_vmin, scale_vmax, parameter = scale_parameters_video(arg.initial_time_step_index,
                                                               arg.final_time_step_index,
                                                               out_folder,
                                                               arg.dataset,
                                                               voxel_size,
                                                               default_flow_axis[flow_axis],
                                                               arg.centerlines_parameter,
                                                               arg.flow_parameter,
                                                               arg.geometry_parameter,
                                                               arg.static_parameter,
                                                               arg.deposition_parameter,
                                                               arg.erosion_parameter,
                                                               arg.precipitation_parameter,
                                                               arg.video_type,
                                                               arg.scale_plot)

    # Define scale maximum and minimum values
    if arg.scale_plot == 'set':
        # Define scale values
        scale_vmin = arg.vmin
        scale_vmax = arg.vmax

    # Define video and folder names
    if arg.video_type == '3d_plot':
        # Define video name
        video_name = f'video_{parameter}_{arg.initial_time_step_index}_\
{arg.final_time_step_index}.mp4'

        # Define snapshots folder name
        snapshots_folder_name = f'{parameter}_video_snapshots'

    elif arg.video_type == 'histogram':
        # Define video name
        video_name = f'video_dist_{parameter}_{arg.initial_time_step_index}_\
{arg.final_time_step_index}.mp4'

        # Define snapshots folder name
        snapshots_folder_name = f'{parameter}_dist_video_snapshots'

    # Create folder for video snapshots
    try:
        os.mkdir(os.path.join(out_folder, snapshots_folder_name))
    except FileExistsError:
        pass

    # Calculate frames per second according to video_duration input
    frames_per_second = arg.final_time_step_index / arg.video_duration

    # Initial time
    video_time_step_index = arg.initial_time_step_index

    # Start scale calculation
    print(f'Creating {arg.final_time_step_index} frames.')

    # Iterative process
    while video_time_step_index <= arg.final_time_step_index:

        # Create variable datasets
        (parameter,
         parameter_array,
         parameter_array_init) = create_datasets(arg.dataset,
                                                 out_folder,
                                                 voxel_size,
                                                 default_flow_axis[flow_axis],
                                                 arg.initial_time_step_index,
                                                 video_time_step_index,
                                                 arg.centerlines_parameter,
                                                 arg.flow_parameter,
                                                 arg.geometry_parameter,
                                                 arg.static_parameter,
                                                 arg.deposition_parameter,
                                                 arg.erosion_parameter,
                                                 arg.precipitation_parameter,
                                                 arg.video_type,
                                                 arg.scale_plot,
                                                 False)

        # Create 3d_plots
        if arg.video_type == '3d_plot':
            plot_3d_parameter(parameter,
                              parameter_array,
                              parameter_array_init,
                              scale_vmin,
                              scale_vmax,
                              arg.final_time_step_index,
                              arg.video_size,
                              arg.mayavi3d_type,
                              x, y, z,
                              X, Y, Z,
                              U,
                              default_flow_axis[flow_axis],
                              out_folder,
                              False,
                              arg.scale_plot,
                              False,
                              arg.transparency,
                              video_time_step_index)

        # Create histogram
        elif arg.video_type == 'histogram':
            plot_histogram(parameter,
                           parameter_array,
                           parameter_array_init,
                           arg.final_time_step_index,
                           default_flow_axis[flow_axis],
                           out_folder,
                           False,
                           arg.scale_plot,
                           False,
                           number_of_links,
                           video_time_step_index)

        # Print video frame time step
        print(f'Time step {video_time_step_index} of {arg.final_time_step_index}')

        # Next iteration
        video_time_step_index += 1

    # Count files
    _, _, files = next(os.walk(os.path.join(out_folder, snapshots_folder_name, '')))

    # Video mp4 format
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    video = cv2.VideoWriter(os.path.join(out_folder, video_name),
                            fourcc,
                            frames_per_second,
                            arg.video_size,
                            isColor=True)

    # Create video using 3d plots
    for img in sorted(glob.glob(os.path.join(out_folder, snapshots_folder_name, '*.png'))):
        img = cv2.imread(img)
        img = cv2.resize(img, arg.video_size)
        video.write(img)

    # Release video
    video.release()

    # Remove snapshots folder
    if arg.save_snapshots is False:
        os.remove(os.path.join(out_folder, snapshots_folder_name))
