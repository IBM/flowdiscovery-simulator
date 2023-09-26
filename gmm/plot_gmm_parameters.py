#!/usr/bin/env python3


"""
This program creates 2d and 3d plots of geometry, fluid flow, and pore-scale
process parameters during the geometry evolution of porous media.
"""

import argparse
import json
import os

import numpy as np
from modules.dynamic import read_network
from modules.visualization import (create_datasets, plot_2d_parameter, plot_3d_parameter,
                                   plot_histogram)

if __name__ == '__main__':

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Geometry evolution parametric plots',
        epilog="""\
            Returns
            -------
            3d plots : parameter_plot_00000n.png
                3d plot of nodal or link parameters at 'n' time step.
            2d plots : parameter_plot_00000n.png
                2d plot of nodal or link parameters at 'n' time step.
            histograms : parameter_dist_00000n.png
                capillary distribution of nodal or link parameters at 'n' time step.
            fluid flow information : str
                information on flow speed and pressure gradients.""")

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
                        choices=['dissolution', 'erosion', 'flow', 'geometry', 'precipitation',
                                 'static'],
                        metavar='dataset',
                        type=str,
                        required=True,
                        help='dataset of centerlines, flow, geometry, and pore-scale parameters.')
    parser.add_argument('--final_time_step_index',
                        action='store',
                        metavar='final_time_step_index',
                        type=int,
                        required=True,
                        help='Final time step. Used for parameter evolution analysis.')
    parser.add_argument('in_folder',
                        action='store',
                        metavar='IN_FOLDER',
                        type=str,
                        help='Directory containing input files.')
    parser.add_argument('--plot_type',
                        action='store',
                        choices=['2d_plot', '3d_plot', 'histogram'],
                        metavar='Plot_type',
                        type=str,
                        required=True,
                        help='Plot type.')

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
                                 'erosion_shear_stress',
                                 'erosion_time_scale'],
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
                        help='Mayavi 3d plot type for link parameters.',
                        default='quiver3d')
    parser.add_argument('--plot_size',
                        action='store',
                        metavar='plot_size',
                        type=str,
                        required=False,
                        default=(800, 700),
                        help='Plot size (width, height).')
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
    parser.add_argument('--print_flow_info',
                        action='store_true',
                        default=False,
                        help='Print fluid flow information.')
    parser.add_argument('--ratio',
                        action='store_true',
                        default=False,
                        help='Plot ratio during t=0 and t=time_step_index.')
    parser.add_argument('--scale_plot',
                        action='store',
                        choices=['auto', 'initial_conditions', 'set_plot_scale'],
                        metavar='scale_plot',
                        type=str,
                        required=False,
                        default='auto',
                        help='Scale type of 3d plots.')
    parser.add_argument('--show',
                        action='store_true',
                        default=False,
                        help='Show plot instead of saving to file.')
    parser.add_argument('--static_parameter',
                        action='store',
                        choices=['pressures', 'flow_rate', 'flow_speed', 'permeability'],
                        metavar='static_parameter',
                        type=str,
                        required=False,
                        help='Static parameter.')
    parser.add_argument('--transparency',
                        action='store_true',
                        default=False,
                        help='Enable 3d plot transparency.')
    parser.add_argument('--vmax',
                        action='store',
                        metavar='vmax',
                        type=float,
                        required=False,
                        help='Maximum scale value.',
                        default=None)
    parser.add_argument('--vmin',
                        action='store',
                        metavar='vmin',
                        type=float,
                        required=False,
                        help='Minimum scale value.',
                        default=None)

    arg = parser.parse_args()

    # Default parameters
    default_flow_axis = ['x', 'y', 'z']

    # Print information on geometry evolution due to pore-scale processes
    print('FLOWSIMULATOR::PLOT GMM PARAMETERS SAYS:')

    # Define out folder
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

    # Define axis
    flow_axis = case_data['simulation']['flow']['experiment']['flow_axis']

    # Create variable datasets including parameter array at the initial conditions
    (parameter,
     parameter_array,
     parameter_array_init) = create_datasets(arg.dataset,
                                             out_folder,
                                             voxel_size,
                                             default_flow_axis[flow_axis],
                                             arg.initial_time_step_index,
                                             arg.final_time_step_index,
                                             arg.centerlines_parameter,
                                             arg.flow_parameter,
                                             arg.geometry_parameter,
                                             arg.static_parameter,
                                             arg.erosion_parameter,
                                             arg.precipitation_parameter,
                                             arg.plot_type,
                                             arg.scale_plot,
                                             arg.print_flow_info,
                                             arg.ratio)

    # PREPARE SCALE FOR 3D PLOTS
    # Define auto scale
    if arg.scale_plot == 'auto':
        scale_vmin = None
        scale_vmax = None

    # Set an input scale
    elif arg.scale_plot == 'set_plot_scale':
        scale_vmin = arg.vmin
        scale_vmax = arg.vmax

    # Define scale based on the initial conditions
    elif arg.scale_plot == 'initial_conditions':
        if any([arg.plot_type == '3d_plot', arg.plot_type == 'histogram']):
            scale_vmin = np.min(parameter_array_init)
            scale_vmax = np.max(parameter_array_init)

    # Plot parameters
    if arg.plot_type == '2d_plot':
        plot_2d_parameter(arg.case_filename,
                          arg.case_folder,
                          parameter,
                          parameter_array,
                          arg.initial_time_step_index,
                          arg.final_time_step_index,
                          out_folder,
                          arg.show,
                          arg.scale_plot)

    elif arg.plot_type == '3d_plot':
        plot_3d_parameter(parameter,
                          parameter_array,
                          parameter_array_init,
                          scale_vmin,
                          scale_vmax,
                          arg.final_time_step_index,
                          arg.plot_size,
                          arg.mayavi3d_type,
                          x, y, z,
                          X, Y, Z,
                          U,
                          default_flow_axis[flow_axis],
                          out_folder,
                          arg.show,
                          arg.scale_plot,
                          arg.ratio,
                          arg.transparency)

    elif arg.plot_type == 'histogram':
        plot_histogram(parameter,
                       parameter_array,
                       parameter_array_init,
                       arg.final_time_step_index,
                       default_flow_axis[flow_axis],
                       out_folder,
                       arg.show,
                       arg.scale_plot,
                       arg.ratio,
                       number_of_links)
