#!/usr/bin/env python3


"""
This program creates 2d and 3d plots of geometry and fluid flow parameters.
"""

import argparse
import os

import numpy as np
from modules.dynamic import read_network
from modules.visualization import create_datasets, plot_3d_parameter, plot_histogram

if __name__ == '__main__':

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Geometry and flow parametric plots',
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
                information on flow speed and pressure gradients distribution.""")

    # Required arguments
    parser.add_argument('--voxel',
                        action='store',
                        metavar='VOXEL_SIZE',
                        type=float,
                        required=True,
                        help='Voxel size [m].')
    parser.add_argument('in_folder',
                        action='store',
                        metavar='IN_FOLDER',
                        type=str,
                        help='Directory containing input files.')

    # Optional arguments
    parser.add_argument('--axis',
                        action='store',
                        choices=['x', 'y', 'z'],
                        metavar='FLOW_AXIS',
                        type=str,
                        required=False,
                        help='Flow axis.')
    parser.add_argument('--centerlines_parameter',
                        action='store',
                        choices=['link_diameters', 'node_diameters'],
                        metavar='centerlines_parameter',
                        type=str,
                        required=False,
                        help='Centerlines parameter.')
    parser.add_argument('--dataset',
                        action='store',
                        choices=['centerlines', 'static'],
                        metavar='dataset',
                        type=str,
                        required=False,
                        help='dataset of centerlines and flow parameters.')
    parser.add_argument('--final_time_step_index',
                        action='store',
                        metavar='final_time_step_index',
                        type=int,
                        required=False,
                        help='Final time step. Used for parameter evolution analysis.')
    parser.add_argument('--initial_time_step_index',
                        action='store',
                        metavar='initial_time_step_index',
                        type=int,
                        required=False,
                        help='Initial time step. Default at initial conditions.',
                        default=0)
    parser.add_argument('--mayavi3d_type',
                        action='store',
                        choices=['points3d', 'quiver3d'],
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
    parser.add_argument('--plot_type',
                        action='store',
                        choices=['3d_plot', 'histogram'],
                        metavar='Plot_type',
                        type=str,
                        required=False,
                        help='Plot type.')
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
                        choices=['flow_rate', 'flow_speed', 'permeability', 'pressure',
                                 'pressures'],
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

    # Print information on geometry evolution due to pore-scale processes
    print('FLOWSIMULATOR::PLOT PARAMETERS SAYS:')

    # Define out folder
    out_folder = os.path.join(arg.in_folder)

    # Time evolution condition
    if (arg.final_time_step_index):
        if arg.final_time_step_index < arg.initial_time_step_index:
            print('Final time step must be greater than initial time step.')
            exit()

    # Use read_network(out_folder, voxel_size):
    x = read_network(arg.in_folder, arg.voxel)[7]
    y = read_network(arg.in_folder, arg.voxel)[8]
    z = read_network(arg.in_folder, arg.voxel)[9]

    # Use read_network(out_folder, voxel_size)
    sources = read_network(arg.in_folder, arg.voxel)[0]
    targets = read_network(arg.in_folder, arg.voxel)[1]

    # Number of links
    number_of_links = len(sources)

    # Calculate node origins (X, Y, Z) and node-to-node vectors (U[0], U[1], U[2])
    X = x[sources]
    Y = y[sources]
    Z = z[sources]
    U = np.array([x[targets] - x[sources],
                  y[targets] - y[sources],
                  z[targets] - z[sources]])

    # Create variable datasets including parameter array at the initial conditions
    (parameter,
     parameter_var,
     parameter_var_init) = create_datasets(arg.dataset,
                                           out_folder,
                                           arg.voxel,
                                           arg.axis,
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

    # Prepare scale for 3d plots
    if arg.scale_plot == 'auto':
        scale_vmin = None
        scale_vmax = None

    elif arg.scale_plot == 'set_plot_scale':
        scale_vmin = arg.vmin
        scale_vmax = arg.vmax

    elif arg.scale_plot == 'initial_conditions':
        scale_vmin = np.min(parameter_var_init)
        scale_vmax = np.max(parameter_var_init)

    # Plot parameters
    if arg.plot_type == '3d_plot':
        plot_3d_parameter(parameter,
                          parameter_var,
                          parameter_var_init,
                          scale_vmin,
                          scale_vmax,
                          arg.final_time_step_index,
                          arg.plot_size,
                          arg.mayavi3d_type,
                          x, y, z,
                          X, Y, Z,
                          U,
                          arg.axis,
                          out_folder,
                          arg.show,
                          arg.scale_plot,
                          arg.ratio,
                          arg.transparency)

    elif arg.plot_type == 'histogram':
        plot_histogram(parameter,
                       parameter_var,
                       parameter_var_init,
                       arg.final_time_step_index,
                       arg.axis,
                       out_folder,
                       arg.show,
                       arg.scale_plot,
                       arg.ratio,
                       number_of_links)
