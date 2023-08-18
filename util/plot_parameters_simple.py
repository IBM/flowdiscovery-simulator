#!/usr/bin/env python3


"""
This program creates 2d and 3d plots of geometry and fluid flow parameters.

Parameters
----------
out_folder : str
    directory containing output files.
voxel : float
    size of the voxel in [m].
parameter : str
    geometry of fluid flow parameter to create 2d or 3d plot.
axis: str
    flow axis. Req: flow rate, flow speed, pressures and pressure grads.
filename : str
    name of the input JSON file.
plot_histogram : bool
    plot histogram. Used for distribution analysis.
print_flow_info: bool
    print fluid flow additional information.
plot_size : int
    plot size (width, height).
vmin : float
    scale minimum value [-].
vmax : float
    scale maximum value [-].
show : bool
    show plot instead of saving to file.

Returns
-------
3d plot : parameter_plot.png
    points plot for nodal parameters or quiver 3d plot of capillary link parameters.
2d histogram plot : parameter_dist.png
    capillary distribution of nodal or capillary link parameters.
fluid flow information : str
    number of capillary links with zero or negative flow speed or pressure gradients.
"""

import argparse
import os

import h5py
import matplotlib.pyplot as plt
import numpy as np
from mayavi import mlab
from modules.dynamic import read_network
from modules.pressure_gradients import pressure_gradients_simple
from modules.info_parameters import info_parameters

if __name__ == '__main__':

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Geometry and fluid flow parameters visualization',
        epilog="""\
            Returns
            -------
            3d plots : parameter_plot_00000n.png
                3d plot of nodal or link parameters at 'n' time step.
            2d plots : parameter_plot_00000n.png
                2d plot of nodal or link parameters at 'n' time step.
            distributions : parameter_dist_00000n.png
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
                        metavar='OUT_FOLDER',
                        type=str,
                        help='Directory containing output files.')
    parser.add_argument('--parameter',
                        action='store',
                        choices=['link_diameters', 'node_diameters', 'flow_rate',
                                 'flow_speed', 'pressures', 'pressure_gradients'],
                        metavar='parameter',
                        type=str,
                        required=True,
                        help='Geometry of fluid flow parameter to create 2d or 3d plot.')

    # Optional arguments
    parser.add_argument('--axis',
                        action='store',
                        choices=['x', 'y', 'z'],
                        metavar='FLOW_AXIS',
                        type=str,
                        required=False,
                        help='Flow axis. Req: flow rate, flow speed, pressures and pressure grads.')
    parser.add_argument('--plot_histogram',
                        action='store_true',
                        default=False,
                        help='Plot histogram. Used for distribution analysis.')
    parser.add_argument('--print_flow_info',
                        action='store_true',
                        default=False,
                        help='Print pressure gradients and flow speed additional information.')
    parser.add_argument('--plot_size',
                        action='store',
                        metavar='plot_size',
                        type=str,
                        required=False,
                        default=(800, 700),
                        help='Plot size (width, height).')
    parser.add_argument('--show',
                        action='store_true',
                        default=False,
                        help='Show plot instead of saving to file.')
    parser.add_argument('--transparency',
                        action='store_true',
                        default=False,
                        help='Enable 3d plot transparency.')
    parser.add_argument('--vmin',
                        action='store',
                        metavar='vmin',
                        type=float,
                        required=False,
                        help='Scale minimum value [-].',
                        default=None)
    parser.add_argument('--vmax',
                        action='store',
                        metavar='vmax',
                        type=float,
                        required=False,
                        help='Scale maximum value [-].',
                        default=None)

    arg = parser.parse_args()

    # Import centerlines parameters using read_network
    sources = read_network(arg.in_folder, arg.voxel)[0]
    targets = read_network(arg.in_folder, arg.voxel)[1]

    # Import visualization parameters
    scale = info_parameters[arg.parameter]['3d_plot']['scale']
    units_conversion = info_parameters[arg.parameter]['units_conversion']

    # Set parameter arrays
    # Import node diameters from centerlines_00000n.json
    if arg.parameter == 'node_diameters':
        # Import node diameters using read_network
        parameter_var = read_network(arg.in_folder, arg.voxel)[5] * units_conversion

    # Import link diameters from centerlines_00000n.json
    elif arg.parameter == 'link_diameters':
        # Import link diameters using read_network
        parameter_var = read_network(arg.in_folder, arg.voxel)[4] * units_conversion

    # Import flow rate from static_results_00000n.h5
    elif arg.parameter == 'flow_rate':
        with h5py.File(os.path.join(arg.in_folder, 'static_results.h5'), "r") as file:
            # Load flow rate input file
            PARAMETER_VAR = np.squeeze(np.array(file[f'flow_rate_{arg.axis}'], dtype=np.double))

            # Convert return variable to apropriate units
            parameter_var = PARAMETER_VAR * units_conversion

    # Import flow speed from static_results_00000n.h5
    elif arg.parameter == 'flow_speed':
        with h5py.File(os.path.join(arg.in_folder, 'static_results.h5'), "r") as file:
            # Load flow speed input file
            PARAMETER_VAR = np.squeeze(np.array(file[f'flow_speed_{arg.axis}'], dtype=np.double))

            # Convert return variable to appropriate units
            parameter_var = PARAMETER_VAR * units_conversion

        # Fluid flow analysis
        if (arg.print_flow_info):
            # Print links id with zero velocity
            print(f'Links with V = 0 : {parameter_var.size - np.count_nonzero(parameter_var)}')

            # Print links id with negative velocity
            print(f'Links with V < 0 : {len(list(filter(lambda x: (x < 0), parameter_var)))}')

    # Import pressures from static_results_00000n.h5
    elif arg.parameter == 'pressures':
        with h5py.File(os.path.join(arg.in_folder, 'static_results.h5'), "r") as file:
            # Load pressure input file
            P = np.squeeze(np.array(file['pressures_' + arg.axis], dtype=np.double)) \
                * units_conversion

        # Calculate parameter array
        parameter_var = P - P.min()

    # Compute and import pressure gradients
    elif arg.parameter == 'pressure_gradients':
        # Compute pressure gradients
        parameter_var = pressure_gradients_simple(arg.in_folder, arg.axis, arg.voxel)[0] \
            * units_conversion

        # Fluid flow analysis
        if (arg.print_flow_info):
            # Import condition for zero pressure gradients
            zero_delta_P = pressure_gradients_simple(arg.in_folder, arg.axis, arg.voxel)[1]

            # Print links id with zero and reverse pressure gradients
            print(f'Links delta_P = 0 : {len(list(filter(lambda x: (x == 1), zero_delta_P)))}')
            print(f'Links delta_P < 0 : {len(list(filter(lambda x: (x < 0), parameter_var)))}')

    # Import link coordinates using read_network
    x = read_network(arg.in_folder, arg.voxel)[7]
    y = read_network(arg.in_folder, arg.voxel)[8]
    z = read_network(arg.in_folder, arg.voxel)[9]

    # Calculate node origins (X, Y, Z) and node-to-node vectors (U[0], U[1], U[2])
    X = x[sources]
    Y = y[sources]
    Z = z[sources]
    U = np.array([x[targets] - x[sources],
                  y[targets] - y[sources],
                  z[targets] - z[sources]])

    # Invert origin and components if (q < 0)
    u = U / np.linalg.norm(U, axis=0)

    # # Define tranmsparent background
    if (arg.transparency):
        mlab_bgcolor = (1, 1, 1)
        mlab_fgcolor = (0.2, 0.2, 0.2)

    # Define black backgroung
    else:
        mlab_bgcolor = (0.1, 0.1, 0.1)
        mlab_fgcolor = None

    # 3d plots of nodal parameters
    if any([arg.parameter == 'node_diameters', arg.parameter == 'pressures']):

        # POINTS Creating 3D visualization for node diameters along the flow axis
        mlab.figure(size=arg.plot_size,  bgcolor=mlab_bgcolor, fgcolor=mlab_fgcolor)
        if arg.parameter == 'node_diameters':
            mlab.points3d(x, y, z, parameter_var, colormap='viridis', vmin=arg.vmin, vmax=arg.vmax)

        # POINTS Creating 3D visualization for pressures along the flow axis
        elif arg.parameter == 'pressures':
            mlab.points3d(x, y, z, parameter_var,
                          scale_mode='none', scale_factor=5, colormap='viridis', vmin=arg.vmin,
                          vmax=arg.vmax)

        # SDefine parameters for scalarbar
        mlab.scalarbar(title=f"{info_parameters[arg.parameter]['3d_plot']['title']} \
[{scale if scale !=1 else ''}\
{info_parameters[arg.parameter]['units']}]", orientation='vertical')

    # 3d plots of capillary link parameters
    elif any([arg.parameter == 'link_diameters',
              arg.parameter == 'flow_rate',
              arg.parameter == 'flow_speed',
              arg.parameter == 'pressure_gradients']):

        # Calculate flow rate and flow speed components
        PARAMETER_VAR = u * (np.finfo(np.double).eps + np.abs(parameter_var))

        # QUIVER Creating 3D visualization for capillary link parameters along the flow axis
        mlab.figure(size=arg.plot_size, bgcolor=mlab_bgcolor, fgcolor=mlab_fgcolor)
        mlab.quiver3d(X, Y, Z, PARAMETER_VAR[0], PARAMETER_VAR[1], PARAMETER_VAR[2],
                      line_width=2, scale_mode='scalar', scale_factor=1, colormap='viridis',
                      scalars=np.linalg.norm(U, axis=0), vmin=arg.vmin, vmax=arg.vmax)
        mlab.vectorbar(title=f"{info_parameters[arg.parameter]['3d_plot']['title']} \
[{scale if scale !=1 else ''}\
{info_parameters[arg.parameter]['units']}]", orientation='vertical')

    # Additional 3d plot parameters
    mlab.outline()
    mlab.orientation_axes()

    # Show plots
    if (arg.show):
        mlab.show()
    else:
        mlab.savefig(os.path.join(arg.in_folder, 'snapshots', f'{arg.parameter}_plot_\
{arg.axis}.png'))

    # Plot histogram
    if (arg.plot_histogram):
        plt.figure(dpi=125)
        plt.hist(parameter_var, bins=150)
        plt.xlabel(f"{info_parameters[arg.parameter]['histogram']['xlabel']} \
[{scale if scale !=1 else ''}\
{info_parameters[arg.parameter]['units']}]", fontsize=16)
        plt.ylabel(info_parameters[arg.parameter]['histogram']['ylabel'], fontsize=16)

        if (arg.show):
            plt.show()
        else:
            plt.savefig(os.path.join(arg.in_folder, 'snapshots',
                        f'{arg.parameter}_dist_{arg.axis}.png'), dpi=125, bbox_inches='tight')
