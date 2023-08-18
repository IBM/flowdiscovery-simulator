#!/usr/bin/env python3

import argparse

from modules.dynamic import (create_video, draw_interfaces, get_time_step, list_output_files,
                             read_network)

if __name__ == '__main__':

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Dynamic flow visualization')
    parser.add_argument('folder',
                        action='store',
                        metavar='OUT_FOLDER',
                        type=str,
                        help='Directory containing output files.')
    parser.add_argument('--voxel',
                        action='store',
                        metavar='VOXEL_SIZE',
                        type=float,
                        required=True,
                        help='Voxel size [m].')
    parser.add_argument('--speed',
                        action='store',
                        metavar='SPEED_MULTIPLIER',
                        type=float,
                        default=1.0,
                        required=False,
                        help='Speed up playback by SPEED times.')
    parser.add_argument('--scale',
                        action='store',
                        metavar='SCALING_FACTOR',
                        type=float,
                        default=1.0,
                        required=False,
                        help='Divide the capillary diameter by SCALING_FACTOR when plotting.')
    parser.add_argument('--plane',
                        action='store',
                        metavar='PLANE',
                        type=str,
                        default=None,
                        required=False,
                        choices=['xy', 'yz', 'zx'],
                        help='Specify \'xy\', \'yz\' or \'zx\' for plotting.')
    parser.add_argument('--clean',
                        action='store_true',
                        default=False,
                        help='Clean image files after creating movie.')
    arg = parser.parse_args()

    # Read network geometry data
    _, _, origins, lengths, diameters = read_network(arg.folder, arg.voxel)

    # Read list of snapshot file
    files = list_output_files(arg.folder)

    # TODO(rneumann): Process different plots in parallel using processes
    # TODO(rneumann): Process different snaptshots in a plot in parallel using pools

    # Draw interfaces for each snapshot
    parameters = ((file,
                   origins,
                   lengths,
                   diameters,
                   arg.folder,
                   arg.plane,
                   arg.scale) for file in files)
    [draw_interfaces(parameter) for parameter in parameters]

    # Build movie out of images
    time_step = get_time_step(files)
    create_video(arg.folder, 'dyn_interfaces', time_step, speed=arg.speed, clean=arg.clean)
