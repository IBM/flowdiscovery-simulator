#!/usr/bin/env python3

"""
This program calculates the initial porosity of a digital rock sample using the microCT Method.
"""

import argparse
import json
import os

from modules.rock_porosity import porosity_microCT_method

if __name__ == '__main__':

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Calculate initial porosity using the microCT method.',
        epilog="""\
                Returns:
                ----------------------------------------
                initial porosity
                initial pore volume
                initial rock volume
               """)

    # Required arguments
    parser.add_argument('--case_folder',
                        action='store',
                        metavar='CASE_FOLDER',
                        type=str,
                        required=True,
                        help='Directory containing CASE.json files.')
    parser.add_argument('--case_filename',
                        action='store',
                        metavar='CASE.json',
                        type=str,
                        required=True,
                        help='Name of the input CASE.json file, containing flow parameters.')
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
    arg = parser.parse_args()

    # Start Geometry Modification Module (GMM)
    print('FLOWSIMULATOR::GMM SAYS:')

    # Open centerlines file
    centerlines_filename_input = os.path.join(arg.in_folder,
                                              f'{arg.centerlines_filename}.json')
    with open(centerlines_filename_input, mode='r') as centerlines_file_data:
        centerlines_data = json.load(centerlines_file_data)

    # Load case parameters from case_m.json
    with open(os.path.join(arg.case_folder, f'{arg.case_filename}.json'),
              mode='r') as case_file_data:
        case_data = json.load(case_file_data)

    # Import data from case_data
    voxel_size = case_data['simulation']['sample']['parameters']['voxel_size']

    # Calculate initial porosity
    initial_porosity = porosity_microCT_method(arg.binary_image_filename,
                                               centerlines_data,
                                               arg.in_folder,
                                               voxel_size)

    # Print initial_porosity
    print('microCT_method')
    print(f'Initial porosity = {initial_porosity[0]}')
    print(f'Initial pore volume = {initial_porosity[1]}')
    print(f'Initial rock volume = {initial_porosity[2]}')
