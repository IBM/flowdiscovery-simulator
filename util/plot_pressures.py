#!/usr/bin/env python3

import argparse
import json

import h5py
import numpy as np
from mayavi import mlab

if __name__ == '__main__':

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Centerline visualization')
    parser.add_argument('out_folder',
                        action='store',
                        metavar='OUT_FOLDER',
                        type=str,
                        help='Directory containing output files.')
    parser.add_argument('--axis',
                        action='store',
                        choices=['x', 'y', 'z'],
                        metavar='FLOW_AXIS',
                        type=str,
                        required=True,
                        help='Flow axis.')
    parser.add_argument('--show',
                        action='store_true',
                        default=False,
                        help='Show plot instead of saving to file.')
    arg = parser.parse_args()

    # Load centrelines input file
    with open(arg.out_folder + '/centerlines.json', mode='r', encoding='utf-8') as file:
        data = json.load(file)

    # Extract node geometry arrays from JSON
    nodes = sorted(data['graph']['nodes'], key=lambda node: int(node['id']))
    x = np.array([node['metadata']['node_coordinates']['x'] for node in nodes])
    y = np.array([node['metadata']['node_coordinates']['y'] for node in nodes])
    z = np.array([node['metadata']['node_coordinates']['z'] for node in nodes])

    # Load pressure input files
    with h5py.File(arg.out_folder + '/static_results.h5', "r") as file:
        P = np.squeeze(np.array(file['pressures_' + arg.axis], dtype=np.double))

    # Creating 3D visualization for flow along U
    mlab.figure(size=(800, 700), bgcolor=(0.1, 0.1, 0.1))
    mlab.points3d(x, y, z, P - P.min(), scale_mode='none', scale_factor=5, colormap='viridis')
    mlab.outline()
    mlab.orientation_axes()
    mlab.scalarbar(title='P - min(P) [Pa]', orientation='vertical')
    if (arg.show):
        mlab.show()
    else:
        mlab.savefig(arg.out_folder + '/pressure_plot_' + arg.axis + '.png')
