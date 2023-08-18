#!/usr/bin/env python3

import argparse
import json

import matplotlib.pyplot as plt
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
    parser.add_argument('--voxel',
                        action='store',
                        metavar='VOXEL_SIZE',
                        type=float,
                        required=True,
                        help='Voxel size [m].')
    parser.add_argument('--filename',
                        action='store',
                        metavar='FILENAME.json',
                        type=str,
                        required=False,
                        default='centerlines',
                        help='Name of the input JSON file.')
    parser.add_argument('--show',
                        action='store_true',
                        default=False,
                        help='Show plot instead of saving to file.')
    arg = parser.parse_args()

    # Load centrelines input file
    with open(arg.out_folder + '/' + arg.filename + '.json', mode='rb') as file:
        data = json.load(file, encoding='utf-8')

    # Extract node geometry arrays from JSON
    nodes = sorted(data['graph']['nodes'], key=lambda node: int(node['id']))
    x = np.array([node['metadata']['node_coordinates']['x'] for node in nodes])
    y = np.array([node['metadata']['node_coordinates']['y'] for node in nodes])
    z = np.array([node['metadata']['node_coordinates']['z'] for node in nodes])
    R2 = np.array([node['metadata']['node_squared_radius'] for node in nodes])

    # Creating 3D visualization
    mlab.figure(size=(800, 700), bgcolor=(0.1, 0.1, 0.1))
    mlab.points3d(x, y, z, 2.0 * np.sqrt(R2) * arg.voxel / 1.0e-6, colormap='viridis')
    mlab.outline()
    mlab.orientation_axes()
    mlab.scalarbar(title='D [1e-6 m]', orientation='vertical')
    if (arg.show):
        mlab.show()
    else:
        mlab.savefig(arg.out_folder + '/' + arg.filename + '_plot.png')

    # Plot pore size distribution
    plt.figure(dpi=125)
    plt.hist(2.0 * np.sqrt(R2) * arg.voxel / 1.0e-6, bins=20)
    plt.xlabel(r'Capillary voxel diameter [$\mu\mathrm{m}$]', fontsize=16)
    plt.ylabel(r'Number of capillary voxels', fontsize=16)
    plt.savefig(arg.out_folder + '/' + arg.filename + '_dist.png', dpi=125, bbox_inches='tight')
