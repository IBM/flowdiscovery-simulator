#!/usr/bin/env python3

import argparse
import json

import h5py
import numpy as np
from mayavi import mlab

if __name__ == '__main__':

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Flow visualization')
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

    # Load centerlines input file
    with open(arg.out_folder + '/centerlines.json', mode='r') as file:
        data = json.load(file, encoding='utf-8')

    # Extract node geometry arrays from JSON
    nodes = sorted(data['graph']['nodes'], key=lambda node: int(node['id']))
    x = np.array([node['metadata']['node_coordinates']['x'] for node in nodes])
    y = np.array([node['metadata']['node_coordinates']['y'] for node in nodes])
    z = np.array([node['metadata']['node_coordinates']['z'] for node in nodes])

    # Extract link geometry arrays from JSON
    edges = sorted(data['graph']['edges'], key=lambda edge: int(edge['id']))
    source = np.array([int(edge['source']) for edge in edges])
    target = np.array([int(edge['target']) for edge in edges])
    link_length = np.array([edge['metadata']['link_length'] for edge in edges])

    with h5py.File(arg.out_folder + '/static_results.h5', "r") as file:
        # Load flow rate input file
        Q = np.squeeze(np.array(file['flow_rate_' + arg.axis], dtype=np.double))
        q = Q * 1e12  # convert to nL/s

        # Load flow speed input file
        V = np.squeeze(np.array(file['flow_speed_' + arg.axis], dtype=np.double))
        v = V * 1e3  # convert to mm/s

    # Calculate node origins (X, Y, Z) and node-to-node vectors (U[0], U[1], U[2])
    X = x[source]
    Y = y[source]
    Z = z[source]
    U = np.array([x[target] - x[source],
                  y[target] - y[source],
                  z[target] - z[source]])

    # Invert origin and components if (q < 0)
    X[q < 0] = x[target][q < 0]
    Y[q < 0] = y[target][q < 0]
    Z[q < 0] = z[target][q < 0]
    U[:, q < 0] *= -1
    u = U / np.linalg.norm(U, axis=0)

    # Calculate flow rate and flow speed components
    Q = u * (np.finfo(np.double).eps + np.abs(q))
    V = u * (np.finfo(np.double).eps + np.abs(v))

    # Creating 3D visualization for flow rate along the flow axis
    mlab.figure(size=(800, 700), bgcolor=(0, 0, 0))
    mlab.quiver3d(X, Y, Z, Q[0], Q[1], Q[2], line_width=2, scale_mode='scalar', scale_factor=1,
                  colormap='viridis', scalars=np.linalg.norm(U, axis=0))
    mlab.outline()
    mlab.orientation_axes()
    mlab.vectorbar(title='Q [nL/s]', orientation='vertical')
    if (arg.show):
        mlab.show()
    else:
        mlab.savefig(arg.out_folder + '/flow_rate_plot_' + arg.axis + '.png')

    # Creating 3D visualization for flow speed along the flow axis
    mlab.figure(size=(800, 700), bgcolor=(0, 0, 0))
    mlab.quiver3d(X, Y, Z, V[0], V[1], V[2], line_width=2, scale_mode='scalar', scale_factor=1,
                  colormap='viridis', scalars=np.linalg.norm(U, axis=0))
    mlab.outline()
    mlab.orientation_axes()
    mlab.vectorbar(title='V [mm/s]', orientation='vertical')
    if (arg.show):
        mlab.show()
    else:
        mlab.savefig(arg.out_folder + '/flow_speed_plot_' + arg.axis + '.png')
