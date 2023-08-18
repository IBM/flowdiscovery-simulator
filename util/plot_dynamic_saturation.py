#!/usr/bin/env python3

import argparse
import json
import os

import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spint
from matplotlib.animation import FuncAnimation


def init():
    plt.xlabel('Injected volume [PV]', fontsize=14)
    plt.ylabel('Fluid saturation', fontsize=14)
    plt.legend(loc='center right', fontsize=14)
    plt.grid(True, color='gray', linestyle=':')
    plt.tick_params(labelsize=14)
    return fig,


def update(n):
    h2o_line.set_data(injected_pore_volumes[:n], h2o_saturation[:n])
    co2_line.set_data(injected_pore_volumes[:n], co2_saturation[:n])
    return co2_line, h2o_line


if __name__ == '__main__':

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Dynamic flow visualization')
    parser.add_argument('folder',
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
    parser.add_argument('--voxel',
                        action='store',
                        metavar='VOXEL_SIZE',
                        type=np.double,
                        required=True,
                        help='Voxel size [m].')
    arg = parser.parse_args()

    # Reading centerlines.json
    with open(os.path.join(arg.folder, 'centerlines.json'), mode='rb') as file:
        data = json.load(file, encoding='utf-8')

    # Extract node geometry arrays from JSON
    nodes = sorted(data['graph']['nodes'], key=lambda node: int(node['id']))
    axis_coord = np.array([node['metadata']['node_coordinates'][arg.axis] for node in nodes],
                          dtype=np.double)

    # Extract link geometry arrays from JSON
    edges = sorted(data['graph']['edges'], key=lambda edge: int(edge['id']))
    sources = np.array([int(edge['source']) for edge in edges])
    targets = np.array([int(edge['target']) for edge in edges])
    link_length = np.array([edge['metadata']['link_length'] for edge in edges])
    link_squared_radius = np.array([edge['metadata']['link_squared_radius'] for edge in edges])

    # Calculate capillary volume in SI units
    pore_volume = np.pi * np.sum(link_squared_radius * link_length) * arg.voxel**3

    # Derive inlet capillary indexes
    inlet_nodes = np.nonzero(axis_coord == axis_coord.min())[0]
    flag = np.zeros_like(sources)
    for node in inlet_nodes:
        flag = np.logical_or(flag, np.logical_or(node == sources, node == targets))
    inlet_capillaries = np.nonzero(flag)[0]

    # Reading datasets into Numpy arrays
    folder = os.path.join(arg.folder, 'snapshots')
    files = sorted([os.path.join(folder, f) for f in os.listdir(folder) if f.endswith('.h5')])
    time = np.array([h5py.File(f, 'r')['current_time'][(0, 0)] for f in files])
    co2_flow_rate = np.array([h5py.File(f, 'r')['flow_rate'][(0)] for f in files])
    h2o_saturation = np.array([h5py.File(f, 'r')['fluid_saturation'][(0, 0)] for f in files])
    co2_saturation = np.array([h5py.File(f, 'r')['fluid_saturation'][(0, 1)] for f in files])

    # Calculate inlet CO2 flow rate [PV/s] and cumulative injected volume [PV]
    inlet_co2_flow_rate = co2_flow_rate[:, inlet_capillaries].sum(axis=1) / pore_volume
    injected_pore_volumes = spint.cumtrapz(inlet_co2_flow_rate, time)

    # Create animation
    fig, ax = plt.subplots()
    co2_line, = ax.plot(injected_pore_volumes, co2_saturation[:-1:], '-', lw=1.5, label='CO2')
    h2o_line, = ax.plot(injected_pore_volumes, h2o_saturation[:-1:], '-', lw=1.5, label='H2O')
    animation = FuncAnimation(fig,
                              func=update,
                              frames=len(files),
                              init_func=init,
                              blit=True)
    animation.save(filename=os.path.join(arg.folder, 'dyn_saturation.mp4'),
                   extra_args=['-vcodec', 'libx264'],
                   writer='ffmpeg',
                   dpi=150,
                   fps=24)
