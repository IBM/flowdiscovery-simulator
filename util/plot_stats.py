#!/usr/bin/env python3

import argparse

import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Centerline statistics.')
    parser.add_argument('odir',
                        action='store',
                        metavar='OUT_FOLDER',
                        type=str,
                        help='Directory containing output files.')
    parser.add_argument('--voxel',
                        action='store',
                        metavar='VOXEL_SIZE',
                        type=float,
                        help='Voxel size [m].')
    arg = parser.parse_args()

    # Load and re-format input file
    length, tortuosity, diameter = np.loadtxt(arg.odir + '/centerlines.stat',
                                              delimiter=',',
                                              dtype=float,
                                              usecols=[0, 1, 2],
                                              unpack=True)

    # Filter valid entries only
    v = (length > 0)
    length, tortuosity, diameter = length[v], tortuosity[v], diameter[v]

    # Convert "voxel" units to micrometers
    diameter *= arg.voxel / 1e-6
    length *= arg.voxel / 1e-6

    # Save length vs. tortuosity plot
    plt.figure(dpi=125)
    plt.xscale('log')
    plt.scatter(x=length, y=tortuosity, s=40, c=diameter, cmap=plt.cm.viridis,
                edgecolors='face', linewidths=3, alpha=0.85)
    plt.xlabel(r'Capillary length [$\mu\mathrm{m}$]', fontsize=16)
    plt.ylabel(r'Capillary tortuosity', fontsize=16)
    plt.colorbar().set_label(r'Average capillary diameter [$\mu\mathrm{m}$]',
                             size=16)
    plt.savefig(arg.odir + "/length_tort.png", dpi=125, bbox_inches='tight')

    # Save diameter vs. tortuosity plot
    plt.figure(dpi=125)
    plt.xscale('log')
    plt.scatter(x=diameter, y=tortuosity, s=40, c=length, cmap=plt.cm.viridis,
                edgecolors='face', alpha=0.85)
    plt.xlabel(r'Average capillary diameter [$\mu\mathrm{m}$]', fontsize=16)
    plt.ylabel(r'Capillary tortuosity', fontsize=16)
    plt.colorbar().set_label(r'Capillary length [$\mu\mathrm{m}$]',
                             size=16)
    plt.savefig(arg.odir + "/diam_tort.png", dpi=125, bbox_inches='tight')

    # Save diameter vs. length plot
    plt.figure(dpi=125)
    plt.xscale('log')
    plt.yscale('log')
    plt.scatter(x=diameter, y=length, s=40, c=tortuosity, cmap=plt.cm.viridis,
                edgecolors='face', alpha=0.85)
    plt.xlabel(r'Average capillary diameter [$\mu\mathrm{m}$]', fontsize=16)
    plt.ylabel(r'Capillary length [$\mu\mathrm{m}$]', fontsize=16)
    plt.colorbar().set_label(r'Capillary tortuosity', size=16)
    plt.savefig(arg.odir + "/diam_length.png", dpi=125, bbox_inches='tight')
