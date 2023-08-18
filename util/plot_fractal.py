#!/usr/bin/env python3

import argparse

import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Box-counting fractal dimension plot')
    parser.add_argument('out_folder',
                        action='store',
                        metavar='OUT_FOLDER',
                        type=str,
                        help='Directory containing output files.')
    parser.add_argument('--show',
                        action='store_true',
                        default=False,
                        help='Show plot instead of saving to file.')
    arg = parser.parse_args()

    # Creating graph for the normalised histogram
    plt.figure(dpi=125)
    plt.title("Box-counting fractal dimension", fontsize=18)
    plt.xlabel(r'$\epsilon$', fontsize=16)
    plt.ylabel(r'-$\frac{\partial \ln{N(\epsilon)}}{\partial \ln{\epsilon}}$', fontsize=16)
    plt.ylim(0, 3.1)
    plt.grid()

    # Reading histogram file from output folder and fitting the curve
    style = {'pore': 'rs-', 'surf': 'go-', 'rock': 'bv-'}
    for phase in ["pore", "surf", "rock"]:
        size, number = np.loadtxt(arg.out_folder + "/" + phase + "_frac_plot.dat", unpack=True)
        frac_dim = - np.gradient(np.log(number), size) / np.gradient(np.log(size), size)
        plt.semilogx(size, frac_dim, style[phase], label=phase)
    plt.legend(loc='best', fancybox=True, shadow=True)

    if (arg.show):
        plt.show()
    else:
        plt.savefig(arg.out_folder + "/frac_plot.png", dpi=125, bbox_inches='tight')
