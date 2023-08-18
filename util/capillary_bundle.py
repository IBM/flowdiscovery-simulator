#!/usr/bin/env python3

import argparse
import json

import numpy as np

if __name__ == '__main__':

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Capillary bundle creation.')
    parser.add_argument('out_folder',
                        action='store',
                        metavar='OUT_FOLDER',
                        type=str,
                        help='Directory containing output files.')
    parser.add_argument('--porosity', '-p',
                        action='store',
                        metavar='POROSITY',
                        type=float,
                        required=True,
                        help='Sample porosity.')
    parser.add_argument('--number', '-n',
                        action='store',
                        metavar='NUMBER_OF_CAPILLARIES',
                        type=int,
                        required=False,
                        help='Number of capillaries.')
    arg = parser.parse_args()

    # Reading centerlines.json
    with open(arg.out_folder + '/centerlines.json', mode='r') as file:
        data = json.load(file, encoding='utf-8')

    # Extract link geometry arrays from JSON and calculate derived quantities
    edges = sorted(data['graph']['edges'], key=lambda edge: int(edge['id']))
    link_squared_radius = np.array([edge['metadata']['link_squared_radius'] for edge in edges],
                                   dtype=np.double)

    # Sample random capillaries from file
    N = link_squared_radius.size if (arg.number is None) else arg.number
    idx = np.random.randint(link_squared_radius.size, size=N)
    R2 = link_squared_radius[idx]
    L = np.sqrt(np.pi * R2.sum() / arg.porosity)

    # Calculate node positions
    x = np.array([(i % 2) * L for i in range(2 * N)])
    y = np.random.rand(N) * L
    z = np.random.rand(N) * L

    # Create graph metadata JSON object
    graph_metadata_obj = {'number_of_nodes': 2 * N,
                          'number_of_links': N}

    # Create nodes JSON object
    nodes_obj = [
        {
            'id': str(i),
            'metadata': {
                'node_squared_radius': R2[int(i/2)],
                'node_coordinates': {
                    'x': x[i],
                    'y': y[int(i/2)],
                    'z': z[int(i/2)]
                }
            }
        } for i in range(graph_metadata_obj['number_of_nodes'])]

    # Create edges JSON object
    edges_obj = [
        {
            'id': str(i),
            'source': str(2 * i),
            'target': str(2 * i + 1),
            'metadata': {
                'link_length': L,
                'link_squared_radius': R2[i]
            }
        } for i in range(graph_metadata_obj['number_of_links'])]

    # Create graph JSON object
    graph_obj = {'metadata': graph_metadata_obj,
                 'nodes': nodes_obj,
                 'edges': edges_obj}

    # Create full JSON object
    json_obj = {'graph': graph_obj}

    # Save output centerlines file with capillary bundle
    with open(arg.out_folder + '/capillary_bundle.json', 'w') as file:
        json.dump(json_obj, file)
