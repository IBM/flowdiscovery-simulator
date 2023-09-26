import os

import numpy as np


# Calculate rock porosity using capillary geometry
def porosity_centerlines_method(centerlines_data,
                                link_length,
                                link_squared_radius,
                                voxel_size):
    """
    Import parameters from h5 files for visualization.

    Parameters
    ----------
    link_lengths : np.array
        (N, 3) array containing the source-to-target distance for all 'N' capillaries.
    link_squared_radius : np.array
        array containing link radius of 'n' capillaries [voxels^2].
    voxel_size : float
        size of the voxel in [m].
    void_space_volume_ratio : float
        Ratio between void_space_volume at the initial conditions and the current time step.

    Returns
    -------
    porosity : float
       Updated rock porosity at the current time step [-].
    """

    # Calculate porosity from centerlines.json
    pore_volume = np.sum(np.pi * link_squared_radius * link_length) * voxel_size**3

    # Extract node geometry arrays from JSON of rock sample centerlines
    nodes = sorted(centerlines_data['graph']['nodes'], key=lambda node: int(node['id']))
    x = np.array([node['metadata']['node_coordinates']['x'] for node in nodes], dtype=float)
    y = np.array([node['metadata']['node_coordinates']['y'] for node in nodes], dtype=float)
    z = np.array([node['metadata']['node_coordinates']['z'] for node in nodes], dtype=float)

    # Define maximum and minimum x, y and z values
    xmax = np.amax(x)
    ymax = np.amax(y)
    zmax = np.amax(z)

    # Calculate dimensions
    X = xmax + 1
    Y = ymax + 1
    Z = zmax + 1

    # Calculate rock volume
    rock_volume = (X) * (Y) * (Z) * voxel_size**3

    # Calculate rock porosity at the initial conditions
    porosity = pore_volume / rock_volume

    return porosity, rock_volume, pore_volume


# Calculate rock porosity using capillary geometry
def porosity_microCT_method(binary_image_filename,
                            centerlines_data,
                            in_folder,
                            voxel_size):
    """
    Import parameters from h5 files for visualization.

    Parameters
    ----------
    link_lengths : np.array
        (N, 3) array containing the source-to-target distance for all 'N' capillaries.
    link_squared_radius : np.array
        array containing link radius of 'n' capillaries [voxels^2].
    voxel_size : float
        size of the voxel in [m].
    void_space_volume_ratio : float
        Ratio between void_space_volume at the initial conditions and the current time step.

    Returns
    -------
    porosity : float
       Updated rock porosity at the current time step [-].
    """

    # Extract node geometry arrays from JSON of rock sample centerlines
    nodes = sorted(centerlines_data['graph']['nodes'], key=lambda node: int(node['id']))
    x = np.array([node['metadata']['node_coordinates']['x'] for node in nodes], dtype=float)
    y = np.array([node['metadata']['node_coordinates']['y'] for node in nodes], dtype=float)
    z = np.array([node['metadata']['node_coordinates']['z'] for node in nodes], dtype=float)

    # Define maximum and minimum x, y and z values
    xmax = np.amax(x)
    ymax = np.amax(y)
    zmax = np.amax(z)
    xmin = np.amin(x)
    ymin = np.amin(y)
    zmin = np.amin(z)

    X = (xmax - xmin) + 1
    Y = (ymax - ymin) + 1
    Z = (zmax - zmin) + 1

    image_dimensions = [int(X), int(Y), int(Z)]

    # Calculate rock volume
    rock_volume = np.prod(image_dimensions) * voxel_size**3

    # Define binary_image_path
    binary_image_path = os.path.join(in_folder, f'{binary_image_filename}.raw')

    # Define rock
    rock = digitalRock(binary_image_path, image_dimensions)

    # Calculate binary-image-based porosity
    porosity_CT = (1 - np.mean(rock.volume))

    # Calculate pore ratio
    pore_volume = rock_volume * porosity_CT

    # Calculate porosity
    porosity = pore_volume / rock_volume

    return porosity, rock_volume, pore_volume


# Define digital rock class
class digitalRock():
    """
    The digital rock class contains the algorithm used for contrast filtering of the rock samples.
    An instance to this class takes as input the path of the raw NxNxN cube file and the number of
    pixels on each size of the cube, N.

    It is important to note, however, that this code won't work as is for the full paralellpided
    files, as it expects cubic volumes. This can be easily circonvented, but may be unusable due to
    memory requirements.

    The equalize function carries out the contrast enhancement filter. As default, the filter cuts
    off the histogram at the 99.8 percentile, but this can be changed in the "clip" variable.
    """

    # Define rock (method N)
    def __init__(self, path, image_dimensions):
        self.image_dimensions = image_dimensions
        self.volume = np.fromfile(path, dtype=np.uint8).reshape(image_dimensions)

    def multp(self, A):
        return np.uint8((A*self.scale)*255)

    def equalize(self, clip=99.8):
        self.filter_thr = int(np.percentile(self.volume, clip))
        clipped_volume = np.clip(self.volume, 0, self.filter_thr)
        self.minInt = clipped_volume.min()
        self.maxInt = clipped_volume.max()
        self.scale = 1/(self.maxInt-self.minInt)

        clipped_volume -= self.minInt
        filtered = np.array(list(map(self.multp, clipped_volume)))

        return filtered
