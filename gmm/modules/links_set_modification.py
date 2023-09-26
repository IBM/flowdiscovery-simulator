import numpy as np


def links_set_modification(initial_radius,
                           modified_location):
    """
    Modifies a set of link radius within a capillary nwtwork
    """

    return np.array([initial_radius[modified_location[link]]
                    for link in range(modified_location.size)])
