def squared_link_radius_clogging(minimum_link_radius,
                                 initial_link_squared_radius,
                                 modified_link_locations):
    """
    Modifies link radius of specific set of links (e.g., inlet face, mineralized, eroded links).
    """

    # Calculate link radius
    modified_link_squared_radius = initial_link_squared_radius

    # Calculate initial_link_radius
    for modified_link in range(modified_link_locations.size):

        # Identify link_radius_index
        link_radius_index = modified_link_locations[modified_link]

        # Calculate initial_link_squared_radius
        initial_link_squared_radius[link_radius_index] = minimum_link_radius**2

    return modified_link_squared_radius
