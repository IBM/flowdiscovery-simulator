import os

import h5py
import numpy as np
from modules.pressure_gradients import pressure_gradients_simple
from modules.rock_porosity import porosity_centerlines_method, porosity_microCT_method


def calculate_flow_geometry_parameters(binary_image_filename,
                                       case_data,
                                       centerlines_data,
                                       default_extract_axis,
                                       geometry_results_filename,
                                       geometry_results_folder,
                                       in_folder,
                                       initial_link_radius,
                                       initial_void_space_volume,
                                       link_length,
                                       link_radius,
                                       out_folder,
                                       previous_void_space_volume,
                                       print_progress,
                                       static_results_filename,
                                       time_step_index):

    """
    Calculates flow and geometry parameters.
    """

    # Import flow experiment parameters from case.json
    case_experiment = case_data['simulation']['flow']['experiment']
    flow_axis = case_experiment['flow_axis']

    # Import sample parameters
    case_sample = case_data['simulation']['sample']
    voxel_size = case_sample['parameters']['voxel_size']

    # Import phasic properites from case.json
    case_liquid = case_data['simulation']['phases']['liquid']
    liquid_density = case_liquid['properties']['liquid_density']
    liquid_dynamic_viscosity = case_liquid['properties']['liquid_dynamic_viscosity']

    # Define flow_results folder
    static_results_folder = os.path.join(out_folder, 'static_results')

    # Calculate pressure gradients
    pressure_gradients = pressure_gradients_simple(default_extract_axis[flow_axis],
                                                   centerlines_data,
                                                   static_results_folder,
                                                   f'{static_results_filename}_\
{time_step_index:06}.h5')[0]

    # Calculate reactive area
    reactive_area = 2 * np.pi * link_radius * link_length

    # Calculate void space volume
    void_space_volume = np.pi * link_radius**2 * link_length * voxel_size**3

    # Calculate void space volume ratio
    void_space_volume_ratio = np.array(np.sum(void_space_volume) /
                                       np.sum(initial_void_space_volume))

    # ACCUMULATED VOLUME WITH RESPECT TO THE INITIAL CONDITIONS
    # Calculate accumulated volume (PER CAPILLARY) between the current and initial time steps
    accumulated_volume_simulation = initial_void_space_volume - void_space_volume

    # Calculate total accumulated volume (SUM) between the current and initial time steps
    total_accumulated_volume_simulation = np.array(np.sum(accumulated_volume_simulation))

    # CALCULATE INITIAL POROSITY
    # Extract sample data
    case_sample = case_data['simulation']['sample']

    # Extract initial porosity
    if case_sample['parameters']['initial_porosity']['porosity_method'] == 'input':
        initial_porosity = case_sample['parameters']['initial_porosity']['value']

    # Calculate initial porosity
    elif case_sample['parameters']['initial_porosity']['porosity_method'] == 'calculation':
        # Calculate initial porosity using binary_image.raw
        if os.path.exists(os.path.join(in_folder, f'{binary_image_filename}.raw')) is True:

            # Calculate porosity using binary_image.raw
            (initial_porosity, initial_rock_volume, _) = \
                porosity_microCT_method(binary_image_filename,
                                        centerlines_data,
                                        in_folder,
                                        voxel_size)
            if (print_progress):
                print('Porosity calculated with the porosity_microCT method')

        # Calculate initial porosity capillary network void space volume
        else:
            # Calculate porosity using centerlines.json
            (initial_porosity, initial_rock_volume, _) = \
                porosity_centerlines_method(centerlines_data,
                                            link_length,
                                            initial_link_radius**2,
                                            voxel_size)

            if (print_progress):
                print('Porosity calculated with the centerlines approximation method')

        if (print_progress):
            # Print initial porosity
            print(f'Initial porosity : {initial_porosity}')

    # CALCULATE UPDATED POROSITY
    # Calculate updated  porosity at the current time step
    porosity = initial_porosity * void_space_volume_ratio

    # Calulate pore_volume
    pore_volume = initial_rock_volume * porosity

    # CALCULATE ACCUMULATED VOLUME AND POROSITY WITH RESPECT TO THE PREVIOUS ITERATION
    if time_step_index == 0:
        # Define previous_void_space_volume
        previous_void_space_volume = initial_void_space_volume

        # Define initial_porosity
        previous_porosity = initial_porosity

    else:
        # Extract previous time step accumulated volume
        with h5py.File(os.path.join(geometry_results_folder, f'{geometry_results_filename}_\
{time_step_index - 1:06}.h5'), 'r') as geometry_file_data:
            # Load void_space_volume input file
            previous_array_void_space_volume = np.array(geometry_file_data['void_space_volume'],
                                                        dtype=np.double)

            # Extract previous_void_space_volume
            previous_void_space_volume = np.squeeze(previous_array_void_space_volume)

            # Load porosity input file
            previous_array_porosity = np.array(geometry_file_data['porosity'], dtype=np.double)

            # Extract previous porosity
            previous_porosity = np.squeeze(previous_array_porosity)

    # Calculate current and previous accumulated volume (PER CAPILLARY)
    previous_acumulated_volume = initial_void_space_volume - previous_void_space_volume
    curent_acumulated_volume = initial_void_space_volume - void_space_volume

    # Calculate total accumulated volume (SUM) between the current and initial time steps
    accumulated_volume_time_step = curent_acumulated_volume - previous_acumulated_volume
    total_accumulated_volume_time_step = np.array(np.sum(accumulated_volume_time_step))

    # Calculate porosity variation during current time step
    porosity_time_step = previous_porosity - porosity

    # Calculate aspect ratio
    aspect_ratio = 2 * link_radius / link_length

    # Calculate wall shear stress
    wall_shear_stress = abs(pressure_gradients) * link_radius / (2 * link_length)

    # CALCULATE REYNOLDS NUMBER
    # Extract flow speed
    with h5py.File(os.path.join(static_results_folder, f'{static_results_filename}_\
{time_step_index:06}.h5'), 'r') as static_file_data:
        # Load flow speed input file
        V = np.squeeze(np.array(static_file_data[f'flow_speed_\
{default_extract_axis[flow_axis]}'], dtype=np.double))

    # Calculate Reynolds number
    reynolds_number = liquid_density * V * 2 * (link_radius * voxel_size) \
        / liquid_dynamic_viscosity

    # CALCULATE MAX FLOW RATE
    # Extract flow rate
    with h5py.File(os.path.join(static_results_folder, f'{static_results_filename}_\
{time_step_index:06}.h5'), 'r') as static_file_data:
        # Load flow rate input file
        Q = np.squeeze(np.array(static_file_data[f'flow_rate_\
{default_extract_axis[flow_axis]}'], dtype=np.double))

    # Calculate max flow rate
    maximum_flow_rate = np.max(Q)

    return (accumulated_volume_simulation,
            accumulated_volume_time_step,
            aspect_ratio,
            maximum_flow_rate,
            pore_volume,
            porosity,
            porosity_time_step,
            pressure_gradients,
            reactive_area,
            reynolds_number,
            total_accumulated_volume_simulation,
            total_accumulated_volume_time_step,
            void_space_volume,
            void_space_volume_ratio,
            wall_shear_stress)


def save_flow_geometry_parameters(accumulated_volume_simulation,
                                  accumulated_volume_time_step,
                                  aspect_ratio,
                                  inlet_link_radius,
                                  link_length,
                                  link_radius,
                                  maximum_flow_rate,
                                  out_folder,
                                  pore_volume,
                                  porosity,
                                  porosity_time_step,
                                  pressure_gradients,
                                  print_progress,
                                  reactive_area,
                                  reynolds_number,
                                  save_flow_parameters,
                                  save_geometry_parameters,
                                  time_step_index,
                                  total_accumulated_volume_simulation,
                                  total_accumulated_volume_time_step,
                                  void_space_volume,
                                  void_space_volume_ratio,
                                  wall_shear_stress,):

    """
    Saves flow and geometry parameters to h5 files.
    """

    # Create HDF5 files for flow results and initialize them
    if (save_flow_parameters):
        # Define flow_results folder
        flow_results_folder = os.path.join(out_folder, 'flow_results')

        # Create geometry_results folder
        try:
            os.mkdir(flow_results_folder)
        except FileExistsError:
            pass

        # Define h5 file
        flow_results = h5py.File(os.path.join(flow_results_folder, f'flow_results_\
{time_step_index:06}.h5'), 'w')

        # Create pressure gradients 1D datasets for flow_results_00000n.h5
        flow_results.create_dataset('maximum_flow_rate', data=maximum_flow_rate)

        # Create pressure gradients datasets for flow_results_00000n.h5
        flow_results.create_dataset('pressure_gradients', data=pressure_gradients)

        # Calculate max Reynolds number
        maximum_reynolds_number = np.max(reynolds_number)

        # Create Max Reynolds number 1D datasets for flow_results_00000n.h5
        flow_results.create_dataset('maximum_reynolds_number', data=maximum_reynolds_number)

        # Create Reynolds number datasets for flow_results_00000n.h5
        flow_results.create_dataset('reynolds_number', data=reynolds_number)

        # Create wall shear stress datasets for flow_results_00000n.h5
        flow_results.create_dataset('wall_shear_stress', data=wall_shear_stress)

        # Close flow_results_00000n.h5 file
        flow_results.close()

    # Create HDF5 files for geometry results and initialize them
    if (save_geometry_parameters):
        # Define geometry_results folder
        geometry_results_folder = os.path.join(out_folder, 'geometry_results')

        # Create geometry_results folder
        try:
            os.mkdir(geometry_results_folder)
        except FileExistsError:
            pass

        # Define h5 file
        geometry_results = h5py.File(os.path.join(geometry_results_folder,
                                                  f'geometry_results_\
{time_step_index:06}.h5'), 'w')

        # Create aspect ratio datasets for geometry_results_00000n.h5
        geometry_results.create_dataset('aspect_ratio', data=aspect_ratio)

        # Create link length datasets for geometry_results_00000n.h5
        geometry_results.create_dataset('link_length', data=link_length)

        # Create link radius for geometry_results_00000n.h5
        geometry_results.create_dataset('link_radius', data=link_radius)

        # Create inlet link radius for geometry_results_00000n.h5
        geometry_results.create_dataset('inlet_link_radius', data=inlet_link_radius)

        # Create porosity datasets for geometry_results_00000n.h5
        geometry_results.create_dataset('porosity', data=porosity)

        # Create porosity datasets for geometry_results_00000n.h5
        geometry_results.create_dataset('porosity_time_step', data=porosity_time_step)

        # Create porosity datasets for geometry_results_00000n.h5
        geometry_results.create_dataset('pore_volume', data=porosity)

        # Print porosity (current time step)
        if (print_progress):
            print(f'Current porosity: {porosity} %')

        # Create reactive areadatasets for geometry_results_00000n.h5
        geometry_results.create_dataset('reactive_area', data=reactive_area)

        # ACCUMULATED VOLUME WITH RESPECT TO THE PREVIOUS TIME STEP
        # Create accumulated volume datasets (PER CAPILLARY) for geometry_results_00000n.h5
        geometry_results.create_dataset('accumulated_volume_time_step',
                                        data=accumulated_volume_time_step)

        # Create total accumulated volume (SUM) datasets for geometry_results_00000n.h5
        geometry_results.create_dataset('total_accumulated_volume_time_step',
                                        data=total_accumulated_volume_time_step)

        # ACCUMULATED VOLUME WITH RESPECT TO THE INITIAL CONDITIONS
        # Create accumulated volume datasets (PER CAPILLARY) for geometry_results_00000n.h5
        geometry_results.create_dataset('accumulated_volume',
                                        data=accumulated_volume_simulation)

        # Create total accumulated volume (SUM) datasets for geometry_results_00000n.h5
        geometry_results.create_dataset('total_accumulated_volume_simulation',
                                        data=total_accumulated_volume_simulation)

        # Create void-space volume datasets for geometry_results_00000n.h5
        geometry_results.create_dataset('void_space_volume', data=void_space_volume)

        # Create void-space volume ratio datasets for geometry_results_00000n.h5
        geometry_results.create_dataset('void_space_volume_ratio',
                                        data=void_space_volume_ratio)

        # Close geometry_results_00000n.h5 file
        geometry_results.close()

    return
