import json
import os

import h5py
import matplotlib.pyplot as plt
import numpy as np
from mayavi import mlab
from modules.dynamic import read_network
from modules.info_parameters import info_parameters


# Import parameters and create datasets
def create_datasets(dataset,
                    out_folder,
                    voxel,
                    axis,
                    initial_time_step_index,
                    time_step_index,
                    centerlines_parameter,
                    flow_parameter,
                    geometry_parameter,
                    static_parameter,
                    deposition_parameter,
                    erosion_parameter,
                    precipitation_parameter,
                    plot_type,
                    scale_plot,
                    print_flow_info,
                    ratio=False):

    """
    Create parameter datasets for plots and videos.

    Parameters
    ----------
    dataset : str
        dataset of centerlines, erosion, flow, geometry, static, and precipitation parameters.
    out_folder : str
        folder containing '/gmm_case_n/'.
    voxel : float
        size of the voxel in [m].
    axis: str
        flow axis. Req: flow rate, flow speed, pressures and pressure grads.
    initial_time_step_index: int
        initial time step. Default at initial conditions.
    time_step_index: int
        final time step. Used for parameter evolution analysis when comparing initial
        conditions against specific time steps.
    flow_parameter : str
        name of the input FLOW_RESULTS.H5 file parameter:
        pressure_gradients,
        reynolds_number,
        wall_shear_stress.
    geometry_parameter : str
        name of the input GEOMETRY_RESULTS.H5 file parameter:
        accumulated_volume_simulation,
        accumulated_volume_time_step
        aspect_ratio,
        porosity,
        reactive_area,
        total_accumulated_volume_simulation,
        total_accumulated_volume_time_step
        void_space_volume,
        void_space_volume_ratio
    static_parameter : str
        name of the input STATIC_RESULTS.H5 file parameter:
        permeability,
        pressure,
        flow_rate,
        flow_speed.
    erosion_parameter : str
        name of the input EROSION_RESULTS.H5 file parameter:
        erosion_onset,
        erosion_rate,
        erosion_shear_stress,
        erosion_time_scale.
    precipitation_parameter : str
        name of the input PRECIPITATION_RESULTS.H5 file parameter:
        precipitation_clogged_capillaries,
        precipitation_clogging_evolution,
        precipitation_clogging_time_step,
        precipitation_number_clogged_links,
        precipitation_number_clogged_inlet_links.
    scale_plot : bool
        Scale type of 3d plots.
    print_flow_info : bool
        Print flow additional information.
    ratio : bool
        Plot ratio during t=0 and t=time_step_index.

    Returns
    -------
    parameter : str
        centerlines, erosion, flow, geometry, static, and precipitation parameters.
    parameter_array : np.array
        array containing parameter of 'n' capillaries [m].
    parameter_array_init: np.array
        array containing parameter of 'n' capillaries [m] at the initial conditions.
    """

    # Obtain parameters from centerlines dataset
    if dataset == 'centerlines':
        # Define parameter
        parameter = centerlines_parameter

        # Convert units
        units_conversion = info_parameters[parameter]['units_conversion']

        if parameter == 'node_diameters':
            # Import node diameters using read_network
            parameter_array = read_network(out_folder, voxel)[5] * units_conversion

    # Obtain parameters from deposition process dataset
    if dataset == 'deposition':
        # Define parameter
        parameter = deposition_parameter

        # Convert units
        units_conversion = info_parameters[parameter]['units_conversion']

        # Import pore-scale process parameters
        if any([parameter == 'deposition_onset',
                parameter == 'deposition_rate',
                parameter == 'deposition_shear_stress']):

            # Import h5 file with gmm_results
            results_filename = os.path.join(out_folder,
                                            f'{dataset}_results_{time_step_index:06}.h5')

            # Import parameter
            with h5py.File(results_filename, 'r') as deposition_file_data:
                parameter_array = (np.squeeze(np.array(deposition_file_data[parameter],
                                                       dtype=np.double))) * units_conversion

            # Define parameter value at initial time step
            if any([(ratio), scale_plot == 'initial_conditions', scale_plot == 'video',
                    scale_plot == 'set']):
                # Import h5 file with gmm_results
                results_filename = os.path.join(out_folder, f'{dataset}_results_\
{initial_time_step_index:06}.h5')

                # Import parameter
                with h5py.File(results_filename, 'r') as deposition_file_data_init:
                    deposition_array = np.array(deposition_file_data_init[parameter],
                                                dtype=np.double)
                    parameter_array_init = (np.squeeze(deposition_array)) * units_conversion

    # Obtain parameters from erosion process dataset
    elif dataset == 'erosion':
        # Define parameter
        parameter = erosion_parameter

        print(parameter)

        # Convert units
        units_conversion = info_parameters[parameter]['units_conversion']
        print(units_conversion)

        # Import pore-scale process parameters
        if any([parameter == 'erosion_onset',
                parameter == 'erosion_rate',
                parameter == 'erosion_shear_stress',
                parameter == 'erosion_time_scale']):

            # Import h5 file with gmm_results
            results_filename = os.path.join(out_folder,
                                            f'{dataset}_results_{time_step_index:06}.h5')

            # Import parameter
            with h5py.File(results_filename, 'r') as erosion_file_data:
                parameter_array = (np.squeeze(np.array(erosion_file_data[parameter],
                                                       dtype=np.double))) * units_conversion

            # Define parameter value at initial time step
            if any([(ratio), scale_plot == 'initial_conditions', scale_plot == 'video',
                    scale_plot == 'set']):
                # Import h5 file with gmm_results
                results_filename = os.path.join(out_folder, f'{dataset}_results_\
{initial_time_step_index:06}.h5')

                # Import parameter
                with h5py.File(results_filename, 'r') as erosion_file_data_init:
                    parameter_array_init = (np.squeeze(np.array(erosion_file_data_init[parameter],
                                                                dtype=np.double))) * \
                                                                    units_conversion

    # Obtain parameters from flow dataset
    elif dataset == 'flow':
        # Define parameter
        parameter = flow_parameter

        # Convert units
        units_conversion = info_parameters[parameter]['units_conversion']

        # Import parameters from flow results files
        if any([parameter == 'pressure_gradients',
                parameter == 'reynolds_number',
                parameter == 'wall_shear_stress']):
            # Import h5 file with gmm_results
            results_filename = os.path.join(out_folder,
                                            f'{dataset}_results_{time_step_index:06}.h5')

            # Import parameter
            with h5py.File(results_filename, 'r') as flow_file_data:
                parameter_array = (np.squeeze(np.array(flow_file_data[parameter],
                                                       dtype=np.double))) * \
                    units_conversion

            # Define parameter value at initial time step
            if any([(ratio), scale_plot == 'initial_conditions', scale_plot == 'video',
                    scale_plot == 'set']):
                # Import h5 file with gmm_results
                results_filename = os.path.join(out_folder, f'{dataset}_results_\
{initial_time_step_index:06}.h5')

                # Import parameter
                with h5py.File(results_filename, 'r') as flow_file_data_init:
                    parameter_array_init = (np.squeeze(np.array(flow_file_data_init[parameter],
                                                                dtype=np.double))) * \
                                                                    units_conversion

    # Obtain parameters from geometry dataset
    elif dataset == 'geometry':
        # Define parameter
        parameter = geometry_parameter

        # Convert units
        units_conversion = info_parameters[parameter]['units_conversion']

        # Import parameters from flow, geometry, or pore-scale processes results files
        if any([parameter == 'accumulated_volume_simulation',
                parameter == 'accumulated_volume_time_step',
                parameter == 'aspect_ratio',
                parameter == 'link_diameters',
                parameter == 'link_length',
                parameter == 'pore_volume',
                parameter == 'porosity',
                parameter == 'reactive_area',
                parameter == 'total_accumulated_volume_simulation',
                parameter == 'total_accumulated_volume_time_step',
                parameter == 'void_space_volume',
                parameter == 'void_space_volume_ratio']):

            # Import 2d parameters
            if any([parameter == 'pore_volume',
                    parameter == 'porosity',
                    parameter == 'void_space_volume_ratio',
                    parameter == 'total_accumulated_volume_simulation',
                    parameter == 'total_accumulated_volume_time_step']):
                # Create array containing parameter value, considering initial conditions
                parameter_array = np.zeros((time_step_index - initial_time_step_index) + 1,
                                           dtype=float)

                # Initial time step
                axis_time_step_index = initial_time_step_index

                # Create dataset for 2d plot horizontal axis
                while axis_time_step_index <= time_step_index:
                    # Import h5 file with gmm_results
                    results_filename = os.path.join(out_folder, f'{dataset}_results_\
{axis_time_step_index:06}.h5')

                    # Import parameter
                    with h5py.File(results_filename, 'r') as geometry_file_data:
                        # Define np.squeeze array argument
                        np_squeeze_array = geometry_file_data[f'{parameter}']

                        # Load void_space_volume_ratio input file
                        parameter_array[axis_time_step_index] = np.squeeze(np_squeeze_array) * \
                            units_conversion

                    # Next iteration
                    axis_time_step_index += 1

            elif parameter == 'link_diameters':
                # Import h5 file with gmm_results
                results_filename = os.path.join(out_folder,
                                                f'{dataset}_results_{time_step_index:06}.h5')

                # Import parameter
                with h5py.File(results_filename, 'r') as geometry_file_data:
                    parameter_array = (np.squeeze(np.array(geometry_file_data['link_radius'],
                                                           dtype=np.double))) * units_conversion

                # Define parameter value at initial time step
                if any([(ratio), scale_plot == 'initial_conditions', scale_plot == 'video',
                        scale_plot == 'set']):
                    # Import h5 file with gmm_results
                    results_filename = os.path.join(out_folder, f'{dataset}_results_\
{initial_time_step_index:06}.h5')

                    # Import parameter
                    with h5py.File(results_filename, 'r') as geometry_file_data_init:
                        # Define np_array
                        np_array = np.array(geometry_file_data_init['link_radius'],
                                            dtype=np.double)

                        # Define parameter_array_init
                        parameter_array_init = (np.squeeze(np_array)) * units_conversion

            else:
                # Import h5 file with gmm_results
                results_filename = os.path.join(out_folder,
                                                f'{dataset}_results_{time_step_index:06}.h5')

                # Import parameter
                with h5py.File(results_filename, 'r') as geometry_file_data:
                    parameter_array = (np.squeeze(np.array(geometry_file_data[parameter],
                                                           dtype=np.double))) * units_conversion

                # Define parameter value at initial time step
                if any([(ratio), scale_plot == 'initial_conditions', scale_plot == 'video',
                        scale_plot == 'set']):
                    # Import h5 file with gmm_results
                    results_filename = os.path.join(out_folder, f'{dataset}_results_\
{initial_time_step_index:06}.h5')

                    # Import parameter
                    with h5py.File(results_filename, 'r') as geometry_file_data_init:
                        # Define numpy arrayarray
                        np_parameter_array = np.array(geometry_file_data_init[parameter],
                                                      dtype=np.double)

                        # Define parameter_array_init
                        parameter_array_init = (np.squeeze(np_parameter_array)) * units_conversion

    # Obtain parameters from precipitation process dataset
    elif dataset == 'precipitation':
        # Define parameter
        parameter = precipitation_parameter

        # Convert units
        units_conversion = info_parameters[parameter]['units_conversion']

        # Import pore-scale process parameters
        if any([parameter == 'precipitation_number_clogged_links',
                parameter == 'precipitation_number_clogged_inlet_links',
                parameter == 'precipitation_clogging_evolution',
                parameter == 'precipitation_rate',
                parameter == 'maximum_precipitation_rate']):

            # Import pore-scale process parameters
            if any([parameter == 'precipitation_number_clogged_links',
                    parameter == 'precipitation_number_clogged_inlet_links',
                    parameter == 'precipitation_clogging_evolution',
                    parameter == 'maximum_precipitation_rate']):

                # Create array containing parameter value
                parameter_array = np.zeros((time_step_index - initial_time_step_index) + 1,
                                           dtype=float)

                # Initial time step
                axis_time_step_index = initial_time_step_index + 1

                # Create dataset for 2d plot horizontal axis
                while axis_time_step_index <= time_step_index:
                    process_results_filename = os.path.join(out_folder, f'{dataset}_results_\
{axis_time_step_index:06}.h5')

                    # Load input file
                    with h5py.File(process_results_filename, 'r') as precipitation_file_data:
                        # Define np_squeeze
                        np_squeeze = np.squeeze(precipitation_file_data[f'{parameter}'])

                        # Define parameter array
                        parameter_array[axis_time_step_index] = np_squeeze * units_conversion

                    # Next iteration
                    axis_time_step_index += 1

            elif parameter == 'precipitation_rate':
                # Import h5 file with gmm_results
                results_filename = os.path.join(out_folder,
                                                f'{dataset}_results_{time_step_index:06}.h5')

                # Import parameter
                with h5py.File(results_filename, 'r') as precipitation_file_data:
                    parameter_array = (np.squeeze(np.array(precipitation_file_data[parameter],
                                                           dtype=np.double))) * units_conversion

                # Define parameter value at initial time step
                if any([(ratio), scale_plot == 'initial_conditions', scale_plot == 'video',
                        scale_plot == 'set']):
                    # Import h5 file with gmm_results
                    results_filename = os.path.join(out_folder, f'{dataset}_results_\
{initial_time_step_index:06}.h5')

                    # Import parameter
                    with h5py.File(results_filename, 'r') as precipitation_file_data_init:
                        # Define numpy arrayarray
                        np_parameter_array = np.array(precipitation_file_data_init[parameter],
                                                      dtype=np.double)

                        # Define parameter_array_init
                        parameter_array_init = (np.squeeze(np_parameter_array)) * units_conversion

    # Obtain parameters from static dataset
    elif dataset == 'static':
        # Define parameter
        parameter = static_parameter

        # Convert units
        units_conversion = info_parameters[parameter]['units_conversion']

        # Import flow rate from static_results_00000n.h5
        if parameter == 'flow_rate':
            # Import h5 file with gmm_results
            results_filename = os.path.join(out_folder,
                                            f'{dataset}_results_{time_step_index:06}.h5')

            with h5py.File(results_filename, 'r') as static_file_data:
                # Load flow rate input file
                PARAMETER_ARRAY = np.squeeze(np.array(static_file_data[f'flow_rate_{axis}'],
                                                      dtype=np.double))

                # Define array containing parameter value
                parameter_array = PARAMETER_ARRAY * units_conversion

            # Define parameter value at initial time step
            if any([(ratio), scale_plot == 'initial_conditions', scale_plot == 'video',
                    scale_plot == 'set']):
                # Import h5 file with gmm_results
                results_filename = os.path.join(out_folder, f'{dataset}_results_\
{initial_time_step_index:06}.h5')

                with h5py.File(results_filename, 'r') as static_file_data_init:
                    # Define np_array
                    np_array = static_file_data_init[f'flow_rate_{axis}']

                    # Load flow rate input file
                    PARAMETER_ARRAY = np.squeeze(np.array(np_array, dtype=np.double))

                    # Define array containing parameter value
                    parameter_array_init = PARAMETER_ARRAY * units_conversion

        # Import flow speed from static_results_00000n.h5
        elif parameter == 'flow_speed':
            # Import h5 file with gmm_results
            results_filename = os.path.join(out_folder,
                                            f'{dataset}_results_{time_step_index:06}.h5')

            with h5py.File(results_filename, 'r') as static_file_data:
                # Load flow speed input file
                PARAMETER_ARRAY = np.squeeze(np.array(static_file_data[f'flow_speed_{axis}'],
                                                      dtype=np.double))

                # Convert return variable to appropriate units
                parameter_array = PARAMETER_ARRAY * units_conversion

            # Define parameter value at initial time step
            if any([(ratio), scale_plot == 'initial_conditions', scale_plot == 'video',
                    scale_plot == 'set']):
                # Import h5 file with gmm_results
                results_filename = os.path.join(out_folder, f'{dataset}_results_\
{initial_time_step_index:06}.h5')

                with h5py.File(results_filename, 'r') as static_file_data_init:
                    # Load flow speed input file
                    PARAMETER_ARRAY = np.squeeze(
                        np.array(static_file_data_init[f'flow_speed{axis}'], dtype=np.double))

                    # Convert return variable to appropriate units
                    parameter_array_init = PARAMETER_ARRAY * units_conversion

            # Fluid flow analysis
            if (print_flow_info):
                # Print links id with zero and negative speed
                print(f'Links with V = 0 : \
{parameter_array.size - np.count_nonzero(parameter_array)}')
                print(f'Links with V < 0 : \
{len(list(filter(lambda x: (x < 0), parameter_array)))}')

# Import permeability from static_results_00000n.h5
        elif parameter == 'permeability':
            # Define array containing parameter value
            parameter_array = np.zeros(time_step_index, dtype=float)

            # Initial time step
            axis_time_step_index = initial_time_step_index

            # Create dataset for 2d plot horizontal axis
            while axis_time_step_index <= time_step_index:
                results_filename = os.path.join(out_folder,
                                                f'{dataset}_results_{time_step_index:06}.h5')

                with h5py.File(results_filename, 'r') as static_file_data:
                    # Define np_squeeze_array
                    np_squeeze_array = np.squeeze(static_file_data[f'permeability_{axis}'],
                                                  dtype=np.double)

                    # Load permeability input file
                    parameter_array[axis_time_step_index] = np_squeeze_array * units_conversion

                # Next iteration
                axis_time_step_index += 1

        # Import pressures from static_results_00000n.h5
        elif parameter == 'pressures':
            # Load pressure input file
            results_filename = os.path.join(out_folder,
                                            f'{dataset}_results_{time_step_index:06}.h5')

            with h5py.File(results_filename, 'r') as static_file_data:
                P = np.squeeze(np.array(static_file_data[f'pressures_{axis}'],
                                        dtype=np.double)) * units_conversion

                # Define array containing parameter value
                parameter_array = P - P.min()

            # Define parameter value at initial time step
            if (ratio):
                results_filename = os.path.join(out_folder, f'{dataset}_results_\
{initial_time_step_index:06}.h5')

                with h5py.File(results_filename, 'r') as static_file_data_init:
                    P = np.squeeze(np.array(static_file_data_init[f'pressures_{axis}'],
                                            dtype=np.double)) * units_conversion

                    # Define array containing parameter value
                    parameter_array_init = P - P.min()

    # Return parameters of scalar and array values
    if any([plot_type == '2d_plot', ratio is False]):
        parameter_array_init = None

    # Return array values
    return parameter, parameter_array, parameter_array_init


# Calculate scale values for videos
def scale_parameters_video(initial_time_step_index,
                           time_step_index,
                           out_folder,
                           dataset,
                           voxel,
                           axis,
                           centerlines_parameter,
                           flow_parameter,
                           geometry_parameter,
                           static_parameter,
                           deposition_parameter,
                           erosion_parameter,
                           precipitation_parameter,
                           plot_type,
                           scale_plot):

    """
    Calculate maximum and minimum parameter values for plot scales.

    Parameters
    ----------
    initial_time_step_index: int
        initial time step. Default at initial conditions.
    time_step_index: int
        final time step. Used for parameter evolution analysis when comparing initial
        conditions against specific time steps.
    out_folder : str
        folder containing '/gmm_case_n/'.
    dataset : str
        dataset of centerlines, erosion, flow, geometry, static, and precipitation parameters.
    voxel : float
        size of the voxel in [m].
    axis: str
        flow axis. Req: flow rate, flow speed, pressures and pressure grads.
    flow_parameter : str
        name of the input FLOW_RESULTS.H5 file parameter: pressure_gradients, reynolds_number,
        wall_shear_stress.
    geometry_parameter : str
        name of the input GEOMETRY_RESULTS.H5 file parameter: aspect_ratio, porosity,
        reactive_area, void_space_volume, void_space_volume_ratio, accumulated_volume.
    static_parameter : str
        name of the input STATIC_RESULTS.H5 file parameter:
        permeability,
        ressure,
        flow_rate,
        flow_speed.
    erosion_parameter : str
        name of the input EROSION_RESULTS.H5 file parameter:
        erosion_onset,
        erosion_rate,
        erosion_shear_stress,
        erosion_time_scale.
    precipitation_parameter : str
        name of the input PRECIPITATION_RESULTS.H5 file parameter:
        precipitation_clogged_capillaries,
        precipitation_clogging_evolution,
    scale_plot : bool
        scale type of 3d plots.

    Returns
    -------
    scale_vmin : float
        minimum value of the specific parameter along the simulation time.
    scale_vmax : float
        maximum value of the specific parameter along the simulation time.
    parameter : str
        centerlines, erosion, flow, geometry, static, and precipitation parameters.
    """

    # Arrays containing maximum and minimum values of video scales during the evolution
    maximum_values = np.zeros(time_step_index - initial_time_step_index + 1, dtype=float)
    min_values = np.zeros(time_step_index - initial_time_step_index + 1, dtype=float)

    # Iterative process for maximum and minimum scale values computation
    scale_time_step_index = initial_time_step_index

    # Iterative process
    while scale_time_step_index <= time_step_index:

        # Create variable datasets
        (parameter, parameter_array, _) = create_datasets(dataset,
                                                          out_folder,
                                                          voxel,
                                                          axis,
                                                          initial_time_step_index,
                                                          scale_time_step_index,
                                                          centerlines_parameter,
                                                          flow_parameter,
                                                          geometry_parameter,
                                                          static_parameter,
                                                          deposition_parameter,
                                                          erosion_parameter,
                                                          precipitation_parameter,
                                                          plot_type,
                                                          scale_plot,
                                                          False,
                                                          False)

        # Scale max and min values
        maximum_values[scale_time_step_index] = max(parameter_array)
        min_values[scale_time_step_index] = min(parameter_array)

        # Next iteration
        scale_time_step_index += 1

    # Max and min scale values
    scale_vmin = min(maximum_values)
    scale_vmax = max(maximum_values)

    return scale_vmin, scale_vmax, parameter


# Plot parameters
def plot_3d_parameter(parameter,
                      parameter_array,
                      parameter_array_init,
                      scale_vmin,
                      scale_vmax,
                      final_time_step_index,
                      plot_size,
                      mayavi3d_type,
                      x, y, z,
                      X, Y, Z,
                      U,
                      axis,
                      out_folder,
                      show,
                      scale_plot,
                      ratio,
                      transparency,
                      video_time_step_index=0):

    """
    Calculate maximum and minimum parameter values for plot scales.

    Parameters
    ----------
    parameter : str
        centerlines, erosion, flow, geometry, static, and precipitation parameters.
    parameter_array : np.array
        array containing parameter of 'n' capillaries [m].
    parameter_array_init: np.array
        array containing parameter of 'n' capillaries [m] at the initial conditions.
    scale_vmin : float
        minimum value of the specific parameter along the simulation time.
    scale_vmax : float
        maximum value of the specific parameter along the simulation time.
    final_time_step_index : int
        final time step. Used for parameter evolution analysis when comparing initial
        conditions against specific time steps.
    plot_size : str
        plot size (width, height).
    mayavi3d_type : str
        mayavi 3d plot type for link parameters.
    x : np.array
        (N,) array containing the node coordinates for all 'N' capillaries.
    y : np.array
        (N,) array containing the node coordinates for all 'N' capillaries.
    z : np.array
        (N,) array containing the node coordinates for all 'N' capillaries.
    X : np.array
        (N,) array containing node origins
    Y : np.array
        (N,) array containing node origins
    Z : np.array
        (N,) array containing node origins
    U : np.array
        (N,) array containing node-to-node vectors
    axis: str
        flow axis. Req: flow rate, flow speed, pressures and pressure grads.
    out_folder : str
        folder containing '/gmm_case_n/'.
    show : bool
        show plot instead of saving to file.
    scale_plot : str
        scale type of 3d plots.
    ratio : bool
        plot ratio during t=0 and t=time_step_index.
    video_time_step_index : str
        time step for frame creation.

    Returns
    -------
    3d plots : parameter_plot_00000n.png
        3d plot of nodal or link parameters at 'n' time step.
    """

    # Invert origin and components if (q < 0)
    u = U / np.linalg.norm(U, axis=0)

    # Define parameter for plot
    if any([parameter == 'aspect_ratio',
            parameter == 'erosion_onset',
            parameter == 'erosion_rate',
            parameter == 'erosion_shear_stress',
            parameter == 'erosion_time_scale',
            parameter == 'flow_rate',
            parameter == 'flow_speed',
            parameter == 'link_diameters',
            parameter == 'link_length',
            parameter == 'node_diameters',
            parameter == 'precipitation_clogged_capillaries',
            parameter == 'precipitation_clogging_evolution',
            parameter == 'precipitation_rate',
            parameter == 'pressures',
            parameter == 'pressure_gradients',
            parameter == 'reactive_area',
            parameter == 'reynolds_number',
            parameter == 'void_space_volume',
            parameter == 'void_space_volume_ratio',
            parameter == 'wall_shear_stress']):

        # Reduce length of scale line
        scale = info_parameters[parameter]['3d_plot']['scale']

        # 3d plots of nodal parameters
        if any([parameter == 'node_diameters', parameter == 'pressures']):
            if parameter == 'node_diameters':
                PARAMETER_ARRAY = np.abs(parameter_array)

            elif parameter == 'pressures':
                if (ratio):
                    PARAMETER_ARRAY = np.abs(parameter_array / (parameter_array_init +
                                                                np.finfo(np.double).eps))

                else:
                    PARAMETER_ARRAY = np.abs(parameter_array)

            # POINTS visualization for nodal parameters along the flow axis
            mlab.figure(size=plot_size, bgcolor=(0.1, 0.1, 0.1))

            mlab.points3d(x, y, z, PARAMETER_ARRAY,
                          colormap='viridis', vmin=scale_vmin, vmax=scale_vmax)

            scalarbar_title = f"{info_parameters[parameter]['3d_plot']['title']} \
[{scale if scale !=1 else ''}{info_parameters[parameter]['units']}]"

            mlab.scalarbar(title=scalarbar_title, orientation='vertical')

        # 3d plots of link parameters
        elif any([parameter == 'aspect_ratio',
                  parameter == 'erosion_onset',
                  parameter == 'erosion_rate',
                  parameter == 'erosion_shear_stress',
                  parameter == 'erosion_time_scale',
                  parameter == 'flow_rate',
                  parameter == 'flow_speed',
                  parameter == 'link_diameters',
                  parameter == 'link_length',
                  parameter == 'precipitation_clogged_capillaries',
                  parameter == 'precipitation_clogging_evolution',
                  parameter == 'precipitation_rate',
                  parameter == 'pressure_gradients',
                  parameter == 'reactive_area',
                  parameter == 'reynolds_number',
                  parameter == 'void_space_volume',
                  parameter == 'void_space_volume_ratio',
                  parameter == 'wall_shear_stress']):

            # Use Quiver3d (arrows)
            if mayavi3d_type == 'quiver3d':
                if (ratio):
                    # Calculate components
                    PARAMETER_ARRAY = u * (np.finfo(np.double).eps +
                                           np.abs(parameter_array / parameter_array_init))

                    # QUIVER3D visualization for capillary link parameters along the flow axis
                    mlab.figure(size=plot_size, bgcolor=(0, 0, 0))

                    mlab.quiver3d(X, Y, Z,
                                  PARAMETER_ARRAY[0], PARAMETER_ARRAY[1], PARAMETER_ARRAY[2],
                                  line_width=2, scale_mode='scalar', scale_factor=1,
                                  colormap='viridis', scalars=np.linalg.norm(U, axis=0),
                                  vmin=scale_vmin, vmax=scale_vmax)

                    # Define vector bar title
                    mlab.vectorbar(title=f"Ratio of \
{info_parameters[parameter]['3d_plot']['title']}", orientation='vertical')

                else:
                    # Calculate components
                    PARAMETER_ARRAY = u * (np.finfo(np.double).eps + np.abs(parameter_array))

                    # QUIVER3D visualization for capillary link parameters along the flow axis
                    mlab.figure(size=plot_size, bgcolor=(0, 0, 0))

                    mlab.quiver3d(X, Y, Z,
                                  PARAMETER_ARRAY[0], PARAMETER_ARRAY[1], PARAMETER_ARRAY[2],
                                  line_width=2, scale_mode='scalar', scale_factor=1,
                                  colormap='viridis', scalars=np.linalg.norm(U, axis=0),
                                  vmin=scale_vmin, vmax=scale_vmax)

                    # Define vector bar title
                    vectorbar_title = f"{info_parameters[parameter]['3d_plot']['title']} \
[{scale if scale !=1 else ''}{info_parameters[parameter]['units']}]"
                    mlab.vectorbar(title=vectorbar_title, orientation='vertical')

            # Use points3d (spheres)
            elif mayavi3d_type == 'points3d':
                # Define transparent background
                if (transparency):
                    mlab_bgcolor = (1, 1, 1)
                    mlab_fgcolor = (0.2, 0.2, 0.2)

                # Define black backgroung
                else:
                    mlab_bgcolor = (0.1, 0.1, 0.1)
                    mlab_fgcolor = None

                # Plot parameters
                if (ratio):
                    # Calculate components
                    PARAMETER_ARRAY = (np.finfo(np.double).eps + np.abs(parameter_array /
                                       parameter_array_init))
                    mlab.figure(size=plot_size, bgcolor=mlab_bgcolor, fgcolor=mlab_fgcolor)
                    mlab.points3d(X, Y, Z, PARAMETER_ARRAY,
                                  colormap='viridis', vmin=scale_vmin, vmax=scale_vmax)

                    # Define scalar bar title
                    mlab.scalarbar(
                        title=f"Ratio of \
{info_parameters[parameter]['3d_plot']['title']}", orientation='vertical')

                else:
                    # Calculate components
                    PARAMETER_ARRAY = np.finfo(np.double).eps + np.abs(parameter_array)
                    mlab.figure(size=plot_size, bgcolor=mlab_bgcolor, fgcolor=mlab_fgcolor)
                    mlab.points3d(X, Y, Z, PARAMETER_ARRAY,
                                  colormap='viridis', vmin=scale_vmin, vmax=scale_vmax)
                    # Define scalar bar title
                    scalarbar_title = f"{info_parameters[parameter]['3d_plot']['title']} \
[{scale if scale !=1 else ''}{info_parameters[parameter]['units']}]"
                    mlab.scalarbar(title=scalarbar_title, orientation='vertical')

            mlab.outline()
            mlab.orientation_axes()

            # Show plot
            if (show):
                mlab.show()

            # Save plot
            else:
                if scale_plot == 'video':
                    snapshots_filename = f'{parameter}_{video_time_step_index:06}.png'
                    mlab.savefig(os.path.join(out_folder, f'{parameter}_video_snapshots',
                                              snapshots_filename))
                else:
                    figure_filename = f"{parameter}_plot\
{f'_{axis}' if any([axis == 'x', axis == 'y', axis == 'z']) else ''}_\
{final_time_step_index:06}.png"
                    mlab.savefig(os.path.join(out_folder, figure_filename))

    return


def plot_2d_parameter(case_filename,
                      case_folder,
                      parameter,
                      parameter_array,
                      initial_time_step_index,
                      final_time_step_index,
                      out_folder,
                      show,
                      scale_plot,
                      video_time_step_index=0):

    """
    Calculate maximum and minimum parameter values for plot scales.

    Parameters
    ----------
    parameter : str
        centerlines, erosion, flow, geometry, static, and precipitation parameters.
    parameter_array : np.array
        array containing parameter of 'n' capillaries [m].
    initial_time_step_index: int
        initial time step. Default at initial conditions.
    final_time_step_index : int
        final time step. Used for parameter evolution analysis when comparing initial
        conditions against specific time steps.
    out_folder : str
        folder containing '/gmm_case_n/'.
    case : int
        case number [-].
    show : bool
        show plot instead of saving to file.
    scale_plot : str
        scale type of 3d plots.
    video_time_step_index : str
        time step for frame creation.

    Returns
    -------
    2d plots : parameter_plot_00000n.png
        2d plot of nodal or link parameters at 'n' time step.
    """

    # Plot geometry and pore-scale process parameters
    if any([parameter == 'maximum_flow_rate',
            parameter == 'maximum_precipitation_rate',
            parameter == 'precipitation_clogging_evolution',
            parameter == 'precipitation_number_clogged_links',
            parameter == 'precipitation_number_clogged_inlet_links',
            parameter == 'permeability',
            parameter == 'porosity',
            parameter == 'total_accumulated_volume_simulation',
            parameter == 'total_accumulated_volume_time_step',
            parameter == 'void_space_volume_ratio']):

        # Reduce length of scale line
        scale = info_parameters[parameter]['2d_plot']['scale']

        if any([parameter == 'maximum_flow_rate',
                parameter == 'maximum_precipitation_rate',
                parameter == 'precipitation_clogging_evolution',
                parameter == 'precipitation_number_clogged_links',
                parameter == 'precipitation_number_clogged_inlet_links',
                parameter == 'porosity',
                parameter == 'total_accumulated_volume_simulation',
                parameter == 'total_accumulated_volume_time_step',
                parameter == 'void_space_volume_ratio']):

            # Load case parameters from case_m.json
            with open(os.path.join(case_folder, f'{case_filename}.json'),
                      mode='r') as case_file_data:
                case_data = json.load(case_file_data)

            # Import reaction time
            if (case_data['setup']['processes']['precipitation']):
                precipitation_data = case_data['simulation']['processes']['precipitation']
                reaction_time = precipitation_data['parameters']['precipitation_reaction_time']
            elif (case_data['setup']['processes']['erosion']):
                erosion_data = case_data['simulation']['processes']['erosion']
                reaction_time = erosion_data['parameters']['erosion_reaction_time']
            elif (case_data['setup']['processes']['dissolution']):
                dissolution_data = case_data['simulation']['processes']['dissolution']
                reaction_time = dissolution_data['parameters']['dissolution_reaction_time']
            elif (case_data['setup']['processes']['deposition']):
                deposition_data = case_data['simulation']['processes']['deposition']
                reaction_time = deposition_data['parameters']['deposition_reaction_time']

            if any([parameter == 'maximum_flow_rate',
                    parameter == 'porosity',
                    parameter == 'void_space_volume_ratio']):

                # Create array containing time
                time_array = np.linspace(start=initial_time_step_index,
                                         stop=reaction_time * final_time_step_index,
                                         num=len(parameter_array))

            elif any([parameter == 'maximum_precipitation_rate',
                      parameter == 'precipitation_clogging_evolution',
                      parameter == 'precipitation_number_clogged_links',
                      parameter == 'precipitation_number_clogged_inlet_links',
                      parameter == 'total_accumulated_volume_simulation',
                      parameter == 'total_accumulated_volume_time_step']):

                # Create array containing time
                time_array = np.linspace(start=initial_time_step_index + 1,
                                         stop=reaction_time * final_time_step_index,
                                         num=len(parameter_array))

            xlabel = f"{info_parameters[parameter]['2d_plot']['xlabel']} \
[{scale if scale !=1 else ''}{info_parameters[parameter]['units']}]"
            plt.xlabel(xlabel, fontsize=16)
            plt.ylabel(info_parameters[parameter]['2d_plot']['ylabel'], fontsize=16)
            plt.plot(time_array, parameter_array)

            plt.grid(True)
            plt.xlim(reaction_time * initial_time_step_index, reaction_time * final_time_step_index)
            plt.title(info_parameters[parameter]['2d_plot']['title'])
            plt.ticklabel_format(style='sci', axis='x', scilimits=(4, 4))

            # Show plot
            if (show):
                plt.show()

            # Save plots
            else:
                if scale_plot == 'video':
                    snapshots_filename = f'{parameter}_{video_time_step_index:06}.png'
                    mlab.savefig(os.path.join(out_folder, f'{parameter}_video_snapshots',
                                              snapshots_filename))
                else:
                    figure_filename = f'plot_gmm_{case_filename}', f'{parameter}.png'
                    plt.savefig(os.path.join(out_folder, figure_filename), dpi=125,
                                bbox_inches='tight')

        # Plot static results parameter
        elif parameter == 'permeability':
            time_array = np.linspace(start=initial_time_step_index + 1, stop=final_time_step_index,
                                     num=final_time_step_index + 1 - initial_time_step_index)

            xlabel = f"{info_parameters[parameter]['xlabel']} \
[{scale if scale !=1 else ''}{info_parameters[parameter]['units']}]"
            plt.xlabel(xlabel, fontsize=16)
            plt.ylabel(info_parameters[parameter]['ylabel'], fontsize=16)
            plt.plot(time_array, parameter_array)

            # Show plot
            if (show):
                plt.show()

            # Save plot
            else:
                if scale_plot == 'video':
                    snapshots_filename = f'{parameter}_{video_time_step_index:06}.png'
                    mlab.savefig(os.path.join(out_folder, f'{parameter}_video_snapshots',
                                              snapshots_filename))
                else:
                    figure_filename = f'plot_gmm_{case_filename}', f'{parameter}.png'
                    plt.savefig(os.path.join(out_folder, figure_filename), dpi=125,
                                bbox_inches='tight')

        # # Save snapshots
        # if scale_plot == 'video':
        #     snapshots_filename = f'{parameter}_{video_time_step_index:06}.png'
        #     mlab.savefig(os.path.join(out_folder, f'{parameter}_video_snapshots',
        #                               snapshots_filename))

    return


def plot_histogram(parameter,
                   parameter_array,
                   parameter_array_init,
                   final_time_step_index,
                   axis,
                   out_folder,
                   show,
                   scale_plot,
                   ratio,
                   number_of_links,
                   video_time_step_index=0):
    """
    Calculate maximum and minimum parameter values for plot scales.

    Parameters
    ----------
    parameter : str
        centerlines, erosion, flow, geometry, static, and precipitation parameters.
    parameter_array : np.array
        array containing parameter of 'n' capillaries [m].
    parameter_array_init: np.array
        array containing parameter of 'n' capillaries [m] at the initial conditions.
    final_time_step_index : int
        final time step. Used for parameter evolution analysis when comparing initial
        conditions against specific time steps.
    axis: str
        flow axis. Req: flow rate, flow speed, pressures and pressure grads.
    out_folder : str
        folder containing '/gmm_case_n/'.
    show : bool
        show plot instead of saving to file.
    scale_plot : str
        scale type of 3d plots.
    ratio : bool
        plot ratio during t=0 and t=time_step_index.
    number_of_links : int
        Number of links.
    video_time_step_index : str
        time step for frame creation.

    Returns
    -------
    histograms : parameter_dist_plot_00000n.png
        2d distribution plot of nodal or link parameters at 'n' time step.
    """

    # Define parameter for plot
    if any([parameter == 'aspect_ratio',
            parameter == 'erosion_onset',
            parameter == 'erosion_rate',
            parameter == 'erosion_shear_stress',
            parameter == 'erosion_time_scale',
            parameter == 'flow_rate',
            parameter == 'flow_speed',
            parameter == 'link_diameters',
            parameter == 'link_length',
            parameter == 'node_diameters',
            parameter == 'precipitation_clogged_capillaries',
            parameter == 'precipitation_clogging_evolution',
            parameter == 'precipitation_rate',
            parameter == 'pressures',
            parameter == 'pressure_gradients',
            parameter == 'reactive_area',
            parameter == 'reynolds_number',
            parameter == 'void_space_volume',
            parameter == 'void_space_volume_ratio',
            parameter == 'wall_shear_stress']):

        # Reduce length of scale line
        scale = info_parameters[parameter]['3d_plot']['scale']

        # 2d plot histogram
        plt.figure(dpi=125)

        if (ratio):
            plt.hist(np.abs(parameter_array / (parameter_array_init +
                                               np.finfo(np.double).eps)), bins=50)
        else:
            # Plot histogram
            plt.hist(parameter_array)

        # Define plot limits
        plt.ylim(0, number_of_links)

        # Define x label
        if (ratio):
            plt.xlabel(f"Ratio of \
{info_parameters[parameter]['histogram']['xlabel']}", fontsize=16)
        else:
            plot_xlabel = f"{info_parameters[parameter]['histogram']['xlabel']} \
[{scale if scale !=1 else ''}{info_parameters[parameter]['units']}]"
            plt.xlabel(plot_xlabel, fontsize=16)

        # Define y label
        plt.ylabel(info_parameters[parameter]['histogram']['ylabel'], fontsize=16)

        if (show):
            plt.show()
        else:
            if scale_plot == 'video':
                snapshots_filename = f'{parameter}_{video_time_step_index:06}.png'
                plt.savefig(os.path.join(out_folder, f'{parameter}_dist_video_snapshots',
                                         snapshots_filename))
            else:
                figure_filename = f"{parameter}_dist\
{f'_{axis}' if any([axis == 'x', axis == 'y', axis == 'z']) else ''}_\
{final_time_step_index:06}.png"
                plt.savefig(os.path.join(out_folder, figure_filename),
                            dpi=125, bbox_inches='tight')

    return
