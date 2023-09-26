import json
import os
import subprocess
import time

import h5py
import numpy as np
from modules.deposition import deposition_geometry
from modules.dissolution import dissolution_geometry
from modules.erosion import erosion_geometry
from modules.flow_geometry_parameters import (calculate_flow_geometry_parameters,
                                              save_flow_geometry_parameters)
from modules.links_set_modification import links_set_modification
from modules.precipitation import precipitation_geometry
from modules.pressure_gradients import pressure_gradients_simple


def coupled_modified_link_radius(case_data,
                                 default_extract_axis,
                                 initial_Q,
                                 inlet_link_location,
                                 link_length,
                                 link_radius,
                                 number_of_inlet_face_adjacent_links,
                                 out_folder,
                                 pressure_gradients_prev,
                                 print_progress,
                                 prev_static_results_filename,
                                 time_step_index):

    """
    Calculates modified link radius.
    """

    # Calculate link radius variation due to EROSION
    link_radius_variation_erosion = erosion_geometry(case_data,
                                                     link_length,
                                                     link_radius,
                                                     out_folder,
                                                     pressure_gradients_prev,
                                                     print_progress,
                                                     time_step_index)

    # Calculate link radius variation due to DEPOSITION
    link_radius_variation_deposition = deposition_geometry(case_data,
                                                           inlet_link_location,
                                                           link_length,
                                                           link_radius,
                                                           number_of_inlet_face_adjacent_links,
                                                           out_folder,
                                                           pressure_gradients_prev,
                                                           print_progress,
                                                           time_step_index)

    # Calculate link radius variation due to DISSOLUTION
    link_radius_variation_dissolution = dissolution_geometry(case_data,
                                                             link_radius,
                                                             out_folder,
                                                             time_step_index)

    # Calculate link radius variation due to PRECIPITATION
    link_radius_variation_precipitation = \
        precipitation_geometry(case_data,
                               default_extract_axis,
                               initial_Q,
                               inlet_link_location,
                               link_length,
                               link_radius,
                               number_of_inlet_face_adjacent_links,
                               out_folder,
                               prev_static_results_filename,
                               print_progress,
                               time_step_index)

    # Calculate updated modified_link_radius due to combined processes
    modified_link_radius = link_radius + link_radius_variation_erosion + \
        link_radius_variation_deposition + link_radius_variation_dissolution + \
        link_radius_variation_precipitation

    # Calculate updated modified_inlet_link_radius due to combined processes
    modified_inlet_link_radius = links_set_modification(modified_link_radius,
                                                        inlet_link_location)

    return modified_link_radius, modified_inlet_link_radius


def geometry_evolution_single_process(binary_image_filename,
                                      case_data,
                                      centerlines_filename,
                                      centerlines_filename_output,
                                      minimum_link_radius,
                                      config_filename,
                                      config_filename_output,
                                      default_extract_axis,
                                      geometry_results_filename,
                                      geometry_results_folder,
                                      in_folder,
                                      initial_link_radius,
                                      initial_Q,
                                      initial_void_space_volume,
                                      inlet_link_location,
                                      link_length,
                                      number_of_inlet_face_adjacent_links,
                                      out_folder,
                                      previous_void_space_volume,
                                      print_flow_info,
                                      print_progress,
                                      save_flow_parameters,
                                      save_geometry_parameters,
                                      save_modified_centerlines,
                                      static_results_filename,
                                      start_transport_reaction,
                                      time_step_index,
                                      total_inlet_links):

    """
    Calculates  geometry evolution of a capillary network.
    """

    # Import time_steps from case_data
    time_steps = case_data['setup']['time_steps']

    # Start iterations
    time_step_index += 1

    # Calculate geometry evolution during every time step
    while time_step_index <= time_steps:
        if (print_progress):
            # Start Geometry Modification Module (GMM)
            print('FLOWSIMULATOR::GMM SAYS:')
            # Print current time step
            print(f'Time step index {time_step_index} of {time_steps}')

        # Define PREVIOUS filenames during the iteration process
        prev_geometry_results_filename = f'{geometry_results_filename}_{time_step_index - 1:06}.h5'
        prev_static_results_filename = f'{static_results_filename}_{time_step_index - 1:06}.h5'

        # Define geometry_results and static_results folder names
        geometry_results_folder = os.path.join(out_folder, 'geometry_results')
        static_results_folder = os.path.join(out_folder, 'static_results')

        # Import fields from geometry_results_00000(n-1).h5 at previous time step
        with h5py.File(os.path.join(geometry_results_folder,
                                    prev_geometry_results_filename),
                       'r') as geometry_file_data:
            # Load flow speed input file
            link_radius = np.squeeze(np.array(geometry_file_data['link_radius'], dtype=np.double))

            # Load flow speed input file
            inlet_link_radius = np.squeeze(np.array(geometry_file_data['inlet_link_radius'],
                                                    dtype=np.double))

        # Import centerlines file
        with open(centerlines_filename_output, mode='r') as centerlines_file_data:
            centerlines_data = json.load(centerlines_file_data)

        # Import data from centerlines file
        number_of_links = centerlines_data['graph']['metadata']['number_of_links']

        # Import flow experiment parameters from case.json
        case_experiment = case_data['simulation']['flow']['experiment']
        flow_axis = case_experiment['flow_axis']

        # Import fields from static_results_00000(n-1).h5 at previous time step
        with h5py.File(os.path.join(static_results_folder,
                                    prev_static_results_filename),
                       'r') as static_file_data:
            # Load flow speed input file
            V = np.squeeze(np.array(static_file_data[f'flow_speed_\
{default_extract_axis[flow_axis]}'], dtype=np.double))

        # Calculate pressure gradients at time_step_index - 1
        pressure_gradients_prev = pressure_gradients_simple(default_extract_axis[flow_axis],
                                                            centerlines_data,
                                                            static_results_folder,
                                                            prev_static_results_filename)[0]

        # Store for fluid flow analysis
        if (print_flow_info):
            # Import zero pressure gradients at time_step_index - 1
            zero_pressure_gradients_prev = \
                pressure_gradients_simple(default_extract_axis[flow_axis],
                                          centerlines_data,
                                          static_results_folder,
                                          prev_static_results_filename)[1]

            # Print links id with near zero and negative speed
            print(f'Links V = 0 : {V.size - np.count_nonzero(V)}')
            print(f'Links V < 0 : {len(list(filter(lambda x: (x < 0), V)))}')

            # Print links id with zero and reverse pressure gradients
            print(f'Links pressure gradients = 0 : \
{len(list(filter(lambda x: (x == 1), zero_pressure_gradients_prev)))}')
            print(f'Links pressure gradients < 0 : \
{len(list(filter(lambda x: (x < 0), pressure_gradients_prev)))}')

        # Calculate modified link radius
        (modified_link_radius,
         modified_inlet_link_radius) = \
            coupled_modified_link_radius(case_data,
                                         default_extract_axis,
                                         initial_Q,
                                         inlet_link_location,
                                         link_length,
                                         link_radius,
                                         number_of_inlet_face_adjacent_links,
                                         out_folder,
                                         pressure_gradients_prev,
                                         print_progress,
                                         prev_static_results_filename,
                                         time_step_index)

        # Print message
        if (print_progress):
            print(f'MAXIMUM updated link radius (in voxels): {np.max(modified_link_radius)}')
            print(f'MINIMUM updated link radius (in voxels): {np.min(modified_link_radius)}')

        # Update link_squared_radius dictionary in centerlines_data
        for link in range(number_of_links):
            edges_data = centerlines_data['graph']['edges']
            edges_data[link]['metadata']['link_squared_radius'] = (modified_link_radius[link])**2

        # Save modified centerlines.json file
        with open(centerlines_filename_output, mode='w') as updated_centerlines_data_out:
            json.dump(centerlines_data, updated_centerlines_data_out, sort_keys=False, indent=2)

        if save_modified_centerlines is True:
            # Define geometry_results folder
            modified_centerlines_folder = os.path.join(out_folder, 'modified_centerlines')

            # Create geometry_results folder
            try:
                os.mkdir(modified_centerlines_folder)
            except FileExistsError:
                pass

            # Save modified centerlines.json file at each time step
            with open(os.path.join(modified_centerlines_folder, f'{centerlines_filename}_\
{time_step_index:06}.json'), mode='w') as updated_centerlines_data_results:
                json.dump(centerlines_data, updated_centerlines_data_results, sort_keys=False,
                          indent=2)

        # Print message
        if (print_progress):
            print(f'Calculating static_results at time step {time_step_index}:')

        # Start process time
        start_flow_simulator = time.process_time()

        # Run flow simulator at current time step
        subprocess.run([os.path.join('bin', f'flow-simulator.x --run_simulation {out_folder}',
                                     f'{config_filename}.json')], shell=True)

        # Print process time
        if (print_progress):
            print(f'Flow simulation time: {time.process_time() - start_flow_simulator} s')

        # Update name of current static_results
        subprocess.run(os.path.join(f'mv {out_folder}',
                                    f'{static_results_filename}.h5 {static_results_folder}',
                                    f'{static_results_filename}_{time_step_index:06}.h5'),
                       shell=True)

        # Calculate additional parameters
        (accumulated_volume_simulation,
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
         wall_shear_stress) = calculate_flow_geometry_parameters(binary_image_filename,
                                                                 case_data,
                                                                 centerlines_data,
                                                                 default_extract_axis,
                                                                 geometry_results_filename,
                                                                 geometry_results_folder,
                                                                 in_folder,
                                                                 initial_link_radius,
                                                                 initial_void_space_volume,
                                                                 link_length,
                                                                 modified_link_radius,
                                                                 out_folder,
                                                                 previous_void_space_volume,
                                                                 print_progress,
                                                                 static_results_filename,
                                                                 time_step_index)

        # Save geometry and flow parameters to h5 files
        save_flow_geometry_parameters(accumulated_volume_simulation,
                                      accumulated_volume_time_step,
                                      aspect_ratio,
                                      inlet_link_radius,
                                      link_length,
                                      modified_link_radius,
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
                                      wall_shear_stress)

        # Define simulation stop criteria
        if any([case_data['setup']['processes']['deposition'],
               case_data['setup']['processes']['precipitation']]):

            # Import precipitation data
            if case_data['setup']['processes']['precipitation']:
                precipitation_data = case_data['simulation']['processes']['precipitation']
                clogging_data = precipitation_data['clogging']

            # Import deposition data
            elif case_data['setup']['processes']['deposition']:
                deposition_data = case_data['simulation']['processes']['deposition']
                clogging_data = deposition_data['clogging']

            # Define simulation stop criteria
            simulation_stop_criteria = clogging_data['simulation_stop_criteria']

            # Define simulation stop criteria
            if simulation_stop_criteria == "inlet_links":
                if all(x <= minimum_link_radius for x in modified_inlet_link_radius):
                    # Print message
                    print('FLOWSIMULATOR::GMM SAYS: Clogged inlet links ')
                    break

            elif simulation_stop_criteria == "capillary_network":
                # STOP CRITERIA FOR TOTAL CLOGGING
                # Stop iterations if total mineral clogging occurs
                if all(x == modified_link_radius[0] for x in modified_link_radius):
                    # Print message
                    print('FLOWSIMULATOR::GMM SAYS: Clogging due to mineral precipitation')
                    break

            # STOP CRITERIA FOR INLET FLOW RATE AND FLOW SPEED THRESHOLD
            elif any([simulation_stop_criteria == "flow_speed",
                     simulation_stop_criteria == "flow_rate"]):
                # Stop iterations if minimum maximum_flow_rate is reached
                # Import fields from static_results_00000(n-1).h5 at previous time step
                with h5py.File(os.path.join(static_results_folder,
                                            prev_static_results_filename),
                               'r') as static_file_data:
                    # Load flow rate input file
                    Q = np.squeeze(np.array(static_file_data[f'flow_rate_\
    {default_extract_axis[flow_axis]}'], dtype=np.double))

                    # Load flow rate input file
                    V = np.squeeze(np.array(static_file_data[f'flow_speed_\
    {default_extract_axis[flow_axis]}'], dtype=np.double))

                # Define precipitation_parameters
                precipitation_parameters = case_data['simulation']['processes']['precipitation']

                if np.max(Q) < precipitation_parameters['parameters']['minimum_flow_rate']:
                    # Print message
                    print('FLOWSIMULATOR::GMM SAYS: Minimum maximum_flow_rate was reached')
                    break

                elif np.max(V) < precipitation_parameters['parameters']['minimum_flow_velocity']:
                    # Print message
                    print('FLOWSIMULATOR::GMM SAYS: Minimum maximum_flow_velocity was reached')
                    break

        # Print message
        print('Capillary network geometry updated with GMM')

        # Print process time
        if (print_progress):
            print(f'TOTAL transport-reaction simulation time: \
{time.process_time() - start_transport_reaction} s')

        # Continue iterations if minimum maximum_flow_rate is reached
        time_step_index += 1

    # Remove copy of config.json
    os.remove(config_filename_output)

    # Remove copy of centerlines.json at the final simulation time step
    os.remove(centerlines_filename_output)

    return
