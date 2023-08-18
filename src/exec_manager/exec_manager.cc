/**
 * \file src/exec_manager/exec_manager.cc
 * \brief Contains the \c ExecutionManager class methods.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This source file contains the methods from the \c ExecutionManager class that run the mode
 * selected by command-line arguments.
 */

#include "src/exec_manager/exec_manager.h"
#include "src/flow_simulator/flow_simulator.h"


void ExecutionManager::RunSimulation(const std::string &json_file_name) const {
/**
 * The \c RunSimulation() method calls \c FlowSimulator methods involved in performing a fluid flow
 * simulation using the geometry obtained previously by the \c DigitalRock class.
 *
 * \param[in] json_file_name  Name of JSON file containing simulation-related parameters.
 */
  // Parse JSON file into structs
  SimulationConfig simulation_cfg;
  config_reader_.PopulateSetupConfig(simulation_cfg, json_file_name);
  config_reader_.PopulateSimulationConfig(simulation_cfg, json_file_name);

  // Perform simulation steps
  FlowSimulator flow_sim;
  flow_sim.PerformFlowSimulation(simulation_cfg);
}  // end of ExecutionManager::RunSimulation method
