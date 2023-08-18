/**
 * \file src/exec_manager/simulation_config.h
 * \brief Contains the \c SimulationConfig class members and methods.
 *
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This header file contains the \c SimulationConfig class.
 */

#ifndef SRC_EXEC_MANAGER_SIMULATION_CONFIG_H_
#define SRC_EXEC_MANAGER_SIMULATION_CONFIG_H_

#include <armadillo>
#include <vector>
#include <string>
#include <utility>
#include "src/exec_manager/setup_config.h"
#include "src/exec_manager/fluid_json.h"
#include "src/exec_manager/wettability_json.h"
#include "src/exec_manager/fluid_interface_json.h"
#include "src/exec_manager/algorithm_json.h"
#include "src/exec_manager/experiment_json.h"

/**
 * \class SimulationConfig simulation_config.h "src/exec_manager/simulation_config.h"
 * \brief SimulationConfig informations for the Capillary Network algorithm.
 *
 * The \c SimulationConfig class gets its information from \c ExecutionManager.
 * The \c SimulationConfig is the in-memory representation of the JSON input file.
 */
class SimulationConfig : public SetupConfig {
 public:
  /**
   * \brief List of fluids (name and properties) used in the flow simulation
   *
   * This vector contains a list of fluid descriptors, containing the fluid name and a list of
   * physical properties that are relevant for the flow simulation.
   */
  std::vector<FluidJSON> fluids_json;

  /**
   * \brief Wettability of two fluids and one rock.
   *
   * This Object variable stores the wettability properties of the interaction between the fluids
   * and the rock.
   */
  WettabilityJSON wettability_json;

  /**
   * \brief Fluid-fluid interface
   *
   * This Object variable stores the Interface properties of the interaction between the fluids.
   */
  FluidInterfaceJSON fluid_interface_json;

  /**
   * \brief Algorithm parameters
   *
   * This Object variable stores the Algorithm parameters of the simulation.
   */
  AlgorithmJSON algorithm_json;

  /**
   * \brief Experiment parameters
   *
   * This Object variable stores the Experiment parameters of the simulation.
   */
  ExperimentJSON experiment_json;
};

#endif  // SRC_EXEC_MANAGER_SIMULATION_CONFIG_H_
