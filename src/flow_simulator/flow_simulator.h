/**
 * \file src/flow_simulator/flow_simulator.h
 * \brief Contains the \c FlowSimulator class members and methods.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This header file contains the \c FlowSimulator class and \c FlowSimConfig struct.
 */

#ifndef SRC_FLOW_SIMULATOR_FLOW_SIMULATOR_H_
#define SRC_FLOW_SIMULATOR_FLOW_SIMULATOR_H_

#include <glog/logging.h>
#include <map>
#include <string>
#include <memory>
#include <utility>
#include "src/exec_manager/config_reader.h"
#include "src/flow_simulator/algorithms/i_algorithm.h"
#include "src/exec_manager/simulation_config.h"

using IAlgorithm = simulator::IAlgorithm;

/**
 * \class FlowSimulator flow_simulator.h "src/flow_simulator/flow_simulator.h"
 * \brief Class in charge of simulating fluid flow experiments using digital rock samples.
 *
 * This class provides methods that build and execute a simulation of a fluid flow experiment
 * using a digital rock representation as geometry. It uses a given algorithm to solve a given
 * physical model for fluid-fluid and fluid-solid interaction.
 */

class FlowSimulator {
 public:
  /// Simulates a fluid flow experiment using a certain algorithm to solve a physical equation
  void PerformFlowSimulation(SimulationConfig &simulation_cfg);

 private:
  /// Returns an object derived from IAlgorithm according to \c name parameter
  std::unique_ptr<IAlgorithm> GetAlgorithm(const std::string &name);

  /// Coordinates the execution of the simulation
  template<typename Algorithm>
  void Execute(const Algorithm &algorithm, SimulationConfig &simulation_cfg) {
  /**
   * The \c Execute() template method coordinates the execution of internal methods of the
   * \c Algorithm class and its \c Physics and \c Experiment member classes in order to run the
   * simulation.
   *
   * \tparam    Algorithm         Class type returned by \c GetAlgorithm().
   * \param[in] algorithm         Instance of \c Algorithm object returned by \c GetAlgorithm().
   * \param[in] simulation_cfg    Instance of \c SimulationConfig with all input parameters.
   */
    std::pair<int, std::string> rc;

    // Build and configure all resources needed for the algorithm
    rc = algorithm->Initialise(simulation_cfg);
    if (rc.first < 0) LOG(FATAL) << rc.second;

    // Build the required geometric representation
    algorithm->BuildGeometricalRepresentation();

    // Build physical equations related to algorithms
    algorithm->BuildPhysicalEquations();

    // Solve physical equations
    algorithm->SolvePhysicalEquations();

    // Calculate derived quantities
    algorithm->CalculateDerivedQuantities();

    // Save results to disk
    algorithm->SaveResultsToDisk();

    // Close all resources and log summary informations
    algorithm->Finalize();
  }  // end of FlowSimulator::Execute() method
};  // end of class FlowSimulator

#endif  // SRC_FLOW_SIMULATOR_FLOW_SIMULATOR_H_
