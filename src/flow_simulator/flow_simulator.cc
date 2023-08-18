/**
 * \file src/flow_simulator/flow_simulator.cc
 * \brief Contains the implementation of \c FlowSimulator class methods.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This source file contains the implementation of \c FlowSimulator class methods.
 * The \c PerformFlowSimulation() instantiates the appropriate \c Algorithm object, initialises its
 * \c Physics, \c Experiment and \c Geometry member objects based on the contents of
 * \c settings and calls the template method \c Execute().
 */

#include "src/flow_simulator/flow_simulator.h"
#include "src/flow_simulator/algorithms/static_capillary_network/static_capillary_network.h"
#include "src/flow_simulator/algorithms/static_capillary_network/static_capillary_network_context.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network_context.h"
#include "src/exec_manager/simulation_config.h"

using IAlgorithm = simulator::IAlgorithm;
using StaticCapillaryNetworkAlgorithm = simulator::StaticCapillaryNetworkAlgorithm;
using DynamicCapillaryNetworkAlgorithm = simulator::DynamicCapillaryNetworkAlgorithm;

void FlowSimulator::PerformFlowSimulation(SimulationConfig &simulation_cfg) {
/**
 * The \c PerformFlowSimulation() method is responsible for retrieving the user-defined execution
 * parameters, instantiating the appropriate \c Algorithm object, initialising its \c Physics,
 * \c Geometry and \c Experiment member objects and calling the \c Execute method.
 *
 * \param[in]   simulation_cfg  Struct that stores parameters from the "simulation" top-level JSON object.
 */
  // Print output messages if necessary
  std::printf("\nFLOWSIMULATOR::PERFORMFLOWSIMULATION SAYS:\n");

  // Instantiate simulator objects according user settings
  const auto algorithm = GetAlgorithm(simulation_cfg.algorithm_json.GetName());

  // Run simulation
  Execute(algorithm, simulation_cfg);
}  // end of FlowSimulator::PerformFlowSimulation() method



std::unique_ptr<IAlgorithm> FlowSimulator::GetAlgorithm(const std::string &name) {
/**
 * The \c GetAlgorithm() returns an object derived from \c IAlgorithm depending on the
 * \c name parameter.
 *
 * \param[in] name      Name of the \c IAlgorithm object to be instantiated.
 * \retval    algorithm Instance of an \c IAlgorithm object.
 */
  if (name == "static") {
    return std::make_unique<StaticCapillaryNetworkAlgorithm>();
  } else if (name == "dynamic") {
    return std::make_unique<DynamicCapillaryNetworkAlgorithm>();
  } else {
    std::fprintf(stderr, "Please select a valid flow simulation algorithm.\n");
    std::exit(-1);
  }
}  // end of FlowSimulator::GetAlgorithm() method
