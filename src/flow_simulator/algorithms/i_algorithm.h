/**
 * \file src/flow_simulator/algorithms/i_algorithm.h
 * \brief Contains the \c IAlgorithm interface class.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This header file contains the \c IAlgorithm interface class from which all \c Algorithm classes
 * derive.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_I_ALGORITHM_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_I_ALGORITHM_H_

#include <armadillo>
#include <map>
#include <string>
#include <memory>
#include <utility>
#include "src/exec_manager/simulation_config.h"

namespace simulator {

/**
 * \class IAlgorithm i_algorithm.h "src/flow_simulator/algorithms/i_algorithm.h"
 * \brief Interface class from which all \c Algorithm objects derive.
 *
 * This class does nothing but define the virtual methods that all \c Algorithm objects should have.
 */

class IAlgorithm {
 public:
  /// Virtual destructor
  virtual ~IAlgorithm() { }

  /// Build and configure all resources needed for the algorithm
  virtual std::pair<int, std::string> Initialise(SimulationConfig &simulation_cfg) = 0;

  /// Builds the required geometric representation depending on the simulation
  virtual void BuildGeometricalRepresentation(void) = 0;

  /// Builds the physical equations depending on the simulation
  virtual void BuildPhysicalEquations(void) = 0;

  /// Solves the physical equations
  virtual void SolvePhysicalEquations(void) = 0;

  /// Calculate derived quantities
  virtual void CalculateDerivedQuantities(void) = 0;

  /// Saves the resulting quantities to disk
  virtual void SaveResultsToDisk(void) = 0;

  /// Close all resources and log summary informations.
  virtual void Finalize(void) = 0;
};  // end of class IAlgorithm

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_I_ALGORITHM_H_
