/**
 * \file src/flow_simulator/algorithms/static_capillary_network/static_capillary_network.h
 * \brief Contains the \c StaticCapillaryNetworkAlgorithm class.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This header file contains the \c StaticCapillaryNetworkAlgorithm class that derives from \c IAlgorithm.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_STATIC_CAPILLARY_NETWORK_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_STATIC_CAPILLARY_NETWORK_H_

#include <armadillo>
#include <map>
#include <string>
#include <memory>
#include <utility>
#include "src/flow_simulator/algorithms/i_algorithm.h"
#include "src/flow_simulator/algorithms/static_capillary_network/models/model_base.h"
#include "src/flow_simulator/algorithms/static_capillary_network/models/static_capillary_physics_base.h"
#include "src/flow_simulator/algorithms/static_capillary_network/models/static_capillary_geometry_base.h"
#include "src/flow_simulator/algorithms/static_capillary_network/experiments/static_capillary_experiment_base.h"
#include "src/exec_manager/simulation_config.h"

namespace simulator {

/**
 * \class StaticCapillaryNetworkAlgorithm static_capillary_network.h "src/flow_simulator/algorithms/static_capillary_network/static_capillary_network.h"
 * \brief Simulator algorithm for Capillary Network using centerline geometries
 *
 * The Static Capillary Network Algorithm performs fluid flow simulations using a set of connected
 * centerlines as representation of the porous rock geometry.
 */

class StaticCapillaryNetworkAlgorithm : public IAlgorithm {
 public:
  /// Build and configure all resources needed for the algorithm.
  std::pair<int, std::string> Initialise(SimulationConfig &simulation_cfg);

  /// Sets the \c Experiment \c protected member object
  void SetExperiment(const double &absolute_pressure,
                     const std::pair<std::string, double> &boundary_condition);

  /// Builds the required geometric representation depending on the simulation
  void BuildGeometricalRepresentation(void);

  /// Builds the physical equations of the Capillary Network Algorithm
  void BuildPhysicalEquations(void);

  /// Calculates left hand side matrix for Capillary Network equations
  void CalculateLeftHandSideMatrix(void);

  /// Calculates right hand side vector for Capillary Network equations
  void CalculateRightHandSideVector(void);

  /// Solves the physical equation of the Capillary Network Algorithm
  void SolvePhysicalEquations(void);

  /// Calculate derived quantities
  void CalculateDerivedQuantities(void);

  /// Saves the resulting quantities to disk
  void SaveResultsToDisk(void);

  /// Close all resources and log summary informations.
  void Finalize(void);

  /// Getter for \c pressures_
  const arma::Col<double> &GetPressures(void) { return pressures_; }

  /// Getter for \c flow_rate_
  const arma::Col<double> &GetFlowRate(void) { return flow_rate_; }

  /// Getter for \c flow_speed_
  const arma::Col<double> &GetFlowSpeed(void) { return flow_speed_; }

  /// Getter for \c permeability_
  double GetPermeability(void) { return permeability_; }

  /// Setter for \c pressures_
  void SetPressures(const arma::Col<double> &pressures) {  pressures_ = pressures; }

 private:
  /**
   * \brief Unique pointer to \c ModelBase instance defined by \c Initialise().
   *
   * Contains getters for \c StaticCapillaryPhysicsBase and \c StaticCapillaryGeometryBase objects.
   */
  std::unique_ptr<ModelBase> model_;

  /**
   * \brief Shared pointer to \c StaticCapillaryNetworkContext instance defined by \c
   * Initialise().
   *
   * Contains variables needed by both \c StaticCapillaryPhysicsBase and \c StaticCapillaryGeometryBase
   */
  std::shared_ptr<StaticCapillaryNetworkContext> context_;

  /**
   * \brief Unique pointer to \c StaticCapillaryExperimentBase instance defined by \c Initialise().
   *
   * This member variable stores a unique pointer to an object derived from
   * \c StaticCapillaryExperimentBase as set by the \c public method \c Initialise() according to
   * the \c name argument.
   */
  std::unique_ptr<StaticCapillaryExperimentBase> experiment_;

  /**
   * \brief Left hand side matrix \f$ A = C G C^{T} \f$ for the Capillary Network Algorithm
   *
   * This matrix is calculated from the connectivity \f$ C \f$ and geometry \f$ G \f$ matrices and
   * goes into the algorithm's matrix equation \f$ A \vec{x} = \vec{b} \f$.
   */
  arma::SpMat<double> lhs_matrix_;

  /**
   * \brief Right hand side vector \f$ \vec{b} \f$ for the Capillary Network Algorithm
   *
   * This vector is calculated from the pressure boundary conditions along the flow direction and
   * goes into the algorithm's matrix equation \f$ A \vec{x} = \vec{b} \f$.
   */
  arma::Col<double> rhs_vector_;

  /**
   * \brief This vector contains the pressures \f$ P_i \f$ at each node \f$ i \f$.
   *
   * The node pressures (in \c [Pa]) are calculated by solving the model-dependent matrix equation.
   */
  arma::Col<double> pressures_;

  /**
   * \brief This vector contains the flow rates \f$ Q_j \f$ at each capillary \f$ j \f$.
   *
   * The flow rates (in \f$ [\text{m}^3 / \text{s}] \f$) are calculated according to the chosen model.
   */
  arma::Col<double> flow_rate_;

  /**
   * \brief This vector contains the flow speeds \f$ V_j \f$ at each capillary \f$ j \f$.
   *
   * It is calculated by dividing the flow rate \f$ Q_j \f$ by the cross-sectional area
   * \f$ \pi R_{j}^2 \f$.
   */
  arma::Col<double> flow_speed_;

  /**
   * \brief Double with component of the permeability tensor corresponding to the inlet/outlet.
   *
   * The permeability components (in \f$ [\text{m}^2] \f$) are calculated from the bulk flow rate
   * (speed) and the overall pressure drop, according to Darcy's Law:
   *
   * \f[
   *    \kappa = \mu \frac{Q}{A} \frac{\Delta L}{\Delta P}
   * \f]
   *
   * Where \f$ Q = \sum\limits_{j \in outlet} Q_j \f$, \f$ A \f$ is the cross sectional area
   * perpendicular to \c flow_axis_, \f$ \Delta L \f$ is the length of the sample along \c flow_axis,
   * and \f$ \Delta P \f$ is the average difference between inlets and outlets.
   */
  double permeability_;
};  // end of class StaticCapillaryNetworkAlgorithm

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_STATIC_CAPILLARY_NETWORK_H_
