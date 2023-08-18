/**
 * \file src/flow_simulator/algorithms/static_capillary_network/models/static_capillary_physics_base.h
 * \brief Contains the \c StaticCapillaryPhysicsBase base class.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This header file contains the \c StaticCapillaryPhysicsBase base class from which all
 * \c CapillaryPhysics objects derive.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_MODELS_STATIC_CAPILLARY_PHYSICS_BASE_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_MODELS_STATIC_CAPILLARY_PHYSICS_BASE_H_

#include <memory>
#include "src/flow_simulator/algorithms/static_capillary_network/models/static_capillary_geometry_base.h"

namespace simulator {

/**
 * \class StaticCapillaryPhysicsBase static_capillary_physics_base.h "src/flow_simulator/algorithms/static_capillary_network/models/static_capillary_physics_base.h"
 * \brief Base class from which all \c CapillaryPhysics objects derive.
 *
 * The capillary physics consists of a set of equations relating physical properties such as
 * flow rate and pressure drop. These equations are used to calculate geometrical parameters for the
 * capillary network algorithm.
 */

class StaticCapillaryPhysicsBase {
 public:
  /// Virtual destructor
  virtual ~StaticCapillaryPhysicsBase() { }

  /// Parametrised constructor
  explicit StaticCapillaryPhysicsBase(std::shared_ptr<StaticCapillaryNetworkContext> context)
    : context_(context) { }

  /// Calculates Right Hand Side Vector Term
  virtual arma::Col<double> CalculateRightHandSideVectorTerm(void) = 0;

  /// Calculates the flow rate at each capillary
  virtual arma::Col<double> CalculateFlowRate(const arma::Col<double> &pressures) const = 0;

  /// Calculates Flow Speeds
  const arma::Col<double> CalculateFlowSpeed(const arma::Col<double> &flow_rate);

  /// Calculates Permeability
  double CalculatePermeability(const arma::Col<double> &flow_rate,
                               const arma::Col<double> &pressures);

  /// Geter for \c dynamic_viscosity
  virtual double GetFluidViscosity(const double pressure) const = 0;

  /// Return true if the parameters fall between the designed behaviour limits, false otherwise
  virtual double IsViscosityBehaviourDesignedForInjectedFluid(const double) const = 0;

 protected:
  /// Sums all the flow rates through the outlets
  double CalculateFlowRateAtOutlet(const arma::Col<double> &flow_rate) const;

  /// Base method that considers only the \f$ G \Delta P \f$ contribution to the flow rate
  arma::Col<double> CalculateFlowRateBase(
    const std::shared_ptr<StaticCapillaryGeometryBase> &geometry,
    const arma::Col<double> &pressures) const;

  /**
   * \brief \c StaticCapillaryNetworkContext instance defined by \c Initialise().
   *
   * Contains variables needed by the simulation
   */
  std::shared_ptr<StaticCapillaryNetworkContext> context_;
};  // end of class StaticCapillaryPhysicsBase

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_MODELS_STATIC_CAPILLARY_PHYSICS_BASE_H_
