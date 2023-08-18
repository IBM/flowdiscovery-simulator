/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/models/dynamic_capillary_physics_base.h
 * \brief Contains the \c DynamicCapillaryPhysicsBase base class.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This header file contains the \c DynamicCapillaryPhysicsBase base class from which all
 * \c CapillaryPhysics objects derive.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_DYNAMIC_CAPILLARY_PHYSICS_BASE_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_DYNAMIC_CAPILLARY_PHYSICS_BASE_H_

#include <armadillo>
#include <memory>
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network_context.h"

namespace simulator {

/**
 * \class DynamicCapillaryPhysicsBase dynamic_capillary_physics_base.h "src/flow_simulator/algorithms/dynamic_capillary_network/models/dynamic_capillary_physics_base.h"
 * \brief Base class from which all \c DynamicCapillaryPhysics objects derive.
 *
 * The capillary physics consists of a set of equations relating physical properties such as
 * flow rate and pressure drop. These equations are used to calculate geometrical parameters for the
 * dynamic capillary network algorithm.
 */

class DynamicCapillaryPhysicsBase {
 public:
  /// Parametrised constructor
  explicit DynamicCapillaryPhysicsBase(std::shared_ptr<DynamicCapillaryNetworkContext> context)
    : context_(context) { }

  /// Virtual destructor
  virtual ~DynamicCapillaryPhysicsBase() { }

  /// Calculates flow speeds in each capillary
  arma::Col<double> CalculateFlowSpeed(const arma::Col<double> &dx) const;

  /// Calculates overall system saturation of each fluid
  arma::Col<double> CalculateFluidSaturation(const arma::Mat<double> &capillary_saturation) const;

  /// Calculates flow rate in each capillary
  arma::Col<double> CalculateFlowRate(const arma::Col<double> &flow_speed) const;

  /// Calculates phase permeabilities
  virtual arma::Col<double> CalculatePhasePermeability(const arma::Col<double> &flow_rate,
                                                       const arma::Col<double> &pressure) const = 0;

  /// Return true if the parameters fall between the designed behaviour limits, false otherwise
  virtual double IsViscosityBehaviourDesignedForResidentFluid(const double lower_pressure,
                                                             const double upper_pressure) const = 0;

  /// Return true if the parameters fall between the designed behaviour limits, false otherwise
  virtual double IsViscosityBehaviourDesignedForInjectedFluid(const double lower_pressure,
                                                             const double upper_pressure) const = 0;

 protected:
  /// Sums all the flow rates through the outlets
  double CalculateFlowRateAtOutlet(const arma::Col<double> &flow_rate) const;

  /// Sums all the flow rates through the inlets
  arma::Row<double> CalculateFlowRateAtInlet(const arma::Mat<double> &flow_rate) const;

  /// Calculates bulk flow rate from phase saturations
  arma::Row<double> CalculateBulkFlowRate(const arma::Row<double> &inlet_flow,
                                          const arma::Row<double> &saturation) const;

  /**
   * \brief \c DynamicCapillaryNetworkContext instance defined by \c
   * Initialise().
   *
   * Contains variables needed by the simulation
   */
  std::shared_ptr<DynamicCapillaryNetworkContext> context_;
};  // end of class DynamicCapillaryPhysicsBase

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_DYNAMIC_CAPILLARY_PHYSICS_BASE_H_
