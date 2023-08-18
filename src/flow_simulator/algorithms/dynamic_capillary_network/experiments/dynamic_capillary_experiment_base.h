/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/experiments/dynamic_capillary_experiment_base.h
 * \brief Contains the \c DynamicCapillaryExperimentBase base class.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 *
 * This header file contains the \c DynamicCapillaryExperimentBase base class from which all
 * \c CapillaryExperiment classes derive.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_EXPERIMENTS_DYNAMIC_CAPILLARY_EXPERIMENT_BASE_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_EXPERIMENTS_DYNAMIC_CAPILLARY_EXPERIMENT_BASE_H_

#include <memory>
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network_context.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/capillary_interface.h"

namespace simulator {

/**
 * \class DynamicCapillaryExperimentBase dynamic_capillary_experiment_base.h "src/flow_simulator/algorithms/dynamic_capillary_network/experiments/dynamic_capillary_experiment_base.h"
 * \brief Base class from which all \c DynamicCapillaryExperiment objects derive.
 *
 * This base class defines virtual methods common to all \c Experiment objects related to the
 * Dynamic Capillary Network Algorithm.
 */

class DynamicCapillaryExperimentBase {
 public:
  /// Virtual destructor
  virtual ~DynamicCapillaryExperimentBase() { }

  ///  Getter for absolute_pressure_
  double GetAbsolutePressure() const { return absolute_pressure_; }

  /// Parametrised constructor
  DynamicCapillaryExperimentBase(const double absolute_pressure,
                                 const double boundary_condition_value,
                                 std::shared_ptr<DynamicCapillaryNetworkContext> context)
    : absolute_pressure_(absolute_pressure),
      boundary_condition_value_(boundary_condition_value),
      context_(context) { }

  /// Applies boundary conditions to induce flow
  virtual arma::Col<double> ApplyBoundaryConditions(
    std::shared_ptr<CapillaryInterface> interface_) = 0;

 protected:
  /**
   * \brief Absolute fluid pressure imposed to the fluid or its boundaries.
   *
   * This variable stores value of the absolute pressure to which all fluids are subjected.
   * Pressures are expressed in \f$ [\text{Pa}] \f$.
   */
  double absolute_pressure_;

  /**
   * \brief Value of the boundary condition to be imposed (\c "pressure_gradient" or \c "flow_rate")
   *
   * This variable stores value of the boundary condition to be imposed in order to induce flow.
   * Pressure gradients are expressed in \f$ [\text{Pa} / \text{m}] \f$ and flow rates are expressed
   * in \f$ [\text{m}^3 / \text{s}] \f$.
   */
  double boundary_condition_value_;

  /**
   * \brief \c DynamicCapillaryNetworkContext instance defined by \c Initialise().
   *
   * Contains variables needed by the simulation
   */
  std::shared_ptr<DynamicCapillaryNetworkContext> context_;
};  // end of class DynamicCapillaryExperimentBase

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_EXPERIMENTS_DYNAMIC_CAPILLARY_EXPERIMENT_BASE_H_
