/**
 * \file src/flow_simulator/algorithms/static_capillary_network/experiments/static_capillary_experiment_base.h
 * \brief Contains the \c StaticCapillaryExperimentBase base class.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This header file contains the \c StaticCapillaryExperimentBase base class from which all
 * \c CapillaryExperiment classes derive.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_EXPERIMENTS_STATIC_CAPILLARY_EXPERIMENT_BASE_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_EXPERIMENTS_STATIC_CAPILLARY_EXPERIMENT_BASE_H_

#include <memory>
#include "src/flow_simulator/algorithms/static_capillary_network/static_capillary_network_context.h"

namespace simulator {

/**
 * \class StaticCapillaryExperimentBase static_capillary_experiment_base.h "src/flow_simulator/algorithms/static_capillary_network/experiments/static_capillary_experiment_base.h"
 * \brief Base class from which all \c StaticCapillaryExperiment objects derive.
 *
 * This base class defines virtual methods common to all \c Experiment objects related to the
 * Static Capillary Network Algorithm.
 */

class StaticCapillaryExperimentBase {
 public:
  /// Virtual destructor
  virtual ~StaticCapillaryExperimentBase() { }

  /// Parametrised constructor
  StaticCapillaryExperimentBase(const double absolute_pressure,
                                const double boundary_condition_value,
                                std::shared_ptr<StaticCapillaryNetworkContext> context)
    : absolute_pressure_(absolute_pressure),
      boundary_condition_value_(boundary_condition_value),
      context_(context) { }

  ///  Getter for absolute_pressure_
  double GetAbsolutePressure() const { return absolute_pressure_; }

  /// Applies boundary conditions to induce flow
  virtual arma::Col<double> ApplyBoundaryConditions(arma::SpMat<double> &matrix,
                                                    arma::Col<double> &term) = 0;

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
   * \brief The \c StaticCapillaryNetworkContext instance defined by \c Initialise().
   *
   * Contains variables needed by the simulation
   */
  std::shared_ptr<StaticCapillaryNetworkContext> context_;
};  // end of class StaticCapillaryExperimentBase

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_EXPERIMENTS_STATIC_CAPILLARY_EXPERIMENT_BASE_H_
