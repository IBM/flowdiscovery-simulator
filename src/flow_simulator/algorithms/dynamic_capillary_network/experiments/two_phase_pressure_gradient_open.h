/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/experiments/two_phase_pressure_gradient_open.h
 * \brief Contains the \c TwoPhasePressureGradientOpen class.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 *
 * This header file contains the \c TwoPhasePressureGradientOpen class derived from \c DynamicCapillaryExperimentBase.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_EXPERIMENTS_TWO_PHASE_PRESSURE_GRADIENT_OPEN_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_EXPERIMENTS_TWO_PHASE_PRESSURE_GRADIENT_OPEN_H_

#include <armadillo>
#include <memory>
#include "src/flow_simulator/algorithms/dynamic_capillary_network/experiments/dynamic_capillary_experiment_base.h"

namespace simulator {

/**
 * \class TwoPhasePressureGradientOpen two_phase_pressure_gradient_open.h "src/flow_simulator/algorithms/dynamic_capillary_network/experiments/two_phase_pressure_gradient.h"
 * \brief Two-phase flow experiment driven by a pressure gradient.
 *
 * This experiment consists of imposing a pressure gradient along \c flow_axis_ to drive two-phase
 * flow through the capillary network geometry.
 */

class TwoPhasePressureGradientOpen : public DynamicCapillaryExperimentBase {
 public:
  /// Parametrised constructor
  TwoPhasePressureGradientOpen(const double absolute_pressure,
                               const double boundary_condition_value,
                               std::shared_ptr<DynamicCapillaryNetworkContext> context)
    : DynamicCapillaryExperimentBase(absolute_pressure, boundary_condition_value, context) { }

  /// Calculate pressure according to linear function
  arma::Col<double> CalculateLinearPressure(const arma::Col<double> &u_coords,
                                            const double midpoint,
                                            const double voxel_size) const;

  /// Apply pressure boundary conditions along flow axis
  arma::Col<double> ApplyBoundaryConditions(std::shared_ptr<CapillaryInterface> interface_);
};  // end of class TwoPhasePressureGradientOpen

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_EXPERIMENTS_TWO_PHASE_PRESSURE_GRADIENT_OPEN_H_
