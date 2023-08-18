/**
 * \file src/flow_simulator/algorithms/static_capillary_network/experiments/single_phase_pressure_gradient.h
 * \brief Contains the \c SinglePhasePressureGradient class.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This header file contains the \c SinglePhasePressureGradient class derived from \c StaticCapillaryExperimentBase.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_EXPERIMENTS_SINGLE_PHASE_PRESSURE_GRADIENT_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_EXPERIMENTS_SINGLE_PHASE_PRESSURE_GRADIENT_H_

#include <armadillo>
#include <memory>
#include "src/flow_simulator/algorithms/static_capillary_network/experiments/static_capillary_experiment_base.h"

namespace simulator {

/**
 * \class SinglePhasePressureGradient single_phase_pressure_gradient.h "src/flow_simulator/algorithms/static_capillary_network/experiments/single_phase_pressure_gradient.h"
 * \brief Single-phase flow experiment driven by a pressure gradient.
 *
 * This experiment consists of imposing a pressure gradient along \c flow_axis_ to drive single-phase
 * flow through the capillary network geometry.
 */

class SinglePhasePressureGradient : public StaticCapillaryExperimentBase {
 public:
  /// Parametrised constructor
  SinglePhasePressureGradient(const double absolute_pressure,
                              const double boundary_condition_value,
                              std::shared_ptr<StaticCapillaryNetworkContext> context)
    : StaticCapillaryExperimentBase(absolute_pressure, boundary_condition_value, context) { }

  /// Calculate pressure according to linear function
  arma::Col<double> CalculateLinearPressure(const arma::Col<arma::sword> &u_coords,
                                            const double midpoint,
                                            const double voxel_size) const;

  /// Apply pressure boundary conditions along flow axis
  arma::Col<double> ApplyBoundaryConditions(arma::SpMat<double> &matrix, arma::Col<double> &term);
};  // end of class SinglePhasePressureGradient

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_EXPERIMENTS_SINGLE_PHASE_PRESSURE_GRADIENT_H_
