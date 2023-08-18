/**
 * \file src/flow_simulator/algorithms/static_capillary_network/experiments/single_phase_flow_rate_open.h
 * \brief Contains the \c SinglePhaseFlowRateOpen class.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This header file contains the \c SinglePhaseFlowRateOpen class derived from \c StaticCapillaryExperimentBase.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_EXPERIMENTS_SINGLE_PHASE_FLOW_RATE_OPEN_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_EXPERIMENTS_SINGLE_PHASE_FLOW_RATE_OPEN_H_

#include <armadillo>
#include <memory>
#include "src/flow_simulator/algorithms/static_capillary_network/experiments/static_capillary_experiment_base.h"

namespace simulator {

/**
 * \class SinglePhaseFlowRateOpen single_phase_flow_rate_open.h "src/flow_simulator/algorithms/static_capillary_network/experiments/single_phase_flow_rate_open.h"
 * \brief Single-phase flow experiment driven by an imposed flow-rate at the inlets and outlets.
 *
 * This experiment consists of imposing a flow-rate along \c flow_axis_ to drive single-phase
 * flow through the capillary network geometry.
 */

class SinglePhaseFlowRateOpen : public StaticCapillaryExperimentBase {
 public:
  /// Parametrised constructor
  SinglePhaseFlowRateOpen(const double absolute_pressure,
                          const double boundary_condition_value,
                          std::shared_ptr<StaticCapillaryNetworkContext> context)
    : StaticCapillaryExperimentBase(absolute_pressure, boundary_condition_value, context) { }

  /// Calculate local flow rates from total flow rate according to pore diameters
  arma::Col<double> CalculateLocalFlowRate(const double total_flow_rate,
                                           const arma::Col<double> &squared_radius) const;

  /// Apply pressure boundary conditions along flow axis
  arma::Col<double> ApplyBoundaryConditions(arma::SpMat<double> &matrix, arma::Col<double> &term);
};  // end of class SinglePhaseFlowRateOpen

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_EXPERIMENTS_SINGLE_PHASE_FLOW_RATE_OPEN_H_
