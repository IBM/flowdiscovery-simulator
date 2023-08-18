/**
 * \file src/flow_simulator/algorithms/static_capillary_network/experiments/single_phase_flow_rate_closed.cc
 * \brief Contains the implementation of \c SinglePhaseFlowRateClosed class methods.
 *
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 *
 * This source file contains the implementation of \c SinglePhaseFlowRateClosed class methods.
 */

#include "src/flow_simulator/algorithms/static_capillary_network/experiments/single_phase_flow_rate_closed.h"

namespace simulator {

arma::Col<double> SinglePhaseFlowRateClosed::CalculateLocalFlowRate(
  const double total_flow_rate,
  const arma::vec &squared_radius) const {
/**
 * The \c SinglePhaseFlowRateClosed::CalculateLocalFlowRate() method returns a column vector with the
 * local flow rates \f$ q_i \f$ at each boundary pore given a total flow rate \f$ Q \f$.
 * We assume a constant fluid velocity \f$ v = \frac{Q}{\Lambda} \f$ perpendicular to the boundary, where
 * \f$ \Lambda = \sum_i \pi R_i^2 \f$ is the total pore area at the boundary and \f$ R_i^2 \f$ are the
 * squared-radii associated to boundary nodes \f$ i \f$.
 *
 * The local flow rates at each boundary pore are calculated as a fraction of the total flow rate
 * that is proportional to the fraction of the pore area it represents.
 *
 * \f[
 *    q_i = v \pi R_i^2 = Q \frac{R_i^2}{\sum_j R_j^2}
 * \f]
 *
 * \param[in] total_flow_rate Total flow rate \f$ Q \f$ to be imposed to a given boundary.
 * \param[in] squared_radius  Vector containing \f$ R^2 \f$ of boundary nodes with imposed flow rate.
 */
  return total_flow_rate * squared_radius / arma::sum(squared_radius);
}  // end of SinglePhaseFlowRateClosed::CalculateLocalFlowRate() method



arma::Col<double> SinglePhaseFlowRateClosed::ApplyBoundaryConditions(arma::SpMat<double> &matrix,
                                                                     arma::Col<double> &term) {
/**
 * The \c SinglePhaseFlowRateClosed::ApplyBoundaryConditions() method returns a column vector with
 * the flow rate boundary conditions to be used as the right hand side vector and adapts the left
 * hand side matrix to it.
 *
 * \param[in, out] matrix      Left hand side matrix with connectivity and geometrical features.
 * \param[in]      term        First term of right hand side vector.
 * \retval         rhs_vector  Right hand side vector with boundary conditions.
 */
  // Initialise vectors to be used in the calculation of boundary conditions
  arma::vec squared_radius = arma::conv_to<arma::vec>::from(context_->ctrl_voxels_.col(3));

  // Apply flow rate boundary conditions on inlet nodes
  term.elem(context_->inlet_nodes_) +=
    CalculateLocalFlowRate(boundary_condition_value_, squared_radius.elem(context_->inlet_nodes_));

  // Absolute pressure on outlet nodes
  term.elem(context_->outlet_nodes_).fill(absolute_pressure_);

  // Replace matrix rows associated to outlet nodes by those of an identity matrix
  for (const auto &i : context_->outlet_nodes_) matrix.row(i).zeros();
  arma::vec matrix_diagonal(matrix.diag());
  matrix_diagonal.elem(context_->outlet_nodes_).ones();
  matrix.diag() = matrix_diagonal;

  return term;
}  // end of SinglePhaseFlowRateClosed::ApplyBoundaryConditions() method

}  // namespace simulator
