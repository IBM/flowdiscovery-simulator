/**
 * \file src/flow_simulator/algorithms/static_capillary_network/experiments/single_phase_pressure_gradient.cc
 * \brief Contains the implementation of \c SinglePhasePressureGradient class methods.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This source file contains the implementation of \c SinglePhasePressureGradient class methods.
 */

#include "src/flow_simulator/algorithms/static_capillary_network/experiments/single_phase_pressure_gradient.h"
#include <algorithm>

namespace simulator {

arma::Col<double> SinglePhasePressureGradient::CalculateLinearPressure(
  const arma::Col<arma::sword> &u_coords,
  const double midpoint,
  const double voxel_size) const {
/**
 * The \c SinglePhasePressureGradient::CalculateLinearPressure() method returns a pressure value as
 * a function of the position of the pressure node along the \c flow_axis_. It is used for applying
 * boundary conditions to the flow and no-flow boundary nodes.
 *
 * \f[
 *    P(u) = P_0 + \omega \left(M_u - u \right) \nu
 * \f]
 *
 * where \f$ P_0 \f$ is the absolute pressure, \f$ \omega = \frac{\partial P}{\partial u} \f$
 * is the pressure gradient along \c flow_axis_, \f$ u \f$ is the \c flow_axis_ coordinate,
 * \f$ M_u \f$ is the midpoint of the sample along \c flow_axis_ and \f$ \nu \f$ is \c voxel_size_.
 *
 * \param[in] u_coords      Spatial coordinates along \c flow_axis_ (in \c [voxel] units).
 * \param[in] midpoint      Midpoint of the sample along \c flow_axis_ (in \c [voxel] units).
 * \param[in] voxel_size    Size of the voxel (in [m]) for unit conversion.
 * \retval    pressures     Local pressure values at positions \c u_coords along \c flow_axis_.
 */
  // Convert relative position to midpoint
  arma::vec relative_position = midpoint - arma::conv_to<arma::vec>::from(u_coords);

  return absolute_pressure_ + boundary_condition_value_ * relative_position * voxel_size;
}  // end of SinglePhasePressureGradient::CalculateLinearPressure() method



arma::Col<double> SinglePhasePressureGradient::ApplyBoundaryConditions(arma::SpMat<double> &matrix,
                                                                       arma::Col<double> &term) {
/**
 * The \c SinglePhasePressureGradient::ApplyBoundaryConditions() method returns a column vector with
 * the pressure boundary conditions to be used as the right hand side vector and adapts the left
 * hand side matrix to it.
 *
 * \param[in, out] matrix      Left hand side matrix with connectivity and geometrical features.
 * \param[in]      term        First term of right hand side vector.
 * \retval         rhs_vector  Right hand side vector with boundary conditions.
 */
  // Extract node coordinates
  arma::ivec u_coords = arma::conv_to<arma::ivec>::from(
                          context_->ctrl_voxels_.col(context_->flow_axis_));
  arma::uvec flow_boundary_nodes = arma::sort(arma::join_cols(context_->inlet_nodes_,
                                                              context_->outlet_nodes_));
  arma::ivec u_coords_boundaries = u_coords.elem(flow_boundary_nodes);

  // Initialise and populate right hand side vector with pressure boundary conditions
  term.elem(flow_boundary_nodes) = CalculateLinearPressure(u_coords_boundaries,
                                                           context_->midpoint_,
                                                           context_->voxel_size_);

  // Replace matrix rows associated to boundary nodes by those of an identity matrix
  for (const auto &i : flow_boundary_nodes) matrix.row(i).zeros();
  arma::vec matrix_diagonal(matrix.diag());
  matrix_diagonal.elem(flow_boundary_nodes).ones();
  matrix.diag() = matrix_diagonal;

  return term;
}  // end of SinglePhasePressureGradient::ApplyBoundaryConditions() method

}  // namespace simulator
