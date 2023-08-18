/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/experiments/two_phase_pressure_gradient_open.cc
 * \brief Contains the implementation of \c TwoPhasePressureGradientOpen class methods.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 *
 * This source file contains the implementation of \c TwoPhasePressureGradientOpen class methods.
 */

#include "src/flow_simulator/algorithms/dynamic_capillary_network/experiments/two_phase_pressure_gradient_open.h"
#include <glog/logging.h>
#include <algorithm>

namespace simulator {

arma::Col<double> TwoPhasePressureGradientOpen::CalculateLinearPressure(
  const arma::Col<double> &u_coords,
  const double midpoint,
  const double voxel_size) const {
/**
 * The \c TwoPhasePressureGradientOpen::CalculateLinearPressure() method returns a pressure value as
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
  arma::vec relative_position = midpoint - u_coords;

  LOG(INFO) << "absolute_pressure_[" << absolute_pressure_ << "].";
  LOG(INFO) << "boundary_condition_value_[" << boundary_condition_value_ << "].";
  LOG(INFO) << "voxel_size[" << voxel_size << "].";

  // Log verbose information
  VLOG(1) << "relative_position: ";
  for (IndexType i = 0U; i != relative_position.n_elem; ++i) {
    VLOG(1) << i << " : " << relative_position(i);
  }

  return absolute_pressure_ + boundary_condition_value_ * relative_position * voxel_size;
}  // end of TwoPhasePressureGradientOpen::CalculateLinearPressure() method



arma::Col<double> TwoPhasePressureGradientOpen::ApplyBoundaryConditions(
  std::shared_ptr<CapillaryInterface> interface_) {
/**
 * The \c TwoPhasePressureGradientOpen::ApplyBoundaryConditions() calculates a linear pressure
 * gradient and applies the pressure in the inlet and outlet nodes of the sample. In the side
 * boundaries the absolute pressure is applied. Then it injects the invading fluid in the
 * capillaries connected to inlet nodes.
 *
 * \param[in] interface_  Shared pointer to the interface information of the capillary network.
 * \retval    pressures   arma::vec with the pressure in all nodes of the  capillary network.
 */
  // Initialise pressure vector to be partially filled
  arma::vec pressures(context_->number_of_nodes_, arma::fill::zeros);

  // Extract node coordinates
  arma::vec u_coords = context_->ctrl_voxels_.col(context_->flow_axis_);
  arma::uvec flow_boundary_nodes = arma::sort(arma::join_cols(context_->inlet_nodes_,
                                                              context_->outlet_nodes_));

  // Set the pressure in the side boundary as the absolute pressure.
  pressures.elem(context_->side_boundary_nodes_).fill(absolute_pressure_);

  // Initialise and populate right hand side vector with pressure boundary conditions
  arma::vec u_coords_boundaries = u_coords.elem(flow_boundary_nodes);
  pressures.elem(flow_boundary_nodes) = CalculateLinearPressure(u_coords_boundaries,
                                                                context_->midpoint_,
                                                                context_->voxel_size_);

  // inject fluid in inlet nodes
  interface_->InjectFluidOnInletNodes();

  // Log verbose information
  VLOG(1) << "side_boundary_nodes_: ";
  for (IndexType i = 0U; i != context_->side_boundary_nodes_.n_elem; ++i) {
    VLOG(1) << i << " : " << context_->side_boundary_nodes_(i);
  }
  VLOG(1) << "interior_nodes_: ";
  for (IndexType i = 0U; i != context_->interior_nodes_.n_elem; ++i) {
    VLOG(1) << i << " : " << context_->interior_nodes_(i);
  }
  VLOG(1) << "u_coords_boundaries: ";
  for (IndexType i = 0U; i != u_coords_boundaries.n_elem; ++i) {
    VLOG(1) << i << " : " << u_coords_boundaries(i);
  }

  return pressures;
}  // end of TwoPhasePressureGradientOpen::ApplyBoundaryConditions() method

}  // namespace simulator
