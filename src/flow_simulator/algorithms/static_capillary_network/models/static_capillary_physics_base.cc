/**
 * \file src/flow_simulator/algorithms/static_capillary_network/models/static_capillary_physics_base.cc
 * \brief Contains the implementation of \c StaticCapillaryPhysicsBase class methods.
 *
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 *
 * This source file contains the implementation of \c StaticCapillaryPhysicsBase class methods.
 */

#include <memory>
#include "src/flow_simulator/algorithms/static_capillary_network/models/static_capillary_physics_base.h"

namespace simulator {

const arma::Col<double> StaticCapillaryPhysicsBase::CalculateFlowSpeed(
    const arma::Col<double> &flow_rate) {
/**
 * Calculates the flow speeds \f$ V_j \f$ at each capillary by dividing the flow rate  \f$ Q \f$
 * by the cross sectional area  \f$ \pi R^2 \f$.
 *
 * \f[
 *    V_{j} = \frac{Q_j}{\pi R^2}
 * \f]
 *
 * \param[in]  flow_rate   Vector containing the flow rate in each capillary
 * \retval     flow_speed  Vector containing the flow speed in each capillary
 */
  return flow_rate / (arma::datum::pi * context_->link_squared_radius_);
}  // end of StaticCapillaryPhysicsBase::CalculateFlowSpeed() method



double StaticCapillaryPhysicsBase::CalculateFlowRateAtOutlet(const arma::Col<double> &flow_rate)
  const {
/**
 *  Calculates the flow sum \f$ Q_{out} = \sum\limits_{j \in outlet} Q_j \f$.
 *
 *  \param[in]  flow_rate         Vector containing the flow rate for each capillary.
 *  \retval     total_flow_rate   Total flow rate of all capillaries connected to outlet nodes.
 */
  // Initialise variable with node and link indexes
  arma::Row<IndexType> node_index = context_->linked_nodes_.row(0);
  arma::Row<IndexType> link_index = context_->linked_nodes_.row(1);

  // Create mask indicating with nodes are at the outlets
  arma::Row<IndexType> is_outlet_node(node_index.n_elem, arma::fill::zeros);
  for (const auto &node : context_->outlet_nodes_) {
    is_outlet_node = is_outlet_node || (node == node_index);
  }

  // Collect flow rates of all links connected to outlet nodes
  arma::uvec outlet_node_indexes = arma::find(is_outlet_node);
  arma::uvec outlet_links = link_index(outlet_node_indexes);
  arma::vec outlet_flow_rates = - context_->link_direction_(outlet_node_indexes)
                              % flow_rate(outlet_links);

  return arma::sum(outlet_flow_rates);
}  // end of StaticCapillaryPhysicsBase::CalculateFlowRateAtOutlet() method



double StaticCapillaryPhysicsBase::CalculatePermeability(const arma::Col<double> &flow_rate,
                                                         const arma::Col<double> &pressures) {
/**
 * Calculates the permeability \f$ \kappa \f$ of the sample along the flow axis as follows:
 *
 * \f[
 *    \kappa = \mu \frac{Q}{A} \frac{\Delta L}{\Delta P}
 * \f]
 *
 * Where \f$ Q = \sum\limits_{j \in outlet} Q_j \f$, \f$ A \f$ is the cross sectional area
 * perpendicular to \c flow_axis_, \f$ \Delta L \f$ is the length of the sample along the flow axis,
 * and \f$ \Delta P \f$ is the difference between the average inlet and outlet.
 *
 * \param[in]  flow_rate           Vector containing the flow rate for each capillary
 * \param[in]  pressures           Local pressures at positions \c u_coords along \c flow_axis_.
 * \retval     permeability        Permeability of the sample in the flow axis.
 */
  // Calculate individual terms
  double total_flow_rate = CalculateFlowRateAtOutlet(flow_rate);
  double delta_L = context_->voxel_size_
                 * static_cast<double>(context_->shape_[context_->flow_axis_]);
  double area = std::pow(context_->voxel_size_, 3.0)
              * static_cast<double>(context_->shape_[0] * context_->shape_[1] * context_->shape_[2])
              / delta_L;
  double delta_P = arma::mean(pressures(context_->inlet_nodes_))
                 - arma::mean(pressures(context_->outlet_nodes_));
  double midpoint_pressure = (arma::mean(pressures(context_->inlet_nodes_))
                           + arma::mean(pressures(context_->outlet_nodes_))) / 2.0;

  return GetFluidViscosity(midpoint_pressure) * total_flow_rate / area * delta_L / delta_P;
}  // end of StaticCapillaryPhysicsBase::CalculatePermeability() method



arma::Col<double> StaticCapillaryPhysicsBase::CalculateFlowRateBase(
  const std::shared_ptr<StaticCapillaryGeometryBase> &geometry,
  const arma::Col<double> &pressures) const {
/**
 * The \c StaticCapillaryPhysicsBase::CalculateFlowRateBase() method calculates the flow rate
 * \f$ Q_j \f$ for each capillary \f$ j \f$ as
 *
 * \f[
 *    Q_j = G_j \Delta P_j \, ,
 * \f]
 *
 * where \f$ G_j = \frac{\pi R_j^4}{8 \mu L_j} \f$ is the geometry matrix diagonal element.
 *
 * \param[in] geometry    Pointer to StaticCapillaryGeometryBase
 * \param[in] pressures   Pressure vector resulting from the solution of the matrix equation.
 * \retval    flow_rate   Vector containing the flow rate for each capillary.
 */
  // Get G matrix diagonal
  arma::vec geometry_matrix_diagonal(geometry->GetGeometryMatrix().diag());

  // Extract node index vector from linked_nodes_
  arma::Row<IndexType> node_index = context_->linked_nodes_.row(0);
  arma::uvec i_nodes = node_index(arma::regspace<arma::uvec>(0, 2, node_index.n_elem - 1));
  arma::uvec ii_nodes = node_index(arma::regspace<arma::uvec>(1, 2, node_index.n_elem - 1));

  // Calculate base flow rate
  // NOTE(rneumann): assuming link_direction_ == 1.0 for even indexes
  arma::vec flow_rate = geometry_matrix_diagonal % (pressures(i_nodes) - pressures(ii_nodes));

  return flow_rate;
}  // end of StaticCapillaryPhysicsBase::CalculateFlowRateBase() method

}  // namespace simulator
