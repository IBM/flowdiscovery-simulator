/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/models/dynamic_capillary_physics_base.cc
 * \brief Contains the implementation of \c DynamicCapillaryPhysicsBase class methods.
 *
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Adolfo Emmanuel Correa López \<adolfo-correa@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 *
 * This source file contains the implementation of \c DynamicCapillaryPhysicsBase class methods.
 */

#include "src/flow_simulator/algorithms/dynamic_capillary_network/models/dynamic_capillary_physics_base.h"

namespace simulator {

arma::Col<double> DynamicCapillaryPhysicsBase::CalculateFlowSpeed(
  const arma::Col<double> &dx) const {
/**
 * Calculates the flow speeds \f$ V_j \f$ in each capillary \f$ j \f$ by taking the time derivative
 * of the fractional interface position \f$ x_j(t) \f$ and multiplying by the \c link_length_.
 *
 * \f[
 *    V_j = L_j \frac{\text{d}}{\text{d} t} x_j(t)
 * \f]
 *
 * \param[in]  dx          Vector containing the derivative of interface position in each capillary.
 * \retval     flow_speed  Vector containing the flow speed in each capillary
 */
  return dx.head(context_->number_of_links_) % context_->link_length_;
}  // end of DynamicCapillaryPhysicsBase::CalculateFlowSpeed() method



arma::Col<double> DynamicCapillaryPhysicsBase::CalculateFluidSaturation(
  const arma::Mat<double> &capillary_saturation) const {
/**
 * Calculates the overall saturations of the capillary network by calculating the weighted
 * average of the capillary saturations.
 *
 * \f[
 *    S_{\alpha} = \frac{\sum\limits_{j} \lambda_{j} S_{\alpha}^{j}}{\sum\limits_{j} \lambda_{j}} ,
 * \f]
 *
 * where, \f$ \lambda_{j} = \pi R_j^2 \cdot L_j \f$ is the volume of capillary \f$ j \f$.
 *
 * \param[in] capillary_saturation  Matrix containing the phase saturations in each capillary.
 * \retval    fluid_saturation      Column vector containing the overall saturations of each fluid.
 */
  // Omitting arma::datum::pi common factor
  const arma::vec capillary_volume = context_->link_squared_radius_ % context_->link_length_;
  const arma::vec volume_weight = capillary_volume / arma::accu(capillary_volume);

  // Calculate volume-weighted average saturations
  return arma::trans(arma::sum(capillary_saturation.each_col() % volume_weight, 0));
}  // end of DynamicCapillaryPhysicsBase::CalculateFluidSaturation() method

arma::Col<double> DynamicCapillaryPhysicsBase::CalculateFlowRate(
  const arma::Col<double> &flow_speed) const {
/**
 * Calculates the flow rate \f$ Q_{j} \f$ in capillary \f$ j \f$
 * by multiplying the flow speed inside each capillary by its cross-section.
 *
 * \f[
 *    Q_{j} = V_j \left( \pi R_j^2 \right) \, ,
 * \f]
 *
 * where, for a given capillary \f$ j \f$, \f$ V_j \f$ is the flow speed and \f$ \pi R_j^2 \f$
 * is the cross-sectional area.
 *
 * \param[in] flow_speed  Vector containing the flow speed in each capillary.
 * \retval    flow_rate   Column vector containing the flow rates in each capillary.
 */

  return arma::datum::pi * context_->link_squared_radius_ % flow_speed;
}  // end of DynamicCapillaryPhysicsBase::CalculateFlowRate() method


double DynamicCapillaryPhysicsBase::CalculateFlowRateAtOutlet(
  const arma::Col<double> &flow_rate) const {
/**
 * Calculates the outlet flow rate by summing the local flow rates at all capillaries
 * connected to outlet nodes
 *
 * \f[
 *    Q_{\text{out}} = \sum\limits_{j \in \text{outlet}} Q_{j} \, .
 * \f]
 *
 * \param[in]  flow_rate   Column vector containing the flow rates in each capillary.
 * \retval     outlet_flow Scalar that represents the outlet flow rate.
 */
  // Initialise variable with node and link indexes
  const arma::urowvec node_index = context_->linked_nodes_.row(0);
  const arma::urowvec link_index = context_->linked_nodes_.row(1);

  // Create mask indicating with nodes are at the outlets
  arma::urowvec is_outlet_node(node_index.n_elem, arma::fill::zeros);
  for (const auto &node : context_->outlet_nodes_) {
    is_outlet_node = is_outlet_node || (node == node_index);
  }

  // Identify all capillaries connected to outlet nodes
  const arma::uvec outlet_node_indexes = arma::find(is_outlet_node);
  const arma::uvec outlet_link_indexes = link_index.elem(outlet_node_indexes);

  // Filter elements related to outlet capillaries
  const arma::vec outlet_link_direction = context_->link_direction_.elem(outlet_node_indexes);
  const arma::vec outlet_flow_rate = flow_rate(outlet_link_indexes);

  return arma::as_scalar(arma::sum(-outlet_link_direction % outlet_flow_rate, 0));
}  // end of DynamicCapillaryPhysicsBase::CalculateFlowRateAtOutlet() method




arma::Row<double> DynamicCapillaryPhysicsBase::CalculateFlowRateAtInlet(
  const arma::Mat<double> &flow_rate) const {
/**
 * Calculates the inlet flow rate for each fluid phase by summing the local flow rates at all
 * capillaries connected to inlet nodes
 *
 * \f[
 *    Q_{\alpha}^{\text{in}} = \sum\limits_{j \in \text{inlet}} Q_{\alpha}^{j} \, .
 * \f]
 *
 * \param[in]  flow_rate  Matrix with 2 columns containing the phase flow rates in each capillary.
 * \retval     inlet_flow Row vector with 2 values representing the inlet flow rate of each phase.
 */
  // Initialise variable with node and link indexes
  const arma::urowvec node_index = context_->linked_nodes_.row(0);
  const arma::urowvec link_index = context_->linked_nodes_.row(1);

  // Create mask indicating with nodes are at the inlets
  arma::urowvec is_inlet_node(node_index.n_elem, arma::fill::zeros);
  for (const auto &node : context_->inlet_nodes_) {
    is_inlet_node = is_inlet_node || (node == node_index);
  }

  // Identify all capillaries connected to inlet nodes
  const arma::uvec inlet_node_indexes = arma::find(is_inlet_node);
  const arma::uvec inlet_link_indexes = link_index.elem(inlet_node_indexes);

  // Filter elements related to inlet capillaries
  const arma::vec inlet_link_direction = context_->link_direction_.elem(inlet_node_indexes);
  const arma::mat inlet_flow_rate = flow_rate.rows(inlet_link_indexes);

  return arma::sum(inlet_link_direction % inlet_flow_rate.each_col(), 0);
}  // end of DynamicCapillaryPhysicsBase::CalculateFlowRateAtInlet() method



arma::Row<double> DynamicCapillaryPhysicsBase::CalculateBulkFlowRate(
  const arma::Row<double> &inlet_flow,
  const arma::Row<double> &saturation) const {
/**
 * Calculates the bulk flow rate by multiplying the overall fluid saturation by the total inlet
 * flow rate.
 *
 * \f[
 *    Q_{\alpha}^{\text{bulk}} = Q_{\text{total}} S_{\alpha} \, .
 * \f]
 *
 * \param[in] inlet_flow  Inlet flow rate of each phase.
 * \param[in] saturation  System saturations of each fluid.
 * \retval    bulk_flow   Bulk flow rate of each phase.
 */
  // Multiply total inlet flow by phase saturation
  return arma::accu(inlet_flow) * saturation;
}  // end of DynamicCapillaryPhysicsBase::CalculateBulkFlowRate() method

}  // namespace simulator
