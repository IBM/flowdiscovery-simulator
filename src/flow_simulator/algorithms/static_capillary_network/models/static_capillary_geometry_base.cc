/**
 * \file src/flow_simulator/algorithms/static_capillary_network/models/static_capillary_geometry_base.cc
 * \brief Contains the implementation of \c StaticCapillaryGeometryBase class methods.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This source file contains the implementation of \c StaticCapillaryGeometryBase class methods.
 */

#include <cassert>
#include <algorithm>
#include "src/flow_simulator/algorithms/network_reader.h"
#include "src/flow_simulator/algorithms/static_capillary_network/models/static_capillary_geometry_base.h"

namespace simulator {

void StaticCapillaryGeometryBase::ReadNetworkFile(void) {
/**
 * The \c ReadNetworkFile() method saves the contents of the \c centerlines.json file into
 * \c DynamicCapillaryGeometryBase or \c DynamicCapillaryNetworkContext members.
 */
  // Read centreline file from disk containing the capillary network definition
  NetworkReader reader;
  NetworkInformation net_info = reader.GetNetwork(folder_ + "/centerlines.json");
  context_->link_length_ = net_info.link_length;
  context_->ctrl_voxels_ = net_info.ctrl_voxels;
  context_->linked_nodes_ = net_info.linked_nodes;
  context_->link_direction_ = net_info.link_direction;
  context_->number_of_nodes_ = net_info.number_of_nodes;
  context_->number_of_links_ = net_info.number_of_links;
  context_->link_squared_radius_ = net_info.link_squared_radius;

  // Convert to SI units
  context_->link_length_ *= context_->voxel_size_;
  context_->link_squared_radius_ *= std::pow(context_->voxel_size_, 2.0);
}  // end of StaticCapillaryGeometryBase::ReadNetworkFile() method



void StaticCapillaryGeometryBase::BuildConnectivityMatrix(void) {
/**
 * The \c BuildConnectivityMatrix() method creates a sparse matrix containing non-zero elements at
 * positions indicated by \c linked_nodes_ and whose values are given by \c link_direction_.
 *
 * This sparse matrix of size \f$ N \times M \f$, where \f$ N \f$ is the number of nodes and
 * \f$ M \f$ is the number of links between nodes, consists of \f$ \ell = \pm 1 \f$ representing
 * the direction of the link that connects two given nodes under 26-neighbourhood scheme.
 *
 * \f[
 *    C_{ij} =
 *    \begin{cases}
 *      \ell, & \text{if node } i \text{ is connected to some other node by link } j, \\
 *      0,    & \text{otherwise}
 *    \end{cases}
 * \f]
 *
 *  The value of \f$ \ell \f$ is determined by the indexes \f$ i \f$ and \f$ i' \f$ of the
 *  nodes connected by link \f$ j \f$.
 *
 * \f[
 *    \ell =
 *    \begin{cases}
 *      +1, & \text{if } i = \text{min}(i, i'), \\
 *      -1, & \text{if } i = \text{max}(i, i')
 *    \end{cases}
 * \f]
 */
  bool sort_locations = true;
  bool check_for_zeros = false;
  connectivity_matrix_ = arma::SpMat<double>(context_->linked_nodes_,
                                             context_->link_direction_,
                                             context_->ctrl_voxels_.n_rows,
                                             context_->link_direction_.n_elem / 2,
                                             sort_locations,
                                             check_for_zeros);
  }  // end of StaticCapillaryGeometryBase::BuildConnectivityMatrix() method



arma::vec StaticCapillaryGeometryBase::CalculateGeometricalFactor(const arma::vec &R2,
                                                                  const arma::vec &L) const {
/**
 * The \c CalculateGeometricalFactor() method calculates the geometrical factors \f$ k_j \f$ of all
 * links as an \c arma::vec containing
 *
 * \f[
 *    k_j = \frac{\pi R_j^4}{8 L_j} \, ,
 * \f]
 *
 * where \f$ R_j^2 \f$ comes from \c link_squared_radius_ and \f$ L_j \f$ comes from \c link_length_,
 * both in SI units.
 *
 * \param[in] R2        Vector with squared-radii of all links.
 * \param[in] L         Vector with lengths of all links.
 * \retval    k         Vector with geometrical factors of all links.
 */
  return arma::datum::pi * arma::square(R2) / (8.0 * L);
}  // end of StaticCapillaryGeometryBase::CalculateGeometricalFactor() method



void StaticCapillaryGeometryBase::BuildGeometryMatrix(const double viscosity) {
/**
 * The \c BuildGeometryMatrix() method creates a diagonal matrix containing \f$ \mu^{-1} k_j \f$
 * for all links. The geometry matrix element \f$ G_j \f$ is given by
 *
 * \f[
 *    G_j = \mu^{-1} k_j = \frac{\pi R_j^4}{8 \mu L_j} \, ,
 * \f]
 *
 * where \f$ L_j \f$ and \f$ R_j^2 \f$ are given by \c link_length_ (calculated by
 * \c CalculateLinkLength() method) and \c link_squared_radius_ (calculated by
 * \c CalculateLinkSquaredRadius() method), respectively, and \f$ \mu \f$ corresponds to the
 * \c Fluid 's dynamic viscosity.
 *
 * \param[in] viscosity  Fluid viscosity
 */
  // Initialise and populate diagonal of geometry matrix
  link_geometry_matrix_ = arma::speye(context_->link_length_.n_elem, context_->link_length_.n_elem);
  link_geometry_matrix_.diag() = CalculateGeometricalFactor(context_->link_squared_radius_,
                                                            context_->link_length_)
                               / viscosity;
}  // end of StaticCapillaryGeometryBase::BuildGeometryMatrix() method



void StaticCapillaryGeometryBase::LocateBoundaryNodes(void) {
/**
 * The \c LocateBoundaryNodes() method identifies the pressure nodes at the sample boundaries and
 * stores their indexes into Armadillo vectors. The \c inlet_nodes_, \c outlet_nodes_ and
 * \c side_boundary_nodes_ vectors store the indexes of the pressure nodes at the flow and
 * no-flow boundary edges. The inlet and outlet edges corresponds to the (2) edges
 * perpendicular to \c flow_axis_ while the side boundary corresponds to the (4) edges parallel
 * to \c flow_axis_.
 */

  // Locate nodes at the flow boundaries
  const arma::vec flow_axis_coords(context_->ctrl_voxels_.col(context_->flow_axis_));
  const arma::uvec ceil_coords(arma::conv_to<arma::uvec>::from(arma::ceil(flow_axis_coords)));
  const arma::uvec floor_coords(arma::conv_to<arma::uvec>::from(arma::floor(flow_axis_coords)));
  const arma::uvec is_inlet(ceil_coords < context_->boundary_thickness_);
  const arma::uvec is_outlet(floor_coords > context_->shape_(context_->flow_axis_) - 1U
                                           - context_->boundary_thickness_);

  // Locate the no-flow axes
  arma::uvec no_flow_axis = {context_->flow_axis_ != 0U ? context_->flow_axis_ - 1U : 2U,
                             context_->flow_axis_ != 2U ? context_->flow_axis_ + 1U : 0U};
  no_flow_axis = no_flow_axis(arma::find(context_->shape_(no_flow_axis) != 1U));


  // Locate nodes at the no-flow axes
  arma::uvec is_side_boundary(context_->ctrl_voxels_.n_rows, arma::fill::zeros);
  for (const auto &axis : no_flow_axis) {
    const arma::vec axis_coords(context_->ctrl_voxels_.col(axis));
    const arma::uvec ceil_axis_coords(arma::conv_to<arma::uvec>::from(arma::ceil(axis_coords)));
    const arma::uvec floor_axis_coords(arma::conv_to<arma::uvec>::from(arma::floor(axis_coords)));

    is_side_boundary = is_side_boundary
                    || ceil_axis_coords == 0U
                    || floor_axis_coords == context_->shape_(axis) - 1U;
  }

  // List indexes of pressure nodes at the flow and no-flow boundaries
  context_->inlet_nodes_ = arma::find(is_inlet);
  context_->outlet_nodes_ = arma::find(is_outlet);
  context_->side_boundary_nodes_ = arma::find(is_side_boundary &&
                                              is_inlet == 0U &&
                                              is_outlet == 0U);
}  // end of StaticCapillaryGeometryBase::LocateBoundaryNodes() method

}  // namespace simulator
