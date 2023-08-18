/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/models/dynamic_capillary_geometry_base.cc
 * \brief Contains the implementation of \c DynamicCapillaryGeometryBase class methods.
 *
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This source file contains the implementation of \c DynamicCapillaryGeometryBase class methods.
 */

#include "src/flow_simulator/algorithms/dynamic_capillary_network/models/dynamic_capillary_geometry_base.h"
#include <glog/logging.h>

namespace simulator {

void DynamicCapillaryGeometryBase::ReadNetworkFile(void) {
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

  // Log verbose information
  VLOG(1) << "link_length_ before conversion to SI units";
  for (IndexType i = 0U; i != context_->link_length_.n_elem; ++i) {
    VLOG(1) << i << " : " << context_->link_length_(i);
  }
  VLOG(1) << "link_squared_radius_ before conversion to SI units";
  for (IndexType i = 0U; i != context_->link_squared_radius_.n_elem; ++i) {
    VLOG(1) << i << " : " << context_->link_squared_radius_(i);
  }

  // Convert to SI units
  context_->link_length_ *= context_->voxel_size_;
  context_->link_squared_radius_ *= std::pow(context_->voxel_size_, 2.0);

  // Log verbose information
  VLOG(1) << "Voxel size: " << context_->voxel_size_;
  VLOG(1) << "link_length_ after conversion to SI units";
  for (IndexType i = 0U; i != context_->link_length_.n_elem; ++i) {
    VLOG(1) << i << " : " << context_->link_length_(i);
  }
  VLOG(1) << "link_squared_radius_ after conversion to SI units";
  for (IndexType i = 0U; i != context_->link_squared_radius_.n_elem; ++i) {
    VLOG(1) << i << " : " << context_->link_squared_radius_(i);
  }
}  // end of DynamicCapillaryGeometryBase::ReadNetworkFile() method



void DynamicCapillaryGeometryBase::LocateCapillariesConnectedToNodes(void) {
/**
 * The \c LocateCapillariesConnectedToNodes() method creates an \c arma::field of length equal to
 * the number \f$ N \f$ of nodes. Each field component is an \c arma::Col<IndexType> containing
 * the indexes of the capillaries connected to that node.
 */
  // Extracts individual \c linked_nodes_ rows into separate variables
  arma::Row<IndexType> node_indexes = context_->linked_nodes_.row(0);

  LOG(INFO) << "Number of nodes: " << context_->number_of_nodes_;
  LOG(INFO) << "Number of links: " << context_->number_of_links_;

  // Build two arma::field with vector of capillary indexes
  context_->capillaries_whose_source_is_.set_size(context_->number_of_nodes_);
  context_->capillaries_whose_target_is_.set_size(context_->number_of_nodes_);

  // Build lists of capillaries connected to each node as source/target
  for (IndexType link_index = 0U; link_index != context_->number_of_links_; ++link_index) {
    const IndexType position = 2 * link_index;
    const arma::Col<IndexType> capillary = { link_index };

    arma::uvec &source_element = context_->capillaries_whose_source_is_(node_indexes(position));
    source_element.insert_rows(source_element.n_elem, capillary);

    arma::uvec &target_element = context_->capillaries_whose_target_is_(node_indexes(position + 1));
    target_element.insert_rows(target_element.n_elem, capillary);
  }

  // Log verbose information
  VLOG(1) << "capillaries_whose_source_is_: ";
  for (IndexType i = 0U; i != context_->capillaries_whose_source_is_.n_elem; ++i) {
    VLOG(1) << i << " : " << context_->capillaries_whose_source_is_(i);
  }
  VLOG(1) << "capillaries_whose_target_is_: ";
  for (IndexType i = 0U; i != context_->capillaries_whose_target_is_.n_elem; ++i) {
    VLOG(1) << i << " : " << context_->capillaries_whose_target_is_(i);
  }
}  // end of DynamicCapillaryGeometryBase::LocateCapillariesConnectedToNodes() method



void DynamicCapillaryGeometryBase::LocateBoundaryNodes(void) {
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

  // List indexes of interior nodes
  context_->interior_nodes_ = arma::find(is_side_boundary == 0U &&
                                         is_inlet == 0U &&
                                         is_outlet == 0U);

  // Check necessary conditions
  if (context_->inlet_nodes_.n_elem == 0U)
    LOG(FATAL) << "There are no inlet nodes.";

  if (context_->outlet_nodes_.n_elem == 0U)
    LOG(FATAL) << "There are no outlet nodes.";

  // Log verbose information
  VLOG(1) << "context_->shape_(context_->flow_axis_) - 1U: "
             << context_->shape_(context_->flow_axis_) - 1U;
  VLOG(1) << "flow_axis_coords: ";
  for (IndexType i = 0U; i != flow_axis_coords.n_elem; ++i) {
    VLOG(1) << i << " : " << flow_axis_coords(i);
  }
  VLOG(1) << "is_inlet: ";
  for (IndexType i = 0U; i != is_inlet.n_elem; ++i) {
    VLOG(1) << i << " : " << is_inlet(i);
  }
  VLOG(1) << "is_outlet: ";
  for (IndexType i = 0U; i != is_outlet.n_elem; ++i) {
    VLOG(1) << i << " : " << is_outlet(i);
  }
  VLOG(1) << "is_side_boundary: ";
  for (IndexType i = 0U; i != is_side_boundary.n_elem; ++i) {
    VLOG(1) << i << " : " << is_side_boundary(i);
  }
  VLOG(1) << "inlet_nodes_: ";
  for (IndexType i = 0U; i != context_->inlet_nodes_.n_elem; ++i) {
    VLOG(1) << i << " : " << context_->inlet_nodes_(i);
  }
  VLOG(1) << "outlet_nodes_: ";
  for (IndexType i = 0U; i != context_->outlet_nodes_.n_elem; ++i) {
    VLOG(1) << i << " : " << context_->outlet_nodes_(i);
  }
  VLOG(1) << "side_boundary_nodes_: ";
  for (IndexType i = 0U; i != context_->side_boundary_nodes_.n_elem; ++i) {
    VLOG(1) << i << " : " << context_->side_boundary_nodes_(i);
  }
  VLOG(1) << "interior_nodes_: ";
  for (IndexType i = 0U; i != context_->interior_nodes_.n_elem; ++i) {
    VLOG(1) << i << " : " << context_->interior_nodes_(i);
  }
}  // end of DynamicCapillaryGeometryBase::LocateBoundaryNodes() method



void DynamicCapillaryGeometryBase::CalculateMidpointAlongFlowAxis(void) {
/**
 * The \c CalculateMidpointAlongFlowAxis() method saves the coordinate of the midpoint along flow
 * axis in the \c midpoint_ member of the \c context_ object.
 */
  const IndexType min_coordinate = context_->ctrl_voxels_.col(context_->flow_axis_).min();
  context_->midpoint_ = static_cast<double>(2 * min_coordinate + GetLengthAlongFlowAxis()) / 2.0;

  // Log debugging information
  LOG(INFO) << "Midpoint: " << context_->midpoint_;
}  // end of DynamicCapillaryGeometryBase::CalculateMidpointAlongFlowAxis() method

}  // namespace simulator
