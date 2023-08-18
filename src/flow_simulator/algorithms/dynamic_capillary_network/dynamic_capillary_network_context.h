/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network_context.h
 * \brief Contains all variables used by both Geometry and CapillaryInterface
 *
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_DYNAMIC_CAPILLARY_NETWORK_CONTEXT_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_DYNAMIC_CAPILLARY_NETWORK_CONTEXT_H_

#include <armadillo>
#include <vector>
#include <unordered_set>
#include <unordered_map>

namespace simulator {

using IndexType = arma::uword;

/**
 * \class DynamicCapillaryNetworkContext dynamic_capillary_network_context.h "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network_context.h"
 * \brief Class that carries the context necessary for different classes in the Dynamic Capillary Network
 *
 * The Dynamic Capillary Network Context class carries context shared among the classes in the Dynamic Capillary Network
 */

class DynamicCapillaryNetworkContext {
 public:
  /// Parametrised constructor
  DynamicCapillaryNetworkContext(const double voxel_size,
                                 const IndexType flow_axis,
                                 const arma::Col<IndexType>::fixed<3> shape,
                                 const IndexType boundary_thickness)
    : voxel_size_(voxel_size),
      flow_axis_(flow_axis),
      shape_(shape),
      boundary_thickness_(boundary_thickness) { }

  /// Getter for the list of capillaries connected to a given \c node_index
  const arma::Col<IndexType> GetCapillariesConnectedToNode(const IndexType node_index) const {
    return arma::join_cols(capillaries_whose_source_is_(node_index),
                           capillaries_whose_target_is_(node_index));
  }

  /**
   * \brief Indexes of pressure nodes at the side boundary
   *
   * This vector contains the indexes of the pressure nodes located at the (4) edges parallel to the
   * \c flow_axis_ direction.
   */
  arma::Col<IndexType> side_boundary_nodes_;

  /**
   * \brief Indexes of pressure nodes at the outlet flow boundary
   *
   * This vector contains the indexes of the pressure nodes located at the outlet along the
   * \c flow_axis_ direction.
   */
  arma::Col<IndexType> outlet_nodes_;

  /**
   * \brief Indexes of pressure nodes at the inlet flow boundary
   *
   * This vector contains the indexes of the pressure nodes located at the inlet along the
   * \c flow_axis_ direction.
   */
  arma::Col<IndexType> inlet_nodes_;

  /**
   * \brief Indexes of pressure nodes in the interior of the sample (not boundary)
   *
   * This vector contains the indexes of all pressure nodes not located at any of the faces
   * of the sample.
   */
  arma::Col<IndexType> interior_nodes_;

  /**
   * \brief Location and squared-radius (non-dimensional) of each centerline voxel
   *
   * This matrix has 4 columns \f$ (x_i, y_i, z_i, R_i^2) \f$ and \f$ N \f$ (number of centerline
   * voxels) rows. The distance is measured in \c [voxel] units.
   */
  arma::Mat<double> ctrl_voxels_;

  /**
   * \brief Length of links between neighbouring nodes (in meters \c [m])
   *
   * This vector with \f$ M \f$ elements, where \f$ M \f$ is the number of links between nodes,
   * contains the links' length. The length \f$ L_j \f$ is calculated as the straight line distance
   * between connected nodes.
   */
  arma::Col<double> link_length_;

  /**
   * \brief Midpoint of the sample relative length along \c flow_axis_
   *
   * The midpoint is calculated according to
   *
   * \f[
   *    M_u = \min(u) + \frac{\max(u) - \min(u)}{2} = \frac{2 \min(u) + L_u}{2}
   * \f]
   *
   * where \f$ u \f$ is the coordinate along \c flow_axis_ and \f$ L_u \f$ is the length of flow_axis_.
   */
  double midpoint_;

  /**
   * \brief Spatial resolution of the tomographic image (in meters \c [m])
   *
   * This variable stores the size associated with each image voxel, i.e. its spatial resolution.
   */
  double voxel_size_;

  /**
   * \brief Axis (\c x, \c y or \c z) along which fluid flow is simulated
   *
   * This integer number represents the direction along which the fluid flow simulation will be
   * performed, according to the convention \f$ (x = 0, y = 1, z = 2) \f$.
   */
  IndexType flow_axis_;

  /**
   * \brief Squared-radius of links between neighbouring nodes (in meters \c [m])
   *
   * This vector with \f$ M \f$ elements, where \f$ M \f$ is the number of links between nodes,
   * contains the squared-radius \f$ R_j^2 \f$ of the links.
   */
  arma::Col<double> link_squared_radius_;

  /**
   * \brief Voxel dimensions (x,y,z) of the Digital Rock sample
   *
   * This integer array contains the values that define the sample dimensions (x,y,z).
   */
  arma::Col<IndexType>::fixed<3> shape_;

  /**
   * \brief Pairs of indexes of nodes and links
   *
   * This \f$ 2 \times 2M \f$ dense matrix, where \f$ M \f$ is the number of links between nodes,
   * contains the indexes of nodes and associated links in each column. The first (second) row
   * contains the node (link) index.
   */
  arma::Mat<IndexType> linked_nodes_;

  /**
   * \brief Field with index of capillaries whose source is a given node
   *
   * This \c arma::field contains one \c arma::Col<IndexType> for each one of the \f$ N \f$ nodes,
   * listing the indexes of capillaries that have that node as a \c source.
   */
  arma::field<arma::Col<IndexType>> capillaries_whose_source_is_;

  /**
   * \brief Field with index of capillaries whose target is a given node
   *
   * This \c arma::field contains one \c arma::Col<IndexType> for each one of the \f$ N \f$ nodes,
   * listing the indexes of capillaries that have that node as a \c target.
   */
  arma::field<arma::Col<IndexType>> capillaries_whose_target_is_;

  /**
   * \brief Directions of links between neighbouring nodes
   *
   * This vector with \f$ 2M \f$ elements, where \f$ M \f$ is the number of links between nodes,
   * contains the directions \f$ \ell = \pm 1 \f$ of links connecting neighbouring nodes.
   * For instance, \c link_direction_(i) is the direction of the link \c linked_nodes_(1,i) with
   * respect to the node \c linked_nodes_(0,i).
   */
  arma::Col<double> link_direction_;

  /**
   * \brief Number of links of the capillary network
   *
   * Stores the number of links, or capillaries, in the dynamic capillary network.
   */
  IndexType number_of_links_;

  /**
   * \brief Number of nodes of the capillary network
   *
   * Stores the number of nodes in the dynamic capillary network.
   */
  IndexType number_of_nodes_;

  /**
   * \brief Thickness of the boundary region (in voxels)
   *
   * This number determines how distant a node can be from the boundary layer but still be
   * considered a boundary node.
   */
  IndexType boundary_thickness_;
};  // end of class DynamicCapillaryNetworkContext

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_DYNAMIC_CAPILLARY_NETWORK_CONTEXT_H_
