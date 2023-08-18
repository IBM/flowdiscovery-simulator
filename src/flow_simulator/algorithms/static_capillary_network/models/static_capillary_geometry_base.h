/**
 * \file src/flow_simulator/algorithms/static_capillary_network/models/static_capillary_geometry_base.h
 * \brief Contains the \c StaticCapillaryGeometryBase base class.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This header file contains the \c StaticCapillaryGeometryBase base class from which all
 * \c CapillaryGeometry objects derive.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_MODELS_STATIC_CAPILLARY_GEOMETRY_BASE_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_MODELS_STATIC_CAPILLARY_GEOMETRY_BASE_H_

#include <armadillo>
#include <string>
#include <memory>
#include "src/flow_simulator/i_geometry.h"
#include "src/flow_simulator/algorithms/static_capillary_network/static_capillary_network_context.h"

namespace simulator {

/**
 * \class StaticCapillaryGeometryBase static_capillary_geometry_base.h \
 * "src/flow_simulator/algorithms/static_capillary_network/models/static_capillary_geometry_base.h"
 * \brief Base class from which all \c CapillaryGeometry objects derive.
 *
 * The capillary geometry consists of a set of "pressure nodes" and "flow links" connecting those
 * nodes. Pressure values \f$ P_i \f$ are defined at the pressure nodes and flow rates \f$ Q_j \f$
 * are defined at the flow links connecting nodes \f$ i \f$ and \f$ i' \f$. This class defines
 * common methods and members that all derived \c StaticCapillaryGeometry objects inherit.
 */

class StaticCapillaryGeometryBase : public IGeometry {
 public:
  /// Virtual destructor
  virtual ~StaticCapillaryGeometryBase() { }

  /// Parametrised constructor
  StaticCapillaryGeometryBase(const std::string &folder,
                              std::shared_ptr<StaticCapillaryNetworkContext> context)
    : IGeometry(folder), context_(context) { }

  /// Reads capillary network data from \c centerlines.json and saves to \c context_ members
  void ReadNetworkFile(void);

  /// Builds sparse matrix representing connectivity between nodes
  void BuildConnectivityMatrix(void);

  /// Builds diagonal matrix representing link geometrical factors
  void BuildGeometryMatrix(const double viscosity);

  /// Locates pressure nodes at the flow and no-flow boundaries
  void LocateBoundaryNodes(void);

  /// Getter for \c link_geometry_matrix_
  const arma::SpMat<double> &GetGeometryMatrix(void) const { return link_geometry_matrix_; }

  /// Getter for \c connectivity_matrix_
  const arma::SpMat<double> &GetConnectivityMatrix(void) const { return connectivity_matrix_; }

  /// Getter for the number of pressure nodes \f$ N \f$
  IndexType GetNumberOfNodes(void) const { return context_->number_of_nodes_; }

  /// Getter for link_length
  const arma::Col<double> &GetLinkLength(void) { return context_->link_length_;}

  /// Getter for the sample length along \c flow_axis_ (in \c [voxel] units)
  IndexType GetLengthAlongFlowAxis(void) const {
    return (context_->ctrl_voxels_.col(context_->flow_axis_).max() -
            context_->ctrl_voxels_.col(context_->flow_axis_).min());
  }

  /// Getter for the midpoint along \c flow_axis_
  double CalculateMidpointAlongFlowAxis(void) const {
    return static_cast<double>(2 * context_->ctrl_voxels_.col(context_->flow_axis_).min() +
                               GetLengthAlongFlowAxis()) / 2.0;
  }

 protected:
  /**
   * \brief \c StaticCapillaryNetworkContext instance defined by \c Initialise().
   *
   * Contains variables needed by the simulation
   */
  std::shared_ptr<StaticCapillaryNetworkContext> context_;

  /**
   * \brief Connectivity matrix \f$ C \f$
   *
   * This sparse matrix of size \f$ N \times M \f$, where \f$ N \f$ is the number of nodes and
   * \f$ M \f$ is the number of links between nodes, consists of \f$ \ell = \pm 1 \f$ representing
   * the direction of the link that connects two given nodes under 26-neighbourhood scheme.
   */
  arma::SpMat<double> connectivity_matrix_;

  /**
   * \brief Link geometry matrix \f$ G \f$ in \f$ [ \text{m}^3 / (\text{Pa} \cdot \text{s}) ] \f$.
   *
   * This diagonal matrix of size \f$ M \times M \f$, where \f$ M \f$ is the number of links between
   * nodes, consists of diagonal terms \f$ G_j \f$ representing the length, radius and orientation
   * of links.
   */
  arma::SpMat<double> link_geometry_matrix_;

  /// Calculates the geometrical factors of links
  arma::vec CalculateGeometricalFactor(const arma::vec &R2, const arma::vec &L) const;
};  // end of class StaticCapillaryGeometryBase

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_MODELS_STATIC_CAPILLARY_GEOMETRY_BASE_H_
