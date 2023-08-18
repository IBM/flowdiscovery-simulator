/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/models/dynamic_capillary_geometry_base.h
 * \brief Contains the \c DynamicCapillaryGeometryBase base class.
 *
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This header file contains the \c DynamicCapillaryGeometryBase base class from which all
 * \c CapillaryGeometry objects derive.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_DYNAMIC_CAPILLARY_GEOMETRY_BASE_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_DYNAMIC_CAPILLARY_GEOMETRY_BASE_H_

#include <armadillo>
#include <string>
#include <memory>
#include "src/flow_simulator/i_geometry.h"
#include "src/flow_simulator/algorithms/network_reader.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network_context.h"

namespace simulator {

/**
 * \class DynamicCapillaryGeometryBase dynamic_capillary_geometry_base.h \
 * "src/flow_simulator/algorithms/dynamic_capillary_network/models/dynamic_capillary_geometry_base.h"
 * \brief Base class from which all \c CapillaryGeometry objects derive.
 *
 * The capillary geometry consists of a set of "pressure nodes" and "flow links" connecting those
 * nodes. Pressure values \f$ P_i \f$ are defined at the pressure nodes and flow rates \f$ Q_j \f$
 * are defined at the flow links connecting nodes \f$ i \f$ and \f$ i' \f$. This class defines
 * common methods and members that all derived \c DynamicCapillaryGeometry objects inherit.
 */

class DynamicCapillaryGeometryBase : public IGeometry {
 public:
  /// Virtual destructor
  virtual ~DynamicCapillaryGeometryBase() { }

  /// Parametrised constructor
  DynamicCapillaryGeometryBase(const std::string &folder,
                              std::shared_ptr<DynamicCapillaryNetworkContext> context)
    : IGeometry(folder), context_(context) { }

  /// Reads capillary network data from \c centerlines.json and saves to \c context_ members
  void ReadNetworkFile(void);

  /// Identifies the list of capillaries connected to each node
  void LocateCapillariesConnectedToNodes(void);

  /// Locates pressure nodes at the flow and no-flow boundaries
  void LocateBoundaryNodes(void);

  /// Getter for the list of capillaries that leaves a given \c node_index
  const arma::Col<IndexType> &GetCapillariesWhoseSourceIs(const IndexType node_index) const {
    return context_->capillaries_whose_source_is_(node_index);
  }

  /// Getter for the list of capillaries that arrives at a given \c node_index
  const arma::Col<IndexType> &GetCapillariesWhoseTargetIs(const IndexType node_index) const {
    return context_->capillaries_whose_target_is_(node_index);
  }

  /// Getter for the sample length along \c flow_axis_ (in \c [voxel] units)
  IndexType GetLengthAlongFlowAxis(void) const {
    return (context_->ctrl_voxels_.col(context_->flow_axis_).max() -
            context_->ctrl_voxels_.col(context_->flow_axis_).min());
  }

  /// Getter for the midpoint along \c flow_axis_
  void CalculateMidpointAlongFlowAxis(void);

 protected:
  /**
   * \brief \c DynamicCapillaryNetworkContext instance defined by \c
   * Initialise().
   *
   * Contains variables needed by the simulation
   */
  std::shared_ptr<DynamicCapillaryNetworkContext> context_;
};  // end of class DynamicCapillaryGeometryBase

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_DYNAMIC_CAPILLARY_GEOMETRY_BASE_H_
