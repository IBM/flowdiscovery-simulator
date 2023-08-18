/**
 * \file src/flow_simulator/algorithms/static_capillary_network/models/hydrostatic.cc
 * \brief Contains the implementation of \c HydrostaticPhysics and \c HydrostaticGeometry
 * class methods.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 *
 * This source file contains the implementation of \c HydrostaticPhysics and
 * \c HydrostaticGeometry class methods.
 */

#include "src/flow_simulator/algorithms/static_capillary_network/models/hydrostatic.h"

namespace simulator {

arma::Mat<double> HydrostaticGeometry::CalculateInternodeDistanceVector(void) const {
/**
 * The \c CalculateInternodeDistanceVector() method creates a matrix containing the distance vectors
 * \f$ \vec{d}_{j} = \vec{d}_{i \rightarrow i'} \f$ between all pairs of linked nodes expressed
 * in \c SI units.
 *
 * \f[
 *    \vec{d}_{j} = \vec{d}_{i \rightarrow i'}
 *                = \nu \left[ (x_{i'}, y_{i'}, z_{i'}) - (x_i, y_i, z_i) \right] \, ,
 * \f]
 *
 * where \f$ (x_i, y_i, z_i) \f$ are the \c voxel coordinates of node \f$ i \f$ and \f$ \nu \f$ is
 * the \c voxel_size_ length conversion factor.
 *
 * \retval   distance_vectors   Vectors connecting all pairs of linked nodes.
 */
  // Extract (x,y,z) coordinate vector from ctrl_voxels_
  arma::mat xyz_coords = arma::conv_to<arma::mat>::from(context_->ctrl_voxels_.cols(0, 2))
                       * context_->voxel_size_;

  // Extract node index vector from linked_nodes_
  arma::Row<IndexType> node_index = context_->linked_nodes_.row(0);

  // Allocate (x,y,z) matrixes with initial and final points of links
  arma::uvec i_nodes = node_index(arma::regspace<arma::uvec>(0, 2, node_index.n_elem - 1));
  arma::uvec ii_nodes = node_index(arma::regspace<arma::uvec>(1, 2, node_index.n_elem - 1));
  arma::mat distance_vectors = xyz_coords.rows(ii_nodes) - xyz_coords.rows(i_nodes);

  return distance_vectors;
}  // end of HydrostaticGeometry::CalculateInternodeDistanceVector() method



void HydrostaticGeometry::CalculateGammaVector(void) {
/**
 * The vertical direction cosine \f$ \gamma_j = \cos( \theta_j ) \f$ is the cosine of the angle
 * \f$ \theta_j \f$ between the line that defines the capillary and the vertical (z axis).
 *
 * \f[
 *    \gamma_j = \cos( \theta_j ) = \left( \frac{L_j^z}{L_j} \right) \, .
 * \f]
 *
 * \c distances contains a matrix where each row is the distance vector between linked nodes
 */
  arma::mat distances = CalculateInternodeDistanceVector();
  gamma_ = distances.col(2) / context_->link_length_;
}  // end of HydrostaticGeometry::CalculateGammaVector() method



arma::Col<double> HydrostaticPhysics::CalculateRightHandSideVectorTerm(void) {
/**
 * The \c HydrostaticPhysics::CalculateRightHandSideVectorTerm method specialises the
 * \c StaticCapillaryPhysicsBase::CalculateRightHandSideVectorTerm with model-dependent
 * calculations related to \c Hydrostatic.
 *
 * For the Hydrostatic model the boundary conditions are applied on top of the
 * \f$ C G h \f$ vector. Since \f$ C Q \f$ is zero for all non-boundary
 * nodes and depends on \c ApplyBoundaryConditions to assign values to the boundary nodes, it does
 * not need to be passed as an argument to \c ApplyBoundaryConditions.
 *
 * The hydrostatic vector \f$ h \f$ is defined as \f$ h_j = \rho g \gamma_j L_j \f$.
 * For the capillary network, the hydrostatic vector is multiplied by the geometry matrix \f$ G \f$
 * and connectivity matrix \f$ C \f$, yielding the vector \f$ C G h \f$.
 *
 * \retval      CGh         Hydrostatic term for the capillary network algorithm
 */
  // Calculate gamma vector
  geometry_->CalculateGammaVector();

  // Calculate hydrostatic vector
  hydrostatic_vector_ = geometry_->GetGammaVector()
                      % geometry_->GetLinkLength()
                      * gravity_
                      * fluid_->GetProperty("density");

  arma::vec CGh = geometry_->GetConnectivityMatrix()
                * geometry_->GetGeometryMatrix()
                * hydrostatic_vector_;

  return CGh;
}  // end of HydrostaticPhysics::CalculateModelRightHandSideVector() method



arma::Col<double> HydrostaticPhysics::CalculateFlowRate(const arma::Col<double> &pressures) const {
/**
 * The \c HydrostaticPhysics::CalculateFlowRate() method calculates the flow rate \f$ Q_j \f$ for
 * each capillary \f$ j \f$ as
 *
 * \f[
 *    Q_j = G_j \Delta P_j - G_j h_j \, ,
 * \f]
 *
 * where \f$ G_j = \frac{\pi R_j^4}{8 \mu L_j} \f$ is the geometry matrix diagonal element and
 * \f$ h_j = \rho g \gamma_j L_j \f$ is the hydrostatic pressure term.
 *
 * \param[in]   pressures   Pressure vector resulting from the matrix equation solution.
 * \retval      flow_rate   Vector containing the flow rate for each capillary
 */
  // Get G matrix diagonal
  arma::vec geometry_matrix_diagonal(geometry_->GetGeometryMatrix().diag());

  // Calculate flow rate
  arma::vec flow_rate = CalculateFlowRateBase(geometry_, pressures)
                      - geometry_matrix_diagonal % hydrostatic_vector_;

  return flow_rate;
}  // end of HydrostaticPhysics::CalculateFlowRate() method

}  // namespace simulator
