/**
 * \file src/flow_simulator/algorithms/static_capillary_network/models/poiseuille.cc
 * \brief Contains the implementation of \c PoiseuillePhysics and \c PoiseuilleGeometry
 * class methods.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This source file contains the implementation of \c PoiseuillePhysics and
 * \c PoiseuilleGeometry class methods.
 */

#include "src/flow_simulator/algorithms/static_capillary_network/models/poiseuille.h"

namespace simulator {

arma::Col<double> PoiseuillePhysics::CalculateRightHandSideVectorTerm(void) {
/**
 * The \c PoiseuillePhysics::CalculateRightHandSideVectorTerm method specialises the
 * \c StaticCapillaryPhysicsBase::CalculateRightHandSideVectorTerm with model-dependent
 * calculations related to \c Poiseuille.
 *
 * For the Poiseuille model the boundary conditions are applied on top of the \f$ C Q \f$
 * vector. This vector has zero values on all non-boundary nodes due to mass conservation.
 * The non-zero values are attributed as part of \c ApplyBoundaryConditions.
 *
 * \retval      rhs_vector  Right-hand side vector of the matrix equation
 */
  // Initialise Cq vector term
  return arma::vec(geometry_->GetNumberOfNodes(), arma::fill::zeros);
}  // end of PoiseuillePhysics::CalculateModelRightHandSideVector() method



arma::Col<double> PoiseuillePhysics::CalculateFlowRate(const arma::Col<double> &pressures) const {
/**
 * The \c PoiseuillePhysics::CalculateFlowRate() method calculates the flow rate \f$ Q_j \f$ for
 * each capillary \f$ j \f$ according to \c StaticCapillaryPhysicsBase::CalculateFlowRateBase().
 *
 * \param[in] pressures   Pressure vector resulting from the solution of the matrix equation.
 * \retval    flow_rate   Vector containing the flow rate for each capillary.
 */
  return CalculateFlowRateBase(geometry_, pressures);
}  // end of PoiseuillePhysics::CalculateFlowRate() method

}  // namespace simulator
