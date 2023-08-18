/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/models/linear_molecular_kinetics/linear_molecular_kinetics_physics.cc
 * \brief Contains the implementation of \c LinearMolecularKineticsPhysics methods.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 *
 * This source file contains the implementation of \c LinearMolecularKineticsPhysics class methods.
 */

#include "src/flow_simulator/algorithms/dynamic_capillary_network/models/linear_molecular_kinetics/linear_molecular_kinetics_physics.h"

namespace simulator {

arma::Col<double> LinearMolecularKineticsPhysics::CalculatePhasePermeability(
  const arma::Col<double> &flow_rate,
  const arma::Col<double> &pressure) const {
/**
 * Calculates the phase permeability \f$ K_{\alpha} \f$ along the \c flow_axis for the resident and
 * injected fluids according to
 *
 * \f[
 *    K_{\alpha} = \frac{Q_{\text{out}} \mu_{\alpha} L}{A \Delta P} \, ,
 * \f]
 *
 * where \f$ Q_{\text{out}} \f$ is the outlet flow rate, \f$ \mu_{\alpha} \f$ is the
 * dynamic viscosity of phase \f$ \alpha \f$, \f$ L \f$ is the length of the sample along
 * \c flow_axis, \f$ A \f$ is the cross-sectional area of the sample perpendicular
 * to \c flow_axis and \f$ \Delta P \f$ is the pressure difference in the sample along \c flow_axis.
 *
 * \param[in]  flow_rate          Column vector containing the flow rates in each capillary.
 * \param[in]  pressure           Column vector containing pressure values at network nodes.
 * \retval     permeability       Column vector with phase permeabilities of the sample along flow
 *                                axis.
 */
  // Calculate individual terms
  double outlet_flow_rate = CalculateFlowRateAtOutlet(flow_rate);
  double midpoint_pressure = (arma::mean(pressure(context_->inlet_nodes_))
                           + arma::mean(pressure(context_->outlet_nodes_))) / 2;
  arma::rowvec viscosity = { GetResidentFluidViscosity(midpoint_pressure),
                             GetInjectedFluidViscosity(midpoint_pressure) };
  double length = context_->voxel_size_
                * static_cast<double>(context_->shape_[context_->flow_axis_]);
  double area = std::pow(context_->voxel_size_, 3.0)
              * static_cast<double>(context_->shape_[0] * context_->shape_[1] * context_->shape_[2])
              / length;
  double delta_P = arma::mean(pressure(context_->inlet_nodes_))
                 - arma::mean(pressure(context_->outlet_nodes_));

  return (length / area / delta_P) * arma::trans(viscosity * outlet_flow_rate);
}  // end of LinearMolecularKineticsPhysics::CalculatePermeability() method

}  // namespace simulator
