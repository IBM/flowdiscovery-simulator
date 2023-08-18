/**
 * \file src/flow_simulator/algorithms/seawater_viscosity_behaviour.h
 * \brief Contains the \c SeawaterViscosityBehaviour class.
 *
 * \authors Adolfo Emmanuel Correa López \<adolfo-correa@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2021
 *
 * This source file contains the implementation of \c SeawaterViscosityBehaviour class methods
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_SEAWATER_VISCOSITY_BEHAVIOUR_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_SEAWATER_VISCOSITY_BEHAVIOUR_H_

#include <glog/logging.h>

#include "src/flow_simulator/algorithms/viscosity_behaviour.h"

namespace simulator {

/**
 * \class SeawaterViscosityBehaviour seawater_viscosity_behaviour.h
 * "src/flow_simulator/algorithms/seawater_viscosity_behaviour.h"
 * \brief Encapsulates the seawater fluid viscosity property behaviour.
 *
 * \c SeawaterViscosityBehaviour is responsible for return accurate information
 * about physical viscosity property based upon a fixed temperature and a variable pressure.
 *
 * The following expression is computed in order to determine the viscosity of
 * pure water:
 *
 * \f[
 *    \mu_\text{seawater} = a_0 + \sum_{i=1}^3 b_i e^{-c_i T} + P\sum_{i=0}^3 d_i(T - 293.15)^i
 * \f]
 *
 * Source for the equation used in this class:
 * http://pubs.geothermal-library.org/lib/grc/1030394.pdf
 */
class SeawaterViscosityBehaviour : public ViscosityBehaviour {
 public:
  /// Parametrised constructor
  explicit SeawaterViscosityBehaviour(double temperature):temperature_(temperature) {
    // Base contribution
    base_and_temperature_contribution_ = coefficients_(0ul, 0ul);

    base_and_temperature_contribution_ += arma::accu(coefficients_.col(1ul) %
                                                     arma::exp(-coefficients_.col(2ul) *
                                                               temperature));

    base_and_temperature_contribution_ /= 1.0e6;  // Pas.s

    // Setup pressure coefficient constant
    pressure_coefficient_ = 0.0;

    const double t = temperature - 293.15;
    arma::vec temperature_powers = { 1.0, t, t, t };
    pressure_coefficient_ += arma::accu(coefficients_.col(3ul) % arma::cumprod(temperature_powers));

    pressure_coefficient_ /= (1.0e5 * 1e6);  // Pascal to Bar, and uPas.s to Pas.s
  }

  /// Get fluid viscosity in Pa.sec
  double GetViscosity(double pressure_in_pascal) const {
    return base_and_temperature_contribution_ + (pressure_in_pascal * pressure_coefficient_);
  }

  /// Return true if the temperature and pressure fall between equation designed limits
  bool IsDesignedFor(const double lower_pressure,
                     const double upper_pressure) const override {
    // Seawater limitations: 20 ºC <= T <= 105 ºC
    if (temperature_ < 293.15 || temperature_ > 378.15) {
      return false;
    }
    // Seawater limitations: 0.1 MPa <= P <= 60 MPa
    if (upper_pressure > 6.0e7 || lower_pressure < 1.0e5) {
      return false;
    }

    return true;
  }

 private:
  /// The constant value of the temperature in wich the experiment is submitted
  double temperature_;

  /// Contains the base and temperature contributions for the final viscosity
  double base_and_temperature_contribution_;

  /// Pre-alculated coefficient for fast execution of \c GetViscosity
  double pressure_coefficient_;
  /// Equation coefficients for the calculation of the pure water viscosity
  const arma::mat coefficients_ = {
    { 2.27027348e+01,             0.0,            0.0, -2.08076250e-01 },
    {            0.0,  1.09512918e+04, 9.89379527e-03,  2.02116962e-02 },
    {            0.0,  9.59424756e+08, 4.38767404e-02, -2.80171705e-04 },
    {            0.0, -8.60130112e+08, 4.41842469e-02,  1.19302430e-06 } };
};  // end of class SeawaterViscosityBehaviour

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_SEAWATER_VISCOSITY_BEHAVIOUR_H_
