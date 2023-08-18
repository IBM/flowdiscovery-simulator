/**
 * \file src/flow_simulator/algorithms/carbonated_water_viscosity_behaviour.h
 * \brief Contains the \c CarbonatedWaterViscosityBehaviour class.
 *
 * \authors Adolfo Emmanuel Correa Lopez \<adolfo-correa@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2020
 *
 * This source file contains the implementation of \c CarbonatedWaterViscosityBehaviour class methods.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_CARBONATED_WATER_VISCOSITY_BEHAVIOUR_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_CARBONATED_WATER_VISCOSITY_BEHAVIOUR_H_

#include "src/flow_simulator/algorithms/viscosity_behaviour.h"
#include "src/flow_simulator/algorithms/pure_water_viscosity_behaviour.h"

using PureWaterViscosityBehaviour = simulator::PureWaterViscosityBehaviour;

namespace simulator {

/**
 * \class CarbonatedWaterViscosityBehaviour carbonated_water_viscosity_behaviour.h
 * "src/flow_simulator/algorithms/carbonated_water_viscosity_behaviour.h"
 * \brief Encapsulates the carbonated water fluid viscosity property behaviour.
 *
 * This class is responsible for calculating the dynamic viscosity of carbonated water
 * as a function of pressure and temperature.
 *
 * This class computes the following expression in order to determine the viscosity of
 * carbonated water:
 *
 * \f[ \mu_r = 1 + \frac{\sum_{i=1}^2a_i x^i_{CO_2}}{\sum_{i=0}^1 b_iT^i}\f]
 * \f[ \mu = \mu_r \times \mu_{h_2O} \f]
 *
 * Source for the equation used in this class:
 * http://pubs.geothermal-library.org/lib/grc/1030394.pdf
 */
class CarbonatedWaterViscosityBehaviour : public ViscosityBehaviour {
 public:
  /// Parametrised constructor
  explicit CarbonatedWaterViscosityBehaviour(double temperature, double co2_fraction) :
    temperature_(temperature), co2_fraction_(co2_fraction / 100.0), pure_water_(temperature) {
    // \f$ \sum_{i=1}^2a_{i-1} x^i_{CO_2} \f$
    arma::vec co2_fraction_powers = { co2_fraction_, co2_fraction_ };

    // Base contribution
    base_and_temperature_contribution_ =
      arma::accu(coefficients_.col(0ul) % arma::cumprod(co2_fraction_powers));

    // \f$ \sum_{i=0}^1 b_iT^i \f$
    const double t = temperature;
    arma::vec temperature_powers = { 1.0, t };
    base_and_temperature_contribution_ /=
      arma::accu(coefficients_.col(1ul) % arma::cumprod(temperature_powers));

    base_and_temperature_contribution_ += 1.0;

    // base_and_temperature_contribution_ contains the computed value of \mu_r
  }

  /// Get fluid viscosity in Pa.sec
  double GetViscosity(double pressure_in_pascal) const {
    return base_and_temperature_contribution_ * pure_water_.GetViscosity(pressure_in_pascal);
  }

  /// Return true if the temperature and pressure fall between equation designed limits
  bool IsDesignedFor(const double lower_pressure,
                     const double upper_pressure) const override {
    // Pure water limitations: 20 ºC <= T <= 105 ºC
    if (!pure_water_.IsDesignedFor(lower_pressure, upper_pressure)) {
      return false;
    }
    // Carbonated water limitations: 20 ºC <= T <= 50 ºC
    if (temperature_ < 293.15 || temperature_ > 323.15) {
      return false;
    }
    // Carbonated water limitations: 0.1 MPa <= P <= 40 MPa
    if (upper_pressure > 4.0e7 || lower_pressure < 1.0e5) {
      return false;
    }
    // Carbonated water limitations: 0% < xCO2 <= 4.8%
    if (co2_fraction_ <= 0 || co2_fraction_ > .048) {
      return false;
    }

    return true;
  }

 private:
  /// The constant value of the temperature in which the experiment is submitted
  double temperature_;

  /// The constant value of the CO2 mass fraction in which the experiment is submitted
  double co2_fraction_;

  /// Contains the base and temperature contributions for the final viscosity
  double base_and_temperature_contribution_;

  /// Pre-computed coefficient for fast execution of \c GetViscosity
  double pressure_coefficient_;

  /// Equation coefficients for the calculation of the carbonated water viscosity
  const arma::Mat<double> coefficients_ = {{ 7.632609119e+02,  -1.047187396332e+04},
                                           {-9.46077673946e+03, 3.68325597e+01    }};

  PureWaterViscosityBehaviour pure_water_;
};  // end of class CarbonatedWaterViscosityBehaviour

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_CARBONATED_WATER_VISCOSITY_BEHAVIOUR_H_
