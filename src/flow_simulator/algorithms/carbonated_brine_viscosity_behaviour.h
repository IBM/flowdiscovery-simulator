/**
 * \file src/flow_simulator/algorithms/carbonated_brine_viscosity_behaviour.h
 * \brief Contains the \c CarbonatedBrineViscosityBehaviour class.
 *
 * \authors Adolfo Emmanuel Correa Lopez \<adolfo-correa@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2021
 *
 * This source file contains the implementation of \c CarbonatedBrineViscosityBehaviour class methods.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_CARBONATED_BRINE_VISCOSITY_BEHAVIOUR_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_CARBONATED_BRINE_VISCOSITY_BEHAVIOUR_H_

#include "src/flow_simulator/algorithms/viscosity_behaviour.h"
#include "src/flow_simulator/algorithms/brine_viscosity_behaviour.h"

using BrineViscosityBehaviour = simulator::BrineViscosityBehaviour;

namespace simulator {

/**
 * \class CarbonatedBrineViscosityBehaviour carbonated_brine_viscosity_behaviour.h
 * "src/flow_simulator/algorithms/carbonated_brine_viscosity_behaviour.h"
 * \brief Encapsulates the cabonated brine fluid viscosity property behaviour.
 *
 * This class is responsible for calculating the dynamic viscosity of carbonated brine
 * as a function of pressure and temperature.
 *
 * This class computes the following expression in order to determine the viscosity of
 * carbonated brine:
 *
 * \f[
       \mu_\text{H$_2$O + NaCl + CO$_2$} = \mu_\text{H$_2$O + NaCl} \left( 1 + 4.65 x_\text{CO$_2$}^{1.0134} \right)
 * \f]
 *
 * Source for the equation used in this class:
 * https://pubs.acs.org/doi/10.1021/ef3006228
 */
class CarbonatedBrineViscosityBehaviour : public ViscosityBehaviour {
 public:
  /**
   * \brief Parametrised constructor
   * \param[in] temp Temperature of the fluid in kelvin
   * \param[in] salinity Parts per million of NaCl
   * \param[in] co2_fraction CO2 - Carbon Dioxide fraction
   */
  explicit CarbonatedBrineViscosityBehaviour(double temp, double salinity, double co2_fraction) :
    temperature_(temp), co2_fraction_(co2_fraction / 100.0), brine_(temp, salinity) {}

  /**
   * \brief Get fluid viscosity in Pa.sec
   * Using the formula: \f[ \mu_{H_2O+NaCl+CO_2} = \mu_{H_2O+NaCl}(1+4.65x_{CO_2}^{1.0134}) \f]
   * \param[in] pressure_in_pascal Pressure of the fluid in Pascal unit
   * \retval Viscosity in Pa.s
   */
  double GetViscosity(double pressure_in_pascal) const {
    return brine_.GetViscosity(pressure_in_pascal) * (1 + 4.65 * std::pow(co2_fraction_, 1.0134));
  }

  /// Return true if the temperature and pressure fall between equation designed limits
  bool IsDesignedFor(const double lower_pressure,
                     const double upper_pressure) const override {
    // Brine limitations: 20 ºC <= T <= 105 ºC
    if (!brine_.IsDesignedFor(lower_pressure, upper_pressure)) {
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

  /// BrineViscosityBehaviour Object
  BrineViscosityBehaviour brine_;
};  // end of class CarbonatedBrineViscosityBehaviour

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_CARBONATED_BRINE_VISCOSITY_BEHAVIOUR_H_
