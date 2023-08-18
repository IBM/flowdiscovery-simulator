/**
 * \file src/flow_simulator/algorithms/supercritical_co2_viscosity_behaviour.h
 * \brief Contains the \c SupercriticalCO2ViscosityBehaviour class.
 *
 * \authors Adolfo Emmanuel Correa Lopez \<adolfo-correa@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2021
 *
 * This source file contains the implementation of \c SupercriticalCO2ViscosityBehaviour class methods.
 */
#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_SUPERCRITICAL_CO2_VISCOSITY_BEHAVIOUR_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_SUPERCRITICAL_CO2_VISCOSITY_BEHAVIOUR_H_

#include "src/flow_simulator/algorithms/viscosity_behaviour.h"

namespace simulator {

/**
 * \class SupercriticalCO2ViscosityBehaviour supercritical_co2_viscosity_behaviour.h
 * "src/flow_simulator/algorithms/supercritical_co2_viscosity_behaviour.h"
 * \brief Encapsulates the supercritical CO2 fluid viscosity property behaviour.
 *
 * This class is responsible for calculating the dynamic viscosity of supercritical CO2
 * as a function of pressure and temperature.
 *
 * This class computes the following expression in order to determine the viscosity of
 * a supercritical CO2:
 *
 * \f[ \mu = f(P,T) \f]
 *
 * \f[ \mu = \frac{A_1 + A_2p + A_3p^2 + A_4 \ln(T) + A_5(\ln(T))^2 + A_6(ln(T))^3}
 * {1+A_7p+A_8 \ln(T)+A_9(\ln(T))^2} \f]
 *
 * In this correlation the unit of viscosity is in Centipoise (cp) whereas
 * temperature and pressure are expressed in K and bar (0.1 MPa), respectively.
 *
 * Source for the equation used in this class:
 * Viscosity of pure carbon dioxide at supercritical region: Measurement and correlation approach
 * https://www.sciencedirect.com/science/article/pii/S0896844610005127
 */
class SupercriticalCO2ViscosityBehaviour : public ViscosityBehaviour {
 public:
  /// Parametrised constructor
  explicit SupercriticalCO2ViscosityBehaviour(double temperature) :temperature_(temperature) {
    numerator_ = coefficients_[0]
               + coefficients_[3] * log(temperature_)
               + coefficients_[4] * pow(log(temperature_), 2.0)
               + coefficients_[5] * pow(log(temperature_), 3.0);
    denominator_ = 1.0
                 + coefficients_[7] * log(temperature_)
                 + coefficients_[8] * pow(log(temperature_), 2.0);
  }

  /// Get fluid viscosity in Pa.sec
  inline double GetViscosity(double pressure_in_pascal) const {
    // Pascal to Bar
    double P = pressure_in_pascal * 1.0e-5;

    double numerator = numerator_
                     + coefficients_[1] * P
                     + coefficients_[2] * pow(P, 2.0);

    double denominator = denominator_ + coefficients_[6] * P;

    // 1 centipoise -> 10ˆ-3 Pa-s
    return numerator / denominator / 1.0e3;
  }

  /// Return true if the temperature and pressure fall between equation designed limits
  bool IsDesignedFor(const double lower_pressure,
                     const double upper_pressure) const override {
    // ScCO2 limitations: 310k <= T <= 900k
    if (temperature_ < 310.0 || temperature_ > 900.0) {
      return false;
    }
    // ScCO2 limitations: 7.5 MPa <= P <= 101.4 MPa
    if (lower_pressure < 7.5e6 || upper_pressure > 101.4e6) {
      return false;
    }
    return true;
  }

 private:
  /// The constant value of the temperature in which the experiment is submitted
  double temperature_;

  /// Numerator value of the temperature of the viscosity formula
  double numerator_;

  /// Denominator value of the temperature of the viscosity formula
  double denominator_;

  /// Equation coefficients for the calculation of the carbonated water viscosity
  const arma::vec coefficients_ = {-1.146067e-01,
                                    6.978380e-07,
                                    3.976765e-10,
                                    6.336120e-02,
                                   -1.166119e-02,
                                    7.142596e-04,
                                    6.519333e-06,
                                   -3.567559e-01,
                                    3.180473e-02};
};  // end of class SupercriticalCO2ViscosityBehaviour

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_SUPERCRITICAL_CO2_VISCOSITY_BEHAVIOUR_H_
