/**
 * \file src/flow_simulator/algorithms/brine_viscosity_behaviour.h
 * \brief Contains the \c BrineViscosityBehaviour class.
 *
 * \authors Adolfo Emmanuel Correa Lopez \<adolfo-correa@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2021
 *
 * This source file contains the implementation of \c BrineViscosityBehaviour class methods.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_BRINE_VISCOSITY_BEHAVIOUR_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_BRINE_VISCOSITY_BEHAVIOUR_H_

#include "src/flow_simulator/algorithms/viscosity_behaviour.h"

namespace simulator {

/**
 * \class BrineViscosityBehaviour brine_viscosity_behaviour.h
 * "src/flow_simulator/algorithms/brine_viscosity_behaviour.h"
 * \brief Encapsulates the brine fluid viscosity property behaviour.
 *
 * This class is responsible for calculating the dynamic viscosity of brine
 * as a function of pressure and temperature.
 *
 * This class computes the following expression in order to determine the viscosity of
 * brine:
 *
 * \f[ \ln \eta_r = A m + B m^2 + C m^3 \f]
 * \f[ A = a_0 + a_1 T + a_2 T^2 \f]
 * \f[ B = b_0 + b_1 T + b_2 T^2 \f]
 * \f[ C = c_0 + c_1T \f]
 * \f[ \eta_r = \frac{\eta_{sol}}{\eta_{H_2O}} \f]
 *
 * where: <br>
 *  \f$\eta \f$ is the viscosity,
 *  \f$m \f$ is molality (mol/kg) of salts
 *
 * Source for the equation used in this class:<br>
 * <a href="https://pubs.acs.org/doi/10.1021/ef3006228">
 * Viscosity Models and Effects of Dissolved CO2</a> <br>
 * <a href="https://link.springer.com/article/10.1007/s10765-009-0646-7">
 * The Viscosity of Aqueous Alkali-Chloride Solutions up to 623 K, 1,000 bar,
 * and High Ionic Strength</a> <br>
 */
class BrineViscosityBehaviour : public ViscosityBehaviour {
 public:
  /**
   * \brief Parametrised constructor
   * \param temperature Temperature of the fluid in kelvin
   * \param ppm Parts per million
   */
  explicit BrineViscosityBehaviour(double temperature, double ppm) :
    temperature_(temperature), molality_(ppm / 58.44277 * 1.0e-3) {
    const double t = temperature;
    arma::vec temperature_powers({1.0, t, t});

    abc_[0] = arma::accu(coefficients_.row(0u).t() % arma::cumprod(temperature_powers));
    abc_[1] = arma::accu(coefficients_.row(1u).t() % arma::cumprod(temperature_powers));
    abc_[2] = arma::accu(coefficients_.row(2u).t() % arma::cumprod(temperature_powers));

    rel_visc_ = std::exp(abc_[0] * molality_
                       + abc_[1] * std::pow(molality_, 2)
                       + abc_[2] * std::pow(molality_, 3));
  }

  /// Get fluid viscosity in Pa.sec
  inline double GetViscosity(double pressure_in_pascal) const {
    return rel_visc_ * WaterViscosity_(pressure_in_pascal);
  }

  /// Return true if the temperature and pressure fall between equation designed limits
  bool IsDesignedFor(const double lower_pressure,
                     const double upper_pressure) const override {
    // brine limitations: 273K <= T <= 623K
    if (temperature_ < 273.0 || temperature_ > 623.0) {
      return false;
    }
    // brine limitations: 1 bar <= P <= 10^3 bar (1 bar --> 10^5 Pa)
    if (upper_pressure > 1.0e8 || lower_pressure < 1.0e5) {
      return false;
    }
    // brine limitations: 0 mol/Kg < molality <= 6 mol/Kg
    if (molality_ <= 0.0 || molality_ > 6.0) {
      return false;
    }

    return true;
  }

 private:
  /// The constant value of the temperature in which the experiment is submitted
  double temperature_;

  /// The constant value of the molality of the salt present in the fluid meassured in mol/kg
  double molality_;

  /// Relative Viscosity
  double rel_visc_;

  /// Coefficients A, B and C used to compute the relative viscosity
  arma::vec::fixed<3> abc_;

  /// Equation coefficients (a, b, c) for the calculation of the brine viscosity
  const arma::Mat<double> coefficients_ = {{-0.21319213,     0.13651589e-2, -0.12191756e-5},
                                           { 0.69161945e-1, -0.27292263e-3,  0.20852448e-6},
                                           {-0.25988855e-2,  0.77989227e-5,  0.0          }};

  /// Equation coefficients for the calculation of the pure water viscosity
  arma::vec water_coef_ =
    {0.28853170e7, -0.11072577e5,  -0.90834095e1, 0.30925651e-1, -0.27407100e-4};

  /// Coefficients of the density contribution for the pure water viscosity calculation
  arma::vec density_coef_ =
    {-0.19283851e7,  0.56216046e4,  0.13827250e2, -0.47609523e-1,  0.35545041e-4};

  /**
   * \brief Get pure water viscosity in Pa.sec
   *
   * In order to do that, it computes the following expression:
   * \f[\ln \eta_{H_2O} = \sum_{i=1}^5 d_i T^{i-3} + \sum_{i=6}^{10} d_i \rho_{H_2O}T^{i-8}\f]
   * \f[\mu_{H_2O} = e^{\eta_{H_2O}}\f]
   */
  double WaterViscosity_(double pressure_in_pascal) const {
    // Auxiliary variables
    arma::vec temperature_powers = arma::ones<arma::vec>(5) * temperature_;
    temperature_powers(0) = std::pow(temperature_, -2);
    double density = WaterDensity_(pressure_in_pascal);
    double lneta = 0.0;

    // Water Viscosity formula
    lneta = arma::sum((water_coef_ + density_coef_ * density)
            % arma::cumprod(temperature_powers));

    return std::exp(lneta);
  }

  /**
   * @brief Get pure water density in Kg/m^3
   *
   * This computes the following expression:
   * \f[\rho_{H_2O} = a_0 + \sum_{i=1}^3 b_i10^{c_iT} + \sum_{i=1}^2 d_i P^i \f]
   */
  double WaterDensity_(double P)const {
    double T = temperature_;

    // Coeficients of the density formula
    double a = 1.34136579e02;
    arma::vec b({ -4.07743800e03, 1.63192756e04, 1.37091355e03});
    arma::vec c({ -5.56126409e-03, -1.07149234e-02, -5.46294495e-04});
    arma::vec d({ 4.45861703e-01, -4.51029739e-04});

    // Partial Expressions
    arma::vec aux1 = arma::exp10(c*T);
    arma::vec P_powers(2, arma::fill::ones);
    P *= 1.0e-6;    //  Pa to MPa
    P_powers *= P;
    arma::vec aux2 = arma::cumprod(P_powers);

    // Density formula
    double rho = a + arma::sum(b % aux1) + arma::sum(d % aux2);
    return rho * 1.0e-3;
  }
};  // end of class BrineViscosityBehaviour

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_BRINE_VISCOSITY_BEHAVIOUR_H_
