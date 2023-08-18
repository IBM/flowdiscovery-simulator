/**
 * \file src/flow_simulator/algorithms/constant_viscosity_behaviour.h
 * \brief Contains the \c ConstantViscosityBehaviour class.
 *
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2020
 *
 * This source file contains the implementation of \c ConstantViscosityBehaviour class methods
 * which are responsible for return the same informations about physical viscosity property.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_CONSTANT_VISCOSITY_BEHAVIOUR_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_CONSTANT_VISCOSITY_BEHAVIOUR_H_

#include "src/flow_simulator/algorithms/viscosity_behaviour.h"

namespace simulator {

/**
 * \class ConstantViscosityBehaviour constant_viscosity_behaviour.h
 * "src/flow_simulator/algorithms/constant_viscosity_behaviour.h"
 * \brief Encapsulates the constant fluid viscosity property behaviour.
 */
class ConstantViscosityBehaviour : public ViscosityBehaviour {
 public:
  explicit ConstantViscosityBehaviour(double dynamic_viscosity)
    : dynamic_viscosity_(dynamic_viscosity) {}

  /// Return the stored dynamic viscosity
  double GetViscosity(double pressure_in_pascal) const override {
    // `pressure_in_pascal` is not used for this implementation
    static_cast<void>(pressure_in_pascal);

    return dynamic_viscosity_;
  }

  /// Return true if the parameters fall between the designed behaviour limits, false otherwise
  bool IsDesignedFor(const double lower_pressure,
                     const double upper_pressure) const override {
    // `lower_pressure` and `upper_pressure` is not used for this implementation
    static_cast<void>(lower_pressure);
    static_cast<void>(upper_pressure);

    return true;
  }

 private:
  /// Stored dynamic viscosity
  double dynamic_viscosity_;
};  // end of class ConstantViscosityBehaviour

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_CONSTANT_VISCOSITY_BEHAVIOUR_H_
