/**
 * \file src/flow_simulator/algorithms/viscosity_behaviour.h
 * \brief Contains the \c ViscosityBehaviour class.
 *
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2020
 *
 * This header file contains the \c ViscosityBehaviour class.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_VISCOSITY_BEHAVIOUR_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_VISCOSITY_BEHAVIOUR_H_

namespace simulator {

/**
 * \class ViscosityBehaviour viscosity_behaviour.h "src/flow_simulator/algorithms/viscosity_behaviour.h"
 * \brief Encapsulates the fluid viscosity property behaviour
 */
class ViscosityBehaviour {
 public:
  /// Default constructor
  ViscosityBehaviour() = default;

  /// Parametrised constructor
  explicit ViscosityBehaviour(double temperature) : temperature_(temperature) { }

  /// Return the dynamic viscosity of the fluid
  virtual double GetViscosity(double pressure_in_pascal) const = 0;

  /// Return true if the parameters fall between the designed behaviour limits, false otherwise
  virtual bool IsDesignedFor(const double lower_pressure,
                             const double upper_pressure) const = 0;
 private:
  /// The temperature of the experiment
  double temperature_;
};  // end of class ViscosityBehaviour

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_VISCOSITY_BEHAVIOUR_H_
