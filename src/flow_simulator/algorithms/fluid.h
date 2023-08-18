/**
 * \file src/flow_simulator/algorithms/fluid.h
 * \brief Contains the \c Fluid class.
 *
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2019
 *
 * This header file contains the \c Fluid class.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_FLUID_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_FLUID_H_

#include <string>
#include <map>
#include <memory>

#include "src/exec_manager/fluid_json.h"
#include "src/flow_simulator/algorithms/viscosity_behaviour.h"
#include "src/flow_simulator/algorithms/viscosity_behaviour_factory.h"


namespace simulator {

/**
 * \class Fluid fluid.h "src/flow_simulator/algorithms/fluid.h"
 * \brief Encapsulates Fluid properties.
 */

class Fluid {
 public:
  /// Parametrised constructor
  Fluid(FluidJSON &fluid_json, ExperimentJSON &experiment_json)
    : name_(fluid_json.GetName()) {
    viscosity_behaviour_ = ViscosityBehaviourFactory::CreateViscosityBehaviour(fluid_json,
                                                                             experiment_json);
    properties_ = fluid_json.GetProperties();
  }

  /// Adds a fluid property
  void AddProperty(std::string property, double value) { properties_[property] = value; }

  /// Getter for name_
  const std::string &GetName() const { return name_; }

  /// Delegate the request for viscosity to ViscosityBehaviour
  double GetViscosity(double pressure_in_pascal) const {
    return viscosity_behaviour_->GetViscosity(pressure_in_pascal);
  }

  /// Getter for \c properties_[name]
  double GetProperty(std::string name) { return properties_[name]; }

  /// Return true if has the property informed by name, false otherwise
  bool HasProperty(std::string name) { return properties_.find(name) != properties_.end(); }

  /// Getter for \c properties_
  const std::map<std::string, double> &GetProperties() const { return properties_; }

  /// Return true if the parameters fall between the designed behaviour limits, false otherwise
  virtual bool IsViscosityBehaviourDesignedFor(const double lower_pressure,
                                              const double upper_pressure) const {
    return viscosity_behaviour_->IsDesignedFor(lower_pressure, upper_pressure);
  }

 private:
  /**
   * \brief The fluid name
   */
  std::string name_;

  /**
   * \brief The fluid viscosity behaviour
   */
  std::unique_ptr<ViscosityBehaviour> viscosity_behaviour_;

  /**
   * \brief The dictionary with all fluid properties
   */
  std::map<std::string, double> properties_;
};  // end of class Fluid

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_FLUID_H_
