/**
 * \file src/flow_simulator/algorithms/fluid_interface.h
 * \brief Contains the \c FluidInterface class.
 *
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2019
 *
 * This header file contains the \c FluidInterface class.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_FLUID_INTERFACE_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_FLUID_INTERFACE_H_

#include <map>
#include <memory>
#include <vector>
#include <string>
#include "src/flow_simulator/algorithms/fluid.h"
#include "src/exec_manager/fluid_interface_json.h"

using Fluid = simulator::Fluid;

namespace simulator {

/**
 * \class FluidInterface fluid_interface.h "src/flow_simulator/algorithms/fluid_interface.h"
 * \brief Encapsulates FluidInterface properties.
 */

class FluidInterface {
 public:
  /// Parametrised constructor
  explicit FluidInterface(FluidInterfaceJSON fluid_interface_json) {
    properties_ = fluid_interface_json.GetProperties();
  }

  /// Adds a FluidInterface property
  void AddProperty(std::string property, double value) { properties_[property] = value; }

  /// Adds a fluid
  void AddFluid(std::shared_ptr<Fluid> &fluid) { fluids_.push_back(fluid); }

  /// Getter for fluids_
  std::vector<std::shared_ptr<Fluid>> &GetFluids() { return fluids_; }

  /// Getter for \c properties_[name]
  double &GetProperty(std::string name) { return properties_[name]; }

  /// Return true if has the property informed by name, false otherwise
  bool HasProperty(std::string name) { return properties_.find(name) != properties_.end(); }

  /// Getter for \c properties_
  std::map<std::string, double> &GetProperties() { return properties_; }

 private:
  /// A vector with Fluid Objects
  std::vector<std::shared_ptr<Fluid>> fluids_;

  /// A dictionary with all FluidInterface properties
  std::map<std::string, double> properties_;
};  // end of class Interface

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_FLUID_INTERFACE_H_
