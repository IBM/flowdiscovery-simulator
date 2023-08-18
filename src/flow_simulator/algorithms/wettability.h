/**
 * \file src/flow_simulator/algorithms/wettability.h
 * \brief Contains the \c Wettability class.
 *
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2019
 *
 * This header file contains the \c Wettability class.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_WETTABILITY_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_WETTABILITY_H_

#include <map>
#include <memory>
#include <string>
#include <vector>
#include "src/flow_simulator/algorithms/fluid.h"
#include "src/exec_manager/wettability_json.h"

using Fluid = simulator::Fluid;

namespace simulator {

/**
 * \class Wettability wettability.h "src/flow_simulator/algorithms/wettability.h"
 * \brief Encapsulates Wettability properties.
 */

class Wettability {
 public:
  /// Parametrised constructor
  explicit Wettability(WettabilityJSON wettability_json) {
    properties_ = wettability_json.GetProperties();
  }

  /// Adds a wettability property
  void AddProperty(const std::string property, const double value) {
    properties_[property] = value;
  }

  /// Adds a fluid
  void AddFluid(std::shared_ptr<Fluid> &fluid) { fluids_.push_back(fluid); }

  /// Getter for fluids_
  const std::vector<std::shared_ptr<Fluid>> &GetFluids() const { return fluids_; }

  /// Getter for \c properties_[name]
  double GetProperty(const std::string name) { return properties_[name]; }

  /// Return true if has the property informed by name, false otherwise
  bool HasProperty(std::string name) { return properties_.find(name) != properties_.end(); }

  /// Getter for \c properties_
  const std::map<std::string, double> &GetProperties() const { return properties_; }

 private:
  /// A vector with Fluid objects
  std::vector<std::shared_ptr<Fluid>> fluids_;

  /// A dictionary with all Wettability properties
  std::map<std::string, double> properties_;
};  // end of class Wettability

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_WETTABILITY_H_
