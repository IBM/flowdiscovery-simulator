/**
 * \file src/exec_manager/fluid_json.h
 * \brief Contains the \c FluidJSON class.
 *
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2019
 *
 * This header file contains the \c FluidJSON class.
 */

#ifndef SRC_EXEC_MANAGER_FLUID_JSON_H_
#define SRC_EXEC_MANAGER_FLUID_JSON_H_

#include <string>
#include <map>

/**
 * \class FluidJSON fluid_json.h "src/exec_manager/fluid_json.h"
 * \brief The representation, in-memory, of the Fluid defined at the input JSON file.
 */

class FluidJSON {
 public:
  /// Parametrised constructor
  FluidJSON(std::string name, std::string viscosity_behaviour) :
    name_(name),
    viscosity_behaviour_(viscosity_behaviour) { }

  /// Adds a fluid input property
  void AddProperty(std::string property, double value) { properties_[property] = value; }

  /// Getter for fluid input name
  const std::string &GetName() const { return name_; }

  /// Getter for fluid input viscosity behaviour
  const std::string &GetViscosityBehaviourName() const { return viscosity_behaviour_; }

  /// Getter for the property with key equals name
  double GetProperty(std::string name) { return properties_[name]; }

  /// Return true if has the property informed by name, false otherwise
  bool HasProperty(std::string name) { return properties_.find(name) != properties_.end(); }

  /// Getter for the dictionary with all fluid input properties
  const std::map<std::string, double> &GetProperties() const { return properties_; }

 private:
  /**
   * \brief The Fluid, in-memory, input parameter name
   */
  std::string name_;

  /**
   * \brief The Fluid, in-memory, input parameter viscosity_behaviour
   */
  std::string viscosity_behaviour_;

  /**
   * \brief The dictionary with all Fluid, in-memory, input properties
   */
  std::map<std::string, double> properties_;
};  // end of class FluidJSON

#endif  // SRC_EXEC_MANAGER_FLUID_JSON_H_
