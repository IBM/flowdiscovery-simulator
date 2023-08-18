/**
 * \file src/exec_manager/fluid_interface_json.h
 * \brief Contains the \c FluidInterfaceJSON class.
 *
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2019
 *
 * This header file contains the \c FluidInterfaceJSON class.
 */

#ifndef SRC_EXEC_MANAGER_FLUID_INTERFACE_JSON_H_
#define SRC_EXEC_MANAGER_FLUID_INTERFACE_JSON_H_

#include <string>
#include <map>
#include <vector>

/**
 * \class FluidInterfaceJSON fluid_interface_json.h "src/exec_manager/model/fluid_interface_json.h"
 * \brief The representation, in-memory, of the fluid_interface defined at the input JSON file.
 */

class FluidInterfaceJSON {
 public:
  /// Adds a fluid_interface input parameter property
  void AddProperty(std::string property, double value) { properties_[property] = value; }

  /// Getter for fluid_interface input parameter name.
  const std::string &GetName() const { return name_; }

  /// Getter for the fluid_interface input parameter property with key equals name.
  double GetProperty(std::string name) { return properties_[name]; }

  /// Getter for the fluid_interface input parameters properties.
  const std::map<std::string, double> &GetProperties() const { return properties_; }

  /// Setter for fluid_interface input parameter name.
  void SetName(std::string name) { name_ = name; }

 private:
  /**
   * \brief The in-memory representation of the JSON's fluid interface parameter \c name
   *
   * This variable stores the name of the fluid interface as provided in the JSON file.
   */
  std::string name_;

  /**
   * \brief The in-memory representation of the fluid interface physical properties
   *
   * This dictionary stores all the values of the physical properties of the fluid interface
   * provided in the JSON file.
   */
  std::map<std::string, double> properties_;
};  // end of class FluidInterfaceJSON

#endif  // SRC_EXEC_MANAGER_FLUID_INTERFACE_JSON_H_
