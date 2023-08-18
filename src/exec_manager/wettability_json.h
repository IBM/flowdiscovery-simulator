/**
 * \file src/exec_manager/wettability_json.h
 * \brief Contains the \c WettabilityJSON class.
 *
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2019
 *
 * This header file contains the \c WettabilityJSON class.
 */

#ifndef SRC_EXEC_MANAGER_WETTABILITY_JSON_H_
#define SRC_EXEC_MANAGER_WETTABILITY_JSON_H_

#include <string>
#include <map>
#include <vector>

/**
 * \class WettabilityJSON wettability_json.h "src/exec_manager/wettability_json.h"
 * \brief The representation, in-memory, of the wettability defined at the input JSON file.
 */

class WettabilityJSON {
 public:
  /// Adds a wettability input parameter property
  void AddProperty(const std::string property, const double value) {
    properties_[property] = value;
  }

  /// Getter for wettability input parameter name
  const std::string &GetName() const { return name_; }

  /// Getter for the wettability input parameter property with key equals name
  double GetProperty(const std::string name) { return properties_[name]; }

  /// Getter for the dictionary with all wettability input properties
  const std::map<std::string, double> &GetProperties() const { return properties_; }

  /// Setter for wettability input parameter name
  void SetName(std::string name) { name_ = name; }

 private:
  /// The wettability, in-memory, input parameter name
  std::string name_;

  /// The dictionary with all wettability, in-memory, input properties
  std::map<std::string, double> properties_;
};  // end of class WettabilityJSON

#endif  // SRC_EXEC_MANAGER_WETTABILITY_JSON_H_
