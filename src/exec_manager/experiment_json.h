/**
 * \file src/exec_manager/experiment_json.h
 * \brief Contains the \c ExperimentJSON class.
 *
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2019
 *
 * This header file contains the \c ExperimentJSON class.
 */

#ifndef SRC_EXEC_MANAGER_EXPERIMENT_JSON_H_
#define SRC_EXEC_MANAGER_EXPERIMENT_JSON_H_

#include <armadillo>
#include <string>
#include <map>
#include <utility>

/**
 * \class ExperimentJSON experiment_json.h "src/exe_manager/experiment_json.h"
 * \brief The representation, in-memory, of the experiment defined at the input JSON file.
 */

class ExperimentJSON {
 public:
  /// Getter for the experiment input parameter flow_axis
  arma::uword GetFlowAxis() const { return flow_axis_; }

  /// Getter for the experiment input parameter temperature
  double GetTemperature() const { return temperature_; }

  /// Getter for the experiment input parameter absolute_pressure
  double GetAbsolutePressure() const { return absolute_pressure_; }

  /// Getter for the experiment input parameter boundary_thickness
  arma::uword GetBoundaryThickness() const { return boundary_thickness_; }

  /// Getter for the experiment input parameter boundary_condition
  const std::pair<std::string, double> &GetBoundaryCondition() const { return boundary_condition_; }

  /// Setter for the experiment input parameter flow_axis
  void SetFlowAxis(arma::uword flow_axis) { flow_axis_ = flow_axis; }

  /// Setter for the experiment input parameter temperature
  void SetTemperature(double temperature) { temperature_ = temperature; }

  /// Setter for the experiment input parameter absolute_pressure
  void SetAbsolutePressure(double absolute_pressure) { absolute_pressure_ = absolute_pressure; }

  /// Setter for the experiment input parameter boundary_thickness
  void SetBoundaryThickness(arma::uword boundary_thickness) {
    boundary_thickness_ = boundary_thickness;
  }

  /// Setter for the experiment input parameter boundary_condition.
  void SetBoundaryCondition(std::string key, double value) {
    boundary_condition_.first = key;
    boundary_condition_.second = value;
  }

 private:
  /**
   * \brief The experiment, in-memory, input parameter flow_axis
   */
  arma::uword flow_axis_;

  /**
   * \brief The experiment, in-memory, input parameter temperature
   */
  double temperature_;

  /**
   * \brief The experiment, in-memory, input parameter absolute_pressure
   */
  double absolute_pressure_;

  /**
   * \brief The experiment, in-memory, input parameter bondary_thickness
   */
  arma::uword boundary_thickness_;

  /**
   * \brief The experiment, in-memory, input parameter boundary_condition
   */
  std::pair<std::string, double> boundary_condition_;
};  // end of class ExperimentJSON

#endif  // SRC_EXEC_MANAGER_EXPERIMENT_JSON_H_
