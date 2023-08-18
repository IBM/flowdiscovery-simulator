/**
 * \file src/exec_manager/algorithm_json.h
 * \brief Contains the \c AlgorithmJSON class.
 *
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2019
 *
 * This header file contains the \c AlgorithmJSON class.
 */

#ifndef SRC_EXEC_MANAGER_ALGORITHM_JSON_H_
#define SRC_EXEC_MANAGER_ALGORITHM_JSON_H_

#include <armadillo>
#include <string>
#include <map>

/**
 * \class AlgorithmJSON algorithm_json.h "src/exec_manager/algorithm_json.h"
 * \brief The representation, in-memory, of the algorithm defined at the input JSON file.
 */

class AlgorithmJSON {
 public:
  /// Getter for algorithm input parameter name.
  const std::string &GetName() const { return name_; }

  /// Getter for algorithm input parameter model.
  const std::string &GetModel() const { return model_; }

  /// Getter for algorithm input parameter initial_time.
  const double &GetInitialTime() const { return initial_time_; }

  /// Getter for algorithm input parameter final_time.
  const double &GetFinalTime() const { return final_time_; }

  /// Getter for algorithm input parameter time step size.
  const double &GetTimeStepSize() const { return time_step_size_; }

  /// Getter for algorithm input parameter relative_tolerance.
  const double &GetRelativeTolerance() const { return relative_tolerance_; }

  /// Getter for algorithm input parameter absolute_link_tolerance.
  const double &GetAbsoluteLinkTolerance() const { return absolute_link_tolerance_; }

  /// Getter for algorithm input parameter absolute_node_tolerance.
  const double &GetAbsoluteNodeTolerance() const { return absolute_node_tolerance_; }

  /// Getter for algorithm input parameter resume.
  const bool &GetResume() const { return resume_; }

  /// Setter for algorithm input parameter name.
  void SetName(std::string name) { name_ = name; }

  /// Setter for algorithm input parameter model.
  void SetModel(std::string model) { model_ = model; }

  /// Setter for algorithm input parameter initial_time.
  void SetInitialTime(double initial_time) { initial_time_ = initial_time; }

  /// Setter for algorithm input parameter final_time.
  void SetFinalTime(double final_time) { final_time_ = final_time; }

  /// Setter for algorithm input parameter time_step_size.
  void SetTimeStepSize(double time_step_size) { time_step_size_ = time_step_size; }

  /// Setter for algorithm input parameter relative_tolerance.
  void SetRelativeTolerance(double relative_tolerance) {
    relative_tolerance_ = relative_tolerance;
  }

  /// Setter for algorithm input parameter absolute_link_tolerance.
  void SetAbsoluteLinkTolerance(double absolute_link_tolerance) {
    absolute_link_tolerance_ = absolute_link_tolerance;
  }

  /// Setter for algorithm input parameter absolute_node_tolerance.
  void SetAbsoluteNodeTolerance(double absolute_node_tolerance) {
    absolute_node_tolerance_ = absolute_node_tolerance;
  }

  /// Setter for algorithm input parameter resume.
  void SetResume(bool resume) { resume_ = resume; }

 private:
  /**
   * \brief The algorithm, in-memory, input parameter name.
   */
  std::string name_;

  /**
   * \brief The algorithm, in-memory, input parameter model.
   */
  std::string model_;

  /**
   * \brief The algorithm, in-memory, input parameter initial_time.
   */
  double initial_time_;

  /**
   * \brief The algorithm, in-memory, input parameter final_time.
   */
  double final_time_;

  /**
   * \brief The algorithm, in-memory, input parameter time_step_size.
   */
  double time_step_size_;

  /**
   * \brief The algorithm, in-memory, input parameter relative_tolerance.
   */
  double relative_tolerance_;

  /**
   * \brief The algorithm, in-memory, input parameter absolute_link_tolerance.
   */
  double absolute_link_tolerance_;

  /**
   * \brief The algorithm, in-memory, input parameter absolute_node_tolerance.
   */
  double absolute_node_tolerance_;

  /**
   * \brief The algorithm, in-memory, input parameter resume.
   */
  bool resume_;
};  // end of class AlgorithmJSON

#endif  // SRC_EXEC_MANAGER_ALGORITHM_JSON_H_
