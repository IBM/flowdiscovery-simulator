/**
 * \file src/exec_manager/exec_manager.h
 * \brief Contains the \c ExecutionManager class declaration.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This header file contains the \c ExecutionManager class that handles different execution modes.
 */

#ifndef SRC_EXEC_MANAGER_EXEC_MANAGER_H_
#define SRC_EXEC_MANAGER_EXEC_MANAGER_H_

#include <string>
#include "src/exec_manager/config_reader.h"

/**
 * \class ExecutionManager exec_manager.h "src/exec_manager/exec_manager.h"
 * \brief Executes the appropriate routines according to command-line options.
 *
 * This class is responsible for calling the appropriate \c FlowSimulator methods
 * depending on the command-line options passed to \c CmdLineParser.
 */

class ExecutionManager {
 public:
  /// Method that executes the "simulation" mode
  void RunSimulation(const std::string &json_file_name) const;

 private:
  /// Private \c ConfigReader member to parse the configuration JSON file
  ConfigReader config_reader_;
};  // end of class ExecutionManager

#endif  // SRC_EXEC_MANAGER_EXEC_MANAGER_H_
