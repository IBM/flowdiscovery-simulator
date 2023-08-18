/**
 * \file src/exec_manager/config_reader.h
 * \brief Contains the declaration of the \c ConfigReader class.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 *
 * This header file contains the declaration of the \c ConfigReader class that extracts configuration
 * parameters from the JSON file.
 */

#ifndef SRC_EXEC_MANAGER_CONFIG_READER_H_
#define SRC_EXEC_MANAGER_CONFIG_READER_H_

#include <rapidjson/document.h>
#include <armadillo>
#include <utility>
#include <vector>
#include <string>
#include <map>
#include "src/json_parser/json_parser.h"
#include "src/exec_manager/setup_config.h"
#include "src/exec_manager/simulation_config.h"

/**
 * \class ConfigReader config_reader.h "src/exec_manager/config_reader.h"
 * \brief Provides an interface to the JSON parsing library.
 *
 * All library-specific details of the JSON parsing procedure should be hidden
 * inside this class that, in turn, is responsible for exposing the JSON configuration parameters to
 * the outside world.
 *
 * \note The current implementation uses the RapidJSON library \<http://rapidjson.org/\> and performs
 *       validation using a JSON schema.
 */

class ConfigReader : public JSONParser {
 public:
  /// Method that parses and validates the setup-related JSON configuration file
  void PopulateSetupConfig(SetupConfig &setup_cfg, const std::string &json_file_name) const;

  /// Method that parses and validates the simulation-related JSON configuration file
  void PopulateSimulationConfig(SimulationConfig &simulation_cfg,
                                const std::string &json_file_name) const;
};  // end of class ConfigReader

#endif  // SRC_EXEC_MANAGER_CONFIG_READER_H_
