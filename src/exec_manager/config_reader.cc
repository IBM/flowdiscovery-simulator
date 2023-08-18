/**
 * \file src/exec_manager/config_reader.cc
 * \brief Contains the methods of the \c ConfigReader class.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 *
 * This source file contains the parse methods that handle each individual JSON file related to each
 * execution mode.
 */

#include <rapidjson/schema.h>
#include "src/exec_manager/config_reader.h"
#include "src/exec_manager/setup_config.h"
#include "src/exec_manager/simulation_config.h"
#include "src/exec_manager/fluid_json.h"

void ConfigReader::PopulateSetupConfig(SetupConfig &setup_cfg,
                                       const std::string &json_file_name) const {
/**
 * This method parses and validates the JSON file containing setup-related configuration parameters.
 *
 * \param[out] setup_cfg       A object containing all setup-related configuration parameters.
 * \param[in]  json_file_name  The address of the JSON file with respect to project root.
 */
  // Parse file into rapidjson::Document object and validate JSON file with respect to JSON schema
  rapidjson::Document json_doc = ParseIntoJsonDocument(json_file_name);
  std::string schema_file_name = "util/config_schema.json";
  ValidateJsonDocument(json_doc, schema_file_name);

  // Check if document contains the required information
  if (json_doc.HasMember("setup")) {
    auto setup_obj = json_doc["setup"].GetObject();

    // Populate object
    setup_cfg.folder = setup_obj["folder"].GetString();
    setup_cfg.input_file = setup_obj["input_file"].GetString();
    setup_cfg.voxel_size = setup_obj["voxel_size"].GetDouble();
    setup_cfg.shape = {setup_obj["shape"]["x"].GetUint64(),
                       setup_obj["shape"]["y"].GetUint64(),
                       setup_obj["shape"]["z"].GetUint64()};
  } else {
    std::fprintf(stderr, "JSONPARSER::POPULATESETUPCONFIG says: ");
    std::fprintf(stderr, "The '%s' file lacks the 'setup' field.\n", json_file_name.c_str());
    std::exit(-1);
  }
}  // end of ConfigReader::PopulateSetupConfig method


void ConfigReader::PopulateSimulationConfig(SimulationConfig &simulation_cfg,
                                            const std::string &json_file_name) const {
/**
 * This method parses and validates the JSON file containing simulation-related configuration
 * parameters.
 *
 * \param[out] simulation_cfg  A object containing all simulation-related configuration parameters.
 * \param[in]  json_file_name  The address of the JSON file with respect to project root.
 */
  // Parse file into rapidjson::Document object and validate JSON file with respect to JSON schema
  rapidjson::Document json_doc = ParseIntoJsonDocument(json_file_name);
  std::string schema_file_name = "util/config_schema.json";
  ValidateJsonDocument(json_doc, schema_file_name);

  // Check if document contains the required information
  if (json_doc.HasMember("simulation")) {
    auto simulation_obj = json_doc["simulation"].GetObject();
    auto fluid_array = simulation_obj["fluid"].GetArray();

    // Populate List of Fluids
    for (const auto &fluid_obj : fluid_array) {
      std::string name = fluid_obj["name"].GetString();
      std::string viscosity_behaviour = fluid_obj["viscosity_behaviour"].GetString();
      FluidJSON my_fluid(name, viscosity_behaviour);
      if (fluid_obj.HasMember("properties")) {
        for (const auto &property_obj : fluid_obj["properties"].GetObject()) {
          my_fluid.AddProperty(property_obj.name.GetString(), property_obj.value.GetDouble());
        }
      }
      simulation_cfg.fluids_json.push_back(my_fluid);
    }

    // Populate Wettability
    if (simulation_obj.HasMember("wettability")) {
      if (simulation_cfg.fluids_json.size() < 2) {
        std::fprintf(stderr, "CONFIGREADER::GETSIMULATIONCONFIG says: ");
        std::fprintf(stderr, "`wettability` needs the definitions of two fluids\n");
        std::exit(-1);
      }
      auto wettability_obj = simulation_obj["wettability"].GetObject();
      simulation_cfg.wettability_json.SetName(wettability_obj["name"].GetString());
      for (const auto &property_obj : wettability_obj["properties"].GetObject()) {
        simulation_cfg.wettability_json.AddProperty(property_obj.name.GetString(),
                                                    property_obj.value.GetDouble());
      }
    }

    // Populate Interface
    if (simulation_obj.HasMember("interface")) {
      if (simulation_cfg.fluids_json.size() < 2) {
        std::fprintf(stderr, "CONFIGREADER::POPULATESIMULATIONCONFIG says: ");
        std::fprintf(stderr, "`interface` needs the definitions of two fluids\n");
        std::exit(-1);
      }
      auto interface_obj = simulation_obj["interface"].GetObject();
      simulation_cfg.fluid_interface_json.SetName(interface_obj["name"].GetString());
      for (const auto &property_obj : interface_obj["properties"].GetObject()) {
        simulation_cfg.fluid_interface_json.AddProperty(property_obj.name.GetString(),
                                                        property_obj.value.GetDouble());
      }
    }

    // Populate Algorithm
    auto algorithm_obj = simulation_obj["algorithm"].GetObject();
    simulation_cfg.algorithm_json.SetName(algorithm_obj["name"].GetString());
    simulation_cfg.algorithm_json.SetModel(algorithm_obj["model"].GetString());
    if (algorithm_obj.HasMember("initial_time"))
      simulation_cfg.algorithm_json.SetInitialTime(algorithm_obj["initial_time"].GetDouble());
    if (algorithm_obj.HasMember("final_time"))
      simulation_cfg.algorithm_json.SetFinalTime(algorithm_obj["final_time"].GetDouble());
    if (algorithm_obj.HasMember("time_step_size"))
      simulation_cfg.algorithm_json.SetTimeStepSize(
        algorithm_obj["time_step_size"].GetDouble());
    if (algorithm_obj.HasMember("relative_tolerance"))
      simulation_cfg.algorithm_json.SetRelativeTolerance(
        algorithm_obj["relative_tolerance"].GetDouble());
    if (algorithm_obj.HasMember("absolute_link_tolerance"))
      simulation_cfg.algorithm_json.SetAbsoluteLinkTolerance(
        algorithm_obj["absolute_link_tolerance"].GetDouble());
    if (algorithm_obj.HasMember("absolute_node_tolerance"))
      simulation_cfg.algorithm_json.SetAbsoluteNodeTolerance(
        algorithm_obj["absolute_node_tolerance"].GetDouble());
    if (algorithm_obj.HasMember("resume"))
      simulation_cfg.algorithm_json.SetResume(
        algorithm_obj["resume"].GetBool());

    // Populate Experiment
    auto experiment_obj = simulation_obj["experiment"].GetObject();
    simulation_cfg.experiment_json.SetFlowAxis(experiment_obj["flow_axis"].GetUint());
    if (experiment_obj.HasMember("temperature"))
      simulation_cfg.experiment_json.SetTemperature(experiment_obj["temperature"].GetDouble());
    simulation_cfg.experiment_json.SetAbsolutePressure(
      experiment_obj["absolute_pressure"].GetDouble());
    simulation_cfg.experiment_json.SetBoundaryThickness(
      experiment_obj["boundary_thickness"].GetUint());
    auto bc_obj = experiment_obj["boundary_condition"].GetObject();
    simulation_cfg.experiment_json.SetBoundaryCondition(bc_obj["driving_force"].GetString(),
                                                        bc_obj["value"].GetDouble());
  } else {
    std::fprintf(stderr, "JSONPARSER::POPULATESIMULATIONCONFIG says: ");
    std::fprintf(stderr, "The '%s' file lacks the 'simulation' field.\n", json_file_name.c_str());
    std::exit(-1);
  }
}  // end of ConfigReader::PopulateSimulationConfig method
