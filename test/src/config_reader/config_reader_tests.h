/**
 * \file test/src/config_reader/config_reader_tests.h
 * \brief Contains regression tests of \c ConfigReader class methods.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 */

#ifndef TEST_SRC_CONFIG_READER_CONFIG_READER_TESTS_H_
#define TEST_SRC_CONFIG_READER_CONFIG_READER_TESTS_H_

#include <gtest/gtest.h>
#include <rapidjson/document.h>
#include <string>
#include "src/exec_manager/config_reader.h"
#include "src/exec_manager/setup_config.h"
#include "src/exec_manager/simulation_config.h"

TEST(ConfigReader, RapidJSON_Test) {
  std::string json_str = R"(
    {
    "hello": "world",
    "t": true ,
    "f": false,
    "n": null,
    "i": 123,
    "pi": 3.1416,
    "a": [1, 2, 3, 4]
    }
  )";

  rapidjson::Document json_doc;
  json_doc.Parse(json_str.c_str());

  // Test rapidjson::Document
  EXPECT_TRUE(json_doc.IsObject());

  // Test member "hello"
  EXPECT_TRUE(json_doc.HasMember("hello"));
  EXPECT_TRUE(json_doc["hello"].IsString());
  EXPECT_STREQ(json_doc["hello"].GetString(), "world");

  // Test member "t"
  EXPECT_TRUE(json_doc.HasMember("t"));
  EXPECT_TRUE(json_doc["t"].IsBool());
  EXPECT_TRUE(json_doc["t"].GetBool());

  // Test member "f"
  EXPECT_TRUE(json_doc.HasMember("f"));
  EXPECT_TRUE(json_doc["f"].IsBool());
  EXPECT_FALSE(json_doc["f"].GetBool());

  // Test member "n"
  EXPECT_TRUE(json_doc.HasMember("n"));
  EXPECT_TRUE(json_doc["n"].IsNull());

  // Test member "i"
  EXPECT_TRUE(json_doc.HasMember("i"));
  EXPECT_TRUE(json_doc["i"].IsNumber());
  EXPECT_TRUE(json_doc["i"].IsInt());
  EXPECT_FALSE(json_doc["i"].IsDouble());

  // Test member "pi"
  EXPECT_TRUE(json_doc.HasMember("pi"));
  EXPECT_TRUE(json_doc["pi"].IsNumber());
  EXPECT_FALSE(json_doc["pi"].IsInt());
  EXPECT_TRUE(json_doc["pi"].IsDouble());

  // Test member "a"
  EXPECT_TRUE(json_doc.HasMember("a"));
  EXPECT_TRUE(json_doc["a"].IsArray());
  uint64_t i = 0U;
  for (const auto &v : json_doc["a"].GetArray()) EXPECT_EQ(v.GetUint64(), ++i);

  // Test "missing" member
  EXPECT_FALSE(json_doc.HasMember("missing"));
}



TEST(ConfigReader, GetSetupConfig_Test) {
  ConfigReader config_reader;
  std::string folder = "test/src/config_reader/input/";
  SetupConfig cfg;

  // config_reader
  config_reader.PopulateSetupConfig(cfg, folder + "setup_config.json");

  EXPECT_STREQ(cfg.folder.c_str(), "results");
  EXPECT_STREQ(cfg.input_file.c_str(), "C1.raw");
  EXPECT_EQ(cfg.shape(0), 400U);
  EXPECT_EQ(cfg.shape(1), 400U);
  EXPECT_EQ(cfg.shape(2), 400U);
  EXPECT_DOUBLE_EQ(cfg.voxel_size, 2.85e-6);
}



TEST(ConfigReader, GetSimulationConfig_Test) {
  ConfigReader config_reader;
  std::string folder = "test/src/config_reader/input/";
  SimulationConfig cfg;
  config_reader.PopulateSimulationConfig(cfg, folder + "simulation_config.json");

  EXPECT_STREQ(cfg.algorithm_json.GetName().c_str(), "dynamic");
  EXPECT_STREQ(cfg.algorithm_json.GetModel().c_str(), "linear_molecular_kinetics");
  EXPECT_DOUBLE_EQ(cfg.algorithm_json.GetInitialTime(), 0.0);
  EXPECT_DOUBLE_EQ(cfg.algorithm_json.GetFinalTime(), 19.2);
  EXPECT_DOUBLE_EQ(cfg.algorithm_json.GetTimeStepSize(), 0.59375);
  EXPECT_DOUBLE_EQ(cfg.algorithm_json.GetRelativeTolerance(), 1.0e-4);
  EXPECT_DOUBLE_EQ(cfg.algorithm_json.GetAbsoluteLinkTolerance(), 1.0e-3);
  EXPECT_DOUBLE_EQ(cfg.algorithm_json.GetAbsoluteNodeTolerance(), 1.0e-3);
  EXPECT_FALSE(cfg.algorithm_json.GetResume());

  EXPECT_DOUBLE_EQ(cfg.wettability_json.GetProperty("contact_angle"), 0.0);
  EXPECT_DOUBLE_EQ(cfg.wettability_json.GetProperty("linear_mk"), 108.63);

  EXPECT_DOUBLE_EQ(cfg.fluid_interface_json.GetProperty("interfacial_tension"), 0.06347);

  EXPECT_DOUBLE_EQ(cfg.experiment_json.GetAbsolutePressure(), 101325.0);
  EXPECT_EQ(cfg.experiment_json.GetFlowAxis(), 0);
  EXPECT_DOUBLE_EQ(cfg.experiment_json.GetTemperature(), 340.0);
  EXPECT_EQ(cfg.experiment_json.GetBoundaryThickness(), 1);
  EXPECT_STREQ(cfg.experiment_json.GetBoundaryCondition().first.c_str(),
               "pressure_gradient_closed");
  EXPECT_DOUBLE_EQ(cfg.experiment_json.GetBoundaryCondition().second, 10132.5);

  EXPECT_STREQ(cfg.fluids_json[0].GetName().c_str(), "water");
  EXPECT_STREQ(cfg.fluids_json[0].GetViscosityBehaviourName().c_str(), "pure_water");
  EXPECT_STREQ(cfg.fluids_json[1].GetName().c_str(), "oil");
  EXPECT_DOUBLE_EQ(cfg.fluids_json[1].GetProperty("dynamic_viscosity"), 1.002e-4);
}

#endif  // TEST_SRC_CONFIG_READER_CONFIG_READER_TESTS_H_
