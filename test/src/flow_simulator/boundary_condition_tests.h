/**
 * \file test/src/flow_simulator/boundary_condition_tests.h
 * \brief Contains regression tests of \c StaticCapillaryNetworkAlgorithm class methods.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 */
#ifndef TEST_SRC_FLOW_SIMULATOR_BOUNDARY_CONDITION_TESTS_H_
#define TEST_SRC_FLOW_SIMULATOR_BOUNDARY_CONDITION_TESTS_H_

#include <gtest/gtest.h>
#include <map>
#include <string>
#include <utility>
#include "src/exec_manager/fluid_json.h"
#include "test/src/utils/flow_simulator_test_utils.h"
#include "src/flow_simulator/algorithms/static_capillary_network/static_capillary_network.h"

using StaticCapillaryNetworkAlgorithm = simulator::StaticCapillaryNetworkAlgorithm;

class StaticCapillaryNetwork_BoundaryCondition: public ::testing::Test {
 public:
  StaticCapillaryNetwork_BoundaryCondition() {
    // initialization code here
  }

  void SetUp() {
    // code here will execute just before the test ensues

    // Define configuration parameters

    // Defining directory where the will be the JSON file centerline.json
    simulation_cfg.folder = "test/results";

    // Inputs for the algorithm execution
    simulation_cfg.shape = {1, 1, 3};  // Default shape
    simulation_cfg.voxel_size = 1.0e-1;
    simulation_cfg.experiment_json.SetFlowAxis(2U);
    simulation_cfg.experiment_json.SetBoundaryThickness(1U);

    simulation_cfg.algorithm_json.SetName("static");
    simulation_cfg.algorithm_json.SetModel("hydrostatic");

    simulation_cfg.experiment_json.SetAbsolutePressure(100.0);
    simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient", 50.0);  // Default bc
    simulation_cfg.experiment_json.SetTemperature(340.0);

    // Defining fluids
    FluidJSON fluid1("water", "constant");

    fluid1.AddProperty("dynamic_viscosity", 1.0);
    fluid1.AddProperty("density", 1000.0);

    simulation_cfg.fluids_json.push_back(fluid1);
  }

  void TearDown() {
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be
  }

  ~StaticCapillaryNetwork_BoundaryCondition()  {
    // cleanup any pending stuff, but no exceptions allowed
  }

  // put in any custom data members that you need
  double gravity = 9.80665;
  double pressure_gradient = 50.0;

  SimulationConfig simulation_cfg;
};

TEST_F(StaticCapillaryNetwork_BoundaryCondition, HydrostaticTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/2_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {1, 1, 3};

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  StaticCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  algorithm.CalculateDerivedQuantities();
  arma::vec result = algorithm.GetPressures();

  // Assert
  EXPECT_DOUBLE_EQ(result(0), p_abs + 1.0 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(result(1), p_abs + 0.0 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(result(2), p_abs - 1.0 * bc.second * voxel_size);
}

TEST_F(StaticCapillaryNetwork_BoundaryCondition, PoiseuilleTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/2_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {1, 1, 3};

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  StaticCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  algorithm.CalculateDerivedQuantities();
  arma::vec result = algorithm.GetPressures();

  // Assert
  EXPECT_DOUBLE_EQ(result(0), p_abs + 1.0 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(result(1), p_abs + 0.0 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(result(2), p_abs - 1.0 * bc.second * voxel_size);
}

TEST_F(StaticCapillaryNetwork_BoundaryCondition, TwoCapillariesThreeVoxelsAwayFromEdgesTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/2_capillaries_3_voxels_away_from_edges.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {1, 1, 9};
  simulation_cfg.experiment_json.SetBoundaryThickness(4U);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  StaticCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  algorithm.CalculateDerivedQuantities();
  arma::vec result = algorithm.GetPressures();

  // Assert
  EXPECT_DOUBLE_EQ(result(0), p_abs + 1.0 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(result(1), p_abs + 0.0 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(result(2), p_abs - 1.0 * bc.second * voxel_size);
}

#endif  // TEST_SRC_FLOW_SIMULATOR_BOUNDARY_CONDITION_TESTS_H_
