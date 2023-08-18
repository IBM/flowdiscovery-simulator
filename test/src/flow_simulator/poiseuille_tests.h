/**
 * \file test/src/flow_simulator/poiseuille_tests.h
 * \brief Contains regression tests of \c StaticCapillaryNetworkAlgorithm class methods using
 * \c PoiseuilleModel, \c PoiseuilleGeometry and \c Fluid.
 *
 * \copyright © IBM Corp.
 * \date 2016

 * \authors Alexandre Ashade Lassance Cunha \<aashade@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 */

#ifndef TEST_SRC_FLOW_SIMULATOR_POISEUILLE_TESTS_H_
#define TEST_SRC_FLOW_SIMULATOR_POISEUILLE_TESTS_H_

#include <gtest/gtest.h>
#include <armadillo>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include "src/flow_simulator/algorithms/static_capillary_network/static_capillary_network.h"
#include "test/src/utils/flow_simulator_test_utils.h"
#include "src/exec_manager/simulation_config.h"
#include "src/exec_manager/fluid_json.h"

using StaticCapillaryNetworkAlgorithm = simulator::StaticCapillaryNetworkAlgorithm;

class StaticCapillaryNetwork_Poiseuille: public ::testing::Test {
 public:
  StaticCapillaryNetwork_Poiseuille() {
    // initialization code here
  }

  void SetUp() {
    // code here will execute just before the test ensues

    // Define configuration parameters

    // Defining directory where the will be the JSON file centerline.json
    simulation_cfg.folder = "test/results";

    // Inputs for the algorithm execution
    simulation_cfg.shape = {1, 1, 3};  // Default shape
    simulation_cfg.voxel_size = 1.0e-6;
    simulation_cfg.experiment_json.SetFlowAxis(2U);
    simulation_cfg.experiment_json.SetBoundaryThickness(1U);

    simulation_cfg.algorithm_json.SetName("static");
    simulation_cfg.algorithm_json.SetModel("poiseuille");

    simulation_cfg.experiment_json.SetAbsolutePressure(100.0);
    // Default bc
    simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient", 50.0);
    simulation_cfg.experiment_json.SetTemperature(340.0);

    // Defining fluids
    std::vector<FluidJSON> my_fluids;
    FluidJSON fluid1("water", "constant");

    fluid1.AddProperty("dynamic_viscosity", 1.0e-3);

    simulation_cfg.fluids_json.push_back(fluid1);
  }

  void TearDown() {
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be
  }

  ~StaticCapillaryNetwork_Poiseuille()  {
    // cleanup any pending stuff, but no exceptions allowed
  }

  // put in any custom data members that you need
  SimulationConfig simulation_cfg;
};

TEST_F(StaticCapillaryNetwork_Poiseuille, DerivedQuantitiesForTwoCapillariesEqualRadiusTest) {
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
  arma::vec result = {p_abs + 1.0 * voxel_size * bc.second,
                      p_abs + 0.0 * voxel_size * bc.second,
                      p_abs - 1.0 * voxel_size * bc.second};
  algorithm.SetPressures(result);
  algorithm.CalculateDerivedQuantities();

  // Assert
  double tolerance = arma::datum::eps;
  arma::vec flow_rate = algorithm.GetFlowRate();
  arma::vec flow_speed = algorithm.GetFlowSpeed();
  double magic = 19634.954084936206 * std::pow(voxel_size, 4.0);
  EXPECT_NEAR(flow_rate(0), magic, tolerance);
  EXPECT_NEAR(flow_rate(1), magic, tolerance);
  EXPECT_DOUBLE_EQ(flow_speed(0), flow_rate(0) / arma::datum::pi / std::pow(voxel_size, 2.0));
  EXPECT_DOUBLE_EQ(flow_speed(1), flow_rate(1) / arma::datum::pi / std::pow(voxel_size, 2.0));
  EXPECT_DOUBLE_EQ(algorithm.GetPermeability(), 5.8904862254808643e-13);
}

TEST_F(StaticCapillaryNetwork_Poiseuille, TShapedCapillariesEqualRadiusTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/t-shaped_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {1, 4, 3};

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
  arma::vec pressure = algorithm.GetPressures();
  arma::vec flow_rate = algorithm.GetFlowRate();

  // Assert
  double tolerance = arma::datum::eps;
  EXPECT_DOUBLE_EQ(pressure(0), p_abs + 0.0 * voxel_size * bc.second);
  EXPECT_NEAR(pressure(1), p_abs + 0.0 * voxel_size * bc.second, tolerance);
  EXPECT_DOUBLE_EQ(pressure(2), p_abs + 0.0 * voxel_size * bc.second);
  EXPECT_DOUBLE_EQ(pressure(3), p_abs + 0.0 * voxel_size * bc.second);
  EXPECT_DOUBLE_EQ(pressure(4), p_abs + 1.0 * voxel_size * bc.second);
  EXPECT_DOUBLE_EQ(pressure(5), p_abs - 1.0 * voxel_size * bc.second);
  EXPECT_NEAR(flow_rate(0), 0.0, tolerance);
  EXPECT_NEAR(flow_rate(1), 0.0, tolerance);
  EXPECT_NEAR(flow_rate(2), 0.0, tolerance);
  EXPECT_NEAR(flow_rate(3), -flow_rate(4), tolerance);
}

TEST_F(StaticCapillaryNetwork_Poiseuille, ThreeCapillariesEqualRadiusTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/3_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {1, 1, 4};

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  // Define configuration parameters
  double p_grad = 50.0;

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
  EXPECT_DOUBLE_EQ(result(0), p_abs + 3.0 * voxel_size * p_grad / 2.0);
  EXPECT_DOUBLE_EQ(result(1), p_abs + 1.0 * voxel_size * p_grad / 2.0);
  EXPECT_DOUBLE_EQ(result(2), p_abs - 1.0 * voxel_size * p_grad / 2.0);
  EXPECT_DOUBLE_EQ(result(3), p_abs - 3.0 * voxel_size * p_grad / 2.0);
}

TEST_F(StaticCapillaryNetwork_Poiseuille, YShapedCapillariesEqualRadiusTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/y-shaped_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {3, 1, 3};

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
  double Pin =  p_abs + 1.0 * voxel_size * bc.second;
  double Pout = p_abs - 1.0 * voxel_size * bc.second;
  double expected_p = (2 * Pout + std::sqrt(2.0) * Pin) / (2 + std::sqrt(2.0));
  EXPECT_DOUBLE_EQ(result(0), Pin);
  EXPECT_DOUBLE_EQ(result(1), expected_p);
  EXPECT_DOUBLE_EQ(result(2), Pout);
  EXPECT_DOUBLE_EQ(result(3), Pout);
}

TEST_F(StaticCapillaryNetwork_Poiseuille, ParallelCapillariesEqualRadiusTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/parallel_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {1, 4, 3};

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  // Define configuration parameters
  double p_grad = 50.0;

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
  EXPECT_DOUBLE_EQ(result(0), p_abs + 1.0 * p_grad * voxel_size);
  EXPECT_DOUBLE_EQ(result(1), p_abs + 0.0 * p_grad * voxel_size);
  EXPECT_DOUBLE_EQ(result(2), p_abs - 1.0 * p_grad * voxel_size);

  EXPECT_DOUBLE_EQ(result(3), p_abs + 1.0 * p_grad * voxel_size);
  EXPECT_DOUBLE_EQ(result(4), p_abs + 0.0 * p_grad * voxel_size);
  EXPECT_DOUBLE_EQ(result(5), p_abs - 1.0 * p_grad * voxel_size);
}

TEST_F(StaticCapillaryNetwork_Poiseuille, TwoCapillariesDifferentRadiusTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/2_capillaries_different_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {1, 1, 3};

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  // Define configuration parameters
  double p_grad = 50.0;

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
  double Pin  = p_abs + 1.0 * p_grad * voxel_size;
  double Pout = p_abs - 1.0 * p_grad * voxel_size;
  double r1_4 = std::pow(1.0, 2);
  double r2_4 = std::pow(1.3719886811400708, 2);
  double expected_p = (r1_4 * Pin + r2_4 * Pout) / (r1_4 + r2_4);
  EXPECT_DOUBLE_EQ(result(0), Pin);
  EXPECT_DOUBLE_EQ(result(1), expected_p);
  EXPECT_DOUBLE_EQ(result(2), Pout);
}

TEST_F(StaticCapillaryNetwork_Poiseuille, FlowRateClosedPressureTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/triple_cross_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {7, 7, 7};
  simulation_cfg.experiment_json.SetAbsolutePressure(0.0);
  simulation_cfg.experiment_json.SetBoundaryCondition("flow_rate_closed", 1.0e-12);

  // Local variables
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
  arma::vec pressures = algorithm.GetPressures();
  arma::vec flow_rate = algorithm.GetFlowRate();

  // Assert
  double tolerance = arma::datum::eps;
  double central_p = 7639.4372684109785;

  // Pressures along flow axis
  for (auto i = 0U; i != 7U; ++i) EXPECT_DOUBLE_EQ(pressures(i), central_p * (6 - i) / 3.0);

  // Pressures everywhere else
  arma::vec constant_p = pressures.subvec(7, 18);
  for (const auto &p : constant_p) EXPECT_DOUBLE_EQ(p, central_p);

  // Flow rates along flow axis
  arma::vec head = flow_rate.head(6);
  for (const auto &q : head) EXPECT_DOUBLE_EQ(q, bc.second);

  // Flow rates everywhere else
  arma::vec tail = flow_rate.tail(12);
  for (const auto &q : tail) EXPECT_NEAR(q, 0.0, tolerance);
}

TEST_F(StaticCapillaryNetwork_Poiseuille, FlowRateClosedContinuityTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/t-shaped_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {1, 4, 3};
  simulation_cfg.experiment_json.SetAbsolutePressure(0.0);
  simulation_cfg.experiment_json.SetBoundaryCondition("flow_rate_closed", 1.0e-12);
  simulation_cfg.experiment_json.SetFlowAxis(1U);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();

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
  arma::vec pressure = algorithm.GetPressures();
  arma::vec flow_rate = algorithm.GetFlowRate();

  // Assert
  EXPECT_TRUE(pressure(0) > pressure(1));
  EXPECT_TRUE(pressure(1) > pressure(2));
  EXPECT_TRUE(pressure(2) > pressure(3));
  EXPECT_DOUBLE_EQ(pressure(3), p_abs);
  EXPECT_DOUBLE_EQ(pressure(4), p_abs);
  EXPECT_DOUBLE_EQ(pressure(5), p_abs);
  EXPECT_DOUBLE_EQ(flow_rate(0), bc.second);
  EXPECT_DOUBLE_EQ(flow_rate(1), bc.second);
  EXPECT_DOUBLE_EQ(flow_rate(2), bc.second);
  EXPECT_DOUBLE_EQ(flow_rate(3), 0.0);
  EXPECT_DOUBLE_EQ(flow_rate(4), 0.0);
}

TEST_F(StaticCapillaryNetwork_Poiseuille, CalculateLocalFlowRateTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/y-shaped_capillaries_different_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {3, 1, 3};
  simulation_cfg.experiment_json.SetAbsolutePressure(0.0);
  simulation_cfg.experiment_json.SetBoundaryCondition("flow_rate_closed", 1.0e-12);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();

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
  arma::vec pressure = algorithm.GetPressures();
  arma::vec flow_rate = algorithm.GetFlowRate();
  arma::vec flow_speed = algorithm.GetFlowSpeed();

  // Assert

  // Inlet pressures must be greater than outlet pressure
  EXPECT_TRUE(pressure(0) > pressure(3));
  EXPECT_TRUE(pressure(1) > pressure(3));
  EXPECT_DOUBLE_EQ(pressure(3), p_abs);
  EXPECT_DOUBLE_EQ(flow_rate(2), bc.second);
  EXPECT_DOUBLE_EQ(flow_rate(0) + flow_rate(1), bc.second);

  // Inlet capillaries have different diameters
  EXPECT_TRUE(pressure(0) > pressure(1));
  EXPECT_DOUBLE_EQ(flow_rate(0), 3.0 * flow_rate(1));
  EXPECT_DOUBLE_EQ(flow_speed(0) * 1.3416407864998738 / 3.0, flow_speed(1));
}

TEST_F(StaticCapillaryNetwork_Poiseuille, FlowRateOpenContinuityTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/triple_cross_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {7, 7, 7};
  simulation_cfg.experiment_json.SetAbsolutePressure(0.0);
  simulation_cfg.experiment_json.SetBoundaryCondition("flow_rate_open", 1.0e-12);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();

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
  arma::vec pressures = algorithm.GetPressures();
  arma::vec flow_rate = algorithm.GetFlowRate();

  // Assert
  double tolerance = arma::datum::eps;

  // Inlet pressure must be highest than all other pressures
  EXPECT_TRUE(arma::all(pressures <= pressures(0)));

  // Pressure at outlet nodes must be equal to absolute pressure
  EXPECT_DOUBLE_EQ(pressures(6), p_abs);
  EXPECT_DOUBLE_EQ(pressures(7), p_abs);
  EXPECT_DOUBLE_EQ(pressures(12), p_abs);
  EXPECT_DOUBLE_EQ(pressures(13), p_abs);
  EXPECT_DOUBLE_EQ(pressures(18), p_abs);

  // Pressure variation from the centre to the sides must be symmetric
  EXPECT_DOUBLE_EQ(pressures(8), pressures(11));
  EXPECT_DOUBLE_EQ(pressures(8), pressures(14));
  EXPECT_DOUBLE_EQ(pressures(8), pressures(17));
  EXPECT_DOUBLE_EQ(pressures(9), pressures(10));
  EXPECT_DOUBLE_EQ(pressures(9), pressures(15));
  EXPECT_DOUBLE_EQ(pressures(9), pressures(16));

  // All inlet capillaries flow rate must be equal
  arma::vec inlet_flow_rate = flow_rate.subvec(0, 2);
  for (const auto &q : inlet_flow_rate) EXPECT_NEAR(std::abs(q), bc.second, tolerance);

  // All outlet capillaries flow rate must be equal
  arma::vec outlet_flow_rate = flow_rate.subvec(3, 17);
  for (const auto &q : outlet_flow_rate) EXPECT_NEAR(std::abs(q), bc.second / 5.0, tolerance);

  // Sum of outlet capillary flow must be equal to inlet flow
  EXPECT_NEAR(std::abs(flow_rate(5)) +
              std::abs(flow_rate(6)) +
              std::abs(flow_rate(11)) +
              std::abs(flow_rate(12)) +
              std::abs(flow_rate(17)),
              bc.second,
              tolerance);
}

TEST_F(StaticCapillaryNetwork_Poiseuille, SaveResultsToDiskTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/2_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {1, 1, 3};

  // Prepare algorithm
  StaticCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  algorithm.CalculateDerivedQuantities();
  algorithm.SaveResultsToDisk();

  constexpr double tolerance = 1.0e-10;

  // Assert
  const std::string expected_results_file = simulation_cfg.folder + "/static_results.h5";

  EXPECT_TRUE(flow_simulator_test_utils::FileExists(expected_results_file));

  arma::vec pressures;
  const arma::vec pressures_expected = {100.00005, 100.0, 99.99995};
  EXPECT_TRUE(pressures.load(arma::hdf5_name(expected_results_file, "pressures_z")));
  EXPECT_TRUE(arma::approx_equal(pressures, pressures_expected, "absdiff", tolerance));

  arma::vec flow_rate;
  const arma::vec flow_rate_expected = {1.9634954085588019e-20, 1.9634954085588019e-20};
  EXPECT_TRUE(flow_rate.load(arma::hdf5_name(expected_results_file, "flow_rate_z")));
  EXPECT_TRUE(arma::approx_equal(flow_rate, flow_rate_expected, "absdiff", tolerance));

  arma::vec flow_speed;
  const arma::vec flow_speed_expected = {6.2500000002074781e-09, 6.2500000002074781e-09};
  EXPECT_TRUE(flow_speed.load(arma::hdf5_name(expected_results_file, "flow_speed_z")));
  EXPECT_TRUE(arma::approx_equal(flow_speed, flow_speed_expected, "absdiff", tolerance));

  arma::vec permeability;
  const arma::vec permeability_expected = {5.8904862254808633e-13};
  EXPECT_TRUE(permeability.load(arma::hdf5_name(expected_results_file, "permeability_z")));
  EXPECT_TRUE(arma::approx_equal(permeability, permeability_expected, "absdiff", tolerance));

  // Clean up
  std::remove((simulation_cfg.folder + "/static_results.h5").c_str());
}

#endif  // TEST_SRC_FLOW_SIMULATOR_POISEUILLE_TESTS_H_
