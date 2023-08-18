/**
 * \file test/src/flow_simulator/hydrostatic_tests.h
 * \brief Contains regression tests of \c StaticCapillaryNetworkAlgorithm class methods using
 * \c Hydrostaticmodel, \c HydrostaticGeometry and \c Fluid.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 */

#ifndef TEST_SRC_FLOW_SIMULATOR_HYDROSTATIC_TESTS_H_
#define TEST_SRC_FLOW_SIMULATOR_HYDROSTATIC_TESTS_H_

#include <gtest/gtest.h>
#include <map>
#include <string>
#include <utility>
#include "src/flow_simulator/algorithms/static_capillary_network/static_capillary_network.h"
#include "test/src/utils/flow_simulator_test_utils.h"
#include "src/exec_manager/simulation_config.h"
#include "src/exec_manager/fluid_json.h"

using StaticCapillaryNetworkAlgorithm = simulator::StaticCapillaryNetworkAlgorithm;

class StaticCapillaryNetwork_Hydrostatic: public ::testing::Test {
 public:
  StaticCapillaryNetwork_Hydrostatic() {
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
    simulation_cfg.experiment_json.SetFlowAxis(0U);
    simulation_cfg.experiment_json.SetBoundaryThickness(1U);

    simulation_cfg.algorithm_json.SetName("static");
    simulation_cfg.algorithm_json.SetModel("hydrostatic");
    simulation_cfg.experiment_json.SetAbsolutePressure(100.0);
    simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient", 50.0);  // Default bc
    simulation_cfg.experiment_json.SetTemperature(340.0);

    // Defining fluids
    FluidJSON fluid1("water", "constant");

    fluid1.AddProperty("dynamic_viscosity", 1.0e-3);
    fluid1.AddProperty("density", 1000.0);

    simulation_cfg.fluids_json.push_back(fluid1);
  }

  void TearDown() {
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be
  }

  ~StaticCapillaryNetwork_Hydrostatic()  {
    // cleanup any pending stuff, but no exceptions allowed
  }

  // put in any custom data members that you need
  double gravity = 9.80665;
  double pressure_gradient = 50.0;

  SimulationConfig simulation_cfg;
};

TEST_F(StaticCapillaryNetwork_Hydrostatic, TripleCrossEqualRadius_PressureGradientTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/triple_cross_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {7, 7, 7};
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient", 50.0);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  double imposed_flow_rate = 19634.954084935922 * std::pow(voxel_size, 4.0);

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
  constexpr double tolerance = 1.0e-10;
  const double pressure_drop = pressure_gradient * voxel_size;
  const double hydrostatic_pressure = simulation_cfg.fluids_json[0].GetProperty("density")
                                    * gravity
                                    * voxel_size;

  // Along X
  double Q_x = imposed_flow_rate;
  double V_x = Q_x / arma::datum::pi / std::pow(voxel_size, 2.0);
  arma::vec flow_rate_x = flow_rate.subvec(6, 11);
  arma::vec flow_speed_x = flow_speed.subvec(6, 11);
  EXPECT_NEAR(pressure(7),  p_abs + 3.0 * pressure_drop, tolerance);
  EXPECT_NEAR(pressure(8),  p_abs + 2.0 * pressure_drop, tolerance);
  EXPECT_NEAR(pressure(9),  p_abs + 1.0 * pressure_drop, tolerance);
  EXPECT_NEAR(pressure(3),  p_abs + 0.0 * pressure_drop, tolerance);
  EXPECT_NEAR(pressure(10), p_abs - 1.0 * pressure_drop, tolerance);
  EXPECT_NEAR(pressure(11), p_abs - 2.0 * pressure_drop, tolerance);
  EXPECT_NEAR(pressure(12), p_abs - 3.0 * pressure_drop, tolerance);
  for (const auto &q_x : flow_rate_x) EXPECT_NEAR(std::abs(q_x), Q_x, tolerance);
  for (const auto &v_x : flow_speed_x) EXPECT_NEAR(std::abs(v_x), V_x, tolerance);

  // Along Y
  double Q_y = 0.0;
  double V_y = 0.0;
  arma::vec pressures_y = pressure.subvec(13, 18);
  arma::vec flow_rate_y = flow_rate.subvec(12, 17);
  arma::vec flow_speed_y = flow_speed.subvec(12, 17);
  for (const auto &p_y : pressures_y) EXPECT_NEAR(p_y, p_abs, tolerance);
  for (const auto &q_y : flow_rate_y) EXPECT_NEAR(std::abs(q_y), Q_y, tolerance);
  for (const auto &v_y : flow_speed_y) EXPECT_NEAR(std::abs(v_y), V_y, tolerance);


  // Along Z
  // TODO(rneumann): check scaling with voxel_size
  double magic = 2.2737367544323206e-12;
  double Q_z = magic * std::pow(voxel_size, 4.0);
  double V_z = Q_z / arma::datum::pi / std::pow(voxel_size, 2.0);
  arma::vec flow_rate_z = flow_rate.subvec(0, 5);
  arma::vec flow_speed_z = flow_speed.subvec(0, 5);
  EXPECT_NEAR(pressure(0), p_abs + 3.0 * hydrostatic_pressure, tolerance);
  EXPECT_NEAR(pressure(1), p_abs + 2.0 * hydrostatic_pressure, tolerance);
  EXPECT_NEAR(pressure(2), p_abs + 1.0 * hydrostatic_pressure, tolerance);
  EXPECT_NEAR(pressure(3), p_abs + 0.0 * hydrostatic_pressure, tolerance);
  EXPECT_NEAR(pressure(4), p_abs - 1.0 * hydrostatic_pressure, tolerance);
  EXPECT_NEAR(pressure(5), p_abs - 2.0 * hydrostatic_pressure, tolerance);
  EXPECT_NEAR(pressure(6), p_abs - 3.0 * hydrostatic_pressure, tolerance);
  for (const auto &q_z : flow_rate_z) EXPECT_NEAR(std::abs(q_z), Q_z, tolerance);
  for (const auto &v_z : flow_speed_z) EXPECT_NEAR(std::abs(v_z), V_z, tolerance);
}

TEST_F(StaticCapillaryNetwork_Hydrostatic, TripleCrossEqualRadius_FlowRateClosedTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/triple_cross_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  double imposed_flow_rate = 19634.954084935922 * std::pow(simulation_cfg.voxel_size, 4.0);

  simulation_cfg.shape = {7, 7, 7};
  simulation_cfg.experiment_json.SetBoundaryCondition("flow_rate_closed", imposed_flow_rate);

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
  arma::vec flow_speed = algorithm.GetFlowSpeed();

  // Assert
  constexpr double tolerance = 1.0e-10;
  const double pressure_drop = pressure_gradient * voxel_size;
  const double hydrostatic_pressure = simulation_cfg.fluids_json[0].GetProperty("density")
                                    * gravity
                                    * voxel_size;

  // Along X
  double Q_x = imposed_flow_rate;
  double V_x = Q_x / arma::datum::pi / std::pow(voxel_size, 2.0);
  arma::vec flow_rate_x = flow_rate.subvec(6, 11);
  arma::vec flow_speed_x = flow_speed.subvec(6, 11);
  EXPECT_NEAR(pressure(7),  p_abs + 6.0 * pressure_drop, tolerance);
  EXPECT_NEAR(pressure(8),  p_abs + 5.0 * pressure_drop, tolerance);
  EXPECT_NEAR(pressure(9),  p_abs + 4.0 * pressure_drop, tolerance);
  EXPECT_NEAR(pressure(3),  p_abs + 3.0 * pressure_drop, tolerance);
  EXPECT_NEAR(pressure(10), p_abs + 2.0 * pressure_drop, tolerance);
  EXPECT_NEAR(pressure(11), p_abs + 1.0 * pressure_drop, tolerance);
  EXPECT_NEAR(pressure(12), p_abs + 0.0 * pressure_drop, tolerance);
  for (const auto &q_x : flow_rate_x) EXPECT_NEAR(std::abs(q_x), Q_x, tolerance);
  for (const auto &v_x : flow_speed_x) EXPECT_NEAR(std::abs(v_x), V_x, tolerance);

  // Along Y
  double Q_y = 0.0;
  double V_y = 0.0;
  arma::vec pressures_y = pressure.subvec(13, 18);
  arma::vec flow_rate_y = flow_rate.subvec(12, 17);
  arma::vec flow_speed_y = flow_speed.subvec(12, 17);
  for (const auto &p_y : pressures_y) EXPECT_NEAR(p_y, p_abs + 3.0 * pressure_drop, tolerance);
  for (const auto &q_y : flow_rate_y) EXPECT_NEAR(std::abs(q_y), Q_y, tolerance);
  for (const auto &v_y : flow_speed_y) EXPECT_NEAR(std::abs(v_y), V_y, tolerance);

  // Along Z
  double Q_z = 0.0;
  double V_z = 0.0;
  arma::vec flow_rate_z = flow_rate.subvec(0, 5);
  arma::vec flow_speed_z = flow_speed.subvec(0, 5);
  EXPECT_NEAR(pressure(0), p_abs + 3.0 * pressure_drop + 3.0 * hydrostatic_pressure, tolerance);
  EXPECT_NEAR(pressure(1), p_abs + 3.0 * pressure_drop + 2.0 * hydrostatic_pressure, tolerance);
  EXPECT_NEAR(pressure(2), p_abs + 3.0 * pressure_drop + 1.0 * hydrostatic_pressure, tolerance);
  EXPECT_NEAR(pressure(3), p_abs + 3.0 * pressure_drop + 0.0 * hydrostatic_pressure, tolerance);
  EXPECT_NEAR(pressure(4), p_abs + 3.0 * pressure_drop - 1.0 * hydrostatic_pressure, tolerance);
  EXPECT_NEAR(pressure(5), p_abs + 3.0 * pressure_drop - 2.0 * hydrostatic_pressure, tolerance);
  EXPECT_NEAR(pressure(6), p_abs + 3.0 * pressure_drop - 3.0 * hydrostatic_pressure, tolerance);
  for (const auto &q_z : flow_rate_z) EXPECT_NEAR(std::abs(q_z), Q_z, tolerance);
  for (const auto &v_z : flow_speed_z) EXPECT_NEAR(std::abs(v_z), V_z, tolerance);
}

TEST_F(StaticCapillaryNetwork_Hydrostatic, TripleCrossEqualRadius_FlowRateOpenTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/triple_cross_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  double imposed_flow_rate = 19634.954084935922 * std::pow(simulation_cfg.voxel_size, 4.0);

  simulation_cfg.shape = {7, 7, 7};
  simulation_cfg.experiment_json.SetBoundaryCondition("flow_rate_open", imposed_flow_rate);

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
  arma::vec flow_speed = algorithm.GetFlowSpeed();

  // Assert
  constexpr double tolerance = 1.0e-10;
  const double gradient_per_outlet = pressure_gradient / 5.0;
  const double flow_rate_per_outlet = imposed_flow_rate / 5.0;
  const double pressure_drop = pressure_gradient * voxel_size;
  const double pressure_drop_per_outlet = gradient_per_outlet * voxel_size;
  const double hydrostatic_pressure = simulation_cfg.fluids_json[0].GetProperty("density")
                                    * gravity
                                    * voxel_size;

  // Inlet along X
  double Q_inlet = imposed_flow_rate;
  double V_inlet = Q_inlet / arma::datum::pi / std::pow(voxel_size, 2.0);
  arma::vec flow_rate_inlet = flow_rate.subvec(6, 8);
  arma::vec flow_speed_inlet = flow_speed.subvec(6, 8);
  EXPECT_NEAR(pressure(7), p_abs + 3.0 * pressure_drop_per_outlet + 3.0 * pressure_drop, tolerance);
  EXPECT_NEAR(pressure(8), p_abs + 3.0 * pressure_drop_per_outlet + 2.0 * pressure_drop, tolerance);
  EXPECT_NEAR(pressure(9), p_abs + 3.0 * pressure_drop_per_outlet + 1.0 * pressure_drop, tolerance);
  for (const auto &q_inlet : flow_rate_inlet) EXPECT_NEAR(std::abs(q_inlet), Q_inlet, tolerance);
  for (const auto &v_inlet : flow_speed_inlet) EXPECT_NEAR(std::abs(v_inlet), V_inlet, tolerance);

  // Outlets in the XY plane
  double Q_plane = flow_rate_per_outlet;
  double V_plane = Q_plane / arma::datum::pi / std::pow(voxel_size, 2.0);
  arma::vec flow_rate_plane = flow_rate.subvec(9, 17);
  arma::vec flow_speed_plane = flow_speed.subvec(9, 17);
  EXPECT_NEAR(pressure(3),  p_abs + 3.0 * pressure_drop_per_outlet, tolerance);  // centre
  EXPECT_NEAR(pressure(10), p_abs + 2.0 * pressure_drop_per_outlet, tolerance);  // centre -> outlet
  EXPECT_NEAR(pressure(11), p_abs + 1.0 * pressure_drop_per_outlet, tolerance);  // centre -> outlet
  EXPECT_NEAR(pressure(12), p_abs + 0.0 * pressure_drop_per_outlet, tolerance);  // outlet
  EXPECT_NEAR(pressure(13), p_abs + 0.0 * pressure_drop_per_outlet, tolerance);  // outlet
  EXPECT_NEAR(pressure(14), p_abs + 1.0 * pressure_drop_per_outlet, tolerance);  // centre -> outlet
  EXPECT_NEAR(pressure(15), p_abs + 2.0 * pressure_drop_per_outlet, tolerance);  // centre -> outlet
  EXPECT_NEAR(pressure(3),  p_abs + 3.0 * pressure_drop_per_outlet, tolerance);  // centre
  EXPECT_NEAR(pressure(16), p_abs + 2.0 * pressure_drop_per_outlet, tolerance);  // centre -> outlet
  EXPECT_NEAR(pressure(17), p_abs + 1.0 * pressure_drop_per_outlet, tolerance);  // centre -> outlet
  EXPECT_NEAR(pressure(18), p_abs + 0.0 * pressure_drop_per_outlet, tolerance);  // outlet
  for (const auto &q_plane : flow_rate_plane) EXPECT_NEAR(std::abs(q_plane), Q_plane, tolerance);
  for (const auto &v_plane : flow_speed_plane) EXPECT_NEAR(std::abs(v_plane), V_plane, tolerance);

  // Downward outlet along Z
  double Q_down = arma::datum::pi / 8.0 * hydrostatic_pressure + flow_rate_per_outlet;
  double V_down = Q_down / arma::datum::pi / std::pow(voxel_size, 2.0);
  arma::vec flow_rate_down = flow_rate.subvec(0, 2);
  arma::vec flow_speed_down = flow_speed.subvec(0, 2);
  EXPECT_NEAR(pressure(0), p_abs + 0.0 * pressure_drop_per_outlet, tolerance);  // outlet
  EXPECT_NEAR(pressure(1), p_abs + 1.0 * pressure_drop_per_outlet, tolerance);  // centre -> outlet
  EXPECT_NEAR(pressure(2), p_abs + 2.0 * pressure_drop_per_outlet, tolerance);  // centre -> outlet
  EXPECT_NEAR(pressure(3), p_abs + 3.0 * pressure_drop_per_outlet, tolerance);  // centre
  EXPECT_NEAR(pressure(4), p_abs + 2.0 * pressure_drop_per_outlet, tolerance);  // centre -> outlet
  EXPECT_NEAR(pressure(5), p_abs + 1.0 * pressure_drop_per_outlet, tolerance);  // centre -> outlet
  EXPECT_NEAR(pressure(6), p_abs + 0.0 * pressure_drop_per_outlet, tolerance);  // outlet
  for (const auto &q_down : flow_rate_down) EXPECT_NEAR(std::abs(q_down), Q_down, tolerance);
  for (const auto &v_down : flow_speed_down) EXPECT_NEAR(std::abs(v_down), V_down, tolerance);

  // Upward outlet along Z
  double Q_up = arma::datum::pi / 8.0 * hydrostatic_pressure - flow_rate_per_outlet;
  double V_up = Q_up / arma::datum::pi / std::pow(voxel_size, 2.0);
  arma::vec flow_rate_up = flow_rate.subvec(3, 5);
  arma::vec flow_speed_up = flow_speed.subvec(3, 5);
  EXPECT_NEAR(pressure(0), p_abs + 0.0 * pressure_drop_per_outlet, tolerance);  // outlet
  EXPECT_NEAR(pressure(1), p_abs + 1.0 * pressure_drop_per_outlet, tolerance);  // centre -> outlet
  EXPECT_NEAR(pressure(2), p_abs + 2.0 * pressure_drop_per_outlet, tolerance);  // centre -> outlet
  EXPECT_NEAR(pressure(3), p_abs + 3.0 * pressure_drop_per_outlet, tolerance);  // centre
  EXPECT_NEAR(pressure(4), p_abs + 2.0 * pressure_drop_per_outlet, tolerance);  // centre -> outlet
  EXPECT_NEAR(pressure(5), p_abs + 1.0 * pressure_drop_per_outlet, tolerance);  // centre -> outlet
  EXPECT_NEAR(pressure(6), p_abs + 0.0 * pressure_drop_per_outlet, tolerance);  // outlet
  for (const auto &q_up : flow_rate_up) EXPECT_NEAR(std::abs(q_up), Q_up, tolerance);
  for (const auto &v_up : flow_speed_up) EXPECT_NEAR(std::abs(v_up), V_up, tolerance);
}

#endif  // TEST_SRC_FLOW_SIMULATOR_HYDROSTATIC_TESTS_H_
