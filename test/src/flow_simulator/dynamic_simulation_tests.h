/**
 * \file test/src/flow_simulator/dynamic_simulation_tests.h
 * \brief Contains regression tests of \c DynamicCapillaryNetworkAlgorithm class methods using
 * \c DynamicPoiseuillemodel, \c DynamicPoiseuilleGeometry and \c Fluid.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 */

#ifndef TEST_SRC_FLOW_SIMULATOR_DYNAMIC_SIMULATION_TESTS_H_
#define TEST_SRC_FLOW_SIMULATOR_DYNAMIC_SIMULATION_TESTS_H_

#include <gtest/gtest.h>
#include <armadillo>
#include <map>
#include <string>
#include <memory>
#include <limits>
#include <utility>
#include <algorithm>
#include <vector>
#include <filesystem>
#include "src/exec_manager/fluid_json.h"
#include "test/src/utils/flow_simulator_test_utils.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network.h"

using CapillaryInterface = simulator::CapillaryInterface;
using DynamicCapillaryNetworkAlgorithm = simulator::DynamicCapillaryNetworkAlgorithm;
using Side = simulator::Side;
using DynamicCapillaryNetworkContext = simulator::DynamicCapillaryNetworkContext;
using IndexType = simulator::IndexType;

class DynamicCapillaryNetwork: public ::testing::Test {
 public:
  DynamicCapillaryNetwork() {
    // initialization code here
  }

  void SetUp() {
    // code here will execute just before the test ensues

    // Define configuration parameters

    // Defining directory where the will be the JSON file centerline.json
    simulation_cfg.folder = "test/results";

    // Inputs for the algorithm execution
    simulation_cfg.shape = {1, 1, 3};  // Default shape
    simulation_cfg.voxel_size = 1.0;
    simulation_cfg.experiment_json.SetFlowAxis(2U);
    simulation_cfg.experiment_json.SetBoundaryThickness(1U);

    simulation_cfg.algorithm_json.SetName("dynamic");
    simulation_cfg.algorithm_json.SetModel("linear_molecular_kinetics");
    simulation_cfg.algorithm_json.SetInitialTime(0.0);
    simulation_cfg.algorithm_json.SetFinalTime(19.0);
    simulation_cfg.algorithm_json.SetTimeStepSize(0.59375);
    simulation_cfg.algorithm_json.SetRelativeTolerance(1.0e-4);
    simulation_cfg.algorithm_json.SetAbsoluteLinkTolerance(1.0e-3);
    simulation_cfg.algorithm_json.SetAbsoluteNodeTolerance(1.0e-3);
    simulation_cfg.algorithm_json.SetResume(false);

    simulation_cfg.wettability_json.AddProperty("contact_angle", 0.0);
    simulation_cfg.wettability_json.AddProperty("linear_mk", 108.63);

    simulation_cfg.fluid_interface_json.AddProperty("interfacial_tension", 0.06347);

    simulation_cfg.experiment_json.SetAbsolutePressure(32.5);
    // Default bc
    simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_open", 10);

    simulation_cfg.experiment_json.SetTemperature(340);

    // Defining fluids
    FluidJSON fluid1("water", "constant");
    FluidJSON fluid2("oil", "constant");

    fluid1.AddProperty("dynamic_viscosity", 1.0);
    fluid2.AddProperty("dynamic_viscosity", 1.0111);

    simulation_cfg.fluids_json.push_back(fluid1);
    simulation_cfg.fluids_json.push_back(fluid2);
  }

  void TearDown() {
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be

    // recreate file .gitignore in snapshots folder (deleted by the simulations)
    std::ofstream out(simulation_cfg.folder + "/snapshots/.gitignore");
    out << "*\n";
  }

  ~DynamicCapillaryNetwork()  {
    // cleanup any pending stuff, but no exceptions allowed
  }

  // put in any custom data members that you need
  SimulationConfig simulation_cfg;
};



bool NearlyEqual(double tested_value, double expected_value, double epsilon) {
/**
 * Compares the tested value to the expected value.
 * Returns true if \c tested_value is approximately equal to \c expected_value following the
 * implementation explained at: https://floating-point-gui.de/errors/comparison/
 *
 * \param[in] tested_value    The value to be tested.
 * \param[in] expected_value  The expected value.
 * \param[in] epsilon         The acceptable tolerance.
 * \retval    equal           A boolean value representing whether test_value ~= expected_value.
 */
  const double abs_difference = std::abs(tested_value - expected_value);

  bool equal;

  if (tested_value == expected_value) {
    // Shortcut to handle infinities
    equal = true;
  } else if (tested_value == 0.0 || expected_value == 0.0 || abs_difference <= arma::datum::eps) {
    // When tested_value and expected_value are both zero (or very close), use absolute error
    equal = abs_difference <= arma::datum::eps;

    if (!equal) {
      std::printf(" --- NearlyEqual --- \n");
      std::printf("The difference between tested_value and expected_value is ");
      std::printf("%.17g, which exceeds tolerance, where\n", abs_difference);
      std::printf("tested_value evaluates to %.17g,\n", tested_value);
      std::printf("expected_value evaluates to %.17g, and\n", expected_value);
      std::printf("tolerance evaluates to %.17g.\n", arma::datum::eps);
    }
  } else {
    // In the general case, use relative error
    equal =  abs_difference / std::abs(expected_value) <= epsilon;

    if (!equal) {
      std::printf(" --- NearlyEqual --- \n");
      std::printf("The relative difference between tested_value and expected_value is ");
      std::printf("%.17g, which exceeds tolerance, where\n", abs_difference);
      std::printf("tested_value evaluates to %.17g,\n", tested_value);
      std::printf("expected_value evaluates to %.17g, and\n", expected_value);
      std::printf("tolerance evaluates to %.17g.\n", epsilon);
    }
  }

  return equal;
}



TEST_F(DynamicCapillaryNetwork, TwoCapillariesDifferentRadiusTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/2_capillaries_different_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {1, 1, 3};
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_open", 10);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;  // dynamic
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  constexpr double tolerance = 1.0e-10;

  // Assert
  // These magic numbers are related to the parameters hardcoded by SetUp().
  double magic = 8.89777343751761407020;

  EXPECT_DOUBLE_EQ(pressures(0), p_abs + 1.0 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(1), p_abs + 0.0 * bc.second * voxel_size + magic);
  EXPECT_DOUBLE_EQ(pressures(2), p_abs - 1.0 * bc.second * voxel_size);

  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 1);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(0), 1.0);
  EXPECT_NEAR(interface->CalculateEffectiveInterfacePositionAtCapillary(1),
              0.75425841099261703,
              tolerance);
  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder
      + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, CheckResumeFile) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/honeycomb_N16L16.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {1, 5, 1};
  simulation_cfg.voxel_size = 1.0e-6;
  simulation_cfg.algorithm_json.SetFinalTime(0.35);
  simulation_cfg.algorithm_json.SetTimeStepSize(0.05);

  // Experiment config
  simulation_cfg.experiment_json.SetFlowAxis(1U);
  simulation_cfg.experiment_json.SetTemperature(340.0);
  simulation_cfg.experiment_json.SetAbsolutePressure(101325.0);
  simulation_cfg.experiment_json.SetBoundaryThickness(1U);
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_open", 10132.5);

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;  // dynamic
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  const std::string resume_file = simulation_cfg.folder + "/resume.h5";

  // Assert
  ASSERT_TRUE(std::filesystem::exists(resume_file + ".tgz"));

  // Extracting resume file
  std::string extract_command = "tar -xzf " + resume_file + ".tgz -C " + simulation_cfg.folder +
                                " &>/dev/null";
  system(extract_command.c_str());

  const std::string resume_file_expected = source_folder + "/resume_honeycomb_N16L16.h5";

  ASSERT_TRUE(std::filesystem::exists(resume_file));

  arma::vec time_info, time_info_expected;
  time_info.load(arma::hdf5_name(resume_file, "time_info"));
  time_info_expected.load(arma::hdf5_name(resume_file_expected, "time_info"));
  EXPECT_TRUE(arma::all(time_info == time_info_expected));

  arma::vec x_vec, x_vec_expected;
  x_vec.load(arma::hdf5_name(resume_file, "x_vec"));
  x_vec_expected.load(arma::hdf5_name(resume_file_expected, "x_vec"));
  EXPECT_TRUE(arma::all(x_vec == x_vec_expected));

  arma::vec dx_vec, dx_vec_expected;
  dx_vec.load(arma::hdf5_name(resume_file, "dx_vec"));
  dx_vec_expected.load(arma::hdf5_name(resume_file_expected, "dx_vec"));
  EXPECT_TRUE(arma::all(dx_vec == dx_vec_expected));

  arma::vec A_vec, A_vec_expected;
  A_vec.load(arma::hdf5_name(resume_file, "A_vec"));
  A_vec_expected.load(arma::hdf5_name(resume_file_expected, "A_vec"));
  EXPECT_TRUE(arma::all(A_vec == A_vec_expected));

  arma::vec previous_x_vec, previous_x_vec_expected;
  previous_x_vec.load(arma::hdf5_name(resume_file, "previous_x"));
  previous_x_vec_expected.load(arma::hdf5_name(resume_file_expected, "previous_x"));
  EXPECT_TRUE(arma::all(previous_x_vec == previous_x_vec_expected));

  arma::vec previous_dx_vec, previous_dx_vec_expected;
  previous_dx_vec.load(arma::hdf5_name(resume_file, "previous_dx"));
  previous_dx_vec_expected.load(arma::hdf5_name(resume_file_expected, "previous_dx"));
  EXPECT_TRUE(arma::all(previous_dx_vec == previous_dx_vec_expected));

  arma::vec previous_A_vec, previous_A_vec_expected;
  previous_A_vec.load(arma::hdf5_name(resume_file, "previous_A"));
  previous_A_vec_expected.load(arma::hdf5_name(resume_file_expected, "previous_A"));
  EXPECT_TRUE(arma::all(previous_A_vec == previous_A_vec_expected));

  const IndexType n_links = 16;
  arma::Col<double> all_interfaces, all_interfaces_expected;
  for (IndexType i = 0U; i < n_links; ++i) {
    all_interfaces.load(arma::hdf5_name(resume_file, "all_interfaces_" + std::to_string(i)));
    all_interfaces_expected.load(arma::hdf5_name(resume_file_expected,
                                                 "all_interfaces_" + std::to_string(i)));
    EXPECT_TRUE(arma::all(all_interfaces == all_interfaces_expected));
  }

  arma::uvec last_jump_info_size, last_jump_info_size_expected;
  last_jump_info_size.load(arma::hdf5_name(resume_file, "last_jump_info_size"));
  last_jump_info_size_expected.load(arma::hdf5_name(resume_file_expected, "last_jump_info_size"));
  EXPECT_TRUE(arma::all(last_jump_info_size == last_jump_info_size_expected));

  arma::umat last_jump_info_data, last_jump_info_data_expected;
  last_jump_info_data.load(arma::hdf5_name(resume_file, "last_jump_info_data"));
  last_jump_info_data_expected.load(arma::hdf5_name(resume_file_expected, "last_jump_info_data"));
  EXPECT_TRUE(arma::all(
    arma::vectorise(last_jump_info_data) == arma::vectorise(last_jump_info_data_expected)));

  arma::Col<IndexType> fluid_at_source, fluid_at_source_expected;
  fluid_at_source.load(arma::hdf5_name(resume_file, "fluid_at_source"));
  fluid_at_source_expected.load(arma::hdf5_name(resume_file_expected, "fluid_at_source"));
  EXPECT_TRUE(arma::all(fluid_at_source == fluid_at_source_expected));

  arma::Col<double> deltas, deltas_expected;
  deltas.load(arma::hdf5_name(resume_file, "deltas"));
  deltas_expected.load(arma::hdf5_name(resume_file_expected, "deltas"));
  EXPECT_TRUE(arma::all(deltas == deltas_expected));

  arma::Col<IndexType> plugs, plugs_expected;
  plugs.load(arma::hdf5_name(resume_file, "plugs"));
  plugs_expected.load(arma::hdf5_name(resume_file, "plugs"));
  EXPECT_TRUE(arma::all(plugs == plugs_expected));

  arma::Col<double> pressures, pressures_expected;
  pressures.load(arma::hdf5_name(resume_file, "pressures"));
  pressures_expected.load(arma::hdf5_name(resume_file_expected, "pressures"));
  EXPECT_TRUE(arma::all(pressures == pressures_expected));

  arma::Col<IndexType> last_plug_info, last_plug_info_expected;
  last_plug_info.load(arma::hdf5_name(resume_file, "last_plug_info"));
  last_plug_info_expected.load(arma::hdf5_name(resume_file_expected, "last_plug_info"));
  EXPECT_TRUE(arma::all(last_plug_info == last_plug_info_expected));

  arma::Col<IndexType> last_plug_inlet_info, last_plug_inlet_info_expected;
  last_plug_inlet_info.load(arma::hdf5_name(resume_file, "last_plug_inlet_info"));
  last_plug_inlet_info_expected.load(arma::hdf5_name(resume_file_expected, "last_plug_inlet_info"));
  EXPECT_TRUE(arma::all(last_plug_inlet_info == last_plug_inlet_info_expected));

  arma::umat destinations_data, destinations_data_expected;
  destinations_data.load(arma::hdf5_name(resume_file, "destinations_data_0"));
  destinations_data_expected.load(arma::hdf5_name(resume_file_expected, "destinations_data_0"));
  EXPECT_TRUE(arma::all(
    arma::vectorise(destinations_data) == arma::vectorise(destinations_data_expected)));

  arma::Col<IndexType> capillaries_with_interfaces, capillaries_with_interfaces_expected;
  capillaries_with_interfaces.load(arma::hdf5_name(resume_file, "capillaries_with_interfaces"));
  capillaries_with_interfaces_expected.load(arma::hdf5_name(resume_file_expected,
                                                            "capillaries_with_interfaces"));
  EXPECT_TRUE(arma::all(capillaries_with_interfaces == capillaries_with_interfaces_expected));

  arma::Col<IndexType> caps_with_interfaces_not_plugged_,
                       caps_with_interfaces_not_plugged_expected_;
  caps_with_interfaces_not_plugged_.load(arma::hdf5_name(resume_file,
                                                         "caps_with_interfaces_not_plugged"));
  caps_with_interfaces_not_plugged_expected_.load(arma::hdf5_name(
                                                              resume_file_expected,
                                                              "caps_with_interfaces_not_plugged"));
  EXPECT_TRUE(arma::all(
    caps_with_interfaces_not_plugged_ == caps_with_interfaces_not_plugged_expected_));

  arma::Col<double> caps_residual_limits_, caps_residual_limits_expected_;
  caps_residual_limits_.load(arma::hdf5_name(resume_file, "caps_residual_limits"));
  caps_residual_limits_expected_.load(arma::hdf5_name(resume_file_expected,
                                                      "caps_residual_limits"));
  EXPECT_TRUE(arma::all(caps_residual_limits_ == caps_residual_limits_expected_));

  arma::Col<int> roots_found_, roots_found_expected_;
  roots_found_.load(arma::hdf5_name(resume_file, "roots_found"));
  roots_found_expected_.load(arma::hdf5_name(resume_file_expected, "roots_found"));
  EXPECT_TRUE(arma::all(roots_found_ == roots_found_expected_));

  std::remove((simulation_cfg.folder + "/resume.h5").c_str());
  std::remove((simulation_cfg.folder + "/resume.h5.tgz").c_str());
  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, ResumeTwoCapillariesDifferentRadiusTest) {
  // OBS: Due to the low number of capillaries, the resume file used in the test was obtained
  // using resume_interval_ = 1
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/2_capillaries_different_radius.json";
  const std::string target_file = "/centerlines.json";
  const std::string resume_file = "/resume_2_capillaries_different_radius.h5.tgz";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);
  flow_simulator_test_utils::CopyFile(source_folder + resume_file,
                                      simulation_cfg.folder + "/resume.h5.tgz");

  // Define local configuration parameters
  simulation_cfg.shape = {1, 1, 3};
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_open", 10);
  simulation_cfg.algorithm_json.SetResume(true);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;  // dynamic
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  constexpr double tolerance = 1.0e-10;

  // Assert
  // These magic numbers are related to the parameters hardcoded by SetUp().
  double magic = 8.89777343751761407020;

  EXPECT_DOUBLE_EQ(pressures(0), p_abs + 1.0 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(1), p_abs + 0.0 * bc.second * voxel_size + magic);
  EXPECT_DOUBLE_EQ(pressures(2), p_abs - 1.0 * bc.second * voxel_size);

  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 1);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(0), 1.0);
  EXPECT_NEAR(interface->CalculateEffectiveInterfacePositionAtCapillary(1),
              0.75425841099261703,
              tolerance);

  // Ensures that the resume files has been loaded by checking if no snapshot from an simulation
  // that would start from the beginning exists.
  EXPECT_FALSE(std::filesystem::exists(simulation_cfg.folder + algorithm.snapshots_folder_ +
                                       "/simres_output_time_000018.h5"));

  std::remove((simulation_cfg.folder + "/resume.h5.tgz").c_str());
  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, CorruptedResumeTwoCapillariesDifferentRadiusTest) {
  // OBS: Due to the low number of capillaries, the resume file used in the test was obtained
  // using resume_interval_ = 1
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/2_capillaries_different_radius.json";
  const std::string target_file = "/centerlines.json";
  const std::string resume_file = "/corrupted_resume_2_capillaries_different_radius.h5.tgz";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);
  flow_simulator_test_utils::CopyFile(source_folder + resume_file,
                                      simulation_cfg.folder + "/resume.h5.tgz");

  // Define local configuration parameters
  simulation_cfg.shape = {1, 1, 3};
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_open", 10);
  simulation_cfg.algorithm_json.SetResume(true);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;  // dynamic
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  constexpr double tolerance = 1.0e-10;

  // Assert
  // These magic numbers are related to the parameters hardcoded by SetUp().
  double magic = 8.89777343751761407020;

  EXPECT_DOUBLE_EQ(pressures(0), p_abs + 1.0 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(1), p_abs + 0.0 * bc.second * voxel_size + magic);
  EXPECT_DOUBLE_EQ(pressures(2), p_abs - 1.0 * bc.second * voxel_size);

  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 1);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(0), 1.0);
  EXPECT_NEAR(interface->CalculateEffectiveInterfacePositionAtCapillary(1),
              0.75425841099261703,
              tolerance);

  // Ensures that the snapshots' files have started from the beggining since the resume file is
  // corrupted.
  EXPECT_TRUE(std::filesystem::exists(simulation_cfg.folder + algorithm.snapshots_folder_ +
                                       "/simres_output_time_000000.h5"));

  std::remove((simulation_cfg.folder + "/resume.h5.tgz").c_str());
  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, ThreeCapillariesEqualRadiusTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/3_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {1, 1, 4};
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_open", 10);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  // These magic numbers are related to the parameters hardcoded by SetUp().
  double magic = 8.99080544625631716826;

  EXPECT_DOUBLE_EQ(pressures(0), p_abs + 1.5 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(1), p_abs + 0.5 * bc.second * voxel_size + 1.0 * magic);
  EXPECT_DOUBLE_EQ(pressures(2), p_abs - 0.5 * bc.second * voxel_size + 2.0 * magic);
  EXPECT_DOUBLE_EQ(pressures(3), p_abs - 1.5 * bc.second * voxel_size);

  EXPECT_EQ(interface->GetNumberOfInterfaces(), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2), 1);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(0), 1.0);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(1), 1.0);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(2),
                   0.37138053565884643);

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, FourCapillariesEqualRadiusTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/4_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {1, 1, 5};
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_open", 10);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  // These magic numbers are related to the parameters hardcoded by SetUp().
  double magic = 8.69924691688713380699;

  EXPECT_DOUBLE_EQ(pressures(0), p_abs + 2.0 * bc.second * voxel_size + 0.0 * magic);
  EXPECT_DOUBLE_EQ(pressures(1), p_abs + 1.0 * bc.second * voxel_size + 1.0 * magic);
  EXPECT_DOUBLE_EQ(pressures(2), p_abs + 0.0 * bc.second * voxel_size + 2.0 * magic);
  EXPECT_DOUBLE_EQ(pressures(3), p_abs - 1.0 * bc.second * voxel_size + 3.0 * magic);
  EXPECT_DOUBLE_EQ(pressures(4), p_abs - 2.0 * bc.second * voxel_size + 0.0 * magic);

  EXPECT_EQ(interface->GetNumberOfInterfaces(), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(3), 1);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(0), 1.0);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(1), 1.0);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(2), 1.0);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(3),
                   0.056777188077616295);

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, YShapedCapillariesEqualRadiusTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/y-shaped_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {5, 1, 3};
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_open", 10);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  // These magic numbers are related to the parameters hardcoded by SetUp().
  double magic = 8.66994572868723878400;

  EXPECT_DOUBLE_EQ(pressures(0), p_abs + 1.0 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(1), p_abs + 0.0 * bc.second * voxel_size + magic);
  EXPECT_DOUBLE_EQ(pressures(2), p_abs - 1.0 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(3), p_abs - 1.0 * bc.second * voxel_size);

  EXPECT_EQ(interface->GetNumberOfInterfaces(), 2);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2), 1);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(0), 1.0);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(1),
                   0.43742068273639401);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(2),
                   0.43742068273639401);

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, ResumeYShapedCapillariesEqualRadiusTest) {
  // OBS: Due to the low number of capillaries, the resume file used in the test was obtained
  // using resume_interval_ = 1
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/y-shaped_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  const std::string resume_file = "/resume_y-shaped_capillaries_equal_radius.h5.tgz";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                    simulation_cfg.folder + target_file);
  flow_simulator_test_utils::CopyFile(source_folder + resume_file,
                                    simulation_cfg.folder + "/resume.h5.tgz");

  // Define local configuration parameters
  simulation_cfg.shape = {5, 1, 3};
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_open", 10);
  simulation_cfg.algorithm_json.SetResume(true);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  // These magic numbers are related to the parameters hardcoded by SetUp().
  double magic = 8.66994572868723878400;

  EXPECT_DOUBLE_EQ(pressures(0), p_abs + 1.0 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(1), p_abs + 0.0 * bc.second * voxel_size + magic);
  EXPECT_DOUBLE_EQ(pressures(2), p_abs - 1.0 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(3), p_abs - 1.0 * bc.second * voxel_size);

  EXPECT_EQ(interface->GetNumberOfInterfaces(), 2);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2), 1);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(0), 1.0);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(1),
                   0.43742068273639401);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(2),
                   0.43742068273639401);

  // Ensures that the resume file has been loaded by checking if no snapshot from an simulation
  // that would start from the beginning exists.
  EXPECT_FALSE(std::filesystem::exists(simulation_cfg.folder + algorithm.snapshots_folder_ +
                                       "/simres_output_time_000018.h5"));

  std::remove((simulation_cfg.folder + "/resume.h5.tgz").c_str());
  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, BigNetworkOneExitClosed) {
// Network Representation
/*                          6
                    ▗▄▛▘▘▌▌
         2       ▗▄▛▀
        ▌▌▄▄▄▄▄▄▖▌▌ 4
 0  ▗▟▀▘         ▀▜▄▄
 ▗▄▛▀               ▝▀▙    5
 ▌▌                    ▝▘▌▌
  ▀▜▄▖
     ▀█▗  1
        ▌▌
         ▝▀▙▄▖    3
             ▀▜▄ ▌▌

  */
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/big_network_one_exit.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {5, 1, 4};
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_closed", 11);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  // These magic numbers are related to the parameters hardcoded by SetUp().
  double rel_tolerance = 1.0e-12;
  double magic1 = 11.126939999999996;
  double magic2 = 8.16826931688562041245;
  double magic3 = 17.16593334836142759059;

  EXPECT_DOUBLE_EQ(pressures(0), p_abs + 1.5 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(1), p_abs + 0.5 * bc.second * voxel_size + magic1);
  EXPECT_DOUBLE_EQ(pressures(2), p_abs + 0.5 * bc.second * voxel_size + magic2);
  EXPECT_DOUBLE_EQ(pressures(3), p_abs + 0.5 * bc.second * voxel_size + magic1);
  EXPECT_DOUBLE_EQ(pressures(4), p_abs - 0.5 * bc.second * voxel_size + magic3);
  EXPECT_DOUBLE_EQ(pressures(5), p_abs - 1.5 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(6), p_abs - 1.5 * bc.second * voxel_size);

  EXPECT_EQ(interface->GetNumberOfInterfaces(), 3);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(3), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(4), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(5), 1);

  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(0),
                          0.0,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(4),
                          0.11703760615178528,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(5),
                          0.11703760615178528,
                          rel_tolerance));

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, BigNetworkOneExitOpen) {
// Network Representation
/*                          6
                    ▗▄▛▘▘▌▌
         2       ▗▄▛▀
        ▌▌▄▄▄▄▄▄▖▌▌ 4
 0  ▗▟▀▘         ▀▜▄▄
 ▗▄▛▀               ▝▀▙    5
 ▌▌                    ▝▘▌▌
  ▀▜▄▖
     ▀█▗  1
        ▌▌
         ▝▀▙▄▖    3
             ▀▜▄ ▌▌

  */
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/big_network_one_exit.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {5, 1, 4};
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_open", 11);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  // These magic numbers are related to the parameters hardcoded by SetUp().
  double rel_tolerance = 1.0e-12;
  double magic1 = -4.71621265673643108585;
  double magic2 = 8.16826931688561330702;
  double magic3 = -5.5;
  double magic4 = 17.16593334836142759059;

  EXPECT_DOUBLE_EQ(pressures(0), p_abs + 1.5 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(1), p_abs + 0.5 * bc.second * voxel_size + magic1);
  EXPECT_DOUBLE_EQ(pressures(2), p_abs + 0.5 * bc.second * voxel_size + magic2);
  EXPECT_DOUBLE_EQ(pressures(3), p_abs + 0.5 * bc.second * voxel_size + magic3);
  EXPECT_DOUBLE_EQ(pressures(4), p_abs - 0.5 * bc.second * voxel_size + magic4);
  EXPECT_DOUBLE_EQ(pressures(5), p_abs - 1.5 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(6), p_abs - 1.5 * bc.second * voxel_size);

  EXPECT_EQ(interface->GetNumberOfInterfaces(), 3);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(3), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(4), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(5), 1);

  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(0),
                          0.93092983079741298,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(4),
                          0.11703760615178686,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(5),
                          0.11703760615178686,
                          rel_tolerance));

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, BigNetworkOneExitSaveTest) {
// Network Representation
/*                          6
                    ▗▄▛▘▘▌▌
         2       ▗▄▛▀
        ▌▌▄▄▄▄▄▄▖▌▌ 4
 0  ▗▟▀▘         ▀▜▄▄
 ▗▄▛▀               ▝▀▙    5
 ▌▌                    ▝▘▌▌
  ▀▜▄▖
     ▀█▗  1
        ▌▌
         ▝▀▙▄▖    3
             ▀▜▄ ▌▌

  */
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/big_network_one_exit.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {5, 1, 4};

  // Local variables
  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_closed", 11);

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  double rel_tolerance = 1.0e-12;

  std::string folder = algorithm.GetFolder() + algorithm.snapshots_folder_;
  /// Load last time_output properties
  std::string t_output_fname = "/simres_output_time_000031.h5";


  /// Get output time
  // arma::vec current_time_at_output;
  arma::mat current_time_at_output;
  current_time_at_output.load(arma::hdf5_name(folder + t_output_fname, "current_time"));
  EXPECT_DOUBLE_EQ(current_time_at_output(0), 19.0);

  /// Get interface positions
  arma::vec interface_positions;
  interface_positions.load(arma::hdf5_name(folder + t_output_fname, "interface_positions"));
  EXPECT_TRUE(NearlyEqual(interface_positions(0), 0.0, rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface_positions(1), 0.11703760615178528, rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface_positions(2), 0.11703760615178528, rel_tolerance));

  /// Get event interface offsets
  arma::uvec interface_offsets;
  interface_offsets.load(arma::hdf5_name(folder + t_output_fname, "interface_offsets"));
  arma::uvec interface_offsets_tester = {0, 1, 1, 1, 1, 2, 3};
  for (arma::uword i = 0; i != 7; ++i) {
    EXPECT_EQ(interface_offsets(i), interface_offsets_tester(i));
  }

  /// Get event fluid at source
  arma::uvec fluid_at_source;
  fluid_at_source.load(arma::hdf5_name(folder + t_output_fname, "fluid_at_source"));
  arma::uvec fluid_at_source_tester = {1, 1, 0, 1, 1, 1};
  for (arma::uword i = 0; i != 6; ++i) {
    EXPECT_EQ(fluid_at_source(i), fluid_at_source_tester(i));
  }

  // Get pressures
  arma::vec loaded_pressures;
  loaded_pressures.load(arma::hdf5_name(folder + t_output_fname, "pressure"));
  for (arma::uword i = 0; i != pressures.n_elem; ++i) {
    EXPECT_DOUBLE_EQ(loaded_pressures(i), pressures(i));
  }

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, BigNetworkTwoExitsClosed) {
  // Network Representation
  /*                        7
                    ▗▄▛▘▘▌▌
         2       ▗▄▛▀
        ▌▌▄▄▄▄▄▄▖▌▌ 4
 0  ▗▟▀▘         ▀▜▄▄
 ▗▄▛▀               ▝▀▙    6
 ▌▌                    ▝▘▌▌
  ▀▜▄▖
     ▀█▗  1
        ▌▌
         ▝▀▙▄▖    3         5
             ▀▜▄ ▌▌▗▄▄▄▄▄ ▌▌

  */
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/big_network_two_exits.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {5, 1, 4};
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_closed", 16.7);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  // These magic numbers are related to the parameters hardcoded by SetUp().
  double rel_tolerance = 1.0e-12;
  double magic1 = -1.8068429321582737401;
  double magic2 = 12.40826883168562666526;
  double magic3 = 26.07355661954086656351;

  EXPECT_DOUBLE_EQ(pressures(0), p_abs + 1.5 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(1), p_abs + 0.5 * bc.second * voxel_size + 1.0 * magic1);
  EXPECT_DOUBLE_EQ(pressures(2), p_abs + 0.5 * bc.second * voxel_size + 1.0 * magic2);
  EXPECT_DOUBLE_EQ(pressures(3), p_abs - 0.5 * bc.second * voxel_size + 2.0 * magic1);
  EXPECT_DOUBLE_EQ(pressures(4), p_abs - 0.5 * bc.second * voxel_size + 1.0 * magic3);
  EXPECT_DOUBLE_EQ(pressures(5), p_abs - 1.5 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(6), p_abs - 1.5 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(7), p_abs - 1.5 * bc.second * voxel_size);

  EXPECT_EQ(interface->GetNumberOfInterfaces(), 2);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(3), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(4), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(5), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(6), 1);

  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(0),
                          1.0,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(1),
                          1.0,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(2),
                          1.0,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(3),
                          1.0,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(4),
                          1.0,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(5),
                          0.97518753575360095,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(6),
                          0.97518753575360095,
                          rel_tolerance));

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, BigNetworkTwoExitsOpen) {
  // Network Representation
  /*                        7
                    ▗▄▛▘▘▌▌
         2       ▗▄▛▀
        ▌▌▄▄▄▄▄▄▖▌▌ 4
 0  ▗▟▀▘         ▀▜▄▄
 ▗▄▛▀               ▝▀▙    6
 ▌▌                    ▝▘▌▌
  ▀▜▄▖
     ▀█▗  1
        ▌▌
         ▝▀▙▄▖    3         5
             ▀▜▄ ▌▌▗▄▄▄▄▄ ▌▌

  */
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/big_network_two_exits.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {5, 1, 4};
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_open", 16.7);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  // These magic numbers are related to the parameters hardcoded by SetUp().
  double rel_tolerance = 1.0e-12;
  double magic1 = 15.50029562932264504127;
  double magic2 = 12.40826882044462031729;
  double magic3 = 8.3499999999999996447;
  double magic4 = 26.07355660035128153140;

  EXPECT_DOUBLE_EQ(pressures(0), p_abs + 1.5 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(1), p_abs + 0.5 * bc.second * voxel_size + 1.0 * magic1);
  EXPECT_DOUBLE_EQ(pressures(2), p_abs + 0.5 * bc.second * voxel_size + 1.0 * magic2);
  EXPECT_DOUBLE_EQ(pressures(3), p_abs - 0.5 * bc.second * voxel_size + 1.0 * magic3);
  EXPECT_DOUBLE_EQ(pressures(4), p_abs - 0.5 * bc.second * voxel_size + 1.0 * magic4);
  EXPECT_DOUBLE_EQ(pressures(5), p_abs - 1.5 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(6), p_abs - 1.5 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(7), p_abs - 1.5 * bc.second * voxel_size);

  EXPECT_EQ(interface->GetNumberOfInterfaces(), 3);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(3), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(4), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(5), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(6), 1);

  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(0),
                          1.0,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(1),
                          1.0,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(2),
                          0.40946692436749665,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(3),
                          1.0,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(4),
                          1.0,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(5),
                          0.97518824657068937,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(6),
                          0.97518824657068937,
                          rel_tolerance));

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, ResumeBigNetworkTwoExitsOpen) {
  // OBS: Due to the low number of capillaries, the resume file used in the test was obtained
  // using resume_interval_ = 1
  // Network Representation
  /*                        7
                    ▗▄▛▘▘▌▌
         2       ▗▄▛▀
        ▌▌▄▄▄▄▄▄▖▌▌ 4
 0  ▗▟▀▘         ▀▜▄▄
 ▗▄▛▀               ▝▀▙    6
 ▌▌                    ▝▘▌▌
  ▀▜▄▖
     ▀█▗  1
        ▌▌
         ▝▀▙▄▖    3         5
             ▀▜▄ ▌▌▗▄▄▄▄▄ ▌▌

  */
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/big_network_two_exits.json";
  const std::string target_file = "/centerlines.json";
  const std::string resume_file = "/resume_big_network_two_exits.h5.tgz";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);
  flow_simulator_test_utils::CopyFile(source_folder + resume_file,
                                      simulation_cfg.folder + "/resume.h5.tgz");

  // Define local configuration parameters
  simulation_cfg.shape = {5, 1, 4};
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_open", 16.7);
  simulation_cfg.algorithm_json.SetResume(true);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  // These magic numbers are related to the parameters hardcoded by SetUp().
  double rel_tolerance = 1.0e-12;
  double magic1 = 15.50029562932264504127;
  double magic2 = 12.40826882044462031729;
  double magic3 = 8.3499999999999996447;
  double magic4 = 26.07355660035128153140;

  EXPECT_DOUBLE_EQ(pressures(0), p_abs + 1.5 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(1), p_abs + 0.5 * bc.second * voxel_size + 1.0 * magic1);
  EXPECT_DOUBLE_EQ(pressures(2), p_abs + 0.5 * bc.second * voxel_size + 1.0 * magic2);
  EXPECT_DOUBLE_EQ(pressures(3), p_abs - 0.5 * bc.second * voxel_size + 1.0 * magic3);
  EXPECT_DOUBLE_EQ(pressures(4), p_abs - 0.5 * bc.second * voxel_size + 1.0 * magic4);
  EXPECT_DOUBLE_EQ(pressures(5), p_abs - 1.5 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(6), p_abs - 1.5 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(7), p_abs - 1.5 * bc.second * voxel_size);

  EXPECT_EQ(interface->GetNumberOfInterfaces(), 3);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(3), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(4), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(5), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(6), 1);

  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(0),
                          1.0,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(1),
                          1.0,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(2),
                          0.40946692436749665,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(3),
                          1.0,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(4),
                          1.0,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(5),
                          0.97518824657068937,
                          rel_tolerance));
  EXPECT_TRUE(NearlyEqual(interface->CalculateEffectiveInterfacePositionAtCapillary(6),
                          0.97518824657068937,
                          rel_tolerance));

  // Ensures that the resume file has been loaded by checking if no snapshot from an simulation
  // that would start from the beginning exists.
  EXPECT_FALSE(std::filesystem::exists(simulation_cfg.folder + algorithm.snapshots_folder_ +
                                       "/simres_output_time_000021.h5"));

  std::remove((simulation_cfg.folder + "/resume.h5.tgz").c_str());
  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, YShapedCapillariesDifferentRadiusPlugTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/y-shaped_capillaries_different_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {3, 1, 3};
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_closed", 8);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  EXPECT_DOUBLE_EQ(pressures(1), p_abs + bc.second * voxel_size);

  EXPECT_TRUE(interface->IsCapillaryPlugged(0));
  EXPECT_FALSE(interface->IsCapillaryPlugged(1));
  EXPECT_FALSE(interface->IsCapillaryPlugged(2));
  EXPECT_EQ(interface->GetNumberOfInterfaces(), 2);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2), 0);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(0), 1.0);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(1),
                   0.87020952768520166);

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, YShapedCapillariesDifferentRadiusUnplugTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/y-shaped_capillaries_different_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {3, 1, 3};
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_closed", 10);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  constexpr double tolerance = 1.0e-10;

  // Assert
  EXPECT_FALSE(interface->IsCapillaryPlugged(0));
  EXPECT_FALSE(interface->IsCapillaryPlugged(1));
  EXPECT_FALSE(interface->IsCapillaryPlugged(2));
  EXPECT_EQ(interface->GetNumberOfInterfaces(), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2), 1);
  EXPECT_DOUBLE_EQ(pressures(1), p_abs + bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(0), 1.0);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(1), 1.0);
  EXPECT_NEAR(interface->CalculateEffectiveInterfacePositionAtCapillary(2),
              0.14190553914266757,
              tolerance);
  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, YShapedCapillariesDifferentRadiusJumpTogetherTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/y-shaped_capillaries_different_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {3, 1, 3};

  // Local variables
  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_closed", 10.0);

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  EXPECT_FALSE(interface->IsCapillaryPlugged(0));
  EXPECT_FALSE(interface->IsCapillaryPlugged(1));
  EXPECT_FALSE(interface->IsCapillaryPlugged(2));
  EXPECT_EQ(interface->GetNumberOfInterfaces(), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2), 1);

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, InitialTimeBiggerThanFinalTimeTest) {
  // Define local configuration parameters
  simulation_cfg.algorithm_json.SetInitialTime(20.0);
  simulation_cfg.algorithm_json.SetFinalTime(10.0);

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;

  // Assert
  std::pair<int, std::string> rc;
  rc = algorithm.Initialise(simulation_cfg);
  EXPECT_EQ(rc.first, -1);

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, LessThanTwoFluids) {
  // Define local configuration parameters
  std::vector<FluidJSON> my_fluids;

  FluidJSON fluid1("water", "constant");

  fluid1.AddProperty("dynamic_viscosity", 1.0);
  my_fluids.push_back(fluid1);
  simulation_cfg.fluids_json = my_fluids;

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;

  // Assert
  std::pair<int, std::string> rc;
  rc = algorithm.Initialise(simulation_cfg);
  EXPECT_EQ(rc.first, -2);

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, AllFluidsMustHaveDynamicViscosityProperty) {
  // Define local configuration parameters
  std::vector<FluidJSON> my_fluids;

  FluidJSON fluid1("water", "constant");
  FluidJSON fluid2("oil", "constant");

  fluid1.AddProperty("dynamic_viscosity", 1.0);
  my_fluids.push_back(fluid1);
  my_fluids.push_back(fluid2);
  simulation_cfg.fluids_json = my_fluids;

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;

  // Assert
  std::pair<int, std::string> rc;
  rc = algorithm.Initialise(simulation_cfg);
  EXPECT_EQ(rc.first, -3);

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, InterfaceJumpThroughInletNodePluggingCapillary) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/1_capillary.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {1, 1, 2};
  simulation_cfg.fluid_interface_json.AddProperty("interfacial_tension", 10.0);
  simulation_cfg.wettability_json.AddProperty("contact_angle", 140.0);  // degrees
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_closed", 1.0e0);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  algorithm.GetInterface()->MoveInterfacePositionsAtCapillaryByDelta(0, 0.5);

  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  double interface_position = algorithm.GetInterface()->GetInterfacePositionAtCapillaryNearSide(
    0, Side::source);
  EXPECT_TRUE(interface_position > 0.0);
  EXPECT_DOUBLE_EQ(pressures(0), p_abs + (bc.second * voxel_size) / 2);
  EXPECT_DOUBLE_EQ(pressures(1), p_abs - (bc.second * voxel_size) / 2);
  EXPECT_TRUE(interface->IsCapillaryPlugged(0));
  EXPECT_EQ(interface->GetNumberOfInterfaces(), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 1);

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, TwoInterfacesJumpThroughInletNodesPluggingCapillary) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/inverted_y-shaped_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {3, 1, 3};

  simulation_cfg.fluid_interface_json.AddProperty("interfacial_tension", 10.0);
  simulation_cfg.wettability_json.AddProperty("contact_angle", 140.0);  // degrees
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_closed", 1.0e0);

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  algorithm.GetInterface()->MoveInterfacePositionsAtCapillaryByDelta(0U, 0.5);
  algorithm.GetInterface()->MoveInterfacePositionsAtCapillaryByDelta(2U, -0.5);

  std::shared_ptr<DynamicCapillaryNetworkContext> myContext = algorithm.GetContext();
  arma::Col<IndexType> &myInletNodes = myContext->inlet_nodes_;
  arma::Col<IndexType> &myOutletNodes = myContext->outlet_nodes_;

  myInletNodes.shed_row(arma::as_scalar(arma::find(myInletNodes == 3U)));
  myOutletNodes.insert_rows(myOutletNodes.n_elem, arma::Row<IndexType>({ 3U }));

  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  EXPECT_TRUE(interface->IsCapillaryPlugged(0));
  EXPECT_FALSE(interface->IsCapillaryPlugged(2));
  EXPECT_TRUE(interface->GetCapillaryPluggedSide(0) == Side::source);
  EXPECT_FALSE(interface->GetCapillaryPluggedSide(2) == Side::target);
  EXPECT_EQ(interface->GetNumberOfInterfaces(), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2), 0);

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, TwoCapillariesThreeVoxelsAwayFromEdgesTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/2_capillaries_3_voxels_away_from_edges.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {1, 1, 9};
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_open", 10);
  simulation_cfg.experiment_json.SetBoundaryThickness(4U);

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;  // dynamic
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();

  std::shared_ptr<DynamicCapillaryNetworkContext> myContext = algorithm.GetContext();
  arma::Col<IndexType> &myInletNodes = myContext->inlet_nodes_;
  arma::Col<IndexType> &myOutletNodes = myContext->outlet_nodes_;

  // algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Assert
  EXPECT_EQ(myInletNodes.n_elem, 1U);
  EXPECT_EQ(myInletNodes(0), 0U);
  EXPECT_EQ(myOutletNodes.n_elem, 1U);
  EXPECT_EQ(myOutletNodes(0), 2U);
}

TEST_F(DynamicCapillaryNetwork, BubbleJumpThroughTarget) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/inverted_y-shaped_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {3, 1, 3};
  simulation_cfg.fluid_interface_json.AddProperty("interfacial_tension", 10.0);

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  algorithm.GetInterface()->AddNewInterface(0, Side::source);
  algorithm.GetInterface()->RemoveInterface(2, Side::target);

  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  EXPECT_EQ(interface->GetNumberOfInterfaces(), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2), 0);

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, BubbleJumpThroughSource) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/inverted_y-shaped_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {3, 1, 3};
  simulation_cfg.fluid_interface_json.AddProperty("interfacial_tension", 10.0);

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  algorithm.GetInterface()->AddNewInterface(2, Side::target);
  algorithm.GetInterface()->RemoveInterface(0, Side::source);

  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  EXPECT_EQ(interface->GetNumberOfInterfaces(), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2), 0);

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, BubblePlusInterfaceJumpThroughTarget) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/inverted_y-shaped_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {3U, 1U, 3U};
  simulation_cfg.fluid_interface_json.AddProperty("interfacial_tension", 5.0);

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  algorithm.GetInterface()->AddNewInterface(0U, Side::source);
  algorithm.GetInterface()->AddNewInterface(0U, Side::source);
  algorithm.GetInterface()->RemoveInterface(2U, Side::target);

  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  EXPECT_EQ(interface->GetNumberOfInterfaces(), 2U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0U), 0U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1U), 1U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2U), 1U);

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, BubblePlusInterfaceJumpThroughSource) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/inverted_y-shaped_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {3U, 1U, 3U};
  simulation_cfg.fluid_interface_json.AddProperty("interfacial_tension", 5.0);

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  algorithm.GetInterface()->AddNewInterface(2U, Side::target);
  algorithm.GetInterface()->AddNewInterface(2U, Side::target);
  algorithm.GetInterface()->RemoveInterface(0U, Side::source);
  algorithm.GetInterface()->FlipFluidAtSource(0U);

  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  EXPECT_EQ(interface->GetNumberOfInterfaces(), 2);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2), 0);

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(DynamicCapillaryNetwork, ViscosityBehaviourTest) {
  // Copy input file
  const std::string source_folder = "test/src/flow_simulator/input";
  const std::string source_file = "/y-shaped_capillaries_equal_radius.json";
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = {5, 1, 3};
  simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_closed", 1.0e1);
  simulation_cfg.experiment_json.SetAbsolutePressure(1.0e6);

  // Redefining fluids
  simulation_cfg.fluids_json.clear();

  FluidJSON fluid1("water", "pure_water");
  FluidJSON fluid2("oil", "constant");

  fluid2.AddProperty("dynamic_viscosity", 1.0111);

  simulation_cfg.fluids_json.push_back(fluid1);
  simulation_cfg.fluids_json.push_back(fluid2);

  // Local variables
  double p_abs = simulation_cfg.experiment_json.GetAbsolutePressure();
  double voxel_size = simulation_cfg.voxel_size;

  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  std::pair<int, std::string> rc = algorithm.Initialise(simulation_cfg);
  ASSERT_GE(rc.first, 0);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();
  arma::vec pressures = algorithm.GetPressures();
  std::shared_ptr<CapillaryInterface> interface = algorithm.GetInterface();

  // Assert
  // These magic numbers are related to the parameters hardcoded by SetUp().
  double magic = 8.63739847473334521055;

  EXPECT_DOUBLE_EQ(pressures(0), p_abs + 1.0 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(1), p_abs + 0.0 * bc.second * voxel_size + magic);
  EXPECT_DOUBLE_EQ(pressures(2), p_abs - 1.0 * bc.second * voxel_size);
  EXPECT_DOUBLE_EQ(pressures(3), p_abs - 1.0 * bc.second * voxel_size);

  EXPECT_EQ(interface->GetNumberOfInterfaces(), 2);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0), 0);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1), 1);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2), 1);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(0), 1.0);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(1),
                   0.48240693765832338);
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(2),
                   0.48240693765832338);

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

#endif  // TEST_SRC_FLOW_SIMULATOR_DYNAMIC_SIMULATION_TESTS_H_
