/**
 * \file test/src/flow_simulator/dynamic_simulation_stress_tests.h
 * \brief It contains the simulator stress tests that use considerably large input files.
 *
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 */

#ifndef TEST_SRC_FLOW_SIMULATOR_DYNAMIC_SIMULATION_STRESS_TESTS_H_
#define TEST_SRC_FLOW_SIMULATOR_DYNAMIC_SIMULATION_STRESS_TESTS_H_

#include <gtest/gtest.h>
#include <map>
#include <string>
#include <memory>
#include <utility>
#include <vector>
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network.h"
#include "src/flow_simulator/algorithms/network_reader.h"
#include "test/src/utils/flow_simulator_test_utils.h"
#include "src/exec_manager/simulation_config.h"
#include "src/json_parser/json_parser.h"
#include "src/exec_manager/fluid_json.h"

using NetworkReader = simulator::NetworkReader;
using NetworkInformation = simulator::NetworkInformation;
using DynamicCapillaryNetworkAlgorithm = simulator::DynamicCapillaryNetworkAlgorithm;

class Stress: public ::testing::Test {
 public:
  Stress() {
    // initialization code here
  }

  void SetUp() {
    // code here will execute just before the test ensues

    // Define configuration parameters

    // Defining directory where the will be the JSON file centerline.json
    simulation_cfg.folder = "test/results";

    // Inputs for the algorithm execution

    // Simulation config
    simulation_cfg.shape = { 1, 1, 1 };  // Default shape
    simulation_cfg.voxel_size = 1.0e-6;

    // Algorithm config
    simulation_cfg.algorithm_json.SetName("dynamic");
    simulation_cfg.algorithm_json.SetModel("linear_molecular_kinetics");
    simulation_cfg.algorithm_json.SetInitialTime(0.0);
    simulation_cfg.algorithm_json.SetFinalTime(0.35);
    simulation_cfg.algorithm_json.SetTimeStepSize(0.05);
    simulation_cfg.algorithm_json.SetRelativeTolerance(1.0e-4);
    simulation_cfg.algorithm_json.SetAbsoluteLinkTolerance(1.0e-3);
    simulation_cfg.algorithm_json.SetAbsoluteNodeTolerance(1.0e-3);

    // Experiment config
    simulation_cfg.experiment_json.SetFlowAxis(1U);
    simulation_cfg.experiment_json.SetTemperature(340.0);
    simulation_cfg.experiment_json.SetAbsolutePressure(101325.0);
    simulation_cfg.experiment_json.SetBoundaryThickness(1U);
    simulation_cfg.experiment_json.SetBoundaryCondition("pressure_gradient_open", 10132.5);

    // Wettability config
    simulation_cfg.wettability_json.SetName("water/rock/oil");
    simulation_cfg.wettability_json.AddProperty("contact_angle", 0.0);
    simulation_cfg.wettability_json.AddProperty("linear_mk", 0.0);

    // Fluit interface config
    simulation_cfg.fluid_interface_json.SetName("water/oil");
    simulation_cfg.fluid_interface_json.AddProperty("interfacial_tension", 0.06347);

    // Defining fluids
    FluidJSON fluid1("water", "constant");
    FluidJSON fluid2("oil", "constant");

    fluid1.AddProperty("dynamic_viscosity", 1.0e-3);
    fluid2.AddProperty("dynamic_viscosity", 1.002e-4);

    simulation_cfg.fluids_json.push_back(fluid1);
    simulation_cfg.fluids_json.push_back(fluid2);
  }

  void TearDown() {
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be
  }

  ~Stress()  {
    // cleanup any pending stuff, but no exceptions allowed
  }

  // put in any custom data members that you need
  SimulationConfig simulation_cfg;
};



arma::uvec::fixed<3> getShape(const std::string json_file_name) {
  NetworkReader reader;
  NetworkInformation net_info = reader.GetNetwork(json_file_name);

  return { static_cast<arma::uword>(floor(net_info.ctrl_voxels.col(0).max())) + 1U,
           static_cast<arma::uword>(floor(net_info.ctrl_voxels.col(1).max())) + 1U,
           static_cast<arma::uword>(floor(net_info.ctrl_voxels.col(2).max())) + 1U };
}



void workUnit(const std::string source_folder,
              const std::string source_file,
              SimulationConfig simulation_cfg) {
  // Copy input file
  const std::string target_file = "/centerlines.json";
  flow_simulator_test_utils::CopyFile(source_folder + source_file,
                                      simulation_cfg.folder + target_file);

  // Define local configuration parameters
  simulation_cfg.shape = getShape(simulation_cfg.folder + "/centerlines.json");  // Default shape

  // Prepare algorithm
  DynamicCapillaryNetworkAlgorithm algorithm;
  algorithm.Initialise(simulation_cfg);
  algorithm.BuildGeometricalRepresentation();
  algorithm.BuildPhysicalEquations();
  std::remove((simulation_cfg.folder + target_file).c_str());

  // Act
  algorithm.SolvePhysicalEquations();

  flow_simulator_test_utils::CleanFolder((simulation_cfg.folder + algorithm.snapshots_folder_));
}

TEST_F(Stress, honeycomb_N1222L1763) {
  struct stat buffer;
  if ((stat("./test/src/flow_simulator/input/honeycomb_N1222L1763.json", &buffer) == 0)) {
    workUnit("test/src/flow_simulator/input", "/honeycomb_N1222L1763.json", simulation_cfg);
  } else {
    std::cerr << "Input file not found!" << std::endl
        << "Please execute the following command:" << std::endl
        << "wget https://raw.github.ibm.com/gist/rneumann/f01e6aa35f378e84637be777c6359ecf/raw"
        << "/b839bce6b7266c2927d60b957f3d411d1412d077/honeycomb_N1222L1763.json?"
        << "token=AACZS7O5SIZJ4RK4RYHUC5LBXRU2A"
        << " -O ./test/src/flow_simulator/input/honeycomb_N1222L1763.json" << std::endl;
    FAIL();
  }
}

#endif  // TEST_SRC_FLOW_SIMULATOR_DYNAMIC_SIMULATION_STRESS_TESTS_H_
