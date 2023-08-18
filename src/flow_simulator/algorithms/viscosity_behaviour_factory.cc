/**
 * \file src/flow_simulator/algorithms/viscosity_behaviour_factory.cc
 * \brief Contains the implementation of \c ViscosityBehaviourFactory class methods.
 *
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2020
 *
 * This source file contains the implementation of \c ViscosityBehaviourFactory class methods.
 * The \c CreateViscosityBehaviour() instantiates the appropriate \c ViscosityBehaviour object, and
 * initialise it.
 */

#include "src/flow_simulator/algorithms/viscosity_behaviour_factory.h"

#include <glog/logging.h>

#include <memory>
#include <string>

#include "src/exec_manager/fluid_json.h"
#include "src/exec_manager/experiment_json.h"

#include "src/flow_simulator/algorithms/viscosity_behaviour.h"
#include "src/flow_simulator/algorithms/constant_viscosity_behaviour.h"
#include "src/flow_simulator/algorithms/pure_water_viscosity_behaviour.h"
#include "src/flow_simulator/algorithms/carbonated_water_viscosity_behaviour.h"
#include "src/flow_simulator/algorithms/supercritical_co2_viscosity_behaviour.h"
#include "src/flow_simulator/algorithms/brine_viscosity_behaviour.h"
#include "src/flow_simulator/algorithms/carbonated_brine_viscosity_behaviour.h"
#include "src/flow_simulator/algorithms/seawater_viscosity_behaviour.h"

using ConstantViscosityBehaviour = simulator::ConstantViscosityBehaviour;
using PureWaterViscosityBehaviour = simulator::PureWaterViscosityBehaviour;
using CarbonatedWaterViscosityBehaviour = simulator::CarbonatedWaterViscosityBehaviour;
using SupercriticalCO2ViscosityBehaviour = simulator::SupercriticalCO2ViscosityBehaviour;
using BrineViscosityBehaviour = simulator::BrineViscosityBehaviour;
using CarbonatedBrineViscosityBehaviour = simulator::CarbonatedBrineViscosityBehaviour;
using SeawaterViscosityBehaviour = simulator::SeawaterViscosityBehaviour;
using ViscosityBehaviour = simulator::ViscosityBehaviour;
using ViscosityBehaviourFactory = simulator::ViscosityBehaviourFactory;

std::unique_ptr<ViscosityBehaviour> ViscosityBehaviourFactory::CreateViscosityBehaviour(
  FluidJSON fluid_json, ExperimentJSON experiment_json) {
/**
 * The \c CreateViscosityBehaviour() method is responsible for retrieving the the appropriate
 * \c ViscosityBehaviour object, and initialise it.
 *
 * \param[in]  fluid_json       Has the name of the selected viscosity behaviour for the fluid, and,
 *                              for some implementations, `dynamic_viscosity` parameter
 * \param[in]  experiment_json  Has the temperature, necessary for some implementations of
 *                              `ViscosityBehaviour`
 */
  std::string name = fluid_json.GetViscosityBehaviourName();

  if (name == "constant") {
    return std::make_unique<ConstantViscosityBehaviour>(
      fluid_json.GetProperty("dynamic_viscosity"));
  } else if (name == "pure_water") {
    return std::make_unique<PureWaterViscosityBehaviour>(experiment_json.GetTemperature());
  } else if (name == "carbonated_water") {
    return std::make_unique<CarbonatedWaterViscosityBehaviour>(
      experiment_json.GetTemperature(),
      fluid_json.GetProperty("CO2_weight_percentage"));
  } else if (name == "supercritical_CO2") {
    return std::make_unique<SupercriticalCO2ViscosityBehaviour>(experiment_json.GetTemperature());
  } else if (name == "brine") {
    return std::make_unique<BrineViscosityBehaviour>(experiment_json.GetTemperature(),
            fluid_json.GetProperty("salinity_ppm"));
  } else if (name == "carbonated_brine") {
    return std::make_unique<CarbonatedBrineViscosityBehaviour>(
            experiment_json.GetTemperature(),
            fluid_json.GetProperty("salinity_ppm"),
            fluid_json.GetProperty("CO2_weight_percentage"));
  } else if (name == "seawater") {
    return std::make_unique<SeawaterViscosityBehaviour>(experiment_json.GetTemperature());
  } else {
    LOG(FATAL) << "Please select a valid viscosity behaviour object.";
  }
}  // end of ViscosityBehaviourFactory::CreateViscosityBehaviour() method
