/**
 * \file src/flow_simulator/algorithms/viscosity_behaviour_factory.h
 * \brief Contains the \c ViscosityBehaviourFactory class.
 *
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2020
 *
 * This header file contains the \c ViscosityBehaviourFactory class.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_VISCOSITY_BEHAVIOUR_FACTORY_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_VISCOSITY_BEHAVIOUR_FACTORY_H_

#include <memory>

#include "src/exec_manager/experiment_json.h"
#include "src/exec_manager/fluid_json.h"

#include "src/flow_simulator/algorithms/pure_water_viscosity_behaviour.h"
#include "src/flow_simulator/algorithms/viscosity_behaviour.h"

namespace simulator {

/**
 * \class ViscosityBehaviourFactory viscosity_behaviour_factory.h
 * "src/flow_simulator/algorithms/viscosity_behaviour_factory.h"
 * \brief Create `ViscosityBehaviour` objects
 */
class ViscosityBehaviourFactory {
 public:
  static std::unique_ptr<ViscosityBehaviour> CreateViscosityBehaviour(FluidJSON, ExperimentJSON);
};  // end of class ViscosityBehaviourFactory

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_VISCOSITY_BEHAVIOUR_FACTORY_H_
