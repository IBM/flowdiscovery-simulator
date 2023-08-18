/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/models/linear_molecular_kinetics/linear_molecular_kinetics_model.h
 * \brief Contains the \c LinearMolecularKineticsModel class.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 *
 * This header file contains the the \c LinearMolecularKineticsModel class that derives from
 * \c DynamicModelBase.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_LINEAR_MOLECULAR_KINETICS_LINEAR_MOLECULAR_KINETICS_MODEL_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_LINEAR_MOLECULAR_KINETICS_LINEAR_MOLECULAR_KINETICS_MODEL_H_

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "src/flow_simulator/algorithms/dynamic_capillary_network/models/dynamic_model_base.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network_context.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/models/linear_molecular_kinetics/linear_molecular_kinetics_physics.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/models/linear_molecular_kinetics/linear_molecular_kinetics_geometry.h"
#include "src/exec_manager/fluid_json.h"
#include "src/exec_manager/wettability_json.h"
#include "src/exec_manager/fluid_interface_json.h"

namespace simulator {

/**
 * \class LinearMolecularKineticsModel linear_molecular_kinetics_model.h "src/flow_simulator/algorithms/dynamic_capillary_network/models/linear_molecular_kinetics/linear_molecular_kinetics_model.h"
 * \brief Linear Molecular Kinetics model for the Dynamic Capillary Network algorithm.
 *
 * Container-like class that hosts the appropriate \c Physics and \c Geometry objects for the
 * Linear Molecular Kinetics model.
 */
class LinearMolecularKineticsModel : public DynamicModelBase {
 public:
  LinearMolecularKineticsModel(const std::string &folder,
                               std::shared_ptr<DynamicCapillaryNetworkContext> context,
                               std::vector<FluidJSON> &fluids_json,
                               WettabilityJSON &wettability_json,
                               FluidInterfaceJSON &fluid_interface_json,
                               ExperimentJSON &experiment_json) {
    // Create and populate resident fluid
    std::shared_ptr<Fluid> resident_fluid = std::make_shared<Fluid>(fluids_json[0],
                                                                    experiment_json);

    // Create and populate injected fluid
    std::shared_ptr<Fluid> injected_fluid = std::make_shared<Fluid>(fluids_json[1],
                                                                    experiment_json);

    // Create and populate wettability
    std::unique_ptr<Wettability> my_wettability = std::make_unique<Wettability>(wettability_json);
    my_wettability->AddFluid(resident_fluid);
    my_wettability->AddFluid(injected_fluid);
    // Adding a derivative property contact angle (in radians)
    my_wettability->AddProperty("contact_angle_in_radians",
                                wettability_json.GetProperty("contact_angle") * arma::datum::pi
                                                                              / 180.0);

    // Create and populate fluid interface
    std::unique_ptr<FluidInterface> my_fluid_interface = std::make_unique<FluidInterface>(
      fluid_interface_json);
    my_fluid_interface->AddFluid(resident_fluid);
    my_fluid_interface->AddFluid(injected_fluid);

    // Create geometry
    geometry_ = std::make_shared<LinearMolecularKineticsGeometry>(folder, context);

    // Create physics
    physics_ = std::make_shared<LinearMolecularKineticsPhysics>(context,
                                                                resident_fluid,
                                                                injected_fluid,
                                                                std::move(my_wettability),
                                                                std::move(my_fluid_interface));
  }

  /// Getter for shared pointer to DynamicCapillaryGeometryBase
  std::shared_ptr<DynamicCapillaryGeometryBase> GetGeometry(void) { return geometry_; }

  /// Getter for shared pointer to DynamicCapillaryPhysicsBase
  std::shared_ptr<DynamicCapillaryPhysicsBase> GetPhysics(void) { return physics_; }

 protected:
  /**
   * \brief Shared pointer to \c LinearMolecularKineticsGeometry instance.
   *
   * This member variable stores a shared pointer to the \c LinearMolecularKineticsGeometry
   * object set by the \c LinearMolecularKineticsModel constructor.
   */
  std::shared_ptr<LinearMolecularKineticsGeometry> geometry_;

  /**
   * \brief Shared pointer to \c LinearMolecularKineticsPhysics instance.
   *
   * This member variable stores a shared pointer to the \c LinearMolecularKineticsPhysics
   * object set by the \c LinearMolecularKineticsModel constructor.
   */
  std::shared_ptr<LinearMolecularKineticsPhysics> physics_;
};  // end of class LinearMolecularKineticsModel

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_LINEAR_MOLECULAR_KINETICS_LINEAR_MOLECULAR_KINETICS_MODEL_H_
