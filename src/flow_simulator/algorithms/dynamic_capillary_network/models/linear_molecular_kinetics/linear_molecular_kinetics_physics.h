/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/models/linear_molecular_kinetics/linear_molecular_kinetics_physics.h
 * \brief Contains the \c LinearMolecularKineticsPhysics class.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 *
 * This header file contains the the \c LinearMolecularKineticsPhysics class that derives from
 * \c DynamicCapillaryPhysicsBase.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_LINEAR_MOLECULAR_KINETICS_LINEAR_MOLECULAR_KINETICS_PHYSICS_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_LINEAR_MOLECULAR_KINETICS_LINEAR_MOLECULAR_KINETICS_PHYSICS_H_

#include <armadillo>
#include <memory>
#include <string>
#include <utility>
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network_context.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/models/dynamic_capillary_physics_base.h"
#include "src/flow_simulator/algorithms/fluid.h"
#include "src/flow_simulator/algorithms/wettability.h"
#include "src/flow_simulator/algorithms/fluid_interface.h"

using Fluid = simulator::Fluid;
using Wettability = simulator::Wettability;
using FluidInterface = simulator::FluidInterface;

namespace simulator {

/**
 * \class LinearMolecularKineticsPhysics linear_molecular_kinetics_physics.h "src/flow_simulator/algorithms/dynamic_capillary_network/models/linear_molecular_kinetics/linear_molecular_kinetics_physics.h"
 * \brief Linear Molecular Kinetics physics for the Dynamic Capillary Network algorithm.
 *
 * ...
 */

class LinearMolecularKineticsPhysics : public DynamicCapillaryPhysicsBase {
 public:
  LinearMolecularKineticsPhysics(std::shared_ptr<DynamicCapillaryNetworkContext> context,
                                 std::shared_ptr<Fluid> resident_fluid,
                                 std::shared_ptr<Fluid> injected_fluid,
                                 std::unique_ptr<Wettability> &&wettability,
                                 std::unique_ptr<FluidInterface> &&fluid_interface)
    : DynamicCapillaryPhysicsBase(context),
      resident_fluid_(resident_fluid),
      injected_fluid_(injected_fluid),
      wettability_(std::move(wettability)),
      fluid_interface_(std::move(fluid_interface)) { }

  /// Calculate phase permeabilities
  arma::Col<double> CalculatePhasePermeability(const arma::Col<double> &flow_rate,
                                               const arma::Col<double> &pressure) const override;

  /// Getter for the resident fluid name
  std::string GetResidentFluidName(void) const {
    return resident_fluid_->GetName();
  }

  /// Getter for the injected fluid name
  std::string GetInjectedFluidName(void) const {
    return injected_fluid_->GetName();
  }

  /// Getter for resident fluid viscosity
  double GetResidentFluidViscosity(double pressure_in_pascal) const {
    return resident_fluid_->GetViscosity(pressure_in_pascal);
  }

  /// Getter for injected fluid viscosity
  double GetInjectedFluidViscosity(double pressure_in_pascal) const {
    return injected_fluid_->GetViscosity(pressure_in_pascal);
  }

  /// Return true if the parameters fall between the designed behaviour limits, false otherwise
  double IsViscosityBehaviourDesignedForResidentFluid(const double lower_pressure,
                                                     const double upper_pressure) const override {
    return resident_fluid_->IsViscosityBehaviourDesignedFor(lower_pressure, upper_pressure);
  }

  /// Return true if the parameters fall between the designed behaviour limits, false otherwise
  double IsViscosityBehaviourDesignedForInjectedFluid(const double lower_pressure,
                                                     const double upper_pressure) const override {
    return injected_fluid_->IsViscosityBehaviourDesignedFor(lower_pressure, upper_pressure);
  }

  /// Getter for wettability contact angle (in degrees)
  double GetWettabilityContactAngle(void) const {
    return wettability_->GetProperty("contact_angle");
  }

  /// Getter for wettability contact angle (in radians)
  double GetWettabilityContactAngleInRadians(void) const {
    return wettability_->GetProperty("contact_angle_in_radians");
  }

  /// Getter for wettability linear mk
  double GetWettabilityLinearMK(void) const {
    return wettability_->GetProperty("linear_mk");
  }

  /// Getter for fluid interface interfacial tension
  double GetFluidInterfaceInterfacialTension(void) const {
    return fluid_interface_->GetProperty("interfacial_tension");
  }

 private:
  /// Pointer to \c Fluid instance for the resident fluid
  std::shared_ptr<Fluid> resident_fluid_;

  /// Pointer to \c Fluid instance for the injected fluid
  std::shared_ptr<Fluid> injected_fluid_;

  /**
   * \brief Pointer to \c Wettability instance for the interaction between the fluids and the rock
   *
   * We expect that this \c Wettability instance to have the following properties:
   *
   * contact_angle - this floating-point variable stores the value (in degrees) of the contact angle.
   *
   * contact_angle_in_radians - this floating-point variable stores the value (in radians) of the
   * contact angle.
   *
   * linear_mk - this only exists for the Linear Molecular Kinetics model, the floating-point
   * variable stores the value of the linear molecular kinetics parameter, which is the product of 2
   * constants (C1 and C2) and has viscosity units (Pa*s).
   *
   * \f[
   *   \cos{\theta_D}\left(\dot{x} \right) = \cos{\theta_S}- \frac{C_2}{\sigma} \sinh^{-1}
   *     \left(C_1 L \dot{x} \right)
   * \f]
   */
  std::unique_ptr<Wettability> wettability_;

  /**
   * \brief Pointer to \c FluidInterface instance for the interaction between the fluids
   *
   * We expect that this \c FluidInterface instance to have the following property:
   *
   * interfacial_tension - this floating-point variable stores the value (in N/m) of the interfacial
   * tension with which the fluid flow simulation is performed.
   */
  std::unique_ptr<FluidInterface> fluid_interface_;
};  // end of class LinearMolecularKineticsPhysics

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_LINEAR_MOLECULAR_KINETICS_LINEAR_MOLECULAR_KINETICS_PHYSICS_H_
