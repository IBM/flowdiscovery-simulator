/**
 * \file src/flow_simulator/algorithms/static_capillary_network/models/poiseuille.h
 * \brief Contains the \c PoiseuillePhysics, \c PoiseuilleGeometry classes.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This header file contains the \c PoiseuillePhysics class that derives from
 * \c StaticCapillaryPhysicsBase, the \c PoiseuilleGeometry class that derives from
 * \c StaticCapillaryGeometryBase, and the \c PoiseuilleModel class that derives from \c ModelBase.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_MODELS_POISEUILLE_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_MODELS_POISEUILLE_H_

#include <armadillo>
#include <map>
#include <string>
#include <memory>
#include <utility>
#include <vector>
#include "src/flow_simulator/algorithms/static_capillary_network/models/static_capillary_physics_base.h"
#include "src/flow_simulator/algorithms/static_capillary_network/models/static_capillary_geometry_base.h"
#include "src/flow_simulator/algorithms/static_capillary_network/models/model_base.h"
#include "src/flow_simulator/algorithms/fluid.h"
#include "src/exec_manager/experiment_json.h"
#include "src/exec_manager/fluid_json.h"

using Fluid = simulator::Fluid;

namespace simulator {

/**
 * \class PoiseuilleGeometry poiseuille.h "src/flow_simulator/algorithms/static_capillary_network/models/poiseuille.h"
 * \brief Poiseuille geometry for the Capillary Network algorithm using centerlines.
 *
 * The Poiseuille geometry employs only link lengths and radii to calculate the geometrical factors.
 */

class PoiseuilleGeometry : public StaticCapillaryGeometryBase {
 public:
  /// Parametrised constructor
  PoiseuilleGeometry(const std::string &folder,
                     std::shared_ptr<StaticCapillaryNetworkContext> context)
    : StaticCapillaryGeometryBase(folder, context) {}
};  // end of class PoiseuilleGeometry



/**
 * \class PoiseuillePhysics poiseuille.h "src/flow_simulator/algorithms/static_capillary_network/models/poiseuille.h"
 * \brief Poiseuille Physics for the Static Capillary Network algorithm using centerlines.
 *
 * The Hagen-Poiseuille equation is a physical law that relates the pressure drop \f$ \Delta P \f$
 * and the flow rate \f$ Q \f$ of an incompressible and Newtonian fluid of dynamic viscosity
 * \f$ \mu \f$, while in laminar flow through a cylindrical capillary of length \f$ L \f$ and
 * constant radius \f$ R \f$.
 *
 * \f[
 *    Q = \frac{\pi R^4}{8 \mu L} \Delta P = \mu^{-1} k \Delta P = G \Delta P
 * \f]
 *
 * where \f$ k = \frac{\pi R^4}{8 L} \f$ is the geometrical factor and \f$ G = \mu^{-1} k \f$ is the
 * geometry matrix element.
 */

class PoiseuillePhysics : public StaticCapillaryPhysicsBase {
 public:
  PoiseuillePhysics(std::unique_ptr<Fluid> &&fluid,
                    const std::shared_ptr<PoiseuilleGeometry> geometry,
                    std::shared_ptr<StaticCapillaryNetworkContext> context)
    : StaticCapillaryPhysicsBase(context),
      fluid_(std::move(fluid)),
      geometry_(geometry) { }

  /// Specialises \c CalculateRightHandSideVectorTerm for Poiseuille
  arma::Col<double> CalculateRightHandSideVectorTerm(void);

  /// Calculates the flow rate at each capillary
  arma::Col<double> CalculateFlowRate(const arma::Col<double> &pressures) const;

  /// Getter for fluid viscosity
  double GetFluidViscosity(const double pressure) const {
    return fluid_->GetViscosity(pressure);
  }

  /// Return true if the parameters fall between the designed behaviour limits, false otherwise
  double IsViscosityBehaviourDesignedForInjectedFluid(const double pressure) const override {
    return fluid_->IsViscosityBehaviourDesignedFor(pressure, pressure);
  }

 private:
  /**
   * \brief Poiseuille fluid instance
   *
   * The \c Fluid requires the definition of dynamic viscosity.
   */
  std::unique_ptr<Fluid> fluid_;

  /**
   * \brief Pointer to \c PoiseuilleGeometry
   *
   * Contains a shared pointer to \c PoiseuilleGeometry
   */
  const std::shared_ptr<PoiseuilleGeometry> geometry_;
};  // end of class PoiseuillePhysics



/**
 * \class PoiseuilleModel poiseuille.h "src/flow_simulator/algorithms/static_capillary_network/models/poiseuille.h"
 * \brief Poiseuille model for the Static Capillary Network algorithm.
 *
 * Container-like class that hosts the appropriate \c Physics and \c Geometry objects for the
 * Poiseuille model.
 */
class PoiseuilleModel : public ModelBase {
 public:
  PoiseuilleModel(const std::string &folder,
                  std::shared_ptr<StaticCapillaryNetworkContext> context,
                  std::vector<FluidJSON> &fluids_json,
                  ExperimentJSON &experiment_json) {
    // Create and populate resident fluid
    std::unique_ptr<Fluid> my_fluid = std::make_unique<Fluid>(fluids_json[0], experiment_json);

    // Create geometry
    geometry_ = std::make_shared<PoiseuilleGeometry>(folder, context);

    // Create physics
    physics_ = std::make_shared<PoiseuillePhysics>(std::move(my_fluid), geometry_, context);
  }

  /// Getter for shared pointer to StaticCapillaryGeometryBase
  std::shared_ptr<StaticCapillaryGeometryBase> GetGeometry(void) { return geometry_; }

  /// Getter for shared pointer to StaticCapillaryPhysicsBase
  std::shared_ptr<StaticCapillaryPhysicsBase> GetPhysics(void) { return physics_; }

 protected:
  /**
   * \brief Shared pointer to \c PoiseuilleGeometry instance.
   *
   * This member variable stores a shared pointer to the \c PoiseuilleGeometry
   * object set by the \c PoiseuilleModel constructor.
   */
  std::shared_ptr<PoiseuilleGeometry> geometry_;

  /**
   * \brief Shared pointer to \c PoiseuillePhysics instance.
   *
   * This member variable stores a shared pointer to the \c PoiseuillePhysics
   * object set by the \c PoiseuilleModel constructor.
   */
  std::shared_ptr<PoiseuillePhysics> physics_;
};  // end of class PoiseuilleModel

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_MODELS_POISEUILLE_H_
