/**
 * \file src/flow_simulator/algorithms/static_capillary_network/models/hydrostatic.h
 * \brief Contains the \c HydrostaticPhysics, \c HydrostaticGeometry, and \c HydrostaticModel
 * classes.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 *
 * This header file contains the \c HydrostaticPhysics class that derives from
 * \c StaticCapillaryPhysicsBase, the \c HydrostaticGeometry class that derives from
 * \c StaticCapillaryGeometryBase, and the \c HydrostaticModel class that
 * derives from \c ModelBase.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_MODELS_HYDROSTATIC_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_MODELS_HYDROSTATIC_H_

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

namespace simulator {

/**
 * \class HydrostaticGeometry hydrostatic.h "src/flow_simulator/algorithms/static_capillary_network/models/hydrostatic.h"
 * \brief Hydrostatic geometry for the Capillary Network algorithm using centerlines.
 *
 * The Hydrostatic geometry employs link lengths, radii and orientations to calculate the
 * geometrical factors.
 */

class HydrostaticGeometry : public StaticCapillaryGeometryBase {
 public:
  /// Parametrised constructor
  HydrostaticGeometry(const std::string &folder,
                      std::shared_ptr<StaticCapillaryNetworkContext> context)
    : StaticCapillaryGeometryBase(folder, context) {}

  /// Getter for gamma vector
  const arma::Col<double> &GetGammaVector(void) { return gamma_; }

  /// Calculates gamma vector
  void CalculateGammaVector(void);

 private:
  /// Calculates distance vectors between linked nodes
  arma::Mat<double> CalculateInternodeDistanceVector(void) const;

  /**
   * \brief Direction cosine with respect to vertical
   *
   * Each \f$ \gamma_j \f$ is defined as the cosine of the angle \f$ \theta_j \f$ between the line
   * that defines the capillary and the vertical (z axis). It is calculated by
   * \c HydrostaticGeometry::CalculateGammaVector.
   */
  arma::Col<double> gamma_;
};  // end of class HydrostaticGeometry



/**
 * \class HydrostaticPhysics hydrostatic.h "src/flow_simulator/algorithms/static_capillary_network/models/hydrostatic.h"
 * \brief Hydrostatic Physics for the Static Capillary Network algorithm using centerlines.
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
 *
 * Considering the action of gravity, we add the hydrostatic contribution to the pressure drop
 *
 * \f[
 *    h_j = \rho g \gamma_j L_j \, ,
 * \f]
 *
 * where \f$ \rho \f$ is the fluid density, \f$ g \f$ is the gravitational acceleration, \f$ L_j \f$
 * is the capillary length and \f$ \gamma_j \f$ is the vertical direction cosine
 * \f$ \gamma_j = cos(\theta_j) = \left( \frac{L_j^z}{L_j} \right) \f$.
 *
 * The final equation for the hydrostatic physics model reads
 *
 * \f[
 *    Q + G h = G \Delta P
 * \f]
 */

class HydrostaticPhysics : public StaticCapillaryPhysicsBase {
 public:
  HydrostaticPhysics(std::unique_ptr<Fluid> &&fluid,
                     const std::shared_ptr<HydrostaticGeometry> geometry,
                     std::shared_ptr<StaticCapillaryNetworkContext> context)
    : StaticCapillaryPhysicsBase(context),
      fluid_(std::move(fluid)),
      geometry_(geometry) { }

  /// Specialises \c CalculateRightHandSideVectorTerm for Hydrostatic
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
   * \brief Hydrostatic fluid instance
   *
   * The \c Fluid requires the definition of dynamic viscosity and density
   */
  std::unique_ptr<Fluid> fluid_;

  /**
   * \brief Pointer to \c HydrostaticGeometry
   *
   * Contains a shared pointer to \c HydrostaticGeometry
   */
  const std::shared_ptr<HydrostaticGeometry> geometry_;

  /**
   * \brief Acceleration of gravity
   *
   * The value is expressed in SI units \f$ \left[ \frac{\text{m}}{\text{s}^2} \right] \f$.
   */
  const double gravity_ = 9.80665;

  /**
   * \brief Vector with the hydrostatic pressure contribution
   *
   * The hydrostatic vector contains the hydrostatic pressure contribution (in \c [Pa]) for
   * each capillary. It is calculated by \c HydrostaticPhysics::CalculateRightHandSideVectorTerm().
   */
  arma::Col<double> hydrostatic_vector_;
};  // end of class HydrostaticPhysics



/**
 * \class HydrostaticModel hydrostatic.h "src/flow_simulator/algorithms/static_capillary_network/models/hydrostatic.h"
 * \brief Hydrostatic model for the Static Capillary Network algorithm.
 *
 * Container-like class that hosts the appropriate \c Physics and \c Geometry objects for the
 * Hydrostatic model.
 */
class HydrostaticModel : public ModelBase {
 public:
  HydrostaticModel(const std::string &folder,
                   std::shared_ptr<StaticCapillaryNetworkContext> context,
                   std::vector<FluidJSON> &fluids_json,
                   ExperimentJSON &experiment_json) {
    // Create and populate resident fluid
    std::unique_ptr<Fluid> my_fluid = std::make_unique<Fluid>(fluids_json[0], experiment_json);

    // Create geometry
    geometry_ = std::make_shared<HydrostaticGeometry>(folder, context);

    // Create physics
    physics_ = std::make_shared<HydrostaticPhysics>(std::move(my_fluid), geometry_, context);
  }

  /// Getter for shared pointer to StaticCapillaryGeometryBase
  std::shared_ptr<StaticCapillaryGeometryBase> GetGeometry(void) { return geometry_; }

  /// Getter for shared pointer to StaticCapillaryPhysicsBase
  std::shared_ptr<StaticCapillaryPhysicsBase> GetPhysics(void) { return physics_; }

 protected:
  /**
   * \brief Shared pointer to \c HydrostaticGeometry instance.
   *
   * This member variable stores a shared pointer to the \c HydrostaticGeometry
   * object set by the \c HydrostaticModel constructor.
   */
  std::shared_ptr<HydrostaticGeometry> geometry_;

  /**
   * \brief Shared pointer to \c HydrostaticPhysics instance.
   *
   * This member variable stores a shared pointer to the \c HydrostaticPhysics
   * object set by the \c HydrostaticModel constructor.
   */
  std::shared_ptr<HydrostaticPhysics> physics_;
};  // end of class HydrostaticModel

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_MODELS_HYDROSTATIC_H_
