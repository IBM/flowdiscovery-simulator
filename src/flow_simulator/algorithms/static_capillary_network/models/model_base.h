/**
 * \file src/flow_simulator/algorithms/static_capillary_network/models/model_base.h
 * \brief Base \c Model class containing virtual methods definition.
 *
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 *
 * This header file contains the base class for model implementation
 */
#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_MODELS_MODEL_BASE_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_MODELS_MODEL_BASE_H_

#include <memory>
#include "src/flow_simulator/algorithms/static_capillary_network/models/static_capillary_physics_base.h"
#include "src/flow_simulator/algorithms/static_capillary_network/models/static_capillary_geometry_base.h"

namespace simulator {

class ModelBase {
 public:
  /// Virtual destructor
  virtual ~ModelBase() { }

  /// Getter for shared pointer to StaticCapillaryPhysicsBase
  virtual std::shared_ptr<StaticCapillaryPhysicsBase> GetPhysics(void) = 0;

  /// Getter for shared pointer to StaticCapillaryGeometryBase
  virtual std::shared_ptr<StaticCapillaryGeometryBase> GetGeometry(void) = 0;
};

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_STATIC_CAPILLARY_NETWORK_MODELS_MODEL_BASE_H_
