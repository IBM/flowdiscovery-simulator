/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/models/dynamic_model_base.h
 * \brief Base \c Model class containing virtual methods definition.
 *
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 *
 * This header file contains the base class for model implementation
 */
#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_DYNAMIC_MODEL_BASE_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_DYNAMIC_MODEL_BASE_H_

#include <memory>
#include "src/flow_simulator/algorithms/dynamic_capillary_network/models/dynamic_capillary_physics_base.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/models/dynamic_capillary_geometry_base.h"

namespace simulator {

class DynamicModelBase {
 public:
  /// Virtual destructor
  virtual ~DynamicModelBase() { }

  /// Getter for shared pointer to DynamicCapillaryPhysicsBase
  virtual std::shared_ptr<DynamicCapillaryPhysicsBase> GetPhysics(void) = 0;

  /// Getter for shared pointer to DynamicCapillaryGeometryBase
  virtual std::shared_ptr<DynamicCapillaryGeometryBase> GetGeometry(void) = 0;
};

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_DYNAMIC_MODEL_BASE_H_
