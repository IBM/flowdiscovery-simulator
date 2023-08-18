/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/models/linear_molecular_kinetics/linear_molecular_kinetics_geometry.h
 * \brief Contains the \c LinearMolecularKineticsGeometry class.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 *
 * This header file contains the the \c LinearMolecularKineticsGeometry class that derives from
 * \c DynamicCapillaryGeometryBase.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_LINEAR_MOLECULAR_KINETICS_LINEAR_MOLECULAR_KINETICS_GEOMETRY_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_LINEAR_MOLECULAR_KINETICS_LINEAR_MOLECULAR_KINETICS_GEOMETRY_H_

#include <string>
#include <memory>
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network_context.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/models/dynamic_capillary_geometry_base.h"

namespace simulator {

/**
 * \class LinearMolecularKineticsGeometry linear_molecular_kinetics_geometry.h "src/flow_simulator/algorithms/dynamic_capillary_network/models/linear_molecular_kinetics/linear_molecular_kinetics_geometry.h"
 * \brief Linear Molecular Kinetics geometry for the Dynamic Capillary Network algorithm.
 *
 * The Linear Molecular Kinetics geometry employs only link lengths and radii.
 */

class LinearMolecularKineticsGeometry : public DynamicCapillaryGeometryBase {
 public:
  /// Parametrised constructor
  LinearMolecularKineticsGeometry(const std::string &folder,
                                  std::shared_ptr<DynamicCapillaryNetworkContext> context)
    : DynamicCapillaryGeometryBase(folder, context) { }
};  // end of class LinearMolecularKineticsGeometry

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_MODELS_LINEAR_MOLECULAR_KINETICS_LINEAR_MOLECULAR_KINETICS_GEOMETRY_H_
