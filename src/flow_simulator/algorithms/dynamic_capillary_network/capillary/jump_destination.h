/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/capillary/jump_destination.h
 * \brief Contains the \c JumpDestination class.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 *
 * This header file contains the \c JumpDestination class.
 */
#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_CAPILLARY_JUMP_DESTINATION_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_CAPILLARY_JUMP_DESTINATION_H_

#include <vector>
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network_context.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/capillary_interface_enums.h"

namespace simulator {

/**
 * \class JumpDestination jump_destination.h "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/jump_destination.h"
 * \brief Class to store informations about a jump destination in a dynamic capillary network simulation.
 *
 * Stores the destination capillary index, the arriving side and the type of a jump.
 */

class JumpDestination {
 public:
  /// Parametrised Constructor
  JumpDestination(const IndexType capillary, const Side jump_side, const JumpType jump_type)
    : destination_capillary_(capillary),
      destination_side_(jump_side),
      destination_type_(jump_type) { }

  /// Gets the capillary index of the given destination
  IndexType GetCapillary(void) const { return destination_capillary_; }

  /// Gets the arriving side of the jump in a given destination capillary
  Side GetSide(void) const { return destination_side_; }

  /// Gets the effect type of the interface jump in the destination
  JumpType GetType(void) const { return destination_type_; }

 private:
  /**
   * \brief Index of the capillary to which the jump is destined
   *
   * Contains the index of the capillary to which the jump is destined.
   */
  IndexType destination_capillary_;

  /**
   * \brief Side of the destination capillary at which the jump arrives
   *
   * Contains the side of the destination capillary at which the jump arrives, according to \c Side.
   */
  Side destination_side_;

  /**
   * \brief Type of the jump (see \c JumpType)
   *
   * Contains the type of the jump, according to \c JumpType.
   */
  JumpType destination_type_;
};  // end of class JumpDestination

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_CAPILLARY_JUMP_DESTINATION_H_
