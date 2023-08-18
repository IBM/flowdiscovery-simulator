/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/capillary/interface_jump.h
 * \brief Contains the \c InterfaceJump class.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 *
 * This header file contains the \c InterfaceJump class.
 */
#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_CAPILLARY_INTERFACE_JUMP_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_CAPILLARY_INTERFACE_JUMP_H_

#include <vector>
#include "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/jump_destination.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network_context.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/capillary_interface_enums.h"

namespace simulator {

/**
 * \class InterfaceJump interface_jump_event_info.h "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/interface_jump_event_info.h"
 * \brief Class to store and manage informations of a jump in a dynamic capillary network simulation.
 *
 * Stores the jump information, such as its origin and destination properties.
 */

class InterfaceJump {
 public:
  /// Parametrised Constructor
  InterfaceJump(const IndexType capillary,
                const Side jump_side,
                const IndexType node,
                const bool is_cross,
                const bool first_in_cross = false)
    : origin_capillary_(capillary),
      origin_side_(jump_side),
      jump_node_(node),
      cross_jump_(is_cross),
      first_in_cross_(first_in_cross) { }

  /// Custom comparison operator required by \c std::find
  inline bool operator== (const IndexType &other) const noexcept {
    return (origin_capillary_ == other);
  }

  /// Return the capillary index of the origin of this jump
  IndexType GetOrigin(void) const { return origin_capillary_; }

  /// Return the side of this jump
  Side GetSide(void) const { return origin_side_; }

  /// Add a destination with all infos to this jump
  IndexType AddDestination(const IndexType destination_capillary,
                           const Side arriving_side,
                           const JumpType jump_type);

  /// Gets the total number of destinations of this interface jump
  IndexType GetNumberOfDestinations(void) const { return destinations_.size(); }

  /// Gets the capillary index of the given destination
  IndexType GetDestinationCapillary(const IndexType destination_index) const {
    return destinations_.at(destination_index).GetCapillary();
  }

  /// Gets the arriving side of the jump in a given destination capillary
  Side GetDestinationSide(const IndexType destination_index) const {
    return destinations_.at(destination_index).GetSide();
  }

  /// Gets the effect type of the interface jump in the destination
  JumpType GetDestinationType(const IndexType destination_index) const {
    return destinations_.at(destination_index).GetType();
  }

  /// Return true if this jump is a cross jump
  bool IsCrossJump(void) const { return cross_jump_; }

  /// Marks this jump as a cross jump
  void SetCrossed(const bool is_first_in_cross);

  /// Marks this jump as first in the cross.
  bool IsFirstInCross(void) const { return first_in_cross_; }

  /// Return the node index where the jump is occurring
  IndexType GetNode(void) const { return jump_node_; }

 private:
  /**
   * \brief Capillary index of the origin of this jump
   *
   * Contains the index of the capillary where the jump originated from
   */
  IndexType origin_capillary_;

  /**
   * \brief Side of the origin capillary where this jump occurred
   *
   * Contains the side of the capillary where this jump occurred
   */
  Side origin_side_;

  /**
   * \brief Node where the jump is taking effect
   *
   * Contains de index of the jump where the jumping is occurring
   */
  IndexType jump_node_;

  /**
   * \brief Identifies if this jump is crossed with another one at the same node
   *
   * Contains a bool identifying if it is a cross jump
   */
  bool cross_jump_;

  /**
   * \brief Identifies if this jump is the first in the cross
   */
  bool first_in_cross_;

  /**
   * \brief \c Vector of JumpDestination to store all the destinations of a jump
   *
   * Contains the destinations of this jump
   */
  std::vector<JumpDestination> destinations_;
};  // end of class InterfaceJump

}  // namespace simulator
#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_CAPILLARY_INTERFACE_JUMP_H_
