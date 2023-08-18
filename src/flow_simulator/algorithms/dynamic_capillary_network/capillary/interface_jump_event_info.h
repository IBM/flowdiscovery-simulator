/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/capillary/interface_jump_event_info.h
 * \brief Contains the \c InterfaceJumpInfo class.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 *
 * This header file contains the \c InterfaceJumpEventInfo class.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_CAPILLARY_INTERFACE_JUMP_EVENT_INFO_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_CAPILLARY_INTERFACE_JUMP_EVENT_INFO_H_

#include <vector>
#include "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/capillary_interface_enums.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/interface_jump.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network_context.h"

namespace simulator {

using interface_jump_vector_t = std::vector<InterfaceJump>;
using interface_jump_vector_iterator_t = interface_jump_vector_t::iterator;
using interface_jump_vector_const_iterator_t = interface_jump_vector_t::const_iterator;

/**
 * \class InterfaceJumpEventInfo interface_jump_event_info.h "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/interface_jump_event_info.h"
 * \brief Class to store and manage jump informations of a given interface jump event in a dynamic
 *  capillary network simulation.
 */
class InterfaceJumpEventInfo {
 public:
  /// Clear all jump informations stored
  void Clear(void) { jumps_.clear(); }

  /// Insert a Jump that occurred in this jump event and returns its index.
  IndexType InsertJump(const IndexType capillary, const Side jump_side, const IndexType node);

  /// Create and insert a Jump into the jumps_ array (used to resume the analysis).
  void AddJumpToJumpVector(const IndexType capillary,
                           const Side jump_side,
                           const IndexType node,
                           bool is_cross,
                           bool first_in_cross);

  /// Add a destination with all infos to a given jump.
  IndexType AddDestination(const IndexType jump_index,
                           const IndexType destination_capillary,
                           const Side arriving_side,
                           const JumpType jump_type);

  /// Return the index of a jump that was originated from a given capillary.
  IndexType GetJumpIndex(const IndexType capillary) const;

  /// Check if a given capillary is origin of a jump.
  bool IsJumpOrigin(const IndexType capillary) const {
    return std::find(jumps_.begin(), jumps_.end(), capillary) != jumps_.end();
  }

  /// Return the capillary index of the origin of the jump.
  IndexType GetJumpOrigin(const IndexType jump_index) const {
    return jumps_.at(jump_index).GetOrigin();
  }

  /// Return the number of jumps registered in the event.
  IndexType GetNumberOfJumps(void) const { return jumps_.size(); }

  /// Return the number of destinations of the jump of index \c jump_index.
  IndexType GetNumberOfDestinations(const IndexType jump_index) const {
    return jumps_.at(jump_index).GetNumberOfDestinations();
  }

  /// Return the capillary index of a given destination of a given jump.
  IndexType GetDestinationCapillary(const IndexType jump_index,
                                    const IndexType destination_index) const {
    return jumps_.at(jump_index).GetDestinationCapillary(destination_index);
  }

  /// Return the arriving side of a given destination of a given jump.
  Side GetDestinationSide(const IndexType jump_index,
                          const IndexType destination_index) const {
    return jumps_.at(jump_index).GetDestinationSide(destination_index);
  }

  /// Return the jump type of a given destination of a given jump.
  JumpType GetDestinationType(const IndexType jump_index,
                              const IndexType destination_index) const {
    return jumps_.at(jump_index).GetDestinationType(destination_index);
  }

  /// Return true if the jump of index \c jump_index is a cross jump;
  bool IsCrossJump(const IndexType jump_index) const { return jumps_.at(jump_index).IsCrossJump(); }

  /// Return the first item iterator to access the \c InterfaceJump elements
  interface_jump_vector_iterator_t begin(void) { return jumps_.begin(); }

  /// Return the last item iterator to access the \c InterfaceJump elements
  interface_jump_vector_iterator_t end(void) { return jumps_.end(); }

  /// Return the first item iterator to access the \c InterfaceJump elements
  interface_jump_vector_const_iterator_t cbegin(void) { return jumps_.cbegin(); }

  /// Return the last item iterator to access the \c InterfaceJump elements
  interface_jump_vector_const_iterator_t cend(void) { return jumps_.cend(); }

 private:
  /**
   * \brief \c Vector of InterfaceJump to store all the jumps of a jump event
   *
   * Contains the jumps of a jump event
   */
  interface_jump_vector_t jumps_;
};  // end of class InterfaceJumpEventInfo

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_CAPILLARY_INTERFACE_JUMP_EVENT_INFO_H_
