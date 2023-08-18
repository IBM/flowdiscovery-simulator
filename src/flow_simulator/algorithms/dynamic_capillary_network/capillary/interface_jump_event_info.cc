/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/capillary/interface_jump_event_info.cc
 * \brief Contains the \c InterfaceJumpEventInfo class implementation.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 *
 * This file contains the implementation of \c InterfaceJumpEventInfo class.
 */

#include "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/interface_jump_event_info.h"

namespace simulator {

IndexType InterfaceJumpEventInfo::InsertJump(const IndexType capillary,
                                             const Side jump_side,
                                             const IndexType node) {
/**
 * Insert a new jump information in the current event.
 *
 * \param[in] capillary    Index of the target capillary.
 * \param[in] jump_side    Side of the capillary from where the jump is originated.
 * \param[in] node         Index of the node where the jump is ocurring
 * \retval    jump_index   Index of the jump in this jump event.
 */
  bool isCrossed = false;
  bool is_first_cross = true;
  for (auto &jump : jumps_) {
    if (jump.GetNode() == node) {
      isCrossed = true;
      jump.SetCrossed(is_first_cross);
      is_first_cross = false;
    }
  }
  jumps_.insert(jumps_.end(), InterfaceJump(capillary, jump_side, node, isCrossed));

  return static_cast<IndexType>(jumps_.size() - 1U);
}  // end of InterfaceJumpEventInfo::InsertJump() method



void InterfaceJumpEventInfo::AddJumpToJumpVector(const IndexType capillary,
                                        const Side jump_side,
                                        const IndexType node,
                                        const bool is_cross,
                                        const bool first_in_cross) {
/**
 * Create and insert a new jump information in the jumps_ array (resume analysis).
 *
 * \param[in] capillary         Index of the target capillary.
 * \param[in] jump_side         Side of the capillary from where the jump is originated.
 * \param[in] node              Index of the node where the jump is ocurring
 * \param[in] is_cross          Identifies if this jump is crossed with another one at the same node
 * \param[in] first_in_cross    Identifies if this jump is the first in the cross
 */
  InterfaceJump jump(capillary, jump_side, node, is_cross, first_in_cross);
  jumps_.insert(jumps_.end(), jump);
}  // end of InterfaceJumpEventInfo::CreateJump() method



IndexType InterfaceJumpEventInfo::AddDestination(const IndexType jump_index,
                                                 const IndexType destination_capillary,
                                                 const Side arriving_side,
                                                 const JumpType jump_type) {
/**
 * Add a destination of the jump of index \c jump_index with its information.
 *
 * \param[in]  jump_index             Index identifying the origin of the jump
 * \param[in]  destination_capillary  Capillary index of the destination of the jump.
 * \param[in]  arriving_side          Side where the jump arrives in the destination capillary.
 * \param[in]  jump_type              Type of the jump, indicating the effects of the given jump.
 * \retval     destination_index      Index of the new destination in the \c std::vector.
 */
  return jumps_.at(jump_index).AddDestination(destination_capillary, arriving_side, jump_type);
}  // end of InterfaceJumpEventInfo::AddDestination() method



IndexType InterfaceJumpEventInfo::GetJumpIndex(const IndexType capillary) const {
/**
 * Return the index of a jump that was originated from a given capillary.
 * \param[in]  capillary     Index of the target capillary.
 * \retval     jump_index    Index of the jump in this jump event.
 */
  auto it = std::find(jumps_.begin(), jumps_.end(), capillary);
  if (it == jumps_.end()) {
    std::fprintf(stderr, "Invalid access at InterfaceJumpEventInfo::GetJumpIndex.\n");
    std::exit(-1);
  }
  return static_cast<IndexType>(std::distance(jumps_.begin(), it));
}  // end of InterfaceJumpEventInfo::GetJumpIndex() method

}  // namespace simulator
