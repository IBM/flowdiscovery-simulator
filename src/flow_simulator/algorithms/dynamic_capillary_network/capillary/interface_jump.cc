/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/capillary/interface_jump.cc
 * \brief Contains the \c InterfaceJump class implementation.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 *
 * This file contains the implementation of \c InterfaceJump classe.
 */

#include "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/interface_jump.h"

namespace simulator {

IndexType InterfaceJump::AddDestination(const IndexType destination_capillary,
                                        const Side arriving_side,
                                        const JumpType jump_type) {
/**
 * Add a destination of this jump with all its informations.
 *
 * \param[in]  destination_capillary  Capillary index of the destination of the jump
 * \param[in]  arriving_side          Side where the jump arrives in the destination capillary
 * \param[in]  jump_type              Type of the jump, indicating the effects of the given jump
 * \retval     destination_index      Index of the new destination in the \c std::vector
 */
  auto last_pos = destinations_.end();
  destinations_.insert(last_pos, JumpDestination(destination_capillary, arriving_side, jump_type));
  return static_cast<IndexType>(destinations_.size() - 1U);
}  // end of InterfaceJump::AddDestination() method



void InterfaceJump::SetCrossed(const bool is_first_in_cross) {
/**
 * Marks this jump as a cross jump.
 * \param[in]  is_first_in_cross        Boolean identifying if this is the first jump in cross
 */
  cross_jump_ = true;
  first_in_cross_ = is_first_in_cross;
}  // end of InterfaceJump::SetCrossed() method

}  // namespace simulator
