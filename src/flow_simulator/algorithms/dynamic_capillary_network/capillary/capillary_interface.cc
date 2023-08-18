/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/capillary/capillary_interface.cc
 * \brief Contains the \c CapillaryInterface class.
 *
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 *
 * This file contains the \c CapillaryInterface implementation.
 */

#include "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/capillary_interface.h"
#include <glog/logging.h>
#include <cassert>
#include "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/interface_jump.h"

namespace simulator {

void CapillaryInterface::InitialiseCapillaryInterface(void) {
/**
 * The \c InitialiseCapillaryInterface() method correctly initialises all variables at
 * correct sizes.
 */
  SetSizes(context_->number_of_links_);
}  // end of CapillaryInterface::InitialiseCapillaryInterface() method



void CapillaryInterface::SetSizes(const IndexType &n_links) {
/**
 * The \c SetSizes method set the size of class members \c all_interfaces_,
 * \c fluid_at_source_ and \c deltas_.
 *
 * \param[in] n_links   Number of links between pairs of nodes.
 */
  all_interfaces_.set_size(n_links);
  fluid_at_source_.zeros(n_links);
  deltas_.zeros(n_links);
  plugs_.zeros(n_links);
}  // end of CapillaryInterface::SetSizes() method



void CapillaryInterface::InjectFluidOnInletNodes(const arma::uword number_interfaces,
                                                 const double distance_interfaces) {
/**
 * The \c InjectFluidOnInletNodes() look for all capillaries in the dynamic
 * capillary network and inject a fluid in the position of the boundary node.
 *
 * \param[in] number_interfaces    Number of interfaces to inject. Default: 1.
 * \param[in] distance_interfaces  Normalised distance between interfaces. Default: 0.
 */
  // Make sure all interfaces fit in a single capillary
  assert(number_interfaces * distance_interfaces < 1.0);

  for (const auto &inlet_node : context_->inlet_nodes_) {
    // Apply the change to all capillaries whose source is inlet_node
    for (const auto &candidate_capillary : context_->capillaries_whose_source_is_(inlet_node)) {
      // Inject N interfaces at the source node of the capillary
      for (arma::uword n = 0U; n != number_interfaces; ++n) {
        MoveInterfacePositionsAtCapillaryByDelta(candidate_capillary, distance_interfaces);
        AddNewInterface(candidate_capillary, Side::source);
        FlipFluidAtSource(candidate_capillary);
      }
    }
    // Apply the change to all capillaries whose target is inlet_node
    for (const auto &candidate_capillary : context_->capillaries_whose_target_is_(inlet_node)) {
      // Inject N interfaces at the target node of the capillary
      for (arma::uword n = 0U; n != number_interfaces; ++n) {
        MoveInterfacePositionsAtCapillaryByDelta(candidate_capillary, -distance_interfaces);
        AddNewInterface(candidate_capillary, Side::target);
      }
    }
  }
}  // end of CapillaryInterface::InjectFluidOnInletNodes() method



void CapillaryInterface::UpdateInterfaceInformation(void) {
/**
 * Updates positions and interface 'jump' and change the fluid at Source of affected capillaries.
 * Considering all (\f$ n \f$) differential equations are set first in \c x_array,
 * the first \f$ n \f$ elements are related to the \f$ n \f$ capillaries in \c all_interfaces_.
 */
  // TODO(ralves): parallelisation is possible
  // Move all interfaces by their respective delta
  for (auto capillary : capillaries_with_interfaces_) {
    if (!IsCapillaryPlugged(capillary)) {
      all_interfaces_(capillary) += deltas_(capillary);
    }
  }
}  // end of CapillaryInterface::UpdateInterfaceInformation() method



void CapillaryInterface::RevertUpdateInterfaceInformation(void) {
/**
 * Undo changes made by the \c UpdateInterfaceInformation() method.
 */
  // TODO(ralves): parallelisation is possible
  for (auto capillary : capillaries_with_interfaces_) {
    all_interfaces_(capillary) -= deltas_(capillary);
  }
}  // end of CapillaryInterface::RevertUpdateInterfaceInformation() method



void CapillaryInterface::MoveInterfacePositionsAtCapillaryByDelta(const IndexType capillary,
                                                                  const double delta) {
/**
 * Move all interfaces of the capillary by a given delta.
 *
 * \param[in] capillary   Index of the target capillary.
 * \param[in] delta       Double value that will be used to move the interfaces.
 */
  all_interfaces_(capillary) += delta;
}  // end of CapillaryInterface::MoveInterfacePositionsAtCapillaryByDelta() method



arma::Col<IndexType> CapillaryInterface::PlugCapillariesJumpAtInletNode(
  const arma::Col<IndexType> &capillaries,
  const arma::Col<IndexType> &caps_roots_directions) {
/**
 * Plug capillaries whose interfaces reach the inlet nodes. It is crucial to save the
 * information of which capillaries have been plugged because this state can be reverted.
 *
 * \param[in] capillaries            Capillaries that were caught at the jump event.
 * \param[in] caps_roots_directions  Inform the direction of the jump event. 0 for source, 1 for
 *                                   target.
 * \retval caps_plugged   A list of capillaries that have been plugged in this method.
 */
  // Resets last inlet plug info
  last_plug_inlet_info_.reset();

  // Identify all capillaries with an interface jumping through inlet nodes
  arma::Col<IndexType> caps_plugged;
  caps_plugged.set_size(0);
  IndexType i = 0;
  for (const auto &capillary : capillaries) {
    // If caps_roots_directions is positive so the jump is in the target, otherwise source.
    const Side jump_side = caps_roots_directions(i) == 1 ? Side::target : Side::source;

    LOG(INFO) << "Check Plugging of jump through inlet node: [" << capillary << "], jump_side["
        << (jump_side == Side::source ? "Source" : "Target") << "].";

    const IndexType event_node = GetNodeIndexAtCapillaryAtSide(capillary, jump_side);

    // Plug capillary wich jump occurs through inlet node
    bool is_event_node_an_inlet_node = arma::any(context_->inlet_nodes_ == event_node);
    if (is_event_node_an_inlet_node) {
      // Plugs the capillary in the correct end
      PlugCapillary(capillary, jump_side);

      caps_plugged.insert_rows(caps_plugged.n_elem, arma::Col<IndexType>({capillary}));

      // Stores the plug information to be used in the revert step.
      last_plug_inlet_info_.insert_rows(last_plug_inlet_info_.n_elem,
                                        arma::Row<IndexType>({capillary}));
      LOG(INFO) << "Capillary[" << capillary << "] plugged.";

      // Move the interfaces to `x = epsilon` or `x = 1-epsilon`, depending on `jump_side`
      RemoveInterface(capillary, jump_side);
      AddNewInterface(capillary, jump_side);
    } else {
      LOG(INFO) << "Capillary[" << capillary << "] not plugged.";
    }
    ++i;
  }
  return caps_plugged;
}  // end of CapillaryInterface::PluggingCapsWithInletNodes() method



bool CapillaryInterface::CheckPlugging(const arma::Col<IndexType> &capillaries,
                                       const arma::Col<IndexType> &caps_roots_directions) {
/**
 * Check if a capillary will get plugged and plugs it when it does.
 *
 * The condition to a capillary get plugged is an interface jumping back to its origin
 * capillary instantaneously after a jump, meaning that the interface would not be able to jump
 * in the first place. Then it is plugged.
 *
 * \param[in] capillaries            Capillaries that were caught at the jump event.
 * \param[in] caps_roots_directions  Inform the direction of the jump event. 0 for source, 1 for
 *                                   target.
 * \retval plugged                   Boolean stating if a plug occurred.
 */
  // Resets plug info
  last_plug_info_.reset();

  // if no jumps were made, no plugging will occur
  bool plugged = false;
  if (last_jump_info_.GetNumberOfJumps() == 0) {
    LOG(INFO) << "No jumps were made, no plugging will occur";
    return plugged;
  }

  // Checks all capillaries with an interface and check if a jump would occur
  IndexType i = 0;
  for (const auto &capillary : capillaries) {
    // If caps_roots_directions is positive so the jump is in the target, otherwise source.
    const Side jump_side = caps_roots_directions(i) == 1 ? Side::target : Side::source;

    LOG(INFO) << "CheckPlugging:capillary[" << capillary << "], jump_side["
        << (jump_side == Side::source ? "Source" : "Target") << "].";

    // Get the node where the interface is crossing to find the capillaries connected to it
    // and excludes the source capillary
    const IndexType event_node = GetNodeIndexAtCapillaryAtSide(capillary, jump_side);
    arma::uvec other_capillaries(context_->GetCapillariesConnectedToNode(event_node));
    other_capillaries.shed_row(arma::as_scalar(arma::find(other_capillaries == capillary)));

    // TODO(rneumann): Plugging only happens when target capillary is connected to 2+ others
    // Verify all capillaries connected to target capillary
    for (const auto &candidate_capillary : other_capillaries) {
      for (const auto &jump : last_jump_info_) {
        // If the connected capillary was the origin of a previous jump, then it is plugged.
        if (candidate_capillary == jump.GetOrigin()) {
          // If it was a cross jump, plugging does not happens
          if (jump.IsCrossJump()) {
            continue;
          }

          // Plugs the capillary in the correct end
          PlugCapillary(candidate_capillary, GetSideOfNodeAtCapillary(candidate_capillary,
                                                                      event_node));

          // Stores the plug information to be used in the revert step.
          last_plug_info_.insert_rows(last_plug_info_.n_elem,
                                      arma::Row<IndexType>({candidate_capillary}));
          plugged = true;

          LOG(INFO) << "Plug Detected on capillary: " << candidate_capillary;
        }
      }
    }
    ++i;
  }
  return plugged;
}  // end of CapillaryInterface::CheckPlugging() method



void CapillaryInterface::Jump(const arma::Col<IndexType> &capillaries,
                              const arma::Col<IndexType> &caps_roots_directions,
                              arma::Col<IndexType> &changed_caps) {
/**
 * Move interface from one capillary to another based on \c roots_found.
 *
 * \param[in]   capillaries            Capillaries that were caught at the jump event.
 * \param[in]   caps_roots_directions  Inform the direction of the jump event. 0 for source, 1 for
 *                                     target.
 * \param[in, out]  changed_caps       Unique indexes of the capillaries that has been changed. In
 *                                     this method, the unique indexes of capillaries that
 *                                     originates the jumps, and capillaries that receive the
 *                                     interfaces are inserted.
 */
  last_jump_info_.Clear();

  // Checks all capillaries with an interface and executes a jump when necessary
  IndexType j = 0;
  for (const auto &capillary : capillaries) {
    changed_caps.insert_rows(changed_caps.n_elem, arma::Col<IndexType>({ capillary }));

    // An interface in this capillary has reached a node
    // Remove the interface from at the correct end of the source of the jump
    const Side jump_side = (caps_roots_directions(j) == 1) ? Side::target : Side::source;

    // Get the node where the interface is crossing to find the capillaries connected to it
    // and excludes the source capillary
    const IndexType event_node = GetNodeIndexAtCapillaryAtSide(capillary, jump_side);
    arma::uvec other_capillaries(context_->GetCapillariesConnectedToNode(event_node));
    other_capillaries.shed_row(arma::as_scalar(arma::find(other_capillaries == capillary)));

    // Insert this jump to last_jump_info
    const IndexType jump_index = last_jump_info_.InsertJump(capillary, jump_side, event_node);

    // Apply the change to all capillaries connected to target capillary
    for (const auto &candidate_capillary : other_capillaries) {
      // Add or annihilates the interface in the correct end of the capillary
      const Side plugged_side = GetCapillaryPluggedSide(candidate_capillary);
      const Side arriving_side = GetNodeSideAtCapillary(event_node, candidate_capillary);
      const JumpType jump_type = GetJumpType(candidate_capillary, arriving_side, plugged_side);
      last_jump_info_.AddDestination(jump_index, candidate_capillary, arriving_side, jump_type);

      changed_caps.insert_rows(changed_caps.n_elem, arma::Col<IndexType>({ candidate_capillary }));
    }

    ++j;

    // Log relevant information
    std::stringstream sstream;
    sstream << "Jump from " << capillary << " to " << other_capillaries.t();
    LOG(INFO) << "interface:Jump:capillary[" << capillary << "].";
    LOG(INFO) << sstream.str();
  }

  // Applies the jump based on gathered informations
  for (const InterfaceJump &jump : last_jump_info_) {
    const Side jump_side = jump.GetSide();
    const IndexType jump_origin = jump.GetOrigin();

    // Remove the jumping interface from the origin capillary
    RemoveInterface(jump_origin, jump_side);

    // Flip the fluid when jump occurred from the source side
    if (jump_side == Side::source) {
      FlipFluidAtSource(jump_origin);
    }

    // Apply jump to destinations depending if its a cross jump or not
    const IndexType number_of_destinations = jump.GetNumberOfDestinations();
    if (jump.IsCrossJump()) {
      // Only if is the first in the cross node, it applies the jump to non cross jump destinations
      if (jump.IsFirstInCross()) {
        for (IndexType i = 0U; i != number_of_destinations; ++i) {
          IndexType dest_cap = jump.GetDestinationCapillary(i);
          if (!last_jump_info_.IsJumpOrigin(dest_cap)) {
            ApplyJumpToDestination(jump, i);
          }
        }
      }
    } else {
      for (IndexType i = 0U; i != number_of_destinations; ++i) {
        ApplyJumpToDestination(jump, i);
      }
    }
  }
  changed_caps = arma::unique(changed_caps);  // Avoid duplicates
}  // end of CapillaryInterface::Jump() method



void CapillaryInterface::RemoveBubblesFromDestination(const Side side,
                                                      const IndexType capillary) {
/**
 * Interfaces with distances smaller than \c arma::datum::eps from \c side would generate a tiny
 * bubble with a new interface or have tiny bubble(s) already.
 *
 * This method eliminates those bubbles, which simplifies `RevertJump()`.
 *
 * \param[in] side        Side from the capillary that the interface(s) will be removed.
 * \param[in] capillary   Capillary that the interfaces will be removed.
 */
  switch (side) {
    case Side::source: {
      double value = 0.0 + arma::datum::eps;
      const arma::vec cap_interfaces = all_interfaces_(capillary);
      const arma::uvec query_index = arma::find(cap_interfaces <= value);
      for (auto interface : cap_interfaces.elem(query_index).eval()) {
        LOG(INFO) << "Removing interface from Source " << interface << " from capillary "
            << capillary;
        RemoveInterface(capillary, side);
      }
      break;
    }
    case Side::target: {
      double value = 1.0 - arma::datum::eps;
      const arma::vec cap_interfaces = all_interfaces_(capillary);
      const arma::uvec query_index = arma::find(all_interfaces_(capillary) >= value);
      for (auto interface : cap_interfaces.elem(query_index).eval()) {
        LOG(INFO) << "Removing interface from Target " << interface << " from capillary "
            << capillary;
        RemoveInterface(capillary, side);
      }
      break;
    }
    case Side::none: {
      LOG(FATAL) << "Error in CapillaryInterface::RemoveInterfacesFromDestinationLimbo. "
                 << "Side::none is not a valid option to add interfaces";
    }
  }
}  // end of CapillaryInterface::RemoveBubblesFromDestination() method



void CapillaryInterface::ApplyJumpToDestination(const InterfaceJump &jump,
                                                const IndexType destination_index) {
/**
 * Correctly apply the jump to destination using the jump information collected beforehand.
 *
 * \param[in] jump              Iterator in the position of the target \c InterfaceJump.
 * \param[in] destination_index Destination index where the jump will be applied.
 */
  // Store jump information inside local variables
  const Side destination_side = jump.GetDestinationSide(destination_index);
  const JumpType destination_type = jump.GetDestinationType(destination_index);
  const IndexType destination_capillary = jump.GetDestinationCapillary(destination_index);

  // Applies the jump
  switch (destination_type) {
    case JumpType::creation: {
      AddNewInterface(destination_capillary, destination_side);
      break;
    }
    case JumpType::annihilation: {
      RemoveInterface(destination_capillary, destination_side);
      UnplugCapillary(destination_capillary);
      break;
    }
    case JumpType::bubble_cleaning: {
      // Acts like create, just remove unecessary bubbles
      RemoveBubblesFromDestination(destination_side, destination_capillary);
      AddNewInterface(destination_capillary, destination_side);
      break;
    }
    case JumpType::unplug_annihilation: {
      // Acts like annihilation, just remove the interface(s)
      RemoveBubblesFromDestination(destination_side, destination_capillary);
      break;
    }
    case JumpType::none: {
      std::fprintf(stderr, "%s %s", "Error in CapillaryInterface::ApplyJumpToDestination.",
                                    "JumpType::none is not a valid option!\n");
      std::exit(-1);
    }
  }

  // Flips the fluid when jump occurred in the source side
  if (destination_side == Side::source) {
    FlipFluidAtSource(destination_capillary);
  }
}  // end of CapillaryInterface::ApplyJumpToDestination() method



void CapillaryInterface::RevertPluggedCapillariesJumpAtInletNode(
  arma::Col<IndexType> &changed_caps) {
/**
 * Revert the last plugged capillaries originated from interfaces reaching inlet nodes.
 *
 * \param[in, out]  changed_caps  Unique indexes of the capillaries that has been changed. In this
 *                                method, the unique indexes of the unplugged capillaries are
 *                                inserted.
 */
  const arma::Col<IndexType> last_plug_unique = arma::unique(last_plug_inlet_info_);

  // Checks all the plugs that occurred
  for (const IndexType &capillary : last_plug_unique) {
    changed_caps.insert_rows(changed_caps.n_elem, arma::Col<IndexType>({ capillary }));
    UnplugCapillary(capillary);
  }
  changed_caps = arma::unique(changed_caps);  // Avoid duplicates
}  // end of CapillaryInterface::RevertPluggedCapillariesJumpAtInletNode() method



void CapillaryInterface::RevertJump(arma::Col<IndexType> &changed_caps) {
/**
 * Revert the last jump originated from plugged capillaries.
 * When a capillary is plugged due to an interface jumping back, all the effects of
 * the original jump have to be reverted to the state it was before the jump.
 *
 * \param[in, out]  changed_caps  Unique indexes of the capillaries that has been changed. In this
 *                                method, the unique indexes of the capillaries that originates the
 *                                previous jumps, and capillaries that receive the interfaces are
 *                                inserted.
 */
  const arma::Col<IndexType> last_plug_unique = arma::unique(last_plug_info_);

  // Checks all the plugs that occurred
  for (const IndexType &capillary : last_plug_unique) {
    changed_caps.insert_rows(changed_caps.n_elem, arma::Col<IndexType>({ capillary }));

    // Get the side of the plug in the plugged capillary
    const Side jump_side = GetCapillaryPluggedSide(capillary);
    const IndexType jump_index = last_jump_info_.GetJumpIndex(capillary);
    const IndexType number_of_destinations = last_jump_info_.GetNumberOfDestinations(jump_index);

    // Recover jumped interface
    AddNewInterface(capillary, jump_side);

    // if jump occurred from source side, flips the fluid
    if (jump_side == Side::source) {
      FlipFluidAtSource(capillary);
    }

    // Accumulate output string into stream
    std::stringstream sstream;
    sstream << "Revert in " << capillary << " applied to";

    // Apply the revert to all capillaries connected to target capillary
    for (IndexType i = 0U; i != number_of_destinations; ++i) {
      const auto destination_capillary = last_jump_info_.GetDestinationCapillary(jump_index, i);
      const auto destination_type = last_jump_info_.GetDestinationType(jump_index, i);
      const auto destination_side = last_jump_info_.GetDestinationSide(jump_index, i);

      changed_caps.insert_rows(
        changed_caps.n_elem, arma::Col<IndexType>({ destination_capillary }));

      switch (destination_type) {
        case JumpType::creation: {
          // If jump created an interface, reverts by removing an interface
          RemoveInterface(destination_capillary, destination_side);
          break;
        }
        case JumpType::annihilation: {
          // If jump annihilated an interface, reverts by adding an interface and plugging it back
          AddNewInterface(destination_capillary, destination_side);
          PlugCapillary(destination_capillary, destination_side);
          break;
        }
        case JumpType::bubble_cleaning: {
          // If jump created an interface, reverts by removing an interface
          RemoveInterface(destination_capillary, destination_side);
          break;
        }
        case JumpType::unplug_annihilation: {
          // If jump annihilated an interface, reverts by adding an interface
          AddNewInterface(destination_capillary, destination_side);
          break;
        }
        case JumpType::none: {
          std::fprintf(stderr, "%s %s", "Error in CapillaryInterface::RevertJump.",
                                        "JumpType::none is not a valid option!\n");
          std::exit(-1);
        }
      }

      // if jump occurred to source side, flips the fluid
      if (destination_side == Side::source) {
        FlipFluidAtSource(destination_capillary);
      }

      // Log reversal information about each destination
      sstream << " : " << destination_capillary
              << "[" << static_cast<int>(destination_type)
              << (destination_side == Side::source ? ",Source]" : ",Target]");
    }

    // Log relevant information
    LOG(INFO) << capillary << ":jump_side[" << (jump_side == Side::source ? "Source" : "Target")
        << "], jump_index[" << jump_index << "], number_of_destinations["
        << number_of_destinations << "].";
    LOG(INFO) << sstream.str();
  }
  last_jump_info_.Clear();
  changed_caps = arma::unique(changed_caps);  // Avoid duplicates
}  // end of CapillaryInterface::RevertJump() method



IndexType CapillaryInterface::GetNodeIndexAtCapillaryAtSide(const IndexType capillary,
                                                            const Side side) const {
/**
 * Get the node index of the capillary based on \c side.
 * Return the node with smaller index if \c side == Side::source, return the node with
 * bigger index if \c side == Side::target.
 *
 * \param[in] capillary   Index of the target capillary.
 * \param[in] side        Side of the capillary requested node.
 * \retval    node_index  Index of the node at the \c side tip of capillary \c capillary.
 */
  const IndexType index = static_cast<IndexType>(side) - 1U;
  return context_->linked_nodes_(0U, 2U * capillary + index);
}  // end of CapillaryInterface::GetNodeIndexAtCapillaryAtSide() method



Side CapillaryInterface::GetNodeSideAtCapillary(const IndexType node,
                                                const IndexType capillary) const {
/**
 * Get the side of a \c node in a given \c capillary.
 *
 * \param[in] node        Index of the node.
 * \param[in] capillary   Index of the capillary.
 * \retval    side        Side at which the node is located with respect to the capillary.
 */
  return GetNodeIndexAtCapillaryAtSide(capillary, Side::source) == node ? Side::source
                                                                        : Side::target;
}  // end of CapillaryInterface::GetNodeSideAtCapillary() method



void CapillaryInterface::AddNewInterface(const IndexType capillary, const Side side) {
/**
 * Adds new interface to \c all_interfaces_(link) at \f$ x = 0 (1) \f$ depending on \c side.
 *
 * \param[in] capillary   Index of the target capillary.
 * \param[in] side        The side of the capillary where the interface will be added.
 */
  if (GetNumberOfInterfacesAtCapillary(capillary) == 0U) {
    const arma::Col<IndexType> temp = { capillary };
    const arma::uvec position = arma::find(capillaries_with_interfaces_ > capillary, 1, "first");
    if (position.n_elem == 0) {
      capillaries_with_interfaces_.insert_rows(capillaries_with_interfaces_.n_elem, temp);
    } else {
      capillaries_with_interfaces_.insert_rows(position(0), temp);
    }
  }

  switch (side) {
    case Side::source: {
      all_interfaces_(capillary).insert_rows(0, arma::zeros<arma::rowvec>(1) + arma::datum::eps);
      break;
    }
    case Side::target: {
      all_interfaces_(capillary).insert_rows(all_interfaces_(capillary).n_elem,
                                             arma::ones<arma::rowvec>(1) - arma::datum::eps);
      break;
    }
    case Side::none: {
      std::fprintf(stderr, "%s %s", "Error in CapillaryInterface::AddNewInterface.",
                                    "Side::none is not a valid option!\n");
      std::exit(-1);
    }
  }
}  // end of CapillaryInterface::AddNewInterface() method



void CapillaryInterface::RemoveInterface(const IndexType capillary, const Side side) {
/**
 * Remove interface from \c all_interfaces_(link) closest to \f$ x = 0 (1) \f$ depending on \c side.
 *
 * \param[in] capillary   Index of the target capillary.
 * \param[in] side   The side of the capillary where the closest interface will be removed.
 */
  // Shed row from all_interfaces_ (depending on side)
  arma::uword row_to_shed = (side == Side::source) ? 0U : all_interfaces_(capillary).n_elem - 1U;
  LOG(INFO) << "Removing interface " << all_interfaces_(capillary)(row_to_shed)
    << " from capillary " << capillary;
  all_interfaces_(capillary).shed_row(row_to_shed);

  // Shed row from capillaries_with_interfaces_ (depending on side)
  if (GetNumberOfInterfacesAtCapillary(capillary) == 0U) {
    row_to_shed = arma::as_scalar(arma::find(capillaries_with_interfaces_ == capillary));
    capillaries_with_interfaces_.shed_row(row_to_shed);
  }
}  // end of CapillaryInterface::RemoveInterface() method



arma::uword CapillaryInterface::GetNumberOfInterfaces(void) const {
/**
 * Returns the number of interfaces of all capillaries in the dynamic capillary network.
 *
 * \retval  total_number_interfaces   Total number of interfaces in the network.
 */
  // TODO(rneumann): Avoid counting altogether by tracking number at each creation/destruction
  // Loop over field and accumulate length of each component
  arma::uword total_number_interfaces = 0U;
  for (const auto &interfaces : all_interfaces_) total_number_interfaces += interfaces.n_elem;
  return total_number_interfaces;
}  // end of CapillaryInterface::GetNumberOfInterfaces() method



arma::uword CapillaryInterface::GetNumberOfInterfacesAtCapillary(const IndexType capillary) const {
/**
 * Returns the number of interfaces of a given capillary in the dynamic capillary network.
 *
 * \param[in] capillary                       Index of the target capillary.
 * \retval    number_interfaces_at_capillary  Number of interfaces at target capillary.
 */
  return all_interfaces_(capillary).n_elem;
}  // end of CapillaryInterface::GetNumberOfInterfacesAtCapillary() method



bool CapillaryInterface::IsNumberOfInterfacesAtCapillaryOdd(const IndexType capillary) const {
/**
 * Return whether or not the number of interfaces at a given capillary is odd.
 *
 * \param[in] capillary Index of the target capillary.
 * \retval    is_odd    Boolean representing whether the target capillary has an odd interface number.
 */
  return GetNumberOfInterfacesAtCapillary(capillary) % 2U == 1U;
}  // end of CapillaryInterface::IsNumberOfInterfacesAtCapillaryOdd() method



bool CapillaryInterface::IsCapillaryPlugged(const IndexType capillary) const {
/**
 * Return whether or not the given capillary is plugged.
 *
 * \param[in] capillary   Index of the target capillary.
 * \retval    is_plugged  Boolean representing whether or not the target capillary is plugged.
 */
  return plugs_(capillary) != 0U;
}  // end of CapillaryInterface::IsCapillaryPlugged() method



Side CapillaryInterface::GetCapillaryPluggedSide(const IndexType capillary) const {
/**
 * Return the plugged end of the capillary. Returns Side::none if it is not plugged.
 *
 * \param[in] capillary   Index of the target capillary.
 * \retval    plug_end    Side of the plugged end of the capillary.
 */
  return static_cast<Side>(plugs_(capillary));
}  // end of CapillaryInterface::IsCapillaryPlugged() method



Side CapillaryInterface::GetSideOfNodeAtCapillary(const IndexType capillary,
                                                  const IndexType node) const {
/**
 * Returns the side of the node at the capillary. Return Side::none if the node is not at the
 * capillary.
 *
 * \param[in] capillary   Index of the target capillary.
 * \param[in] node        Index of the node.
 * \retval    node_side   Side of the node at the capillary.
 */
  Side node_side;
  if (node == GetNodeIndexAtCapillaryAtSide(capillary, Side::source)) {
    node_side = Side::source;
  } else if (node == GetNodeIndexAtCapillaryAtSide(capillary, Side::target)) {
    node_side = Side::target;
  } else {
    node_side = Side::none;
  }
  return node_side;
}  // end of CapillaryInterface::GetSideOfNodeAtCapillary() method



JumpType CapillaryInterface::GetJumpType(const IndexType capillary,
                                         const Side arriving_side,
                                         const Side plugged_side) const {
/**
 * Returns the jump type, depending on whether a new interface should be created, an existing
 * interface should be annihilated, tiny bubbles should be eliminated or unplugged interface
 * with possible (or not) tiny bubbles should be annihilated.
 *
 * \param[in] capillary       Index of the target capillary.
 * \param[in] arriving_side   Side at which a new interface is arriving.
 * \param[in] plugged_side    Side at which the capillary is plugged (if any).
 * \retval    jump_type       Jump type based on the incoming and existing sides.
 */
  if (arriving_side == plugged_side) {
    return JumpType::annihilation;
  } else {
    IndexType capillaries_to_be_removed = 0U;
    switch (arriving_side) {
      case Side::source: {
        for (auto interface : all_interfaces_(capillary)) {
          double value = 0.0 + arma::datum::eps;
          if (interface <= value) {
            capillaries_to_be_removed++;
          }
        }
        break;
      }
      case Side::target: {
        for (auto interface : all_interfaces_(capillary)) {
          double value = 1.0 - arma::datum::eps;
          if (interface >= value) {
            capillaries_to_be_removed++;
          }
        }
        break;
      }
      case Side::none: {
        LOG(FATAL) << "Error in CapillaryInterface::GetJumpType. "
                   << "Side::none is not a valid option to add interfaces";
      }
    }

    if (capillaries_to_be_removed == 0U) {
      return JumpType::creation;
    } else {
      if (capillaries_to_be_removed % 2U == 0U) {
        return JumpType::bubble_cleaning;
      } else {
        return JumpType::unplug_annihilation;
      }
    }
  }
}  // end of CapillaryInterface::GetJumpType() method



double CapillaryInterface::CalculateEffectiveInterfacePositionAtCapillary(
  const IndexType capillary) const {
/**
 * Returns the effective interface position at a given capillary.
 *
 * \param[in] capillary     Index of the target capillary.
 * \retval    eff_position  Effective interface position at target capillary.
 */
  // Prepare capillary interface vector
  arma::vec capillary_interfaces = {0.0, 1.0};
  capillary_interfaces.insert_rows(1U, all_interfaces_(capillary));

  // diffs
  const arma::vec diffs = arma::diff(capillary_interfaces);

  // sum
  const arma::uvec even_indexes = arma::regspace<arma::uvec>(0U, 2U, diffs.n_elem - 1U);
  return arma::sum(diffs.elem(even_indexes));
}  // end of CapillaryInterface::CalculateEffectiveInterfacePositionAtCapillary() method



void CapillaryInterface::CalculateEffectiveInterfacePosition(arma::Col<double> &x_vec) const {
/**
 * Calculates the effective interface position at all capillaries and saves it in a vector.
 *
 * \param[out] x_vec    Effective interface positions in all capillaries.
 */
  // TODO(ralves): parallelisation is possible
  for (IndexType i = 0U; i < context_->number_of_links_; ++i) {
    x_vec(i) = CalculateEffectiveInterfacePositionAtCapillary(i);
  }
}  // end of CapillaryInterface::CalculateEffectiveInterfacePosition() method



double CapillaryInterface::GetInterfacePositionAtCapillaryNearSide(const IndexType capillary,
                                                                   const Side side) const {
/**
 * Returns the position of the interface closest to a given \c side at a given \c capillary, or -1
 * when there is no interface.
 *
 * \param[in] capillary               Index of the target capillary.
 * \param[in] side                    Side near which we are querying the interface position.
 * \retval    interface_position      Position of the interface closest to \c side.
 */
  const arma::vec capillary_interfaces(all_interfaces_(capillary));
  const arma::uword index = (side == Side::target) ? capillary_interfaces.n_elem - 1U : 0U;
  return capillary_interfaces.is_empty() ? -1.0 : capillary_interfaces(index);
}  // end of CapillaryInterface::GetInterfacePositionAtCapillaryNearSide() method



const arma::Col<IndexType> &CapillaryInterface::GetIndexOfCapillariesWithInterfaces(void) const {
/**
 * Return a vector with the index of each capillary with at least one interface in
 * the dynamic capillary network.
 *
 * \retval  capillaries_with_interfaces_  Index of capillaries with interfaces.
 */
  return capillaries_with_interfaces_;
}  // end of CapillaryInterface::GetIndexOfCapillariesWithInterfaces() method



arma::uword CapillaryInterface::GetNumberOfCapillariesWithInterfaces(void) const {
/**
 * Return the number of capillaries that have at least one interface in the dynamic capillary network.
 *
 * \retval  number_capillaries_with_interfaces_ Number of capillaries with interfaces.
 */
  return capillaries_with_interfaces_.n_elem;
}  // end of CapillaryInterface::GetNumberOfCapillariesWithInterfaces() method



void CapillaryInterface::FlipFluidAtSource(const IndexType capillary) {
/**
 * Swap the fluid number at source at the target capillary
 *
 * \param[in] capillary   Index of the target capillary.
 */
  LOG(INFO) << "Flip fluid at source of capillary " << capillary;
  fluid_at_source_(capillary) = !fluid_at_source_(capillary);
}  // end of CapillaryInterface::FlipFluidAtSource() method



IndexType CapillaryInterface::GetFluidAtCapillaryAtSide(const IndexType capillary,
                                                        const Side side) const {
/**
 * Get the fluid number at \c side at the target \c capillary.
 *
 * \param[in] capillary         Index of the target capillary.
 * \param[in] side              Side at which we are querying the interface position.
 * \retval    fluid_at_side     Fluid type wetting the \c side node of target \c capillary.
 */
  return (side == Side::source) ? fluid_at_source_(capillary)
                                : IsNumberOfInterfacesAtCapillaryOdd(capillary) ?
                                  !fluid_at_source_(capillary) : fluid_at_source_(capillary);
}  // end of CapillaryInterface::GetFluidAtCapillaryAtSide() method



void CapillaryInterface::PlugCapillary(const IndexType capillary, const Side side) {
/**
 * Plugs the capillary at the specified end.
 *
 * \param[in] capillary    Index of the target capillary.
 * \param[in] side         Specify which side of the capillary will be plugged
 */
  plugs_(capillary) = static_cast<IndexType>(side);
}  // end of CapillaryInterface::PlugCapillary() method



void CapillaryInterface::UnplugCapillary(const IndexType capillary) {
/**
 * Unplugs the capillary.
 *
 * \param[in] capillary    Index of the target capillary.
 */
  plugs_(capillary) = static_cast<IndexType>(Side::none);
}  // end of CapillaryInterface::PlugCapillary() method



void CapillaryInterface::ComputeDeltas(const arma::Col<double> &x_vec) {
/**
 * Computes the delta of interface movement in all capillaries and stores int vector
 * \c deltas_
 *
 * \param[in] x_vec  Double array with the new effective interface position at each capillary.
 */
  for (IndexType i = 0U; i < context_->number_of_links_; ++i) {
    if (IsNumberOfInterfacesAtCapillaryOdd(i)) {
      deltas_(i) = x_vec(i) - CalculateEffectiveInterfacePositionAtCapillary(i);
    } else {
      deltas_(i) = x_vec(i) - GetInterfacePositionAtCapillaryNearSide(i, Side::source);
    }
  }
}  // end of CapillaryInterface::ComputeDeltas() method



void CapillaryInterface::ComputeInterfacePositionsAndOffsets(arma::Col<double> &positions,
                                                             arma::Col<arma::uword> &offsets) {
/**
 * This method takes newly-created (and uninitialised) \c positions and \c offsets Armadillo vectors
 * and populate them together. The \c positions vector contains the coordinate of all the interfaces
 * in the system, while the \c offsets vector indicates the index at which the first interface of
 * each capillary can be found.
 *
 * \par Example:
 * If \c offsets(4) = 10 and \c offsets(5) = 13, that means that one can find the coordinates of the
 * 3 interfaces of capillary number 4 between \c positions(10) and \c positions(12).
 *
 * \param[in,out] positions Positions of all the interfaces in the system.
 * \param[in,out] offsets   Indexes at which interfaces of a capillary are found in \c positions.
 */
  arma::uword counter = 0U;
  const arma::uword number_of_capillaries = all_interfaces_.n_elem;

  // Iterate over all capillaries and accumulate their interfaces
  for (IndexType capillary_idx = 0U; capillary_idx != number_of_capillaries; ++capillary_idx) {
    for (const double interface : all_interfaces_(capillary_idx)) {
      positions(counter++) = interface;
    }
    offsets(capillary_idx + 1) = counter;
  }
}  // end of CapillaryInterface::ComputeInterfacePositionsAndOffsets() method



arma::Mat<double> CapillaryInterface::ComputeCapillarySaturation(void) const {
/**
 * This method calculates the phase saturation in each capillary from the position of the
 * interfaces along the capillaries and the fluid that wets the capillary's \c source node.
 *
 * \f[
 *    S_{\alpha}^{j} =
 *    \begin{cases}
 *      x_j,      & \text{if fluid } \alpha \text{ wets the source node of capillary } j, \\
 *      1 - x_j,  & \text{otherwise.}
 *    \end{cases}
 * \f]
 *
 * \retval    saturation  Matrix with 2 columns containing the phase saturations in each capillary.
 */
  // Create local matrix
  arma::mat phase_saturation(context_->number_of_links_, 2, arma::fill::zeros);

  // Iterate over capillaries
  // TODO(ralves): parallelisation is possible
  for (IndexType i = 0U; i < context_->number_of_links_; ++i) {
    const IndexType fluid_at_source = GetFluidAtCapillaryAtSide(i, Side::source);
    const double effective_interface_position = CalculateEffectiveInterfacePositionAtCapillary(i);
    phase_saturation(i, fluid_at_source) = effective_interface_position;
    phase_saturation(i, !fluid_at_source) = 1.0 - effective_interface_position;
  }

  return phase_saturation;
}  // end of CapillaryInterface::ComputeCapillarySaturation() method

}  // namespace simulator
