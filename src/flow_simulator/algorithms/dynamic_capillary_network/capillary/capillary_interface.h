/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/capillary/capillary_interface.h
 * \brief Contains the \c CapillaryInterface class.
 *
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 *
 * This header file contains the \c CapillaryInterface class.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_CAPILLARY_CAPILLARY_INTERFACE_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_CAPILLARY_CAPILLARY_INTERFACE_H_

#include <armadillo>
#include <memory>
#include <cstdint>
#include <vector>
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network_context.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/capillary_interface_enums.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/interface_jump_event_info.h"

namespace simulator {

/**
 * \class CapillaryInterface capillary_interface.h "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/capillary_interface.h"
 * \brief Class that manages interface information in the dynamic capillary network
 *
 * The Capillary Interface class control the interface informations of
 * all capillaries in the dynamic capillary network.
 */

class CapillaryInterface {
 public:
  /// Parametrised constructor
  explicit CapillaryInterface(std::shared_ptr<DynamicCapillaryNetworkContext> context)
    : context_(context) { }

  /// Initialises the \c CapillaryInterface instance
  void InitialiseCapillaryInterface(void);

  /// Injects a fluid in all inlet nodes, creating a new interface
  void InjectFluidOnInletNodes(const arma::uword number_interfaces = 1U,
                               const double distance_interfaces = 0);

  /// Updates the interface information of all capillaries
  void UpdateInterfaceInformation(void);

  /// Undo the movement of interfaces performed by \c UpdateInterfaceInformation()
  void RevertUpdateInterfaceInformation(void);

  /// Moves the interfaces of a given capillary by a value \c delta
  void MoveInterfacePositionsAtCapillaryByDelta(const IndexType capillary, const double delta);

  /// Plug capillaries witch interfaces reached an inlet node
  arma::Col<IndexType> PlugCapillariesJumpAtInletNode(
    const arma::Col<IndexType> &capillaries,
    const arma::Col<IndexType> &caps_roots_directions);

  /// Check if there is a plugging
  bool CheckPlugging(const arma::Col<IndexType> &capillaries,
                     const arma::Col<IndexType> &caps_roots_directions);

  /// Executes all interface jumps given by \c roots_found
  void Jump(const arma::Col<IndexType> &capillaries,
            const arma::Col<IndexType> &caps_roots_directions,
            arma::Col<IndexType> &changed_caps);

  /// Remove interfaces that are too close from the destination side at a capillary
  void RemoveBubblesFromDestination(const Side side, const IndexType capillary);

  /// Applies the jump to a given destination
  void ApplyJumpToDestination(const InterfaceJump &jump, const IndexType destination_index);

  /// Revert the last plugged capillaries originated from interfaces reaching inlet nodes
  void RevertPluggedCapillariesJumpAtInletNode(arma::Col<IndexType> &changed_caps);

  /// Revert the plugged capillaries based on the last jump info
  void RevertJump(arma::Col<IndexType> &changed_caps);

  /// Adds a new interface to a given capillary
  void AddNewInterface(const IndexType capillary, const Side side);

  /// Removes an interface from a given capillary
  void RemoveInterface(const IndexType capillary, const Side side);

  /// Returns the number of interfaces in all capillaries
  arma::uword GetNumberOfInterfaces(void) const;

  /// Returns the number of interfaces at a given capillary
  arma::uword GetNumberOfInterfacesAtCapillary(const IndexType capillary) const;

  /// Returns whether or not the number of interfaces at a given capillary is odd
  bool IsNumberOfInterfacesAtCapillaryOdd(const IndexType capillary) const;

  /// Returns whether or not the given capillary is plugged
  bool IsCapillaryPlugged(const IndexType capillary) const;

  /// Returns the side where the capillary is plugged
  Side GetCapillaryPluggedSide(const IndexType capillary) const;

  /// Returns the side of the node at the capillary
  Side GetSideOfNodeAtCapillary(const IndexType capillary, const IndexType node) const;

  // Returns the type of jump
  JumpType GetJumpType(const IndexType capillary,
                       const Side arriving_side,
                       const Side plugged_side) const;

  /// Returns the effective interface position at a given capillary
  double CalculateEffectiveInterfacePositionAtCapillary(const IndexType capillary) const;

  /// Returns the effective interface position at all capillaries
  void CalculateEffectiveInterfacePosition(arma::Col<double> &x_vec) const;

  /// Returns the position of the interface closest to a given \c side of a given \c capillary
  double GetInterfacePositionAtCapillaryNearSide(const IndexType capillary, const Side side) const;

  /// Getter for \c capillaries_with_interfaces_
  const arma::Col<IndexType> &GetIndexOfCapillariesWithInterfaces(void) const;

  /// Returns the number of elements in \c capillaries_with_interfaces_
  arma::uword GetNumberOfCapillariesWithInterfaces(void) const;

  /// Swaps the fluid number at "source" at a given capillary
  void FlipFluidAtSource(const IndexType capillary);

  /// Returns a vector of integer identifying the fluid at the "source" node of each capillary
  const arma::Col<IndexType> &GetFluidAtSource(void) const { return fluid_at_source_; }

  /// Returns the fluid number at \c side at a given \c capillary
  IndexType GetFluidAtCapillaryAtSide(const IndexType capillary, const Side side) const;

  /// Plugs Capillary at the specified end
  void PlugCapillary(const IndexType capillary, const Side side);

  // Unplugs the specified capillary
  void UnplugCapillary(const IndexType capillary);

  /// Compute interface \c deltas_ for all capillaries
  void ComputeDeltas(const arma::Col<double> &x_vec);

  /// Compute interface positions and offsets for all capillaries
  void ComputeInterfacePositionsAndOffsets(arma::Col<double> &positions,
                                           arma::Col<arma::uword> &offsets);

  /// Compute the saturations of fluids 0 and 1 in each capillary
  arma::Mat<double> ComputeCapillarySaturation(void) const;

  /// Returns the node index at a particular \c side of a given \c capillary
  IndexType GetNodeIndexAtCapillaryAtSide(const IndexType capillary, const Side side) const;

  /// Returns the \c Side at which a given \c node is located in a \c capillary
  Side GetNodeSideAtCapillary(const IndexType node, const IndexType capillary) const;

  /// Returns arma::field with interface positions at all capillaries
  const arma::field<arma::Col<double>> &GetInterfacePositions(void) const {
    return all_interfaces_;
  }

  /// Returns arma::field with interface positions at a given \c capillary
  const arma::Col<double> &GetInterfacePositionsAtCapillary(const IndexType capillary) const {
    return all_interfaces_(capillary);
  }

  /// Getter for \c last_plug_info_
  const arma::Col<IndexType> &GetLastPlugInfo() { return last_plug_info_; }

  /// Getter for \c last_plug_inlet_info_
  const arma::Col<IndexType> &GetLastPlugInletInfo() { return last_plug_inlet_info_; }

  /// Getter for \c deltas_
  const arma::Col<double> &GetDeltas() { return deltas_; }

  /// Getter for \c plugs_
  const arma::Col<IndexType> &GetPlugs() { return plugs_; }

  /// Getter for \c last_jump_info_
  const InterfaceJumpEventInfo &GetLastJumpInfo() { return last_jump_info_; }

  /// Setter for \c all_interfaces_
  void SetInterfacePositions(const arma::field<arma::Col<double>> &all_interfaces) {
    all_interfaces_ = all_interfaces;
  }

  /// Setter for \c fluid_at_source_
  void SetFluidAtSource(const arma::Col<IndexType> &fluid_at_source) {
    fluid_at_source_ = fluid_at_source;
  }

  /// Setter for \c deltas_
  void SetDeltas(const arma::Col<double> &deltas) {
    deltas_ = deltas;
  }

  /// Setter for \c plugs_
  void SetPlugs(const arma::Col<IndexType> &plugs) {
    plugs_ = plugs;
  }

  /// Setter for \c last_plug_info_
  void SetLastPlugInfo(const arma::Col<IndexType> &last_plug_info) {
    last_plug_info_ = last_plug_info;
  }

  /// Setter for \c last_plug_inlet_info_
  void SetLastPlugInletInfo(const arma::Col<IndexType> &last_plug_inlet_info) {
    last_plug_inlet_info_ = last_plug_inlet_info;
  }

  /// Setter for \c last_jump_info_
  void SetLastJumpInfo(const InterfaceJumpEventInfo &last_jump_info) {
    last_jump_info_ = last_jump_info;
  }

  /// Setter for \c capillaries_with_interfaces_
  void SetCapillariesWithInterfaces(const arma::Col<IndexType> &capillaries_with_interfaces) {
    capillaries_with_interfaces_ = capillaries_with_interfaces;
  }

 private:
  /// Sets the sizes of vectors in this class
  void SetSizes(const IndexType &n_links);

  /**
   * \brief Mapping of the interface positions in all capillaries.
   *
   * This member variable stores the interface position of all capillaries in the
   *  dynamic capillary network.
   */
  arma::field<arma::Col<double>> all_interfaces_;

  /**
   * \brief Vector of fluid number at source of each capillary.
   *
   * Contains the fluid number at source of each capillary in the dynamic capillary
   * network.
   */
  arma::Col<IndexType> fluid_at_source_;

  /**
   * \brief Vector mapping wether a capillary is plugged and the plugged node.
   *
   * Contains 0 for unplugged capillary, 1 for plugged at source and 2 for plugged at target.
   */
  arma::Col<IndexType> plugs_;

  /**
   * \brief \c InterfaceJumpInfo instance.
   *
   * Instance of InterfaceJumpInfo that store where a jump took place in
   * the last jump event, and contains informations about that jump.
   */
  InterfaceJumpEventInfo last_jump_info_;

  /**
   * \brief Vector mapping what capillaries were plugged in the last update.
   *
   * Contains the index of capillaries plugged in the last update.
   */
  arma::Col<IndexType> last_plug_info_;

  /**
   * \brief Vector mapping what capillaries were plugged in the last update because of jumps at
   * inlet nodes.
   *
   * Contains the index of capillaries plugged in the last update because of inlet nodes.
   */
  arma::Col<IndexType> last_plug_inlet_info_;

  /**
   * \brief Shared pointer to \c StaticCapillaryNetworkContext instance
   *
   * Contains variables needed to some of the methods of this class
   */
  std::shared_ptr<DynamicCapillaryNetworkContext> context_;

  /**
   * \brief Vector that stores the index of the capillaries that have at least one interface
   *
   * Contains the index of all capillaries of the dynamic interface network that contains
   * at least one interface.
   */
  arma::Col<IndexType> capillaries_with_interfaces_;

  /**
   * \brief Vector that stores the interface displacement in each capillaries
   *
   * Contains the interface displacement in each capillaries since the last call of \c
   * ComputeDeltas()
   */
  arma::Col<double> deltas_;
};  // end of class CapillaryInterface

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_CAPILLARY_CAPILLARY_INTERFACE_H_
