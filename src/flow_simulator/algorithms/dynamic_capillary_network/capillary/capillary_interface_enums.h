/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/capillary/capillary_interface.h
 * \brief Contains the \c CapillaryInterface class.
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 *
 * This header file contains the \c CapillaryInterface class.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_CAPILLARY_CAPILLARY_INTERFACE_ENUMS_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_CAPILLARY_CAPILLARY_INTERFACE_ENUMS_H_

#include <cstdint>

namespace simulator {

enum class Side : uint8_t {
  none = 0,
  source = 1,
  target = 2
};  // end of class Side

enum class JumpType : uint8_t {
  none = 0,
  annihilation = 1,
  creation = 2,
  bubble_cleaning = 3,
  unplug_annihilation = 4
};  // end of class JumpType

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_CAPILLARY_CAPILLARY_INTERFACE_ENUMS_H_
