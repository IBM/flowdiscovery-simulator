/*
 * \file test/main_test.cc
 * \brief Contains the regression tests.
 *
 * \authors Alexandre Ashade Lassance Cunha \<aashade@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Hugo de Oliveira Barbalho \<hugob@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 */

#include <gtest/gtest.h>
#include <glog/logging.h>
#include "test/src/config_reader/config_reader_tests.h"
#include "test/src/flow_simulator/poiseuille_tests.h"
#include "test/src/flow_simulator/hydrostatic_tests.h"
#include "test/src/flow_simulator/network_reader_tests.h"
#include "test/src/flow_simulator/boundary_condition_tests.h"
#include "test/src/flow_simulator/dynamic_simulation_tests.h"
#include "test/src/flow_simulator/dynamic_capillary_interface_tests.h"
#include "test/src/flow_simulator/jump_event_info_tests.h"
#include "test/src/flow_simulator/viscosity_behaviour_tests.h"

int main(int argc, char** argv) {
  // Initialise Google's logging library
  ::google::InitGoogleLogging(argv[0]);

  // Initialise Google's testing library
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
