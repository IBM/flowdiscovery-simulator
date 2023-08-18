/*
 * \file test/main_stress_test.cc
 * \brief Contains the stress testing.
 *
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 */

#include <gtest/gtest.h>
#include <glog/logging.h>
#include "test/src/flow_simulator/dynamic_simulation_stress_tests.h"

int main(int argc, char** argv) {
  // Initialise Google's logging library
  ::google::InitGoogleLogging(argv[0]);

  // Initialise Google's testing library
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
