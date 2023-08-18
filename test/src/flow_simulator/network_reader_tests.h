/**
 * \file test/src/flow_simulator/network_reader_tests.h
 * \brief Contains regression tests of \c NetworkReader class methods.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 */

#ifndef TEST_SRC_FLOW_SIMULATOR_NETWORK_READER_TESTS_H_
#define TEST_SRC_FLOW_SIMULATOR_NETWORK_READER_TESTS_H_

#include <gtest/gtest.h>
#include <string>
#include "src/flow_simulator/algorithms/network_reader.h"

using NetworkReader = simulator::NetworkReader;
using NetworkInformation = simulator::NetworkInformation;

TEST(NetworkReader, GetNetwork_Test) {
  std::string folder = "test/src/flow_simulator/input";
  NetworkReader reader;
  NetworkInformation net_info = reader.GetNetwork(folder + "/JGF_template.json");

  // Define complex targets
  arma::mat ctrl_voxels = {{0.0, 0.0, 0.0, 4.0}, {1.0, 1.0, 1.0, 9.0}};
  arma::umat linked_nodes = {{0, 1}, {0, 0}};
  arma::vec link_direction = {{+1.0, -1.0}};

  // Check equality
  double tolerance = 1.0e-12;
  EXPECT_EQ(net_info.number_of_nodes, 2);
  EXPECT_EQ(net_info.number_of_links, 1);
  EXPECT_DOUBLE_EQ(arma::as_scalar(net_info.link_length), 1.73205080757);
  EXPECT_DOUBLE_EQ(arma::as_scalar(net_info.link_squared_radius), 5.169298742047715);
  EXPECT_TRUE(arma::approx_equal(net_info.ctrl_voxels, ctrl_voxels, "absdiff", tolerance));
  EXPECT_TRUE(arma::approx_equal(net_info.linked_nodes, linked_nodes, "absdiff", tolerance));
  EXPECT_TRUE(arma::approx_equal(net_info.link_direction, link_direction, "absdiff", tolerance));
}

#endif  // TEST_SRC_FLOW_SIMULATOR_NETWORK_READER_TESTS_H_
