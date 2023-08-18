/**
 * \file test/src/flow_simulator/jump_event_info_tests.h
 * \brief Contains regression tests of \c InterfaceJumpEventInfo class methods
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright IBM Corp. All Rights Reserved
 * \date 2018
 */

#ifndef TEST_SRC_FLOW_SIMULATOR_JUMP_EVENT_INFO_TESTS_H_
#define TEST_SRC_FLOW_SIMULATOR_JUMP_EVENT_INFO_TESTS_H_

#include <gtest/gtest.h>
#include "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/interface_jump_event_info.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/capillary_interface_enums.h"

using Side = simulator::Side;
using JumpType = simulator::JumpType;
using IndexType = simulator::IndexType;
using InterfaceJumpEventInfo = simulator::InterfaceJumpEventInfo;

TEST(JumpEventInfo, InsertJumpWithNoCrossed_Return_Correct_Insertions) {
  // arrange
  InterfaceJumpEventInfo jump_event;

  // act
  IndexType jump1 = jump_event.InsertJump(0, Side::target, 0);
  IndexType jump2 = jump_event.InsertJump(1, Side::target, 1);
  IndexType jump3 = jump_event.InsertJump(2, Side::source, 2);

  // assert
  IndexType source1 = jump_event.GetJumpOrigin(jump1);
  IndexType source2 = jump_event.GetJumpOrigin(jump2);
  IndexType source3 = jump_event.GetJumpOrigin(jump3);
  IndexType number_of_jumps = jump_event.GetNumberOfJumps();

  EXPECT_EQ(source1, 0);
  EXPECT_EQ(source2, 1);
  EXPECT_EQ(source3, 2);
  EXPECT_EQ(number_of_jumps, 3);
}

TEST(JumpEventInfo, GetJumpIndex_ReturnCorrectIndexOfJumps) {
  // arrange
  InterfaceJumpEventInfo jump_event;
  IndexType jump1 = jump_event.InsertJump(4, Side::target, 0);
  IndexType jump2 = jump_event.InsertJump(2, Side::target, 1);
  IndexType jump3 = jump_event.InsertJump(6, Side::source, 2);

  // act
  IndexType jump_index_of_source_6 = jump_event.GetJumpIndex(6);
  IndexType jump_index_of_source_4 = jump_event.GetJumpIndex(4);
  IndexType jump_index_of_source_2 = jump_event.GetJumpIndex(2);

  // assert
  EXPECT_EQ(jump_index_of_source_6, jump3);
  EXPECT_EQ(jump_index_of_source_4, jump1);
  EXPECT_EQ(jump_index_of_source_2, jump2);
}

TEST(JumpEventInfo, IsJumpSource_ReturnTrueWhenSource) {
  // arrange
  InterfaceJumpEventInfo jump_event;
  jump_event.InsertJump(4, Side::target, 0);
  jump_event.InsertJump(2, Side::target, 1);
  jump_event.InsertJump(6, Side::source, 2);

  // act
  bool is_4_origin = jump_event.IsJumpOrigin(4);
  bool is_7_origin = jump_event.IsJumpOrigin(7);
  bool is_2_origin = jump_event.IsJumpOrigin(2);

  // assert

  EXPECT_TRUE(is_4_origin);
  EXPECT_FALSE(is_7_origin);
  EXPECT_TRUE(is_2_origin);
}

TEST(JumpEventInfo, InsertJumpWithWithCrossed_Return_Correct_Crossed_Insertions) {
  // arrange
  InterfaceJumpEventInfo jump_event;

  // act
  IndexType jump1 = jump_event.InsertJump(0, Side::target, 0);
  IndexType jump2 = jump_event.InsertJump(1, Side::target, 1);
  IndexType jump3 = jump_event.InsertJump(2, Side::source, 2);
  IndexType jump4 = jump_event.InsertJump(6, Side::source, 1);
  IndexType jump5 = jump_event.InsertJump(9, Side::source, 0);

  // assert
  bool is_crossed_1 = jump_event.IsCrossJump(jump1);
  bool is_crossed_2 = jump_event.IsCrossJump(jump2);
  bool is_crossed_3 = jump_event.IsCrossJump(jump3);
  bool is_crossed_4 = jump_event.IsCrossJump(jump4);
  bool is_crossed_5 = jump_event.IsCrossJump(jump5);

  EXPECT_TRUE(is_crossed_1);
  EXPECT_TRUE(is_crossed_2);
  EXPECT_FALSE(is_crossed_3);
  EXPECT_TRUE(is_crossed_4);
  EXPECT_TRUE(is_crossed_5);
}

TEST(JumpEventInfo, AddDestination_ReturnCorrectDestinationCapillaries) {
  // arrange
  InterfaceJumpEventInfo jump_event;
  IndexType jump1 = jump_event.InsertJump(0, Side::target, 0);
  IndexType jump2 = jump_event.InsertJump(1, Side::target, 1);

  // act
  IndexType destination_1_1 = jump_event.AddDestination(jump1,
                                                        2,
                                                        Side::target,
                                                        JumpType::creation);
  IndexType destination_1_2 = jump_event.AddDestination(jump1,
                                                        5,
                                                        Side::target,
                                                        JumpType::annihilation);
  IndexType destination_2_1 = jump_event.AddDestination(jump2,
                                                        3,
                                                        Side::target,
                                                        JumpType::creation);
  IndexType destination_2_2 = jump_event.AddDestination(jump2,
                                                        6,
                                                        Side::source,
                                                        JumpType::creation);

  // assert
  IndexType number_of_destinations_Jump_1 = jump_event.GetNumberOfDestinations(jump1);
  IndexType number_of_destinations_Jump_2 = jump_event.GetNumberOfDestinations(jump2);
  IndexType destination_capillary_1_1 = jump_event.GetDestinationCapillary(jump1, destination_1_1);
  IndexType destination_capillary_1_2 = jump_event.GetDestinationCapillary(jump1, destination_1_2);
  IndexType destination_capillary_2_1 = jump_event.GetDestinationCapillary(jump2, destination_2_1);
  IndexType destination_capillary_2_2 = jump_event.GetDestinationCapillary(jump2, destination_2_2);

  EXPECT_EQ(number_of_destinations_Jump_1, 2);
  EXPECT_EQ(number_of_destinations_Jump_2, 2);
  EXPECT_EQ(destination_capillary_1_1, 2);
  EXPECT_EQ(destination_capillary_1_2, 5);
  EXPECT_EQ(destination_capillary_2_1, 3);
  EXPECT_EQ(destination_capillary_2_2, 6);
}

TEST(JumpEventInfo, AddDestination_ReturnCorrectDestinationSides) {
  // arrange
  InterfaceJumpEventInfo jump_event;
  IndexType jump1 = jump_event.InsertJump(0, Side::target, 0);

  // act
  IndexType destination_1_1 = jump_event.AddDestination(jump1,
                                                        2,
                                                        Side::target,
                                                        JumpType::creation);
  IndexType destination_1_2 = jump_event.AddDestination(jump1,
                                                        5,
                                                        Side::source,
                                                        JumpType::annihilation);

  // assert
  Side destination_side_1_1 = jump_event.GetDestinationSide(jump1, destination_1_1);
  Side destination_side_1_2 = jump_event.GetDestinationSide(jump1, destination_1_2);

  EXPECT_TRUE(destination_side_1_1 == Side::target);
  EXPECT_TRUE(destination_side_1_2 == Side::source);
}

TEST(JumpEventInfo, AddDestination_ReturnCorrectDestinationType) {
  // arrange
  InterfaceJumpEventInfo jump_event;
  IndexType jump1 = jump_event.InsertJump(0, Side::target, 0);

  // act
  IndexType destination_1_1 = jump_event.AddDestination(jump1,
                                                        2,
                                                        Side::target,
                                                        JumpType::creation);
  IndexType destination_1_2 = jump_event.AddDestination(jump1,
                                                        5,
                                                        Side::source,
                                                        JumpType::annihilation);
  IndexType destination_1_3 = jump_event.AddDestination(jump1,
                                                        8,
                                                        Side::target,
                                                        JumpType::bubble_cleaning);
  IndexType destination_1_4 = jump_event.AddDestination(jump1,
                                                        11,
                                                        Side::source,
                                                        JumpType::unplug_annihilation);

  // assert
  JumpType destination_type_1_1 = jump_event.GetDestinationType(jump1, destination_1_1);
  JumpType destination_type_1_2 = jump_event.GetDestinationType(jump1, destination_1_2);
  JumpType destination_type_1_3 = jump_event.GetDestinationType(jump1, destination_1_3);
  JumpType destination_type_1_4 = jump_event.GetDestinationType(jump1, destination_1_4);

  EXPECT_TRUE(destination_type_1_1 == JumpType::creation);
  EXPECT_TRUE(destination_type_1_2 == JumpType::annihilation);
  EXPECT_TRUE(destination_type_1_3 == JumpType::bubble_cleaning);
  EXPECT_TRUE(destination_type_1_4 == JumpType::unplug_annihilation);
}

TEST(JumpEventInfo, InsertJumpWithWithCrossed_Return_Correct_First_Cross_Jumps) {
  // arrange
  InterfaceJumpEventInfo jump_event;

  // act
  jump_event.InsertJump(0, Side::target, 0);
  jump_event.InsertJump(1, Side::target, 1);
  jump_event.InsertJump(2, Side::source, 2);
  jump_event.InsertJump(6, Side::source, 1);
  jump_event.InsertJump(9, Side::source, 0);

  // assert
  auto jump_iterator = jump_event.begin();
  bool is_first_in_cross_1 = jump_iterator->IsFirstInCross();
  jump_iterator++;
  bool is_first_in_cross_2 = jump_iterator->IsFirstInCross();
  jump_iterator++;
  bool is_first_in_cross_3 = jump_iterator->IsFirstInCross();
  jump_iterator++;
  bool is_first_in_cross_4 = jump_iterator->IsFirstInCross();
  jump_iterator++;
  bool is_first_in_cross_5 = jump_iterator->IsFirstInCross();

  EXPECT_TRUE(is_first_in_cross_1);
  EXPECT_TRUE(is_first_in_cross_2);
  EXPECT_FALSE(is_first_in_cross_3);
  EXPECT_FALSE(is_first_in_cross_4);
  EXPECT_FALSE(is_first_in_cross_5);
}

#endif  // TEST_SRC_FLOW_SIMULATOR_JUMP_EVENT_INFO_TESTS_H_
