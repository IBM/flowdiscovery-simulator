/**
 * \file test/src/flow_simulator/dynamic_capillary_interface_tests.h
 * \brief Contains regression tests of \c CapillaryInterface class methods
 *
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 */

#ifndef TEST_SRC_FLOW_SIMULATOR_DYNAMIC_CAPILLARY_INTERFACE_TESTS_H_
#define TEST_SRC_FLOW_SIMULATOR_DYNAMIC_CAPILLARY_INTERFACE_TESTS_H_

#include <gtest/gtest.h>
#include <map>
#include <string>
#include <memory>
#include <utility>
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network_context.h"

using Side = simulator::Side;
using IndexType = simulator::IndexType;
using CapillaryInterface = simulator::CapillaryInterface;
using DynamicCapillaryNetworkContext = simulator::DynamicCapillaryNetworkContext;

TEST(CapillaryInterface, GetFluidAtCapillaryAtSource_Return_0_For_Both_Capillaries) {
  // Arrange
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {1U, 1U, 3U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 2U;
  interface->InitialiseCapillaryInterface();

  // Assert
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(0U, Side::source), 0U);
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(1U, Side::source), 0U);
}

TEST(CapillaryInterface, GetFluidAtCapillaryAtTarget_Return_0_For_Both_Capillaries) {
  // Arrange
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {1U, 1U, 3U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 2U;
  interface->InitialiseCapillaryInterface();

  // Assert
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(0U, Side::target), 0U);
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(1U, Side::target), 0U);
}

TEST(CapillaryInterface, FlipFluidAtSource_Return_Flipped_Fluid) {
  // Arrange
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {1U, 1U, 2U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 1U;
  interface->InitialiseCapillaryInterface();

  // Act
  IndexType fluid_before = interface->GetFluidAtCapillaryAtSide(0, Side::source);
  interface->FlipFluidAtSource(0U);
  IndexType fluid_after = interface->GetFluidAtCapillaryAtSide(0, Side::source);

  // Assert
  EXPECT_EQ(fluid_before, 0U);
  EXPECT_EQ(fluid_after, 1U);
}

TEST(CapillaryInterface, GetNumberOfInterfacesAtCapillary_Return_0_For_Both_Capillaries) {
  // Arrange
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {1U, 1U, 3U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 2U;
  interface->InitialiseCapillaryInterface();

  // Assert
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0U), 0U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1U), 0U);
}

TEST(CapillaryInterface, AddNewInterfaceAtSource_Return_0_Interfaces_Before_And_1_After) {
  // Arrange
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {1U, 1U, 2U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 1U;
  interface->InitialiseCapillaryInterface();

  // Act
  arma::uword number_of_interfaces_before = interface->GetNumberOfInterfacesAtCapillary(0U);
  interface->AddNewInterface(0U, Side::source);
  arma::uword number_of_interfaces_after = interface->GetNumberOfInterfacesAtCapillary(0U);

  // Assert
  EXPECT_EQ(number_of_interfaces_before, 0U);
  EXPECT_EQ(number_of_interfaces_after, 1U);
}

TEST(CapillaryInterface, AddNewInterfaceAtTarget_Return_0_Interfaces_Before_And_1_After) {
  // Arrange
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {1U, 1U, 2U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 1U;
  interface->InitialiseCapillaryInterface();

  // Act
  arma::uword number_of_interfaces_before = interface->GetNumberOfInterfacesAtCapillary(0U);
  interface->AddNewInterface(0U, Side::target);
  arma::uword number_of_interfaces_after = interface->GetNumberOfInterfacesAtCapillary(0U);

  // Assert
  EXPECT_EQ(number_of_interfaces_before, 0U);
  EXPECT_EQ(number_of_interfaces_after, 1U);
}

TEST(CapillaryInterface, RemoveInterface_Return_1_Interface_Before_And_0_After) {
  // Arrange
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {1U, 1U, 2U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 1U;
  interface->InitialiseCapillaryInterface();
  interface->AddNewInterface(0U, Side::source);

  // Act
  arma::uword number_of_interfaces_before = interface->GetNumberOfInterfacesAtCapillary(0U);
  interface->RemoveInterface(0U, Side::source);
  arma::uword number_of_interfaces_after = interface->GetNumberOfInterfacesAtCapillary(0U);

  // Assert
  EXPECT_EQ(number_of_interfaces_before, 1U);
  EXPECT_EQ(number_of_interfaces_after, 0U);
}

TEST(CapillaryInterface, GetNumberOfInterfaces_Return_2_Interfaces) {
  // Arrange
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {1U, 1U, 3U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 2U;
  interface->InitialiseCapillaryInterface();
  interface->AddNewInterface(0U, Side::source);
  interface->AddNewInterface(1U, Side::source);

  // Assert
  EXPECT_EQ(interface->GetNumberOfInterfaces(), 2U);
}

TEST(CapillaryInterface,
  IsNumberOfInterfacesAtCapillaryOdd_Return_True_For_First_And_False_For_Second) {
  // Arrange
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {1U, 1U, 3U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 2U;
  interface->InitialiseCapillaryInterface();
  interface->AddNewInterface(0, Side::source);

  // Assert
  EXPECT_TRUE(interface->IsNumberOfInterfacesAtCapillaryOdd(0U));
  EXPECT_FALSE(interface->IsNumberOfInterfacesAtCapillaryOdd(1U));
}

TEST(CapillaryInterface, GetIndexOfCapillariesWithInterfaces_Return_Vector_With_Index_0_And_2) {
  // Arrange
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {1U, 1U, 4U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 3U;
  interface->InitialiseCapillaryInterface();
  interface->AddNewInterface(2U, Side::source);
  interface->AddNewInterface(0U, Side::source);

  // Act
  arma::Col<IndexType> capillaries = interface->GetIndexOfCapillariesWithInterfaces();

  // Assert
  EXPECT_EQ(capillaries(0U), 0U);
  EXPECT_EQ(capillaries(1U), 2U);
}

TEST(CapillaryInterface, GetInterfacePositionAtCapillaryNearTarget_Return_Near_1_At_Both) {
  // Arrange
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {1U, 1U, 3U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 2U;
  interface->InitialiseCapillaryInterface();
  interface->AddNewInterface(0U, Side::target);
  interface->AddNewInterface(1U, Side::source);
  interface->AddNewInterface(1U, Side::target);

  // Assert
  double tolerance = 1.0e-15;
  EXPECT_NEAR(interface->GetInterfacePositionAtCapillaryNearSide(0U, Side::target), 1.0, tolerance);
  EXPECT_NEAR(interface->GetInterfacePositionAtCapillaryNearSide(1U, Side::target), 1.0, tolerance);
}

TEST(CapillaryInterface, GetInterfacePositionAtCapillaryNearSource_Return_Near_0_At_Both) {
  // Arrange
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {1U, 1U, 3U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 2U;
  interface->InitialiseCapillaryInterface();
  interface->AddNewInterface(0U, Side::source);
  interface->AddNewInterface(1U, Side::source);
  interface->AddNewInterface(1U, Side::target);

  // Assert
  double tol = 1.0e-15;
  EXPECT_NEAR(interface->GetInterfacePositionAtCapillaryNearSide(0, Side::source), 0.0, tol);
  EXPECT_NEAR(interface->GetInterfacePositionAtCapillaryNearSide(1, Side::source), 0.0, tol);
}

TEST(CapillaryInterface,
  MoveInterfacePositionsAtCapillaryByDelta_Return_Correct_Value_For_Both_Interfaces) {
  // Arrange
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {1U, 1U, 2U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 1U;
  interface->InitialiseCapillaryInterface();
  interface->AddNewInterface(0U, Side::source);

  // Act
  interface->MoveInterfacePositionsAtCapillaryByDelta(0U, 0.1);
  interface->AddNewInterface(0U, Side::source);
  interface->MoveInterfacePositionsAtCapillaryByDelta(0U, 0.2);

  // Assert
  double tol = 1.0e-15;
  EXPECT_NEAR(interface->GetInterfacePositionAtCapillaryNearSide(0U, Side::target), 0.3, tol);
  EXPECT_NEAR(interface->GetInterfacePositionAtCapillaryNearSide(0U, Side::source), 0.2, tol);
}

TEST(CapillaryInterface,
  CalculateEffectiveInterfacePositionAtCapillary_Return_The_Effective_Interface_Value) {
  // Arrange
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {1U, 1U, 2U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 1U;
  interface->InitialiseCapillaryInterface();
  interface->AddNewInterface(0U, Side::source);
  interface->MoveInterfacePositionsAtCapillaryByDelta(0U, 0.1);
  interface->AddNewInterface(0U, Side::source);
  interface->MoveInterfacePositionsAtCapillaryByDelta(0U, 0.2);

  // Assert
  EXPECT_DOUBLE_EQ(interface->CalculateEffectiveInterfacePositionAtCapillary(0U), 0.9);
}

TEST(CapillaryInterface, GetFluidAtCapillaryAtTarget_Return_1_When_There_Is_One_Interface) {
  // Arrange
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {1U, 1U, 2U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 1U;
  interface->InitialiseCapillaryInterface();
  interface->AddNewInterface(0U, Side::source);

  // Assert
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(0U, Side::target), 1U);
}

TEST(CapillaryInterface, InjectFluidOnInletNodes_Inject_Fluid_Properly) {
  // Arrange
  /*   0
      ▌▌▄▄
         ▝▀▙
            ▝▀▙  1      2
              ▌▌▄▄▄▄▄▄▖▌▌
            ▗▄▛▘
      3   ▗▄▛▘
     ▌▌▗▄▛▀
*/
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {3U, 1U, 3U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 3;
  context->inlet_nodes_ = {0, 3};
  context->link_direction_ = {+1, -1, +1, -1, +1, -1};
  context->linked_nodes_ = {{0U, 1U, 1U, 2U, 1U, 3U},
                            {0U, 0U, 1U, 1U, 2U, 2U}};
  context->capillaries_whose_source_is_ = {{0U}, {1U, 2U}, {}, {}};
  context->capillaries_whose_target_is_ = {{}, {0U}, {1U}, {2U}};

  interface->InitialiseCapillaryInterface();

  // Act
  interface->InjectFluidOnInletNodes();

  // Assert
  EXPECT_EQ(interface->GetNumberOfInterfaces(), 2U);
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(0U, Side::source), 1U);
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(2U, Side::source), 0U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0U), 1U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2U), 1U);
}

TEST(CapillaryInterface, UpdateInterfaceInformation_Jump_From_Second_To_Other_Capillaries) {
  // Arrange
  /*   0
      ▌▌▄▄
         ▝▀▙
            ▝▀▙  1      2
              ▌▌▄▄▄▄▄▄▖▌▌
            ▗▄▛▘
      3   ▗▄▛▘
     ▌▌▗▄▛▀
*/
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {3U, 1U, 3U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  arma::Col<IndexType> changed_caps;

  context->number_of_links_ = 3U;
  context->inlet_nodes_ = {0U, 3U};
  context->link_direction_ = {+1, -1, +1, -1, +1, -1};
  context->linked_nodes_ = {{0U, 1U, 1U, 2U, 1U, 3U},
                            {0U, 0U, 1U, 1U, 2U, 2U}};
  context->capillaries_whose_source_is_ = {{0U}, {1U, 2U}, {}, {}};
  context->capillaries_whose_target_is_ = {{}, {0U}, {1U}, {2U}};

  interface->InitialiseCapillaryInterface();

  arma::Col<IndexType> caps_with_roots = {1U};
  arma::Col<IndexType> caps_roots_directions = {0U};
  arma::vec::fixed<3> x_vec = {0.6, 0.0, 0.2};

  interface->AddNewInterface(0U, Side::source);
  interface->FlipFluidAtSource(0U);
  interface->MoveInterfacePositionsAtCapillaryByDelta(0U, 0.8);
  interface->AddNewInterface(1U, Side::source);
  interface->MoveInterfacePositionsAtCapillaryByDelta(1U, 0.2);
  interface->AddNewInterface(2U, Side::source);
  interface->MoveInterfacePositionsAtCapillaryByDelta(2U, 0.4);

  // Act
  interface->ComputeDeltas(x_vec);
  interface->UpdateInterfaceInformation();
  interface->Jump(caps_with_roots, caps_roots_directions, changed_caps);

  // Assert
  EXPECT_EQ(interface->GetNumberOfInterfaces(), 4U);
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(0U, Side::source), 1U);
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(2U, Side::source), 1U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0U), 2U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2U), 2U);
  EXPECT_DOUBLE_EQ(interface->GetInterfacePositionAtCapillaryNearSide(0U, Side::source), 0.6);
  EXPECT_DOUBLE_EQ(interface->GetInterfacePositionAtCapillaryNearSide(2U, Side::target), 0.2);
  EXPECT_EQ(changed_caps.n_elem, 3U);
  EXPECT_EQ(changed_caps(0U), 0);
  EXPECT_EQ(changed_caps(1U), 1);
  EXPECT_EQ(changed_caps(2U), 2);
}

TEST(CapillaryInterface, PlugCapillary_Return_False_Before_And_True_After) {
  // Arrange
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {1U, 1U, 2U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 1U;
  interface->InitialiseCapillaryInterface();

  // Act
  bool is_capillary_plugged_before = interface->IsCapillaryPlugged(0U);
  interface->PlugCapillary(0U, Side::target);
  bool is_capillary_plugged_after = interface->IsCapillaryPlugged(0U);

  // Assert
  EXPECT_FALSE(is_capillary_plugged_before);
  EXPECT_TRUE(is_capillary_plugged_after);
}

TEST(CapillaryInterface, CheckPlugging_Test_Plug_After_Jump_Return_True_For_Plugged) {
  // Arrange
  /*   0
      ▌▌▄▄
         ▝▀▙
            ▝▀▙  1      2
              ▌▌▄▄▄▄▄▄▖▌▌
            ▗▄▛▘
      3   ▗▄▛▘
     ▌▌▗▄▛▀
*/
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {3U, 1U, 3U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  arma::Col<IndexType> changed_caps;

  context->number_of_links_ = 3U;
  context->inlet_nodes_ = {0U, 3U};
  context->link_direction_ = {+1, -1, +1, -1, +1, -1};
  context->linked_nodes_ = {{0U, 1U, 1U, 2U, 1U, 3U},
                            {0U, 0U, 1U, 1U, 2U, 2U}};
  context->capillaries_whose_source_is_ = {{0U}, {1U, 2U}, {}, {}};
  context->capillaries_whose_target_is_ = {{}, {0U}, {1U}, {2U}};

  interface->InitialiseCapillaryInterface();

  arma::Col<IndexType> caps_with_roots = {0U};
  arma::Col<IndexType> caps_roots_directions = {0U};
  arma::vec::fixed<3> x_vec_initial = {1.0, 0.4, 1.0};

  interface->AddNewInterface(0U, Side::source);
  interface->MoveInterfacePositionsAtCapillaryByDelta(0U, 0.8);
  interface->AddNewInterface(1U, Side::source);
  interface->MoveInterfacePositionsAtCapillaryByDelta(1U, 0.2);

  interface->ComputeDeltas(x_vec_initial);
  interface->UpdateInterfaceInformation();
  interface->Jump(caps_with_roots, caps_roots_directions, changed_caps);

  // Act
  caps_with_roots = {1U};
  arma::vec::fixed<3> x_vec = {1.0, 0.4, 0.0};
  interface->ComputeDeltas(x_vec);

  // Assert
  EXPECT_TRUE(interface->CheckPlugging(caps_with_roots,
                                       caps_roots_directions));
  EXPECT_TRUE(interface->IsCapillaryPlugged(0U));
  EXPECT_FALSE(interface->IsCapillaryPlugged(1U));
  EXPECT_FALSE(interface->IsCapillaryPlugged(2U));
  EXPECT_EQ(interface->GetCapillaryPluggedSide(0U), Side::target);
  EXPECT_EQ(changed_caps.n_elem, 1U);
  EXPECT_EQ(changed_caps(0U), 0);
}

TEST(CapillaryInterface, RevertJump_Test_After_Plug_Return_Reverted_Interfaces) {
  // Arrange
  /*   0
      ▌▌▄▄
         ▝▀▙
            ▝▀▙  1      2
              ▌▌▄▄▄▄▄▄▖▌▌
            ▗▄▛▘
      3   ▗▄▛▘
     ▌▌▗▄▛▀
  */
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {3U, 1U, 3U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  arma::Col<IndexType> changed_caps;

  context->number_of_links_ = 3U;
  context->inlet_nodes_ = {0, 3};
  context->link_direction_ = {+1, -1, +1, -1, +1, -1};
  context->linked_nodes_ = {{0U, 1U, 1U, 2U, 1U, 3U},
                            {0U, 0U, 1U, 1U, 2U, 2U}};
  context->capillaries_whose_source_is_ = {{0U}, {1U, 2U}, {}, {}};
  context->capillaries_whose_target_is_ = {{}, {0U}, {1U}, {2U}};

  interface->InitialiseCapillaryInterface();

  arma::Col<IndexType> caps_with_roots = {0U};
  arma::Col<IndexType> caps_roots_directions = {1U};
  arma::vec::fixed<3> x_vec_initial = {1.0, 0.4, 1.0};

  interface->AddNewInterface(0U, Side::source);
  interface->FlipFluidAtSource(0U);
  interface->MoveInterfacePositionsAtCapillaryByDelta(0U, 0.8);
  interface->AddNewInterface(1U, Side::source);
  interface->MoveInterfacePositionsAtCapillaryByDelta(1U, 0.2);

  interface->ComputeDeltas(x_vec_initial);
  interface->UpdateInterfaceInformation();
  interface->Jump(caps_with_roots, caps_roots_directions, changed_caps);

  caps_with_roots = {1U};
  caps_roots_directions = {0U};
  arma::vec::fixed<3> x_vec = {1.0, 0.4, 0.0};
  interface->ComputeDeltas(x_vec);
  changed_caps.set_size(0U);
  interface->CheckPlugging(caps_with_roots, caps_roots_directions);

  // Act
  interface->RevertJump(changed_caps);

  // Assert
  EXPECT_EQ(interface->GetNumberOfInterfaces(), 2U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0U), 1U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2U), 0U);
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(0U, Side::source), 1U);
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(2U, Side::source), 0U);
  EXPECT_DOUBLE_EQ(interface->GetInterfacePositionAtCapillaryNearSide(0U, Side::source), 1.0);
  EXPECT_DOUBLE_EQ(interface->GetInterfacePositionAtCapillaryNearSide(1U, Side::target), 0.4);
  EXPECT_EQ(changed_caps.n_elem, 3U);
  EXPECT_EQ(changed_caps(0U), 0);
  EXPECT_EQ(changed_caps(1U), 1);
  EXPECT_EQ(changed_caps(2U), 2);
}

TEST(CapillaryInterface, Jump_Bubble_Cleaning) {
  // Arrange
  /*   0
      ▌▌▄▄
         ▝▀▙
            ▝▀▙  1      2
              ▌▌▄▄▄▄▄▄▖▌▌
            ▗▄▛▘
      3   ▗▄▛▘
     ▌▌▗▄▛▀
*/
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {3U, 1U, 3U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  arma::Col<IndexType> changed_caps;

  context->number_of_links_ = 3U;
  context->inlet_nodes_ = {0U, 3U};
  context->link_direction_ = {+1, -1, +1, -1, +1, -1};
  context->linked_nodes_ = {{0U, 1U, 1U, 2U, 1U, 3U},
                            {0U, 0U, 1U, 1U, 2U, 2U}};
  context->capillaries_whose_source_is_ = {{0U}, {1U, 2U}, {}, {}};
  context->capillaries_whose_target_is_ = {{}, {0U}, {1U}, {2U}};

  interface->InitialiseCapillaryInterface();

  interface->AddNewInterface(0U, Side::source);
  interface->FlipFluidAtSource(0U);
  interface->MoveInterfacePositionsAtCapillaryByDelta(0U, 0.8);
  interface->AddNewInterface(1U, Side::source);
  interface->MoveInterfacePositionsAtCapillaryByDelta(1U, 0.2);

  arma::vec::fixed<3> x_vec_initial = {1.0, 0.4, 1.0};

  interface->ComputeDeltas(x_vec_initial);
  interface->UpdateInterfaceInformation();

  // Inserting Bubble
  interface->AddNewInterface(1U, Side::source);
  interface->AddNewInterface(1U, Side::source);

  // Act
  auto number_of_interfaces_before = interface->GetNumberOfInterfaces();

  arma::Col<IndexType> caps_with_roots = {0U};
  arma::Col<IndexType> caps_roots_directions = {1U};

  interface->Jump(caps_with_roots, caps_roots_directions, changed_caps);
  auto number_of_interfaces_after = interface->GetNumberOfInterfaces();

  // Assert
  EXPECT_EQ(number_of_interfaces_before, 4U);
  EXPECT_EQ(number_of_interfaces_after, 3U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0U), 0U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1U), 2U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2U), 1U);
  EXPECT_DOUBLE_EQ(interface->GetInterfacePositionAtCapillaryNearSide(1, Side::target), 0.4);
  EXPECT_EQ(changed_caps.n_elem, 3U);
  EXPECT_EQ(changed_caps(0U), 0);
  EXPECT_EQ(changed_caps(1U), 1);
  EXPECT_EQ(changed_caps(2U), 2);
}

TEST(CapillaryInterface, Jump_Unplug_Annihilation) {
  // Arrange
  /*   0
      ▌▌▄▄
         ▝▀▙
            ▝▀▙  1      2
              ▌▌▄▄▄▄▄▄▖▌▌
            ▗▄▛▘
      3   ▗▄▛▘
     ▌▌▗▄▛▀
*/
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {3U, 1U, 3U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  arma::Col<IndexType> changed_caps;

  context->number_of_links_ = 3U;
  context->inlet_nodes_ = {0U, 3U};
  context->link_direction_ = {+1, -1, +1, -1, +1, -1};
  context->linked_nodes_ = {{0U, 1U, 1U, 2U, 1U, 3U},
                            {0U, 0U, 1U, 1U, 2U, 2U}};
  context->capillaries_whose_source_is_ = {{0U}, {1U, 2U}, {}, {}};
  context->capillaries_whose_target_is_ = {{}, {0U}, {1U}, {2U}};

  interface->InitialiseCapillaryInterface();

  arma::Col<IndexType> caps_with_interfaces_not_plugged_ = {0U};
  arma::Col<IndexType> caps_roots_directions = {1U};
  arma::vec::fixed<3> x_vec_initial = {1.0, 0.4, 1.0};

  interface->AddNewInterface(0U, Side::source);
  interface->FlipFluidAtSource(0U);
  interface->MoveInterfacePositionsAtCapillaryByDelta(0U, 0.8);
  interface->AddNewInterface(1U, Side::source);
  interface->MoveInterfacePositionsAtCapillaryByDelta(1U, 0.2);

  interface->ComputeDeltas(x_vec_initial);
  interface->UpdateInterfaceInformation();

  // Inserting Bubble plus one interface at distance arma::eps
  interface->AddNewInterface(1U, Side::source);
  interface->AddNewInterface(1U, Side::source);
  interface->AddNewInterface(1U, Side::source);

  // Act
  auto number_of_interfaces_before = interface->GetNumberOfInterfaces();
  interface->Jump(caps_with_interfaces_not_plugged_, caps_roots_directions, changed_caps);
  auto number_of_interfaces_after = interface->GetNumberOfInterfaces();

  // Assert
  EXPECT_EQ(number_of_interfaces_before, 5U);
  EXPECT_EQ(number_of_interfaces_after, 2U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0U), 0U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1U), 1U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2U), 1U);
  EXPECT_DOUBLE_EQ(interface->GetInterfacePositionAtCapillaryNearSide(1U, Side::target), 0.4);
  EXPECT_EQ(changed_caps.n_elem, 3U);
  EXPECT_EQ(changed_caps(0U), 0);
  EXPECT_EQ(changed_caps(1U), 1);
  EXPECT_EQ(changed_caps(2U), 2);
}

TEST(CapillaryInterface, RevertJump_Bubble_Cleaning) {
  // Arrange
  /*   0
      ▌▌▄▄
         ▝▀▙
            ▝▀▙  1      2
              ▌▌▄▄▄▄▄▄▖▌▌
            ▗▄▛▘
      3   ▗▄▛▘
     ▌▌▗▄▛▀
  */
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {3U, 1U, 3U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  arma::Col<IndexType> changed_caps;

  context->number_of_links_ = 3U;
  context->inlet_nodes_ = {0U, 3U};
  context->link_direction_ = {+1, -1, +1, -1, +1, -1};
  context->linked_nodes_ = {{0U, 1U, 1U, 2U, 1U, 3U},
                            {0U, 0U, 1U, 1U, 2U, 2U}};
  context->capillaries_whose_source_is_ = {{0U}, {1U, 2U}, {}, {}};
  context->capillaries_whose_target_is_ = {{}, {0U}, {1U}, {2U}};

  interface->InitialiseCapillaryInterface();

  interface->AddNewInterface(0U, Side::source);
  interface->FlipFluidAtSource(0U);
  interface->MoveInterfacePositionsAtCapillaryByDelta(0U, 0.8);
  interface->AddNewInterface(1U, Side::source);
  interface->MoveInterfacePositionsAtCapillaryByDelta(1U, 0.2);

  arma::vec::fixed<3> x_vec_initial = {1.0, 0.4, 1.0};

  interface->ComputeDeltas(x_vec_initial);
  interface->UpdateInterfaceInformation();

  // Inserting Bubble
  interface->AddNewInterface(1U, Side::source);
  interface->AddNewInterface(1U, Side::source);

  arma::Col<IndexType> caps_with_roots = {0U};
  arma::Col<IndexType> caps_roots_directions = {1U};

  interface->Jump(caps_with_roots, caps_roots_directions, changed_caps);

  caps_with_roots = {1U};
  caps_roots_directions = {0U};
  arma::vec::fixed<3> x_vec = {1.0, 0.4, 0.0};
  interface->ComputeDeltas(x_vec);
  interface->CheckPlugging(caps_with_roots, caps_roots_directions);

  // Act
  auto number_of_interfaces_before = interface->GetNumberOfInterfaces();
  changed_caps.set_size(0U);
  interface->RevertJump(changed_caps);
  auto number_of_interfaces_after = interface->GetNumberOfInterfaces();

  // Assert
  EXPECT_EQ(number_of_interfaces_before, 3U);
  EXPECT_EQ(number_of_interfaces_after, 2U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0U), 1U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2U), 0U);
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(0U, Side::source), 1U);
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(2U, Side::source), 0U);
  EXPECT_DOUBLE_EQ(interface->GetInterfacePositionAtCapillaryNearSide(0U, Side::source), 1.0);
  EXPECT_DOUBLE_EQ(interface->GetInterfacePositionAtCapillaryNearSide(1U, Side::target), 0.4);
  EXPECT_EQ(changed_caps.n_elem, 3U);
  EXPECT_EQ(changed_caps(0U), 0);
  EXPECT_EQ(changed_caps(1U), 1);
  EXPECT_EQ(changed_caps(2U), 2);
}

TEST(CapillaryInterface, RevertJump_Unplug_Annihilation) {
  // Arrange
  /*   0
      ▌▌▄▄
         ▝▀▙
            ▝▀▙  1      2
              ▌▌▄▄▄▄▄▄▖▌▌
            ▗▄▛▘
      3   ▗▄▛▘
     ▌▌▗▄▛▀
  */
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {3U, 1U, 3U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  arma::Col<IndexType> changed_caps;

  context->number_of_links_ = 3U;
  context->inlet_nodes_ = {0U, 3U};
  context->link_direction_ = {+1, -1, +1, -1, +1, -1};
  context->linked_nodes_ = {{0U, 1U, 1U, 2U, 1U, 3U},
                            {0U, 0U, 1U, 1U, 2U, 2U}};
  context->capillaries_whose_source_is_ = {{0U}, {1U, 2U}, {}, {}};
  context->capillaries_whose_target_is_ = {{}, {0U}, {1U}, {2U}};

  interface->InitialiseCapillaryInterface();

  arma::Col<IndexType> caps_with_roots = {0U};
  arma::Col<IndexType> caps_roots_directions = {1U};
  arma::vec::fixed<3> x_vec_initial = {1.0, 0.4, 1.0};

  interface->AddNewInterface(0U, Side::source);
  interface->FlipFluidAtSource(0U);
  interface->MoveInterfacePositionsAtCapillaryByDelta(0U, 0.8);
  interface->AddNewInterface(1U, Side::source);
  interface->MoveInterfacePositionsAtCapillaryByDelta(1U, 0.2);

  interface->ComputeDeltas(x_vec_initial);
  interface->UpdateInterfaceInformation();

  // Inserting Bubble plus one interface at distance arma::datum::eps
  interface->AddNewInterface(1U, Side::source);
  interface->AddNewInterface(1U, Side::source);
  interface->AddNewInterface(1U, Side::source);

  interface->Jump(caps_with_roots, caps_roots_directions, changed_caps);

  caps_with_roots = {1U};
  caps_roots_directions = {0U};
  arma::vec::fixed<3> x_vec = {1.0, 0.4, 0.0};
  interface->ComputeDeltas(x_vec);
  interface->CheckPlugging(caps_with_roots, caps_roots_directions);

  // Act
  auto number_of_interfaces_before = interface->GetNumberOfInterfaces();
  changed_caps.set_size(0U);
  interface->RevertJump(changed_caps);
  auto number_of_interfaces_after = interface->GetNumberOfInterfaces();

  // Assert
  EXPECT_EQ(number_of_interfaces_before, 2U);
  EXPECT_EQ(number_of_interfaces_after, 3U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0U), 1U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2U), 0U);
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(0U, Side::source), 1U);
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(2U, Side::source), 0U);
  EXPECT_DOUBLE_EQ(interface->GetInterfacePositionAtCapillaryNearSide(0U, Side::source), 1.0);
  EXPECT_DOUBLE_EQ(interface->GetInterfacePositionAtCapillaryNearSide(1U, Side::target), 0.4);
  EXPECT_EQ(changed_caps.n_elem, 3U);
  EXPECT_EQ(changed_caps(0U), 0);
  EXPECT_EQ(changed_caps(1U), 1);
  EXPECT_EQ(changed_caps(2U), 2);
}

TEST(CapillaryInterface, RevertJumpAndUnplugInletCapillary) {
  // Arrange
  /*   0
      ▌▌▄▄
         ▝▀▙
            ▝▀▙  1      2
              ▌▌▄▄▄▄▄▄▖▌▌
            ▗▄▛▘
      3   ▗▄▛▘
     ▌▌▗▄▛▀
  */
  double eps = arma::datum::eps;
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {3U, 1U, 3U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  arma::Col<IndexType> changed_caps;

  context->number_of_links_ = 3U;
  context->inlet_nodes_ = {0U, 3U};
  context->link_direction_ = {+1, -1, +1, -1, +1, -1};
  context->linked_nodes_ = {{0U, 1U, 1U, 2U, 1U, 3U},
                            {0U, 0U, 1U, 1U, 2U, 2U}};
  context->capillaries_whose_source_is_ = {{0U}, {1U, 2U}, {}, {}};
  context->capillaries_whose_target_is_ = {{}, {0U}, {1U}, {2U}};

  interface->InitialiseCapillaryInterface();

  interface->AddNewInterface(2U, Side::source);
  interface->FlipFluidAtSource(2U);

  // The first act, jump

  arma::vec::fixed<3> x_vec = {0.0, 0.0, -eps};
  arma::Col<IndexType> roots_found = {2U};
  arma::Col<IndexType> caps_roots_directions = {0U};

  interface->ComputeDeltas(x_vec);
  interface->UpdateInterfaceInformation();
  interface->Jump(roots_found, caps_roots_directions, changed_caps);

  // The second act, revert jump while interface at 0 reaches inlet node

  // Add new interface at capillary 0 that will reach the inlet node
  interface->AddNewInterface(0U, Side::source);
  interface->FlipFluidAtSource(0U);

  roots_found = {0U, 1U};
  caps_roots_directions = {0U, 0U};
  x_vec = {-eps, -eps, 0.0};

  interface->ComputeDeltas(x_vec);
  interface->UpdateInterfaceInformation();

  // Remove capillaries plugged because of jumps through inlet nodes
  arma::Col<IndexType> caps_to_remove = interface->PlugCapillariesJumpAtInletNode(
    roots_found,
    caps_roots_directions);
  for (const auto &capillary : caps_to_remove) {
    roots_found.shed_row(arma::as_scalar(arma::find(roots_found == capillary)));
  }

  interface->CheckPlugging(roots_found, caps_roots_directions);

  bool is_zero_plugged_before = interface->IsCapillaryPlugged(0);
  auto number_of_interfaces_before = interface->GetNumberOfInterfaces();
  changed_caps.set_size(0);
  interface->RevertJump(changed_caps);
  interface->RevertPluggedCapillariesJumpAtInletNode(changed_caps);
  bool is_zero_plugged_after = interface->IsCapillaryPlugged(0);
  auto number_of_interfaces_after = interface->GetNumberOfInterfaces();

  // Assert
  EXPECT_TRUE(is_zero_plugged_before);
  EXPECT_FALSE(is_zero_plugged_after);
  EXPECT_EQ(number_of_interfaces_before, 3U);
  EXPECT_EQ(number_of_interfaces_after, 2U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(0U), 1U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(1U), 0U);
  EXPECT_EQ(interface->GetNumberOfInterfacesAtCapillary(2U), 1U);
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(0U, Side::source), 1U);
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(1U, Side::source), 0U);
  EXPECT_EQ(interface->GetFluidAtCapillaryAtSide(2U, Side::source), 1U);
  EXPECT_EQ(interface->GetInterfacePositionAtCapillaryNearSide(0U, Side::source), eps);
  EXPECT_DOUBLE_EQ(interface->GetInterfacePositionAtCapillaryNearSide(2, Side::source), eps);
  EXPECT_EQ(changed_caps.n_elem, 3U);
  EXPECT_EQ(changed_caps(0U), 0);
  EXPECT_EQ(changed_caps(1U), 1);
  EXPECT_EQ(changed_caps(2U), 2);
}

TEST(CapillaryInterface, ComputeDeltas) {
  // Arrange
  /*   0
      ▌▌▄▄
         ▝▀▙
            ▝▀▙  1      2
              ▌▌▄▄▄▄▄▄▖▌▌
            ▗▄▛▘
      3   ▗▄▛▘
     ▌▌▗▄▛▀
  */
  double voxel_size = 1.0;
  arma::uword flow_axis = 2U;
  arma::uvec::fixed<3> shape = {3U, 1U, 3U};
  arma::uword boundary_thickness = 1U;
  auto context = std::make_shared<DynamicCapillaryNetworkContext>(voxel_size,
                                                                  flow_axis,
                                                                  shape,
                                                                  boundary_thickness);
  auto interface = std::make_shared<CapillaryInterface>(context);
  context->number_of_links_ = 3U;
  context->inlet_nodes_ = {0U, 3U};
  context->link_direction_ = {+1, -1, +1, -1, +1, -1};
  context->linked_nodes_ = {{0U, 1U, 1U, 2U, 1U, 3U},
                            {0U, 0U, 1U, 1U, 2U, 2U}};
  context->capillaries_whose_source_is_ = {{0U}, {1U, 2U}, {}, {}};
  context->capillaries_whose_target_is_ = {{}, {0U}, {1U}, {2U}};

  // The first act, compute delta of capillary with one interface

  interface->InitialiseCapillaryInterface();
  interface->AddNewInterface(0U, Side::source);

  // x_vec(0) represents the effective interface position at capillary 0
  arma::vec::fixed<3> x_vec = {0.25, 0.0, 0.0};

  interface->ComputeDeltas(x_vec);  // Compute by effetive interface position
  interface->UpdateInterfaceInformation();

  double interface_position_1 = interface->GetInterfacePositionAtCapillaryNearSide(0, Side::source);

  // The second act, compute delta of capillary with one interface

  interface->AddNewInterface(0U, Side::source);

  // x_vec(0) represents the interface position near of source at capillary 0
  x_vec = {0.5, 0.0, 0.0};

  interface->ComputeDeltas(x_vec);  // Compute by minimal x: interface position near source
  interface->UpdateInterfaceInformation();

  double interface_position_2 = interface->GetInterfacePositionAtCapillaryNearSide(0, Side::source);
  double interface_position_3 = interface->GetInterfacePositionAtCapillaryNearSide(0, Side::target);

  // Assert
  EXPECT_EQ(interface_position_1, 0.25);
  EXPECT_EQ(interface_position_2, 0.5);
  EXPECT_EQ(interface_position_3, 0.75 - arma::datum::eps);
}

#endif  // TEST_SRC_FLOW_SIMULATOR_DYNAMIC_CAPILLARY_INTERFACE_TESTS_H_
