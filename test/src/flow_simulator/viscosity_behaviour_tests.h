/**
 * \file test/src/flow_simulator/viscosity_behaviour_tests.h
 * \brief Contains regression tests of \c StaticViscosityBehaviour, \c PureWaterViscosityBehaviour,
 * \c CarbonatedWaterViscosityBehaviour, \c SupercriticalCO2ViscosityBehaviour,
 * \c BrineViscosityBehaviour and \c CarbonatedBrineViscosityBehaviour
 * objects.
 *
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2020
 */

#ifndef TEST_SRC_FLOW_SIMULATOR_VISCOSITY_BEHAVIOUR_TESTS_H_
#define TEST_SRC_FLOW_SIMULATOR_VISCOSITY_BEHAVIOUR_TESTS_H_

#include <gtest/gtest.h>
#include <vector>

#include "src/flow_simulator/algorithms/constant_viscosity_behaviour.h"
#include "src/flow_simulator/algorithms/pure_water_viscosity_behaviour.h"
#include "src/flow_simulator/algorithms/carbonated_water_viscosity_behaviour.h"
#include "src/flow_simulator/algorithms/supercritical_co2_viscosity_behaviour.h"
#include "src/flow_simulator/algorithms/brine_viscosity_behaviour.h"
#include "src/flow_simulator/algorithms/carbonated_brine_viscosity_behaviour.h"
#include "src/flow_simulator/algorithms/seawater_viscosity_behaviour.h"

using std::vector;
using ConstantViscosityBehaviour = simulator::ConstantViscosityBehaviour;
using PureWaterViscosityBehaviour = simulator::PureWaterViscosityBehaviour;
using CarbonatedWaterViscosityBehaviour = simulator::CarbonatedWaterViscosityBehaviour;
using SupercriticalCO2ViscosityBehaviour = simulator::SupercriticalCO2ViscosityBehaviour;
using BrineViscosityBehaviour = simulator::BrineViscosityBehaviour;
using CarbonatedBrineViscosityBehaviour = simulator::CarbonatedBrineViscosityBehaviour;
using SeawaterViscosityBehaviour = simulator::SeawaterViscosityBehaviour;

TEST(ViscosityBehaviour, ConstantViscosityBehaviour_GetViscosity) {
  // Arrange
  double dynamic_viscosity_constant = 8.90e-4;
  double pressure = 101325.0;  // in pascal

  ConstantViscosityBehaviour constant_viscosity_behaviour(dynamic_viscosity_constant);

  double dynamic_viscosity = constant_viscosity_behaviour.GetViscosity(pressure);  // Pa.s

  // Assert
  EXPECT_DOUBLE_EQ(dynamic_viscosity, dynamic_viscosity_constant);  // Pa.s
}

TEST(ViscosityBehaviour, PureWaterViscosityBehaviour_GetViscosity) {
  // Arrange
  double temperature_25_degree_celsius = 298.15;
  double temperature_30_degree_celsius = 303.15;
  double temperature_35_degree_celsius = 308.15;
  double temperature_40_degree_celsius = 313.15;
  double pressure = 101325.0;  // in pascal

  PureWaterViscosityBehaviour pure_water_viscosity_behaviour1(temperature_25_degree_celsius);
  PureWaterViscosityBehaviour pure_water_viscosity_behaviour2(temperature_30_degree_celsius);
  PureWaterViscosityBehaviour pure_water_viscosity_behaviour3(temperature_35_degree_celsius);
  PureWaterViscosityBehaviour pure_water_viscosity_behaviour4(temperature_40_degree_celsius);

  double dynamic_viscosity1 = pure_water_viscosity_behaviour1.GetViscosity(pressure);  // Pa.s
  double dynamic_viscosity2 = pure_water_viscosity_behaviour2.GetViscosity(pressure);  // Pa.s
  double dynamic_viscosity3 = pure_water_viscosity_behaviour3.GetViscosity(pressure);  // Pa.s
  double dynamic_viscosity4 = pure_water_viscosity_behaviour4.GetViscosity(pressure);  // Pa.s

  // Assert
  EXPECT_DOUBLE_EQ(dynamic_viscosity1, 0.0008900289582347362);  // Pa.s
  EXPECT_DOUBLE_EQ(dynamic_viscosity2, 0.00079713097880318104);  // Pa.s
  EXPECT_DOUBLE_EQ(dynamic_viscosity3, 0.00071853518325418197);  // Pa.s
  EXPECT_DOUBLE_EQ(dynamic_viscosity4, 0.00065162625088585201);  // Pa.s
}

TEST(ViscosityBehaviour, PureWaterViscosityBehaviour_IsDesignedFor) {
  // Arrange
  double temperature_above_limit = 378.15 + 0.01;
  double temperature_below_limit = 293.15 - 0.01;
  double temperature_within_limits = 308.15;

  PureWaterViscosityBehaviour pure_water_viscosity_behaviour1(temperature_above_limit);
  PureWaterViscosityBehaviour pure_water_viscosity_behaviour2(temperature_below_limit);
  PureWaterViscosityBehaviour pure_water_viscosity_behaviour3(temperature_within_limits);

  // Assert
  EXPECT_FALSE(pure_water_viscosity_behaviour1.IsDesignedFor(1.0e6 - 1.0e5, 1.0e6 + 1.0e5));
  EXPECT_FALSE(pure_water_viscosity_behaviour2.IsDesignedFor(1.0e6 - 1.0e5, 1.0e6 + 1.0e5));
  EXPECT_TRUE(pure_water_viscosity_behaviour3.IsDesignedFor(1.0e6 - 1.0e5, 1.0e6 + 1.0e5));
  EXPECT_FALSE(pure_water_viscosity_behaviour3.IsDesignedFor(1.0e6 / 1.0e2, 1.0e6 + 1.0e1));
  EXPECT_FALSE(pure_water_viscosity_behaviour3.IsDesignedFor(1.0e6 - 1.0e5, 1.0e6 * 1.0e3));
}

TEST(ViscosityBehaviour, CarbonatedWaterViscosityBehaviour_GetViscosity) {
  // Arrange
  double temp_25_celsius = 298.15;
  double temp_30_celsius = 303.15;
  double temp_35_celsius = 308.15;
  double temp_40_celsius = 313.15;
  double pressure = 101325.0;  // in pascal
  double co2 = 2.71;  // 0 to 1 dimensionless

  CarbonatedWaterViscosityBehaviour carbonated_water_viscosity_behaviour1(temp_25_celsius, co2);
  CarbonatedWaterViscosityBehaviour carbonated_water_viscosity_behaviour2(temp_30_celsius, co2);
  CarbonatedWaterViscosityBehaviour carbonated_water_viscosity_behaviour3(temp_35_celsius, co2);
  CarbonatedWaterViscosityBehaviour carbonated_water_viscosity_behaviour4(temp_40_celsius, co2);

  double dynamic_viscosity1 = carbonated_water_viscosity_behaviour1.GetViscosity(pressure);  // Pa.s
  double dynamic_viscosity2 = carbonated_water_viscosity_behaviour2.GetViscosity(pressure);  // Pa.s
  double dynamic_viscosity3 = carbonated_water_viscosity_behaviour3.GetViscosity(pressure);  // Pa.s
  double dynamic_viscosity4 = carbonated_water_viscosity_behaviour4.GetViscosity(pressure);  // Pa.s

  // Assert
  EXPECT_DOUBLE_EQ(dynamic_viscosity1, 0.0009140124784973225);  // Pa.s
  EXPECT_DOUBLE_EQ(dynamic_viscosity2, 0.0008129104212821088);  // Pa.s
  EXPECT_DOUBLE_EQ(dynamic_viscosity3, 0.0007297756276792505);  // Pa.s
  EXPECT_DOUBLE_EQ(dynamic_viscosity4, 0.0006600526930944584);  // Pa.s
}

TEST(ViscosityBehaviour, CarbonatedWaterViscosityBehaviour_IsDesignedFor) {
  // This test checks the operation ranges of the viscosity modelling
  // 20C <= T <= 105C
  // 100kPa <= P <= 40MPa
  // 0% < CO2 < 4.8%
  // Unit Convertion
  // 1Bar -> 100kPa

  // Arrange
  double temperature_above_limit = 378.15 + 0.01;  // 105 C + 0.01 C
  double temperature_below_limit = 293.15 - 0.01;  // 20 C - 0.01 C
  double temperature_within_limits = 308.15;  // 35 C
  double co2 = 2.71;  // 0 to 4.8%
  double co2_upper = 4.81;  // 0 to 4.8%
  double co2_lower = 0;  // 0 to 4.8%

  CarbonatedWaterViscosityBehaviour carbonated_water_1(temperature_above_limit, co2);
  CarbonatedWaterViscosityBehaviour carbonated_water_2(temperature_below_limit, co2);
  CarbonatedWaterViscosityBehaviour carbonated_water_3(temperature_within_limits, co2);
  CarbonatedWaterViscosityBehaviour carbonated_water_4(temperature_within_limits, co2_upper);
  CarbonatedWaterViscosityBehaviour carbonated_water_5(temperature_within_limits, co2_lower);

  // Assert
  // Pressure is in the range, but Temperature is out of range
  // 900KPa, 1.1MPa
  EXPECT_FALSE(carbonated_water_1.IsDesignedFor(1.0e6 - 1.0e5, 1.0e6 + 1.0e5));
  // 900kPa, 1.1MPa
  EXPECT_FALSE(carbonated_water_2.IsDesignedFor(1.0e6 - 1.0e5, 1.0e6 + 1.0e5));
  // Temperature is in the range, but Pressure is out of range
  // 1kPa, 1MPa + 10Pa
  EXPECT_FALSE(carbonated_water_3.IsDesignedFor(1.0e6 / 1.0e2, 1.0e6 + 1.0e1));
  // 900kPa, 1GPa
  EXPECT_FALSE(carbonated_water_3.IsDesignedFor(1.0e6 - 1.0e5, 1.0e6 * 1.0e3));
  // Both temperature and Pressure are within the range
  // 900kPa, 1.1MPa
  EXPECT_TRUE(carbonated_water_3.IsDesignedFor(1.0e6 - 1.0e5, 1.0e6 + 1.0e5));
  // Temperature and Pressure are within the range but CO2 is out of range
  // 900kPa, 1.1MPa
  EXPECT_FALSE(carbonated_water_4.IsDesignedFor(1.0e6 - 1.0e5, 1.0e6 + 1.0e5));
  EXPECT_FALSE(carbonated_water_5.IsDesignedFor(1.0e6 - 1.0e5, 1.0e6 + 1.0e5));
}


TEST(ViscosityBehaviour, SupercriticalCO2ViscosityBehaviour_GetViscosity) {
  // Arrange
  double temp[4] = {298.15, 323.15, 473.15, 910.15};  // Kelvin
  double pressure = 10.0e6;  // in pascal

  SupercriticalCO2ViscosityBehaviour scco2_viscosity_behaviour1(temp[0]);
  SupercriticalCO2ViscosityBehaviour scco2_viscosity_behaviour2(temp[1]);
  SupercriticalCO2ViscosityBehaviour scco2_viscosity_behaviour3(temp[2]);
  SupercriticalCO2ViscosityBehaviour scco2_viscosity_behaviour4(temp[3]);

  double viscosity[4] = {scco2_viscosity_behaviour1.GetViscosity(pressure),   // Pa.s
                         scco2_viscosity_behaviour2.GetViscosity(pressure),   // Pa.s
                         scco2_viscosity_behaviour3.GetViscosity(pressure),   // Pa.s
                         scco2_viscosity_behaviour4.GetViscosity(pressure)};  // Pa.s

  // Assert
  ASSERT_NEAR(viscosity[0], 6.366658117001264e-05, 1.0e-17);  // Pa.s
  ASSERT_NEAR(viscosity[1], 3.259271283304327e-05, 3.0e-17);
  ASSERT_NEAR(viscosity[2], 2.3837638205007856e-05, 1.0e-17);
  ASSERT_NEAR(viscosity[3], 3.7690237411647795e-05, 1.0e-17);
}

TEST(ViscosityBehaviour, SupercriticalCO2ViscosityBehaviour_IsDesignedFor) {
  // This test checks the operation ranges of the viscosity modelling
  // 310K <= T <= 900K
  // 7.5MPa <= P <= 101.4MPa
  // Unit Convertion
  // 1Bar -> 100kPa

  // Arrange
  double temperature_above_limit = 900.01;
  double temperature_below_limit = 293.15;
  double temperature_within_limits = 408.15;

  SupercriticalCO2ViscosityBehaviour ScCO2_1(temperature_above_limit);
  SupercriticalCO2ViscosityBehaviour ScCO2_2(temperature_below_limit);
  SupercriticalCO2ViscosityBehaviour ScCO2_3(temperature_within_limits);

  // Assert
  // Pressure is in the range, but Temperature is out of range
  EXPECT_FALSE(ScCO2_1.IsDesignedFor(7.5e6, 101.4e6));
  EXPECT_FALSE(ScCO2_2.IsDesignedFor(7.5e6, 101.4e6));
  // Temperature is in the range, but Pressure is out of range
  EXPECT_FALSE(ScCO2_3.IsDesignedFor(7.49e6, 101.4e6));
  EXPECT_FALSE(ScCO2_3.IsDesignedFor(7.5e6, 101.5e6));
  // Both temperature and Pressure are within the range
  EXPECT_TRUE(ScCO2_3.IsDesignedFor(7.5e6, 101.4e6));
}

TEST(ViscosityBehaviour, BrineViscosityBehaviour_GetViscosity) {
  // Arrange
  arma::vec::fixed<4> temp({298.15, 303.15, 308.15, 313.15});  // Kelvin
  const double pressure = 101325.0;  // in pascal
  const double salinity = 175328.31;  //  parts per million


  BrineViscosityBehaviour viscosity_1(temp(0), salinity);
  BrineViscosityBehaviour viscosity_2(temp(1), salinity);
  BrineViscosityBehaviour viscosity_3(temp(2), salinity);
  BrineViscosityBehaviour viscosity_4(temp(3), salinity);

  arma::vec::fixed<4> viscosity({viscosity_1.GetViscosity(pressure),
                                 viscosity_2.GetViscosity(pressure),
                                 viscosity_3.GetViscosity(pressure),
                                 viscosity_4.GetViscosity(pressure)});  // Pa.s

  // Assert
  ASSERT_NEAR(viscosity(0), 1.2080804272531652E-03, 1.0e-17);  // Pa.s
  ASSERT_NEAR(viscosity(1), 1.0857236328013992E-03, 1.0e-17);
  ASSERT_NEAR(viscosity(2), 9.8237382021643506E-04, 1.0e-17);
  ASSERT_NEAR(viscosity(3), 8.9432761157560243E-04, 1.0e-17);
}

TEST(ViscosityBehaviour, BrineViscosityBehaviour_IsDesignedFor) {
  // This test checks the operation ranges of the viscosity modelling
  // 273K <= T <= 623K
  // 1Bar <= P <= 10^3 Bar
  // 0 mol/Kg < m <= 6 mol/Kg
  // Unit Convertion
  // 1Bar -> 100kPa

  // Arrange
  const double temperature_above_limit = 623.01;
  const double temperature_below_limit = 272.99;
  const double temperature_within_limits = 378.15;
  const double salinity = 4.0e4;  // parts per million

  BrineViscosityBehaviour Obj1(temperature_above_limit, salinity);
  BrineViscosityBehaviour Obj2(temperature_below_limit, salinity);
  BrineViscosityBehaviour Obj3(temperature_within_limits, salinity);

  // Assert
  // Pressure is in the range, but Temperature is out of range
  EXPECT_FALSE(Obj1.IsDesignedFor(0.1e6, 100.0e6));
  EXPECT_FALSE(Obj2.IsDesignedFor(0.1e6, 100.0e6));
  // Temperature is in the range, but Pressure is out of range
  EXPECT_FALSE(Obj3.IsDesignedFor(0.99e5, 100.0e6));
  EXPECT_FALSE(Obj3.IsDesignedFor(0.1e6, 101.5e6));
  // Both temperature and Pressure are within the range
  EXPECT_TRUE(Obj3.IsDesignedFor(0.1e6, 60.0e6));
}

TEST(ViscosityBehaviour, CarbonatedBrineViscosityBehaviour_GetViscosity) {
  // Arrange
  arma::vec::fixed<4> temp({298.15, 303.15, 308.15, 313.15});  // Kelvin
  const double pressure = 101325.0;  // in pascal
  const double salinity = 175328.31;  //  parts per million
  double co2 = 2.71;  // 0 to 4.8%

  vector<CarbonatedBrineViscosityBehaviour> viscosity_obj;
  viscosity_obj.reserve(4);
  vector<double> viscosity(4);

  for (unsigned int i=0; i < 4; i++) {
    viscosity_obj.push_back(CarbonatedBrineViscosityBehaviour(temp(i), salinity, co2));
    viscosity[i] = viscosity_obj[i].GetViscosity(pressure);
  }

  // Assert
  ASSERT_NEAR(viscosity[0], 1.3531311444229149E-03, 1.0e-17);  // Pa.s
  ASSERT_NEAR(viscosity[1], 1.2160833241210121E-03, 1.0e-17);
  ASSERT_NEAR(viscosity[2], 1.1003245989366663E-03, 1.0e-17);
  ASSERT_NEAR(viscosity[3], 1.0017069370884771E-03, 1.0e-17);
}

TEST(ViscosityBehaviour, CarbonatedBrineViscosityBehaviour_IsDesignedFor) {
  // This test checks the operation ranges of the viscosity modelling
  // 273K <= T <= 623K
  // 1Bar <= P <= 10^3 Bar
  // 0 mol/Kg < m <= 6 mol/Kg
  // Unit Convertion
  // 1Bar -> 100kPa

  arma::vec::fixed<3> temp({623.01, 272.99, 318.15});

  enum {out_low, in_low, in_upper, out_upper};

  arma::Mat<arma::uword> tests_pressure({{in_low,  in_upper},
                                         {in_low,  in_upper},
                                         {out_low, in_upper},
                                         {in_low,  out_upper},
                                         {in_low,  in_upper}});

  arma:: vec::fixed<4> pressures({0.99e5, 0.1e6, 100.0e6, 101.5e6});
  const double salinity = 175328.31;
  const double co2 = 2.71;  // 0 to 4.8%
  vector<bool> truth_table(3);
  vector<CarbonatedBrineViscosityBehaviour> viscosity_obj;
  viscosity_obj.reserve(3);

  for (unsigned int i=0; i < 3; i++) {
    viscosity_obj.push_back(CarbonatedBrineViscosityBehaviour(temp(i), salinity, co2));
    truth_table[i] = viscosity_obj[i].IsDesignedFor(pressures(tests_pressure(i, 0)),
                                                    pressures(tests_pressure(i, 1)));
    EXPECT_FALSE(truth_table[i]);
  }

  unsigned int i = 3;
  EXPECT_FALSE(viscosity_obj[2].IsDesignedFor(pressures(tests_pressure(i, 0)),
                                              pressures(tests_pressure(i, 1))));
  i++;

  EXPECT_TRUE(viscosity_obj[2].IsDesignedFor(pressures(tests_pressure(i, 0)), 40.0e6));
}

TEST(ViscosityBehaviour, SeawaterViscosityBehaviour_GetViscosity) {
  // Arrange
  arma::vec::fixed<4> temp({298.15, 303.15, 308.15, 313.15});  // Kelvin
  const double pressure = 101325.0;  // in pascal


  SeawaterViscosityBehaviour viscosity_1(temp(0));
  SeawaterViscosityBehaviour viscosity_2(temp(1));
  SeawaterViscosityBehaviour viscosity_3(temp(2));
  SeawaterViscosityBehaviour viscosity_4(temp(3));

  arma::vec::fixed<4> viscosity({viscosity_1.GetViscosity(pressure),
                                 viscosity_2.GetViscosity(pressure),
                                 viscosity_3.GetViscosity(pressure),
                                 viscosity_4.GetViscosity(pressure)});  // Pa.s

  // Assert
  ASSERT_NEAR(viscosity(0), 9.597086571567677E-04, 1.0e-6);  // Pa.s
  ASSERT_NEAR(viscosity(1), 8.624127049540729E-04, 1.0e-6);
  ASSERT_NEAR(viscosity(2), 7.797683601171808E-04, 1.0e-6);
  ASSERT_NEAR(viscosity(3), 7.091671441330667E-04, 1.0e-6);
}

TEST(ViscosityBehaviour, SeawaterViscosityBehaviour_IsDesignedFor) {
  // This test checks the operation ranges of the viscosity modelling
  // 20C <= T <= 105C
  // 1Bar <= P <= 600 Bar
  // Unit Convertion
  // 1Bar -> 100kPa
  // 0 C -> 273.15K

  // Arrange
  const double temp_max = 378.15;
  const double temp_min = 293.15;
  const double temp_within_limits = 308.15;
  const double p_max = 6.0e7;
  const double p_min = 1.0e5;

  SeawaterViscosityBehaviour Obj1(temp_max + 1.0);
  SeawaterViscosityBehaviour Obj2(temp_min - 1.0);
  SeawaterViscosityBehaviour Obj3(temp_within_limits);

  // Assert
  // Pressure is in the range, but Temperature is out of range
  EXPECT_FALSE(Obj1.IsDesignedFor(p_min, p_max));
  EXPECT_FALSE(Obj2.IsDesignedFor(p_min, p_max));
  // Temperature is in the range, but Pressure is out of range
  EXPECT_FALSE(Obj3.IsDesignedFor(p_min, p_max + 1.0));
  EXPECT_FALSE(Obj3.IsDesignedFor(p_min - 1.0, p_max));
  // Both temperature and Pressure are within the range
  EXPECT_TRUE(Obj3.IsDesignedFor(p_min, p_max));
}
#endif  // TEST_SRC_FLOW_SIMULATOR_VISCOSITY_BEHAVIOUR_TESTS_H_
