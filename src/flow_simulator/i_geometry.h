/**
 * \file src/flow_simulator/i_geometry.h
 * \brief Contains the \c IGeometry interface class.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This header file contains the \c IGeometry interface class from which all \c Geometry objects
 * derive.
 */

#ifndef SRC_FLOW_SIMULATOR_I_GEOMETRY_H_
#define SRC_FLOW_SIMULATOR_I_GEOMETRY_H_

#include <string>

namespace simulator {

/**
 * \class IGeometry i_geometry.h "src/flow_simulator/i_geometry.h"
 * \brief Interface class from which all \c Geometry objects derive.
 *
 * This class does nothing but define the virtual methods that all \c Geometry objects must have.
 */

class IGeometry {
 public:
  /// Virtual destructor
  virtual ~IGeometry() { }

  /**
   * \brief Name of the folder that contains the geometrical data
   *
   * This variable stores the name of the folder from which geometry data should be read.
   */
  std::string folder_;

 protected:
  /// Parametrised constructor
  explicit IGeometry(const std::string &folder) : folder_(folder) { }
};  // end of class IGeometry

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_I_GEOMETRY_H_
