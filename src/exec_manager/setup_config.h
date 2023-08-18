/**
 * \file src/exec_manager/setup_config.h
 * \brief Contains the declaration of the \c SetupConfig class.
 *
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2019
 *
 * This header file contains the declaration of the \c SetupConfig class that stores the
 * configuration parameters related to config from the JSON file.
 */

#ifndef SRC_EXEC_MANAGER_SETUP_CONFIG_H_
#define SRC_EXEC_MANAGER_SETUP_CONFIG_H_

#include <armadillo>
#include <string>

/**
 * \class SetupConfig setup_config.h "src/exec_manager/setup_config.h"
 * \brief Provides a superclass that stores parameters from the "setup" top-level JSON object.
 */
class SetupConfig {
 public:
  /// Virtual destructor
  virtual ~SetupConfig(void) { }

  /**
   * \brief Spatial resolution of the tomographic image (in meters)
   *
   * This floating-point variable stores the size associated with each image voxel.
   */
  double voxel_size;

  /**
   * \brief Name of the folder in which input/output files are stored
   *
   * This string holds the name of the folder storing all input and output files.
   */
  std::string folder;

  /**
   * \brief Name of the input file to be loaded
   *
   * This string holds the name of the input file with the raw sample data to be loaded.
   */
  std::string input_file;

  /**
   * \brief Voxel dimensions (x,y,z) of the Digital Rock sample
   *
   * This integer array contains the values that define the sample dimensions (x,y,z).
   */
  arma::uvec::fixed<3> shape;
};

#endif  // SRC_EXEC_MANAGER_SETUP_CONFIG_H_
