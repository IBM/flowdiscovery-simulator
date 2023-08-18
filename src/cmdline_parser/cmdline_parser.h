/**
 * \file src/cmdline_parser/cmdline_parser.h
 * \brief Contains the class that holds the command-line arguments.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This header file contains the declaration of the \c CmdLineParser class that
 * handles the command-line arguments.
 */

#ifndef SRC_CMDLINE_PARSER_CMDLINE_PARSER_H_
#define SRC_CMDLINE_PARSER_CMDLINE_PARSER_H_

#include <tclap/CmdLine.h>
#include <string>

/**
 * \class CmdLineParser cmdline_parser.h "src/cmdline_parser/cmdline_parser.h"
 * \brief Provides an interface to the command-line argument parsing library.
 *
 * All library-specific details of the argument parsing procedure
 * should be hidden inside this class that, in turn, is responsible for exposing
 * the command-line options to the outside world.
 *
 * \note The current implementation uses the TCLAP library
 *       <http://tclap.sourceforge.net>.
 */

class CmdLineParser {
 public:
  /// Method that parses the command-line options
  void Parse(int argc, const char **argv);

  /// Name of the JSON file to be parsed
  std::string json_file_name_;

  /// Boolean variable that triggers the execution of the setup step
  bool run_setup_;

  /// Boolean variable that triggers the execution of the segmentation step
  bool run_segmentation_;

  /// Boolean variable that triggers the execution of the morphology step
  bool run_morphology_;

  /// Boolean variable that triggers the execution of the simulation step
  bool run_simulation_;
};  // end of class CmdLineParser

#endif  // SRC_CMDLINE_PARSER_CMDLINE_PARSER_H_
