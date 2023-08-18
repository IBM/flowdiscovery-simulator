/**
 * \file src/main.cc
 * \brief Run the flow-simulator.x executable.
 *
 * \author Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \author Alexandre Ashade Lassance Cunha \<aashade@br.ibm.com\>
 * \author Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This source file contains the \c main() function that reads the execution mode from the
 * command-line and calls the appropriate methods.
 */

#include "src/main.h"
#include <glog/logging.h>
#include "src/exec_manager/exec_manager.h"
#include "src/cmdline_parser/cmdline_parser.h"

/**
 * The \c main() function reads the configuration parameters from JSON files,
 * loads the desired input file, applies the required algorithms and writes
 * output data files.
 */
int main(int argc, const char **argv) {
  // Initialise Google's logging library
  ::google::InitGoogleLogging(argv[0]);

  // Instantiate and populate CmdLineParser object
  CmdLineParser cmd;
  cmd.Parse(argc, argv);

  // Run the desired execution mode
  ExecutionManager exec_manager;
  if (cmd.run_simulation_) exec_manager.RunSimulation(cmd.json_file_name_);
  if (!cmd.run_simulation_) {
    std::fprintf(stderr, "FLOWDISCOVERY says: ");
    std::fprintf(stderr, "Please define at least one execution mode. See '--help' flag.\n");
    std::exit(-1);
  }
  return 0;
}  // end of main()
