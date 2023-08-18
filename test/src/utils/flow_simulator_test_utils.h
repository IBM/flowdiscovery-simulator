/**
 * \file test/src/utils/flow_simulator_test_utils.h
 * \brief Contains helper functions for the \c FlowSimulator tests.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 */

#ifndef TEST_SRC_UTILS_FLOW_SIMULATOR_TEST_UTILS_H_
#define TEST_SRC_UTILS_FLOW_SIMULATOR_TEST_UTILS_H_

#include <dirent.h>
#include <string>
#include <fstream>

namespace flow_simulator_test_utils {

bool CopyFile(const std::string &source_file_name, const std::string &destination_file_name) {
/**
 * This function copies a file from \c source_file_name to \c destination_file_name.
 *
 * \param[in] source_file_name        Source file name.
 * \param[in] destination_file_name   Destination file name.
 * \retval    success                 Boolean representing the success of the copy operation.
 */
  // Open file streams
  std::ifstream source(source_file_name, std::ios::binary);
  std::ofstream destination(destination_file_name, std::ios::binary);

  // Write to file stream
  destination << source.rdbuf();

  // Close file streams
  destination.close();
  source.close();

  return destination.good();
}

inline bool FileExists(const std::string &file_name) {
/**
 * This function checks if \c file_name exists.
 *
 * \param[in] file_name File name.
 * \retval    success   Boolean representing the success of the operation.
 */
  std::ifstream file(file_name.c_str());
  return file.good();
}

int CleanFolder(const std::string &folder_name) {
/**
 * This function cleans the contents of the \c folder_name folder.
 *
 * \param[in] folder_name Folder name.
 * \retval    success     Boolean representing the success of the operation.
 */
  // Erases the content of the folder from previous simulation
  DIR *theFolder = opendir(folder_name.c_str());
  struct dirent *next_file;
  std::string path;

  // Cycle through all files in the folder to delete them
  while ( (next_file = readdir(theFolder)) != NULL ) {
    // build the path for each file in the folder
    if (next_file->d_name[0] != '.') {
      path = folder_name + "/" + std::string(next_file->d_name);
      std::remove(path.c_str());
    }
  }
  int status = closedir(theFolder);
  return status;
}

}  // namespace flow_simulator_test_utils

#endif  // TEST_SRC_UTILS_FLOW_SIMULATOR_TEST_UTILS_H_
