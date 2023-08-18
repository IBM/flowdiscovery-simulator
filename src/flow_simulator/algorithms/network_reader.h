/**
 * \file src/flow_simulator/algorithms/network_reader.h
 * \brief Contains the definition of \c NetworkReader class methods.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 *
 * This header file contains the \c NetworkReader class which translates the centerlines file into
 * the information required by \c FlowSimulator.
 */

#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_NETWORK_READER_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_NETWORK_READER_H_

#include <rapidjson/document.h>
#include <armadillo>
#include <string>
#include "src/json_parser/json_parser.h"

namespace simulator {

/**
 * \struct NetworkInformation network_reader.h "src/flow_simulator/algorithms/network_reader.h"
 * \brief Struct that stores the capillary network definition.
 *
 * \note Read the documentation of \c StaticCapillaryNetworkContext and
 * \c DynamicCapillaryNetworkContext members for more information on how those variables are used
 * inside the Capillary Network algorithm.
 */
struct NetworkInformation {
  /// Parametrised struct constructor
  NetworkInformation(const arma::uword &num_nodes, const arma::uword &num_links)
    : number_of_nodes(num_nodes),
      number_of_links(num_links) {
        link_length.set_size(num_links);
        ctrl_voxels.set_size(num_nodes, 4);
        link_direction.set_size(2 * num_links);
        linked_nodes.set_size(2, 2 * num_links);
        link_squared_radius.set_size(num_links);
      }

  /// Number of nodes of the capillary network
  arma::uword number_of_nodes;

  /// Number of links of the capillary network
  arma::uword number_of_links;

  /// Length of links between neighbouring nodes (in voxels)
  arma::Col<double> link_length;

  /// Directions of links between neighbouring nodes
  arma::Col<double> link_direction;

  /// Location and squared-radius of each centreline voxel (in voxels)
  arma::Mat<double> ctrl_voxels;

  /// Pairs of indexes of nodes and links
  arma::Mat<arma::uword> linked_nodes;

  /// Squared-radius of links between neighbouring nodes (in voxels)
  arma::Col<double> link_squared_radius;
};



/**
 * \class NetworkReader network_reader.h "src/flow_simulator/algorithms/network_reader.h"
 * \brief Reads the Capillary Network file and extracts node and link information.
 *
 * The Network Reader parses the \c centerlines.json file into the Armadillo containers used
 * by the Capillary Network Algorithm.
 */

class NetworkReader : public JSONParser {
 public:
  /// Extract network information from JSON file
  NetworkInformation GetNetwork(const std::string &json_file_name) const;
};  // end of class NetworkReader

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_NETWORK_READER_H_
