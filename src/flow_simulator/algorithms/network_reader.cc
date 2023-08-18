/**
 * \file src/flow_simulator/algorithms/network_reader.cc
 * \brief Contains the implementation of \c NetworkReader class methods.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2018
 *
 * This source file contains the implementation of \c NetworkReader class methods which are
 * responsible for reading the network information from the \c centerlines.json file and
 * returning the Armadillo containers used by \c FlowSimulator.
 */

#include "src/flow_simulator/algorithms/network_reader.h"

namespace simulator {

NetworkInformation NetworkReader::GetNetwork(const std::string &json_file_name) const {
/**
 * The \c GetNetwork() method parses the \c centerlines.json file and returns a \c struct with the
 * information.
 *
 * \param[in] json_file_name        Name of the JSON file to be parsed.
 * \retval    network_information   Struct containing the data in the JSON in Armadillo format.
 */
  // Parse file into rapidjson::Document object and validate JSON file with respect to JSON schema
  rapidjson::Document json_doc = ParseIntoJsonDocument(json_file_name);
  std::string schema_file_name = "util/json_graph_schema.json";
  ValidateJsonDocument(json_doc, schema_file_name);

  // Get root "graph" object
  assert(json_doc.HasMember("graph"));
  auto graph_obj = json_doc["graph"].GetObject();

  // Get "metadata" object
  assert(graph_obj.HasMember("metadata"));
  auto graph_metadata_obj = graph_obj["metadata"].GetObject();
  arma::uword number_of_nodes = graph_metadata_obj["number_of_nodes"].GetUint64();
  arma::uword number_of_links = graph_metadata_obj["number_of_links"].GetUint64();

  // Initialise return struct
  NetworkInformation network_information(number_of_nodes, number_of_links);

  // Get "nodes" array
  assert(graph_obj.HasMember("nodes"));
  auto nodes_array = graph_obj["nodes"].GetArray();
  for (const auto &node_obj : nodes_array) {
    // Get metadata objects
    auto node_metadata_obj = node_obj["metadata"].GetObject();
    auto node_coordinates_obj = node_metadata_obj["node_coordinates"].GetObject();

    // Get node information
    arma::uword node_id = std::stoull(node_obj["id"].GetString());
    double node_squared_radius = node_metadata_obj["node_squared_radius"].GetDouble();
    double node_x = node_coordinates_obj["x"].GetDouble();
    double node_y = node_coordinates_obj["y"].GetDouble();
    double node_z = node_coordinates_obj["z"].GetDouble();

    // Create temporary vectors
    arma::rowvec ctrl_voxels_row = {node_x, node_y, node_z, node_squared_radius};

    // Populate struct members
    network_information.ctrl_voxels.row(node_id) = ctrl_voxels_row;
  }

  // Get "edges" array
  assert(graph_obj.HasMember("edges"));
  auto edges_array = graph_obj["edges"].GetArray();
  for (const auto &edge_obj : edges_array) {
    // Get metadata objects
    auto edge_metadata_obj = edge_obj["metadata"].GetObject();

    // Get link information
    arma::uword link_id = std::stoull(edge_obj["id"].GetString());
    arma::uword source = std::stoull(edge_obj["source"].GetString());
    arma::uword target = std::stoull(edge_obj["target"].GetString());
    double link_length = edge_metadata_obj["link_length"].GetDouble();
    double link_squared_radius = edge_metadata_obj["link_squared_radius"].GetDouble();

    // Create temporary vectors
    arma::vec link_direction_rows = {+1.0, -1.0};
    arma::uvec linked_nodes_source_col = {source, link_id};
    arma::uvec linked_nodes_target_col = {target, link_id};

    // Populate struct members
    network_information.link_length(link_id) = link_length;
    network_information.link_squared_radius(link_id) = link_squared_radius;
    network_information.linked_nodes.col(2U * link_id) = linked_nodes_source_col;
    network_information.linked_nodes.col(2U * link_id + 1U) = linked_nodes_target_col;
    network_information.link_direction.rows(2U * link_id, 2U * link_id + 1U) = link_direction_rows;
  }

  return network_information;
}  // end of NetworkReader::GetNetwork() method

}  // namespace simulator
