/**
 * \file src/json_parser/json_parser.cc
 * \brief Contains the methods of the \c JSONParser class.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \author Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This source file contains the parse methods that handle each individual JSON file.s
 */

#include "src/json_parser/json_parser.h"
#include <rapidjson/schema.h>
#include <fstream>
#include <sstream>

rapidjson::Document JSONParser::ParseIntoJsonDocument(const std::string &json_file_name) const {
/**
 * This \c private method parses the content of a given JSON file into a \c rapidjson::Document.
 *
 * \param[in] json_file_name  The address of the JSON file with respect to project root.
 * \retval    json_doc        The parsed document containing data from the JSON file.
 */
  // Read JSON file
  const std::ifstream json_file{json_file_name, std::ifstream::in};
  if (!json_file.is_open()) throw std::ios::failure{"Impossible to read " + json_file_name};

  // Convert content to std::stringstream
  std::stringstream json_stream;
  json_stream << json_file.rdbuf();

  // Parse string stream into Document
  rapidjson::Document json_doc;
  json_doc.Parse(json_stream.str().c_str());

  // Return document if parsing was successful
  if (json_doc.HasParseError()) {
    std::fprintf(stderr, "JSONPARSER::PARSEINTOJSONDOCUMENT says: ");
    std::fprintf(stderr, "The '%s' file is not a valid JSON.\n", json_file_name.c_str());
    std::exit(-1);
  } else {
    return json_doc;
  }
}  // end of JSONParser::ParseIntoJsonDocument method



void JSONParser::ValidateJsonDocument(const rapidjson::Document &json_doc,
                                      const std::string &schema_file_name) const {
/**
 * This \c private method validates the a given \c rapidjson::Document using a JSON schema file.
 *
 * \param[in] json_doc          The parsed document to be validated.
 * \param[in] schema_file_name  The schema against which to validate.
 */
  // Parse schema file as generic JSON document
  rapidjson::Document schema_doc = ParseIntoJsonDocument(schema_file_name.c_str());

  // Create validator from schema JSON file
  rapidjson::SchemaDocument schema(schema_doc);
  rapidjson::SchemaValidator validator(schema);

  if (!json_doc.Accept(validator)) {
    std::fprintf(stderr, "JSONPARSER::VALIDATEJSONDOCUMENT says: ");
    std::fprintf(stderr, "The JSON file failed the validation test ");
    std::fprintf(stderr, "against '%s'.\n", schema_file_name.c_str());
    std::exit(-1);
  }
}  // end of JSONParser::ValidateJsonDocument method
