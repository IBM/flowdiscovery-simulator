/**
 * \file src/json_parser/json_parser.h
 * \brief Contains the declaration of the \c JSONParser class.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \author Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This header file contains the declaration of the \c JSONParser class that loads a JSON file and
 * validates it against a JSON schema.
 */

#ifndef SRC_JSON_PARSER_JSON_PARSER_H_
#define SRC_JSON_PARSER_JSON_PARSER_H_

#include <rapidjson/document.h>
#include <string>

/**
 * \class JSONParser json_parser.h "src/json_parser/json_parser.h"
 * \brief Provides an interface to the JSON parsing library.
 *
 * All library-specific details of the JSON parsing procedure should be hidden
 * inside this class that, in turn, is responsible for exposing the JSON configuration parameters to
 * the outside world.
 *
 * \note The current implementation uses the RapidJSON library \<http://rapidjson.org/\> and performs
 *       validation using a JSON schema.
 */

class JSONParser {
 protected:
  /// Method that parses a generic JSON file and returns a \c rapidjson::Document
  rapidjson::Document ParseIntoJsonDocument(const std::string &json_file_name) const;

  /// Method that validates a parsed \c rapidjson::Document against a JSON schema
  void ValidateJsonDocument(const rapidjson::Document &json_doc,
                            const std::string &schema_file_name) const;
};  // end of class JSONParser

#endif  // SRC_JSON_PARSER_JSON_PARSER_H_
