/**
 * \file src/cluster_label/cluster_label.h
 * \brief Function declaration for cluster labelling algorithms.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This header file contains the function declarations for cluster labelling
 * algorithms that identify clusters of flagged sites.
 */

#ifndef SRC_CLUSTER_LABEL_CLUSTER_LABEL_H_
#define SRC_CLUSTER_LABEL_CLUSTER_LABEL_H_

#include <armadillo>

/// Enhanced Hoshen & Kopelman cluster labelling algorithm
void EnhancedHoshenKopelman(arma::ucube &l);

/// HK proper label
arma::uword HKProper(const arma::uword &label, arma::ivec &n);

/// Previously labelled neighbours up to third order (26-neighbourhood)
void PreviousNeighbours(const arma::uword i, const arma::uword j, const arma::uword k,
                        const arma::ucube &l, arma::ivec &n, arma::uvec &nnl);

/// HK consistency check
void HKCheck(arma::ivec &n, arma::ucube &l);

/// All neighbouring labels up to third order (26-neighbourhood)
void AllNeighbours(const arma::uword i, const arma::uword j, const arma::uword k,
                   const arma::ucube &l, arma::uvec &nnl);

#endif  // SRC_CLUSTER_LABEL_CLUSTER_LABEL_H_
