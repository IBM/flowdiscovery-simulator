/**
 * \file util/tif2raw.cc
 * \brief Converts a set of TIF images (2D slices) to a single binary (3D cube).
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2014
 *
 * This source file contains a function that converts a set of TIF images with
 * 2D slices of a rock sample into a single 3D binary file with all the data.
 * The ROCKETH group (ETH - Zurich) provides several tomography samples using
 * this format at http://www.rockphysics.ethz.ch/downloads.
 */

#include <armadillo>

/**
 * The \c main() function reads the list of input files, loads each one into
 * a slice of an \c arma::Cube and outputs the whole \c arma::Cube to single
 * binary (raw) file.
 */
int main(int argc, const char **argv) {
  arma::uword size = argc-1;
  arma::Mat<uint16_t> M;
  arma::Cube<uint8_t> C(size, size, size);

  for (int i = 1; i != argc; ++i) {
    M.load(argv[i], arma::raw_binary);
    M.reshape(size, size);
    C.slice(i-1) = arma::conv_to< arma::Mat<uint8_t> >::from(M / 257.0);
  }

  C.save("geometry.raw", arma::raw_binary);

  return 0;
}  // end of main()
