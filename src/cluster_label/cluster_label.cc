/**
 * \file src/cluster_label/cluster_label.cc
 * \brief Implements cluster labelling algorithms.
 *
 * \author Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \author Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This source file implements cluster labelling algorithms. They take as input
 * an \c arma::ucube indicating the flagged sites and returns a fully labelled
 * \f$ l(\vec{r}) \f$ and an \c arma::ivec \f$ n(l) \f$ with the number of sites
 * for each cluster, in descending order.
 */

#include "src/cluster_label/cluster_label.h"
#include <algorithm>
#include <cassert>

void EnhancedHoshenKopelman(arma::ucube &l) {
/**
 * This function reads the \c arma::ucube with the flagged sites and performs
 * the Hoshen & Kopelman cluster labelling algorithm (Hoshen and Kopelman, 1976)
 * in place in order to identify the infinite cluster.
 * The "Enhanced" extension of the algorithm (Hoshen et al., 1997) simply adds
 * the calculation of more cluster properties to the algorithm, such as, the
 * cluster boundary, spatial moments, bounding box, radius of giration, squared
 * correlation length, fractal dimension, etc. to identity and better characterise
 * the infinite cluster.
 *
 * \param[in,out] l         Flagged sites as initialised by \c DigitalRock::DetermineMorphology().
 *
 * The Hoshen & Kopelman algorithm is a simple algorithm for labelling clusters
 * on a grid where each site may be either "occupied" or "unoccupied". It is an
 * efficient means of identifying clusters of contiguous sites in a single pass
 * through the grid. It can also be seen as a particular application of the
 * \c Union-Find algorithm, well known to computer scientists.
 *
 * In a single pass, it sequentially looks for flagged (occupied) sites in
 * \f$ l(\vec{r}) \f$ and assigns them a new cluster label, in case there is no
 * neighbouring labelled site, or a pre-existing cluster label, in case there is
 * a neighbouring flagged site. When the site in question connects two or more
 * clusters with different labels, the clusters are merged and the new site is
 * labelled with the proper label of the newly merged cluster. The merging
 * procedure consists of adding the merged number of sites to the proper label
 * and setting up links from the other labels to the proper with the aid of the
 * \c HKProper() function. Other quantities, that are related to the infinite
 * cluster, are also calculated in the first pass through the \c arma::ucube.
 *
 * The \c arma::ucube is swept a second time by the function \c HKCheck()
 * in order to check the correctness of the final \f$ l(\vec{r}) \f$, ensuring
 * that every flagged site has the same label as its flagged neighbours.
 *
 * Then, the infinite cluster(s) is (are) identified using the bounding boxes of
 * the clusters. The sites associated with finite (or "isolated") clusters then
 * have their flags removed, in order to discard them in the next calculations.
 */

  // Create arma:uvec with the required length
  arma::uvec nnl(13);               // Nearest neighbour labels in (i,j,k)
  arma::uvec uniq_non0;             // Unique list of non-zero nnl
  arma::uword largest_label = 0;    // Stores the largest label assigned so far

  // Create arma::ivec for the cluster sizes (or links to proper labels)
  arma::ivec n = arma::ivec(l.n_elem/2, arma::fill::zeros);

  // Allocate arma::umat for the bounding box of each cluster
  arma::umat bbox_min(l.n_elem/2, 3);  // (xmin, ymin, zmin)
  arma::umat bbox_max(l.n_elem/2, 3);  // (xmax, ymax, zmax)
  arma::urowvec::fixed<3> t;           // temporary (i, j, k) vector

  // Sweep arma::Cube optimizing memory access pattern
  for (arma::uword k = 0U; k != l.n_slices; ++k)  // For every slice
  for (arma::uword j = 0U; j != l.n_cols;   ++j)  // For every column
  for (arma::uword i = 0U; i != l.n_rows;   ++i)  // For every row
  if (l(i, j, k)) {                 // Check if the site is flagged (l != 0)
    // Look for labelled sites among nearest neighbours, respecting boundaries
    PreviousNeighbours(i, j, k, l, n, nnl);

    // Consolidate proper labels and links
    nnl.transform( [&n](arma::uword x) { return HKProper(x, n); } );

    // Eliminate non-flagged sites (label = 0) and duplicate labels
    uniq_non0 = arma::unique(nnl.elem(arma::find(nnl)));

    // Store temporary arma::uvec with current site coordinates
    t = {i, j, k};

    // Determine whether site belongs to a new or to one, two or three clusters
    switch (uniq_non0.n_elem) {
      case 1: {  // This site belongs to an existing cluster
        // Assign label to site and increment cluster size
        l(i, j, k) = uniq_non0(0);
        ++n(uniq_non0(0));

        // Update bounding box of the cluster
        bbox_min.row(l(i, j, k)) = arma::min(t, bbox_min.row(l(i, j, k)));
        bbox_max.row(l(i, j, k)) = arma::max(t, bbox_max.row(l(i, j, k)));
        break;
      }
      case 0: {  // This site belongs to a new cluster
        // Create new label and assign it to the site
        l(i, j, k) = ++largest_label;
        ++n(l(i, j, k));

        // Initialise the bounding box for the newly created cluster
        bbox_min.row(l(i, j, k)) = t;
        bbox_max.row(l(i, j, k)) = t;
        break;
      }
      default: {  // This site connects existing clusters that need to be merged
        // Assign lowest label to site, update cluster size and create links
        l(i, j, k) = uniq_non0(0);
        n(uniq_non0(0)) += 1 + arma::sum(n.elem(uniq_non0.tail(uniq_non0.n_elem-1)));
        n(uniq_non0.tail(uniq_non0.n_elem-1)).fill(static_cast<arma::sword>(-uniq_non0(0)));

        // Merge bounding boxes of linked clusters and erase obsolete ones
        bbox_min.row(l(i, j, k)) = arma::min(t,
                                   arma::min(bbox_min.rows(uniq_non0), 0));
        bbox_max.row(l(i, j, k)) = arma::max(t,
                                   arma::max(bbox_max.rows(uniq_non0), 0));
        bbox_min.rows(uniq_non0.tail(uniq_non0.n_elem-1)).zeros();
        bbox_max.rows(uniq_non0.tail(uniq_non0.n_elem-1)).zeros();
        break;
      }
    }
  }

  // Check for consistency within neighbouring labels
  HKCheck(n, l);

  // Search for infinite clusters in descending order of cluster size
  arma::uword infinite_clusters = 0;
  arma::urowvec one = {1, 1, 1};
  arma::uvec idx = arma::sort_index(n, "descend");
  for (auto x = idx.cbegin(); x != idx.cend(); ++x) {
    if (arma::prod(bbox_max.row(*x) - bbox_min.row(*x) + one) == l.n_elem) {
      ++infinite_clusters;
    } else {
      break;
    }
  }

  // Eliminate clusters that are not connected to the infinite clusters
  for (auto x = l.begin(); x != l.end(); ++x) {
    // Check if the proper label corresponds to any of the infinity clusters
    (*x) = arma::any(HKProper(*x, n) == idx.head(infinite_clusters)) ? 1 : 0;
  }

  // Print output messages
  std::printf("Amongst the %llu ", n(arma::find(n > 0)).eval().n_elem);
  std::printf("clusters found, only %llu of those ", infinite_clusters);
  std::printf("percolate(s) the system. ");
}  // end of EnhancedHoshenKopelman function



arma::uword HKProper(const arma::uword &label, arma::ivec &n) {
/**
 * This function takes a label as input and returns the associated proper label.
 *
 * \param[in]     label  Input label (may or not be a 'proper label')
 * \param[in,out] n      Number of sites (or links) assigned to each label.
 * \retval        proper Proper label associated with input label.
 */

  // Initialise static variables to hold the values from previous function call
  static arma::uword prev_label = -1U;
  static arma::uword prev_proper = -1U;

  // Take short-cut in case proper label is already known
  if (prev_label == label) return prev_proper;

  arma::uword proper = label;  // Starts as "input" and ends as "proper" label

  while (n(proper) < 0)   proper = static_cast<arma::uword>(-n(proper));  // Follow link to proper
  if    (proper != label) n(label) = -static_cast<arma::sword>(proper);   // Label points to proper

  // Store values to prepare for next iteration
  prev_label  = label;
  prev_proper = proper;

  return proper;
}  // end of HKProper function



void PreviousNeighbours(const arma::uword i, const arma::uword j, const arma::uword k,
                        const arma::ucube &l, arma::ivec &n, arma::uvec &nnl) {
/**
 * This function populates an \c arma::uvec vector with the proper labels of the
 * previously labelled neighbours of a site, up to third order (26-neighbourhood), given its
 * coordinates.
 *
 * \param[in]     i   Row index.
 * \param[in]     j   Column index.
 * \param[in]     k   Slice index.
 * \param[in]     l   Cube that stores labelled (previous) and flagged (future) sites.
 * \param[in,out] n   Number of sites (or links) assigned to each label.
 * \param[out]    nnl Previously labelled nearest neighbour labels around (i,j,k).
 */

  // Pre-calculate the maximum indexes for each dimension
  const arma::uword imax = l.n_rows - 1U;
  const arma::uword jmax = l.n_cols - 1U;

  // Consider only previously labelled neighbours up to third order
  nnl(0)  = HKProper(i && j && k           ? l(i-1, j-1, k-1) : 0, n);
  nnl(1)  = HKProper(j && k                ? l(i  , j-1, k-1) : 0, n);
  nnl(2)  = HKProper(i-imax && j && k      ? l(i+1, j-1, k-1) : 0, n);
  nnl(3)  = HKProper(i && k                ? l(i-1, j  , k-1) : 0, n);
  nnl(4)  = HKProper(k                     ? l(i  , j  , k-1) : 0, n);
  nnl(5)  = HKProper(i-imax && k           ? l(i+1, j  , k-1) : 0, n);
  nnl(6)  = HKProper(i && j-jmax && k      ? l(i-1, j+1, k-1) : 0, n);
  nnl(7)  = HKProper(j-jmax && k           ? l(i  , j+1, k-1) : 0, n);
  nnl(8)  = HKProper(i-imax && j-jmax && k ? l(i+1, j+1, k-1) : 0, n);
  nnl(9)  = HKProper(i && j                ? l(i-1, j-1, k)   : 0, n);
  nnl(10) = HKProper(j                     ? l(i  , j-1, k)   : 0, n);
  nnl(11) = HKProper(i-imax && j           ? l(i+1, j-1, k)   : 0, n);
  nnl(12) = HKProper(i                     ? l(i-1, j  , k)   : 0, n);
}  // end of PreviousThirdNeighbours function



void HKCheck(arma::ivec &n, arma::ucube &l) {
/**
 * This function simply goes over the \f$ l(\vec{r}) \f$ \c arma::ucube checking,
 * for each labelled site, whether its neighbours have the same label. It is a
 * simple consistency check.
 *
 * \param[in,out]    n  Number of sites (or links) assigned to each label.
 * \param[in]        l  Labelled sites according to \c HKRelabel().
 */
  // 26 nearest neighbour labels in (i,j,k)
  arma::uvec nnl(26);

  // Update each site with its proper label
  l.transform([&n](arma::uword label) { return HKProper(label, n); });

  // Sweep arma::ucube checking if labelled sites are equal to its neighbours
  for (arma::uword k = 0; k != l.n_slices; ++k)  // For every slice
  for (arma::uword j = 0; j != l.n_cols;   ++j)  // For every column
  for (arma::uword i = 0; i != l.n_rows;   ++i)  // For every row
    if (l(i, j, k)) {                            // Check if the site is labelled (l != 0)
      // Make sure all nearest neighbour labels are either -1U, 0U or l(i,j,k)
      AllNeighbours(i, j, k, l, nnl);
      assert(arma::all(nnl == -1U || nnl == 0U || nnl == l(i, j, k)));
    }
}  // end of HKCheck function



void AllNeighbours(const arma::uword i, const arma::uword j, const arma::uword k,
                   const arma::ucube &l, arma::uvec &nnl) {
/**
 * This function populates an \c arma::uvec vector with the labels of the
 * neighbouring sites, up to third order (26-neighbourhood).
 * When a site is located at one of the 6 external faces, the label of the non-existing neighbour
 * is taken as \c -1U.
 *
 * \param[in]   i   Row index.
 * \param[in]   j   Column index.
 * \param[in]   k   Slice index.
 * \param[in]   l   Cube that stores all labelled sites.
 * \param[out]  nnl Nearest neighbour labels around (i,j,k).
 */

  // Pre-calculate the maximum indexes for each dimension
  const arma::uword imax = l.n_rows - 1U;
  const arma::uword jmax = l.n_cols - 1U;
  const arma::uword kmax = l.n_slices - 1U;

  // Consider all neighbours up to third order
  nnl(0)  = (i && j && k                      ? l(i-1, j-1, k-1) : -1U);
  nnl(1)  = (j && k                           ? l(i  , j-1, k-1) : -1U);
  nnl(2)  = ((i-imax) && j && k               ? l(i+1, j-1, k-1) : -1U);
  nnl(3)  = (i && k                           ? l(i-1, j  , k-1) : -1U);
  nnl(4)  = (k                                ? l(i  , j  , k-1) : -1U);
  nnl(5)  = ((i-imax) && k                    ? l(i+1, j  , k-1) : -1U);
  nnl(6)  = (i && (j-jmax) && k               ? l(i-1, j+1, k-1) : -1U);
  nnl(7)  = ((j-jmax) && k                    ? l(i  , j+1, k-1) : -1U);
  nnl(8)  = ((i-imax) && (j-jmax) && k        ? l(i+1, j+1, k-1) : -1U);
  nnl(9)  = (i && j                           ? l(i-1, j-1, k)   : -1U);
  nnl(10) = (j                                ? l(i  , j-1, k)   : -1U);
  nnl(11) = ((i-imax) && j                    ? l(i+1, j-1, k)   : -1U);
  nnl(12) = (i                                ? l(i-1, j  , k)   : -1U);
  nnl(13) = ((i-imax)                         ? l(i+1, j  , k)   : -1U);
  nnl(14) = (i && (j-jmax)                    ? l(i-1, j+1, k)   : -1U);
  nnl(15) = ((j-jmax)                         ? l(i  , j+1, k)   : -1U);
  nnl(16) = ((i-imax) && (j-jmax)             ? l(i+1, j+1, k)   : -1U);
  nnl(17) = (i && j && (k-kmax)               ? l(i-1, j-1, k+1) : -1U);
  nnl(18) = (j && (k-kmax)                    ? l(i  , j-1, k+1) : -1U);
  nnl(19) = ((i-imax) && j && (k-kmax)        ? l(i+1, j-1, k+1) : -1U);
  nnl(20) = (i && (k-kmax)                    ? l(i-1, j  , k+1) : -1U);
  nnl(21) = ((k-kmax)                         ? l(i  , j  , k+1) : -1U);
  nnl(22) = ((i-imax) && (k-kmax)             ? l(i+1, j  , k+1) : -1U);
  nnl(23) = (i && (j-jmax) && (k-kmax)        ? l(i-1, j+1, k+1) : -1U);
  nnl(24) = ((j-jmax) && (k-kmax)             ? l(i  , j+1, k+1) : -1U);
  nnl(25) = ((i-imax) && (j-jmax) && (k-kmax) ? l(i+1, j+1, k+1) : -1U);
}  // end of AllThirdNeighbours function
