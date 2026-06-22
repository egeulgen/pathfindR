#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>

using namespace Rcpp;

// =============================================================================
// subnetwork_scoring.cpp
//
// Hot-path helper shared by the simulated-annealing and genetic-algorithm
// searches.  Both algorithms evaluate a candidate solution (a set of "on"
// nodes) tens of thousands of times, and for each evaluation they need the
// score of every connected component induced by the on-nodes.
//
// =============================================================================

// CSR adjacency:
//   csr_offsets : length N+1, neighbours of node i are
//                 csr_nbrs[csr_offsets[i] .. csr_offsets[i+1]-1]
//   csr_nbrs    : flattened 0-based neighbour ids
//   on          : logical vector aligned to the node order (TRUE = node switched on)
//   z           : per-node z-scores aligned to the node order
//   means, stds : score-context calibration vectors (indexed by size, 0-based)
//
// Returns the component scores in descending order.
// [[Rcpp::export]]
NumericVector component_scores_sorted(LogicalVector on,
                                      IntegerVector csr_offsets,
                                      IntegerVector csr_nbrs,
                                      NumericVector z,
                                      NumericVector means,
                                      NumericVector stds) {
  int N = on.size();

  // Fast raw pointers (avoid Rcpp proxy overhead in the inner loop).
  const int* offsets = INTEGER(csr_offsets);
  const int* nbrs    = INTEGER(csr_nbrs);
  const double* zp   = REAL(z);
  const double* mp   = REAL(means);
  const double* sp   = REAL(stds);
  const int* onp     = LOGICAL(on);

  std::vector<char> visited(N, 0);
  std::vector<int> stack;
  stack.reserve(64);
  std::vector<double> scores;

  for (int s = 0; s < N; ++s) {
    if (!onp[s] || visited[s]) continue;

    // Iterative DFS over the on-nodes reachable from s.
    stack.clear();
    stack.push_back(s);
    visited[s] = 1;
    int size = 0;
    double zsum = 0.0;

    while (!stack.empty()) {
      int u = stack.back();
      stack.pop_back();
      ++size;
      zsum += zp[u];
      for (int e = offsets[u]; e < offsets[u + 1]; ++e) {
        int v = nbrs[e];
        if (onp[v] && !visited[v]) {
          visited[v] = 1;
          stack.push_back(v);
        }
      }
    }

    double sc;
    if (size == 1) {
      sc = 0.0;
    } else {
      sc = (zsum / std::sqrt((double)size) - mp[size - 1]) / sp[size - 1];
    }
    scores.push_back(sc);
  }

  std::sort(scores.begin(), scores.end(), std::greater<double>());
  return wrap(scores);
}

// -----------------------------------------------------------------------------
// Component membership labels for the on-nodes, for materialising the final /
// best solution without igraph.  Returns an integer vector of length N: 0 for
// off-nodes, otherwise a 1-based component id.  R can then split the node names
// by this label to reconstruct the subnetworks.
// [[Rcpp::export]]
IntegerVector component_labels(LogicalVector on,
                               IntegerVector csr_offsets,
                               IntegerVector csr_nbrs) {
  int N = on.size();
  const int* offsets = INTEGER(csr_offsets);
  const int* nbrs    = INTEGER(csr_nbrs);
  const int* onp     = LOGICAL(on);

  IntegerVector label(N, 0);
  std::vector<int> stack;
  stack.reserve(64);
  int comp = 0;

  for (int s = 0; s < N; ++s) {
    if (!onp[s] || label[s] != 0) continue;
    ++comp;
    stack.clear();
    stack.push_back(s);
    label[s] = comp;
    while (!stack.empty()) {
      int u = stack.back();
      stack.pop_back();
      for (int e = offsets[u]; e < offsets[u + 1]; ++e) {
        int v = nbrs[e];
        if (onp[v] && label[v] == 0) {
          label[v] = comp;
          stack.push_back(v);
        }
      }
    }
  }
  return label;
}
