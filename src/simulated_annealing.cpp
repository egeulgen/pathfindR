#include "java_compat.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>

// =============================================================================
// simulated_annealing.cpp
//
// Bit-for-bit C++ port of the Java SimulatedAnnealing (ActiveSubnetworkSearch).
// Java's SA is deterministic (a single seeded java.util.Random, single-threaded),
// so it can be reproduced exactly.
//
// Returns the final accepted "on" mask (length n_nodes, aligned to node order);
// the R side materialises the subnetwork objects once via .find_subnetworks().
// =============================================================================

// Sorted-descending component scores for the current on-mask (CSR adjacency).
static std::vector<double> sa_component_scores(
    const std::vector<char>& on, int N,
    const int* offsets, const int* nbrs,
    const double* z, const double* means, const double* stds,
    std::vector<char>& visited, std::vector<int>& stack) {

  std::fill(visited.begin(), visited.end(), 0);
  std::vector<double> scores;

  for (int s = 0; s < N; ++s) {
    if (!on[s] || visited[s]) continue;
    stack.clear();
    stack.push_back(s);
    visited[s] = 1;
    int size = 0;
    double zsum = 0.0;
    while (!stack.empty()) {
      int u = stack.back();
      stack.pop_back();
      ++size;
      zsum += z[u];
      for (int e = offsets[u]; e < offsets[u + 1]; ++e) {
        int v = nbrs[e];
        if (on[v] && !visited[v]) {
          visited[v] = 1;
          stack.push_back(v);
        }
      }
    }
    double sc = (size == 1) ? 0.0
    : (zsum / std::sqrt((double)size) - means[size - 1]) / stds[size - 1];
    scores.push_back(sc);
  }
  std::sort(scores.begin(), scores.end(), std::greater<double>());
  return scores;
}

// [[Rcpp::export]]
LogicalVector run_simulated_annealing(IntegerVector csr_offsets,
                                      IntegerVector csr_nbrs,
                                      NumericVector z,
                                      NumericVector means,
                                      NumericVector stds,
                                      int n_nodes,
                                      double gene_init_prob,
                                      bool start_with_all_positives,
                                      double sa_initial_temp,
                                      double sa_final_temp,
                                      int sa_iterations,
                                      int seed) {
  const int N = n_nodes;
  const int* offsets = INTEGER(csr_offsets);
  const int* nbrs    = INTEGER(csr_nbrs);
  const double* zp   = REAL(z);
  const double* mp   = REAL(means);
  const double* sp   = REAL(stds);

  JavaRandom rng(seed);

  std::vector<char> on(N, 0);

  // --- initial candidate solution -------------------------------------------
  if (start_with_all_positives) {
    // Java: add nodes with z > 0; no RNG consumed.
    for (int i = 0; i < N; ++i) on[i] = (zp[i] > 0) ? 1 : 0;
  } else {
    // Java: iterate nodesOffSet (== node order) drawing nextDouble per node.
    for (int i = 0; i < N; ++i) on[i] = (rng.nextDouble() < gene_init_prob) ? 1 : 0;
  }

  std::vector<char> visited(N, 0);
  std::vector<int> stack;
  stack.reserve(64);

  std::vector<double> cur = sa_component_scores(on, N, offsets, nbrs, zp, mp, sp, visited, stack);

  const double T0 = sa_initial_temp;
  const double Tf = sa_final_temp;
  const int total = sa_iterations;
  double T = T0;
  const double temp_step = 1.0 - std::pow(Tf / T0, 1.0 / total);

  for (int iter = 0; iter < total; ++iter) {
    int idx = rng.nextInt(N);
    on[idx] = !on[idx];

    std::vector<double> nw = sa_component_scores(on, N, offsets, nbrs, zp, mp, sp, visited, stack);

    bool keep = false, decision = false;
    int m = (int)std::min(cur.size(), nw.size());
    int k = 0;
    while (!decision && k < m) {
      double delta = nw[k] - cur[k];
      if (delta > 0.001) {
        keep = true;
        decision = true;
      } else if (rng.nextDouble() > std::exp(delta / T)) {
        keep = false;
        decision = true;
      }
      ++k;
    }

    if (keep) {
      cur.swap(nw);
    } else {
      on[idx] = !on[idx];   // revert the toggle
    }

    T *= (1.0 - temp_step);
  }

  LogicalVector out(N);
  for (int i = 0; i < N; ++i) out[i] = (bool)on[i];
  return out;
}
