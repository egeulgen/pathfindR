#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

using namespace Rcpp;

// Structure to track mutable state per seed without R environment overhead
struct SearchState {
  std::vector<bool> comp_members;
  std::vector<bool> removable_vec;
  std::vector<int> node2predecessor;
  std::vector<int> node2dep_count;
};

struct ExpandResult {
  bool improved;
  double best_score;
  int size;
  double zsum;
};

// 1. Fast C++ BFS for reachability
std::vector<bool> greedy_init_max_depth_cpp(const List& nbr_idx, int n_nodes, int start_idx, int depth) {
  std::vector<bool> within_vec(n_nodes, false);
  within_vec[start_idx] = true;
  if (depth == 0) return within_vec;

  std::vector<int> frontier;
  frontier.push_back(start_idx);
  std::vector<int> dist(n_nodes, 0);

  while (!frontier.empty()) {
    std::vector<int> nf;
    for (int cur : frontier) {
      int d_next = dist[cur] + 1;
      if (d_next <= depth) {
        IntegerVector nbrs = nbr_idx[cur];
        for (int i = 0; i < nbrs.size(); ++i) {
          int nb = nbrs[i] - 1; // R 1-based to C++ 0-based
          if (!within_vec[nb]) {
            within_vec[nb] = true;
            dist[nb] = d_next;
            nf.push_back(nb);
          }
        }
      }
    }
    frontier = std::move(nf);
  }
  return within_vec;
}

// 2. Fast C++ Recursive Expansion
ExpandResult greedy_expand_cpp(
    const List& nbr_idx, const NumericVector& z_vec,
    const NumericVector& sc_means, const NumericVector& sc_stds,
    const std::vector<bool>& within_vec, bool use_within,
    int search_depth, int depth, int size, double zsum,
    double score, double best_score, SearchState& state, int last_added
) {
  bool improved = false;
  if (score > best_score) {
    depth = search_depth;
    improved = true;
    best_score = score;
  }

  if (depth > 0) {
    bool any_improved = false;
    state.removable_vec[last_added] = false;
    int dep_count = 0;

    IntegerVector nbrs = nbr_idx[last_added];
    for (int i = 0; i < nbrs.size(); ++i) {
      int nb = nbrs[i] - 1;
      bool within_ok = !use_within || within_vec[nb];

      if (within_ok && !state.comp_members[nb]) {
        int new_size = size + 1;
        double new_zsum = zsum + z_vec[nb];
        // sc_means/sc_stds are 1-indexed from R, so new_size maps to new_size - 1
        double new_score = (new_zsum / std::sqrt(new_size) - sc_means[new_size - 1]) / sc_stds[new_size - 1];

        state.comp_members[nb] = true;
        state.removable_vec[nb] = true;

        ExpandResult res = greedy_expand_cpp(
          nbr_idx, z_vec, sc_means, sc_stds, within_vec, use_within,
          search_depth, depth - 1, new_size, new_zsum, new_score, best_score,
          state, nb
        );
        best_score = res.best_score;

        if (!res.improved) {
          state.comp_members[nb] = false;
          state.removable_vec[nb] = false;
        } else {
          size = res.size;
          zsum = res.zsum;
          dep_count++;
          any_improved = true;
          state.node2predecessor[nb] = last_added + 1; // Keep 1-based for R consistency if needed
        }
      }
    }

    improved = improved || any_improved;
    if (dep_count > 0) {
      state.removable_vec[last_added] = false;
      state.node2dep_count[last_added] = dep_count;
    }
  }

  return {improved, best_score, size, zsum};
}

// 3. Fast C++ Removal Pass
double greedy_removal_cpp(
    SearchState& state, const NumericVector& z_vec,
    const NumericVector& sc_means, const NumericVector& sc_stds,
    int& size, double& zsum, double best_score, int n_nodes
) {
  for (int cur = 0; cur < n_nodes; ++cur) {
    if (state.removable_vec[cur]) {
      int new_size = size - 1;
      double new_zsum = zsum - z_vec[cur];
      double new_score = (new_size <= 1) ? 0.0 :
        (new_zsum / std::sqrt(new_size) - sc_means[new_size - 1]) / sc_stds[new_size - 1];

      if (new_score > best_score) {
        best_score = new_score;
        size = new_size;
        zsum = new_zsum;
        state.comp_members[cur] = false;
        int pred = state.node2predecessor[cur] - 1;
        if (pred >= 0) {
          int dep = state.node2dep_count[pred] - 1;
          if (dep == 0) {
            state.removable_vec[pred] = true;
          } else {
            state.node2dep_count[pred] = dep;
          }
        }
      }
    }
  }
  return best_score;
}

// 4. Main Exported Driver Loop over all Seeds
// [[Rcpp::export]]
List run_greedy_search_cpp(
    List nbr_idx, NumericVector z_vec, NumericVector sc_means, NumericVector sc_stds,
    int max_depth, int search_depth, int n_nodes
) {
  List node2best(n_nodes);

  for (int seed_id = 0; seed_id < n_nodes; ++seed_id) {
    std::vector<bool> within_vec;
    bool use_within = (max_depth > 0);
    if (use_within) {
      within_vec = greedy_init_max_depth_cpp(nbr_idx, n_nodes, seed_id, max_depth);
    }

    SearchState state;
    state.comp_members.assign(n_nodes, false);
    state.removable_vec.assign(n_nodes, false);
    state.node2predecessor.assign(n_nodes, 0);
    state.node2dep_count.assign(n_nodes, 0);

    state.comp_members[seed_id] = true;
    state.node2dep_count[seed_id] = 1;

    double seed_zsum = z_vec[seed_id];

    ExpandResult res = greedy_expand_cpp(
      nbr_idx, z_vec, sc_means, sc_stds, within_vec, use_within,
      search_depth, search_depth, 1, seed_zsum, 0.0, -1e9, state, seed_id
    );

    double final_best = greedy_removal_cpp(
      state, z_vec, sc_means, sc_stds, res.size, res.zsum, res.best_score, n_nodes
    );

    // Collect final members
    std::vector<int> final_idx;
    for (int i = 0; i < n_nodes; ++i) {
      if (state.comp_members[i]) {
        final_idx.push_back(i + 1); // Return 1-based indices to R
      }
    }

    if (final_idx.size() >= 2 && final_best > 0) {
      for (int nd : final_idx) {
        int nd_idx = nd - 1;
        bool replace = false;
        if (node2best[nd_idx] == R_NilValue) {
          replace = true;
        } else {
          List ex = node2best[nd_idx];
          double ex_score = ex["score"];
          if (final_best > ex_score) replace = true;
        }

        if (replace) {
          node2best[nd_idx] = List::create(
            Named("idx") = wrap(final_idx),
            Named("score") = final_best
          );
        }
      }
    }
  }

  return node2best;
}
