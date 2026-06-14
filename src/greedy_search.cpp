#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

using namespace Rcpp;

// Reusable allocation structures to completely eliminate heap allocation during loops
struct SearchState {
  std::vector<bool> comp_members;
  std::vector<bool> removable_vec;
  std::vector<int> node2predecessor;
  std::vector<int> node2dep_count;
  std::vector<bool> within_vec;
  std::vector<int> dist;

  void reset(size_t n) {
    std::fill(comp_members.begin(), comp_members.end(), false);
    std::fill(removable_vec.begin(), removable_vec.end(), false);
    std::fill(node2predecessor.begin(), node2predecessor.end(), 0);
    std::fill(node2dep_count.begin(), node2dep_count.end(), 0);
    std::fill(within_vec.begin(), within_vec.end(), false);
    std::fill(dist.begin(), dist.end(), 0);
  }
};

struct ExpandResult {
  bool improved;
  double best_score;
  int size;
  double zsum;
};

struct Subnetwork {
  std::vector<int> idx; // Sorted 1-based indices
  double score;

  // Fast descending sort
  bool operator<(const Subnetwork& other) const {
    if (std::abs(score - other.score) > 1e-9) {
      return score > other.score;
    }
    return idx.size() > other.idx.size();
  }
};

// Fast BFS that modifies state in-place with zero allocations
void greedy_init_max_depth_fast(const List& nbr_idx, int start_idx, int depth, SearchState& state) {
  state.within_vec[start_idx] = true;
  if (depth == 0) return;

  // Reusable frontier buffer to avoid reallocation
  static std::vector<int> frontier;
  static std::vector<int> next_frontier;
  frontier.clear();
  frontier.push_back(start_idx);

  while (!frontier.empty()) {
    next_frontier.clear();
    for (int cur : frontier) {
      int d_next = state.dist[cur] + 1;
      if (d_next <= depth) {
        IntegerVector nbrs = nbr_idx[cur];
        int n_nbrs = nbrs.size();
        for (int i = 0; i < n_nbrs; ++i) {
          int nb = nbrs[i] - 1;
          if (!state.within_vec[nb]) {
            state.within_vec[nb] = true;
            state.dist[nb] = d_next;
            next_frontier.push_back(nb);
          }
        }
      }
    }
    frontier.swap(next_frontier);
  }
}

ExpandResult greedy_expand_fast(
    const List& nbr_idx, const NumericVector& z_vec,
    const NumericVector& sc_means, const NumericVector& sc_stds,
    bool use_within, int search_depth, int depth, int size, double zsum,
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
    int n_nbrs = nbrs.size();
    for (int i = 0; i < n_nbrs; ++i) {
      int nb = nbrs[i] - 1;
      if (use_within && !state.within_vec[nb]) continue;

      if (!state.comp_members[nb]) {
        int new_size = size + 1;
        double new_zsum = zsum + z_vec[nb];
        double new_score = (new_zsum / std::sqrt(new_size) - sc_means[new_size - 1]) / sc_stds[new_size - 1];

        state.comp_members[nb] = true;
        state.removable_vec[nb] = true;

        ExpandResult res = greedy_expand_fast(
          nbr_idx, z_vec, sc_means, sc_stds, use_within,
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
          state.node2predecessor[nb] = last_added + 1;
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

// Linear time Removal Pass
double greedy_removal_fast(
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

// [[Rcpp::export]]
List run_greedy_search(
    List nbr_idx, NumericVector z_vec, NumericVector sc_means, NumericVector sc_stds,
    CharacterVector node_names, int max_depth, int search_depth, int n_nodes,
    double overlap_threshold, int subnetwork_num
) {
  std::vector<Subnetwork> seed_best_subnetworks(n_nodes, {{}, -1e9});
  std::vector<bool> has_valid_subnetwork(n_nodes, false);

  // Instantiate memory structures once outside the loop
  SearchState state;
  state.comp_members.resize(n_nodes);
  state.removable_vec.resize(n_nodes);
  state.node2predecessor.resize(n_nodes);
  state.node2dep_count.resize(n_nodes);
  state.within_vec.resize(n_nodes);
  state.dist.resize(n_nodes);

  bool use_within = (max_depth > 0);

  for (int seed_id = 0; seed_id < n_nodes; ++seed_id) {
    state.reset(n_nodes); // Linear memset style speed instead of reallocation

    if (use_within) {
      greedy_init_max_depth_fast(nbr_idx, seed_id, max_depth, state);
    }

    state.comp_members[seed_id] = true;
    state.node2dep_count[seed_id] = 1;

    double seed_zsum = z_vec[seed_id];

    ExpandResult res = greedy_expand_fast(
      nbr_idx, z_vec, sc_means, sc_stds, use_within,
      search_depth, search_depth, 1, seed_zsum, 0.0, -1e9, state, seed_id
    );

    double final_best = greedy_removal_fast(
      state, z_vec, sc_means, sc_stds, res.size, res.zsum, res.best_score, n_nodes
    );

    if (final_best > 0) {
      std::vector<int> final_idx;
      // Pre-reserve to avoid reallocations
      final_idx.reserve(res.size);
      for (int i = 0; i < n_nodes; ++i) {
        if (state.comp_members[i]) {
          final_idx.push_back(i + 1); // 1-based for R
        }
      }

      if (final_idx.size() >= 2) {
        // Ensure indices are naturally sorted (they are because loop goes 0 to n_nodes)
        for (int nd : final_idx) {
          int nd_idx = nd - 1;
          if (!has_valid_subnetwork[nd_idx] || final_best > seed_best_subnetworks[nd_idx].score) {
            seed_best_subnetworks[nd_idx] = {final_idx, final_best};
            has_valid_subnetwork[nd_idx] = true;
          }
        }
      }
    }
  }

  // Optimization: Unique checks via raw index vector comparison (Skipping std::set overhead)
  std::vector<Subnetwork> unique_candidates;
  unique_candidates.reserve(n_nodes);

  // Sort components internally to make duplication checking trivial
  std::sort(seed_best_subnetworks.begin(), seed_best_subnetworks.end(), [](const Subnetwork& a, const Subnetwork& b){
    return a.idx < b.idx;
  });

  for (int i = 0; i < n_nodes; ++i) {
    if (seed_best_subnetworks[i].idx.empty()) continue;
    if (i == 0 || seed_best_subnetworks[i].idx != seed_best_subnetworks[i-1].idx) {
      unique_candidates.push_back(seed_best_subnetworks[i]);
    }
  }

  // Sort unique modules descending by score
  std::sort(unique_candidates.begin(), unique_candidates.end());

  // Two-pointer Jaccard intersect calculation
  std::vector<Subnetwork> filtered_networks;
  filtered_networks.reserve(subnetwork_num);

  for (const auto& cand : unique_candidates) {
    if ((int)filtered_networks.size() >= subnetwork_num) break;

    bool keep = true;
    for (const auto& accepted : filtered_networks) {
      int intersection_count = 0;
      size_t p1 = 0, p2 = 0;

      // Linear scan because indices are sorted
      while (p1 < cand.idx.size() && p2 < accepted.idx.size()) {
        if (cand.idx[p1] == accepted.idx[p2]) {
          intersection_count++;
          p1++; p2++;
        } else if (cand.idx[p1] < accepted.idx[p2]) {
          p1++;
        } else {
          p2++;
        }
      }

      int union_count = cand.idx.size() + accepted.idx.size() - intersection_count;
      double jaccard_overlap = (double)intersection_count / union_count;

      if (jaccard_overlap > overlap_threshold) {
        keep = false;
        break;
      }
    }

    if (keep) {
      filtered_networks.push_back(cand);
    }
  }

  // Export final allocations back to R
  List result_list(filtered_networks.size());
  for (size_t i = 0; i < filtered_networks.size(); ++i) {
    CharacterVector out_nodes(filtered_networks[i].idx.size());
    for (size_t j = 0; j < filtered_networks[i].idx.size(); ++j) {
      out_nodes[j] = node_names[filtered_networks[i].idx[j] - 1];
    }

    result_list[i] = List::create(
      Named("nodes") = out_nodes,
      Named("score") = filtered_networks[i].score
    );
  }

  return result_list;
}
