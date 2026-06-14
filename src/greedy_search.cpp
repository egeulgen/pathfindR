#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <set>

using namespace Rcpp;

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

// Helper struct to sort and filter final modules
struct Subnetwork {
  std::vector<int> idx; // 1-based indices for R consistency
  double score;

  // Sort descending by score
  bool operator<(const Subnetwork& other) const {
    if (std::abs(score - other.score) > 1e-9) {
      return score > other.score;
    }
    return idx.size() > other.idx.size();
  }
};

// 1. Fast C++ BFS
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
          int nb = nbrs[i] - 1;
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

// 3. Removal Pass
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

// [[Rcpp::export]]
List run_greedy_search(
    List nbr_idx, NumericVector z_vec, NumericVector sc_means, NumericVector sc_stds,
    CharacterVector node_names, int max_depth, int search_depth, int n_nodes,
    double overlap_threshold, int subnetwork_num
) {
  // 1. Core Map Tracker
  std::vector<Subnetwork> seed_best_subnetworks(n_nodes, {{}, -1e9});
  std::vector<bool> has_valid_subnetwork(n_nodes, false);

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

    std::vector<int> final_idx;
    for (int i = 0; i < n_nodes; ++i) {
      if (state.comp_members[i]) {
        final_idx.push_back(i + 1); // 1-based indexing for R
      }
    }

    if (final_idx.size() >= 2 && final_best > 0) {
      // Update mapping for all members belonging to this discovery
      for (int nd : final_idx) {
        int nd_idx = nd - 1;
        if (!has_valid_subnetwork[nd_idx] || final_best > seed_best_subnetworks[nd_idx].score) {
          seed_best_subnetworks[nd_idx] = {final_idx, final_best};
          has_valid_subnetwork[nd_idx] = true;
        }
      }
    }
  }

  // 2. De-duplication using std::set
  std::vector<Subnetwork> unique_candidates;
  std::set<std::vector<int>> seen_sets;

  for (int i = 0; i < n_nodes; ++i) {
    if (has_valid_subnetwork[i]) {
      if (seen_sets.find(seed_best_subnetworks[i].idx) == seen_sets.end()) {
        seen_sets.insert(seed_best_subnetworks[i].idx);
        unique_candidates.push_back(seed_best_subnetworks[i]);
      }
    }
  }

  // 3. Sort candidates descending by score
  std::sort(unique_candidates.begin(), unique_candidates.end());

  // 4. Fast Jaccard Overlap Filtering
  std::vector<Subnetwork> filtered_networks;

  for (const auto& cand : unique_candidates) {
    if ((int)filtered_networks.size() >= subnetwork_num) break;

    bool keep = true;
    // Convert current subnetwork to a fast-lookup set
    std::set<int> cand_set(cand.idx.begin(), cand.idx.end());

    for (const auto& accepted : filtered_networks) {
      int intersection_count = 0;
      for (int id : accepted.idx) {
        if (cand_set.count(id)) {
          intersection_count++;
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

  // 5. Build final R output list structure
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
