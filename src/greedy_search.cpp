#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <random>

using namespace Rcpp;

// Reusable allocation structures
struct SearchState {
  std::vector<bool> comp_members;
  std::vector<bool> removable_vec;
  std::vector<int> node2predecessor;
  std::vector<int> node2dep_count;
  std::vector<bool> within_vec;
  std::vector<int> dist;

  std::vector<int> frontier;
  std::vector<int> next_frontier;

  void reset(size_t n) {
    std::fill(comp_members.begin(), comp_members.end(), false);
    std::fill(removable_vec.begin(), removable_vec.end(), false);
    std::fill(node2predecessor.begin(), node2predecessor.end(), 0);
    std::fill(node2dep_count.begin(), node2dep_count.end(), 0);
    std::fill(within_vec.begin(), within_vec.end(), false);
    std::fill(dist.begin(), dist.end(), 0);
    frontier.clear();
    next_frontier.clear();
  }
};

struct ExpandResult {
  bool improved;
  double best_score;
  int size;
  double zsum;
};

struct Subnetwork {
  std::vector<int> idx;
  double score;

  bool operator<(const Subnetwork& other) const {
    if (std::abs(score - other.score) > 1e-9) {
      return score > other.score;
    }
    return idx.size() > other.idx.size();
  }
};

void greedy_init_max_depth(
    const int* flat_nbrs, const int* nbr_offsets,
    int start_idx, int depth, SearchState& state
) {
  state.within_vec[start_idx] = true;
  if (depth == 0) return;

  state.frontier.push_back(start_idx);

  while (!state.frontier.empty()) {
    state.next_frontier.clear();
    for (int cur : state.frontier) {
      int d_next = state.dist[cur] + 1;
      if (d_next <= depth) {
        int start_offset = nbr_offsets[cur];
        int end_offset = nbr_offsets[cur + 1];

        for (int i = start_offset; i < end_offset; ++i) {
          int nb = flat_nbrs[i];
          if (!state.within_vec[nb]) {
            state.within_vec[nb] = true;
            state.dist[nb] = d_next;
            state.next_frontier.push_back(nb);
          }
        }
      }
    }
    state.frontier.swap(state.next_frontier);
  }
}

ExpandResult greedy_expand(
    const int* flat_nbrs, const int* nbr_offsets,
    const double* z_vec, const double* sc_means, const double* sc_stds,
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

    int start_offset = nbr_offsets[last_added];
    int end_offset = nbr_offsets[last_added + 1];

    std::vector<int> neighbors;
    neighbors.reserve(end_offset - start_offset);
    for (int i = start_offset; i < end_offset; ++i) {
      int nb = flat_nbrs[i];

      if (use_within && !state.within_vec[nb]) continue;

      if (!state.comp_members[nb]) {
        int new_size = size + 1;
        double new_zsum = zsum + z_vec[nb];
        double new_score = (new_zsum / std::sqrt(new_size) - sc_means[new_size - 1]) / sc_stds[new_size - 1];

        state.comp_members[nb] = true;
        state.removable_vec[nb] = true;

        ExpandResult res = greedy_expand(
          flat_nbrs, nbr_offsets, z_vec, sc_means, sc_stds, use_within,
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

double greedy_removal(
    SearchState& state, const double* z_vec,
    const double* sc_means, const double* sc_stds,
    int& size, double& zsum, double best_score, int n_nodes
) {
  std::vector<int> removable_list;
  for (int cur = 0; cur < n_nodes; ++cur) {
    if (state.removable_vec[cur]) {
      removable_list.push_back(cur);
    }
  }

  for (int cur : removable_list) {
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
  return best_score;
}

// [[Rcpp::export]]
List run_greedy_search(
    List nbr_idx, NumericVector z_vec, NumericVector sc_means, NumericVector sc_stds,
    CharacterVector node_names, int max_depth, int search_depth, int n_nodes,
    double overlap_threshold, int subnetwork_num
) {
  // 1. Flatten the R List into contiguous standard vectors
  std::vector<int> flat_nbrs;
  std::vector<int> nbr_offsets(n_nodes + 1, 0);

  // Estimate size to reserve capacity up front
  int total_edges = 0;
  for (int i = 0; i < n_nodes; ++i) {
    IntegerVector nbrs = nbr_idx[i];
    total_edges += nbrs.size();
  }
  flat_nbrs.reserve(total_edges);

  for (int i = 0; i < n_nodes; ++i) {
    IntegerVector nbrs = nbr_idx[i];
    nbr_offsets[i] = flat_nbrs.size();
    int n_nbrs = nbrs.size();
    for (int j = 0; j < n_nbrs; ++j) {
      flat_nbrs.push_back(nbrs[j] - 1); // convert to 0-based index
    }
  }
  nbr_offsets[n_nodes] = flat_nbrs.size();

  // 2. Extract pointers directly to bypass Rcpp wrapper bounds checking inside loop
  const int* p_flat_nbrs = flat_nbrs.data();
  const int* p_nbr_offsets = nbr_offsets.data();
  const double* p_z_vec = REAL(z_vec);
  const double* p_sc_means = REAL(sc_means);
  const double* p_sc_stds = REAL(sc_stds);

  std::vector<Subnetwork> seed_best_subnetworks(n_nodes, {{}, -1e9});
  std::vector<bool> has_valid_subnetwork(n_nodes, false);

  SearchState state;
  state.comp_members.resize(n_nodes);
  state.removable_vec.resize(n_nodes);
  state.node2predecessor.resize(n_nodes);
  state.node2dep_count.resize(n_nodes);
  state.within_vec.resize(n_nodes);
  state.dist.resize(n_nodes);

  bool use_within = (max_depth > 0);

  for (int seed_id = 0; seed_id < n_nodes; ++seed_id) {
    state.reset(n_nodes);

    if (use_within) {
      greedy_init_max_depth(p_flat_nbrs, p_nbr_offsets, seed_id, max_depth, state);
    }

    state.comp_members[seed_id] = true;
    state.node2dep_count[seed_id] = 1;

    double seed_zsum = p_z_vec[seed_id];

    ExpandResult res = greedy_expand(
      p_flat_nbrs, p_nbr_offsets, p_z_vec, p_sc_means, p_sc_stds, use_within,
      search_depth, search_depth, 1, seed_zsum, 0.0, -1e9, state, seed_id
    );

    double final_best = greedy_removal(
      state, p_z_vec, p_sc_means, p_sc_stds, res.size, res.zsum, res.best_score, n_nodes
    );

    if (final_best > 0) {
      std::vector<int> final_idx;
      final_idx.reserve(res.size);
      for (int i = 0; i < n_nodes; ++i) {
        if (state.comp_members[i]) {
          final_idx.push_back(i + 1);
        }
      }

      if (final_idx.size() >= 2) {
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

  // 3. Unique filtering
  std::vector<Subnetwork> unique_candidates;
  unique_candidates.reserve(n_nodes);

  std::sort(seed_best_subnetworks.begin(), seed_best_subnetworks.end(), [](const Subnetwork& a, const Subnetwork& b){
    return a.idx < b.idx;
  });

  for (int i = 0; i < n_nodes; ++i) {
    if (seed_best_subnetworks[i].idx.empty()) continue;
    if (i == 0 || seed_best_subnetworks[i].idx != seed_best_subnetworks[i-1].idx) {
      unique_candidates.push_back(seed_best_subnetworks[i]);
    }
  }

  std::sort(unique_candidates.begin(), unique_candidates.end());

  // 4. Overlap filtering with inverted index.
  // For each kept subnetwork we record which nodes it owns. For each new
  // candidate we only compare against keepers that share at least one node
  // (found via the inverted index), skipping all non-overlapping pairs.
  // idx values are 1-based; node2kept is sized n_nodes+1 to accommodate them.
  std::vector<Subnetwork> filtered_networks;
  filtered_networks.reserve(subnetwork_num);

  std::vector<std::vector<int>> node2kept(n_nodes + 1);

  for (const auto& cand : unique_candidates) {
    if ((int)filtered_networks.size() >= subnetwork_num) break;

    // Collect the set of already-kept subnetworks that share any node
    std::vector<int> candidates_to_check;
    for (int ni : cand.idx) {
      for (int ki : node2kept[ni]) {
        candidates_to_check.push_back(ki);
      }
    }
    // Deduplicate
    std::sort(candidates_to_check.begin(), candidates_to_check.end());
    candidates_to_check.erase(
      std::unique(candidates_to_check.begin(), candidates_to_check.end()),
      candidates_to_check.end()
    );

    bool keep = true;
    for (int ki : candidates_to_check) {
      const auto& accepted = filtered_networks[ki];
      // Merge-based intersection on sorted idx vectors (both already sorted)
      int intersection_count = 0;
      size_t p1 = 0, p2 = 0;
      while (p1 < cand.idx.size() && p2 < accepted.idx.size()) {
        if      (cand.idx[p1] == accepted.idx[p2]) { intersection_count++; p1++; p2++; }
        else if (cand.idx[p1] <  accepted.idx[p2]) { p1++; }
        else                                        { p2++; }
      }

      int min_size = (int)std::min(cand.idx.size(), accepted.idx.size());
      if ((double)intersection_count / min_size > overlap_threshold) {
        keep = false;
        break;
      }
    }

    if (keep) {
      int kept_idx = (int)filtered_networks.size();
      filtered_networks.push_back(cand);
      for (int ni : cand.idx) node2kept[ni].push_back(kept_idx);
    }
  }

  // 5. Construct final R List output
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
