#include "java_compat.h"

// =============================================================================
// java_compat.cpp
//
// Rcpp-exported wrappers over the Java collection-ordering helpers defined in
// java_compat.h.  Keeping these exports in their own translation unit means:
//
//   * Rcpp::compileAttributes() generates exactly one forward declaration in
//     RcppExports.cpp that resolves to exactly one definition here — no ODR
//     issues regardless of how many other .cpp files #include java_compat.h.
//
//   * Any future algorithm file (SA, GA, …) that needs java_node_order or
//     java_neighbour_order can call them as ordinary R-visible functions
//     without duplicating the implementation.
// =============================================================================

// -----------------------------------------------------------------------------
// Reconstruct Java's networkNodeList order from SIF edges.
//
// Java's SIFReader inserts nodes into adjacency (a HashMap) in order of first
// appearance scanning the file top-to-bottom, node1 before node2 per line;
// self-loops are skipped.  networkNodeList = new ArrayList<>(adjacency.keySet())
// then gives the HashMap key iteration order.
//
// `src` and `tgt` must be the raw, upper-cased, in-file-order endpoint columns
// (NOT de-duplicated or reordered by the caller).
//
// Returns the node names in Java networkNodeList order.
// [[Rcpp::export]]
CharacterVector java_node_order(CharacterVector src, CharacterVector tgt) {
  int m = src.size();

  std::vector<std::string> insertion;
  std::vector<int>         insertion_hash;
  std::unordered_map<std::string, char> present;
  present.reserve(m * 2);

  auto add = [&](const std::string& nm) {
    if (present.find(nm) == present.end()) {
      present.emplace(nm, 1);
      insertion.push_back(nm);
      insertion_hash.push_back(java_string_hashcode(nm));
    }
  };

  for (int i = 0; i < m; ++i) {
    std::string a = as<std::string>(src[i]);
    std::string b = as<std::string>(tgt[i]);
    if (a == b) continue;   // self-loop — discarded as in Java
    add(a);
    add(b);
  }

  std::vector<int> order = java_iteration_order(insertion_hash);
  CharacterVector out(order.size());
  for (std::size_t i = 0; i < order.size(); ++i) out[i] = insertion[order[i]];
  return out;
}

// -----------------------------------------------------------------------------
// Reorder one node's neighbour list into Java HashSet iteration order.
//
// `nbr_insertion_ids`   — 1-based node ids in the order edges were first seen
//                         in the SIF for this node
// `nbr_insertion_names` — corresponding node names (same order)
//
// Returns the ids permuted into Java HashSet iteration order, matching the
// order that Java's greedy expansion visits neighbours.
// [[Rcpp::export]]
IntegerVector java_neighbour_order(IntegerVector   nbr_insertion_ids,
                                   CharacterVector nbr_insertion_names) {
  int k = nbr_insertion_ids.size();
  std::vector<int> hashes(k);
  for (int i = 0; i < k; ++i)
    hashes[i] = java_string_hashcode(as<std::string>(nbr_insertion_names[i]));

  std::vector<int> order = java_iteration_order(hashes);
  IntegerVector out(k);
  for (int i = 0; i < k; ++i) out[i] = nbr_insertion_ids[order[i]];
  return out;
}
