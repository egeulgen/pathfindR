# =============================================================================
# Genetic-algorithm active subnetwork search. Each individual is a logical
# genome over the network nodes (TRUE = node switched on). Fitness is the set
# of connected subnetworks induced by the on-nodes, compared rank-by-rank by
# score. The population evolves through rank-based selection, uniform
# crossover, optional mutation, periodic random injection and elitism.
# =============================================================================

#' Build a GA individual from a logical genome
#'
#' Computes the connected subnetworks induced by the on-nodes and stores their
#' scores (sorted descending) for fast fitness comparison. The scores are
#' obtained with the C++ component scorer over the precomputed CSR adjacency,
#' which is the only thing the fitness comparison (\code{.ga_compare} /
#' \code{.ga_sort_desc}) ever reads. The node membership of a genome is not
#' materialized here — it is reconstructed once, for the best individual, at the
#' end of \code{.genetic_algorithm()}.
#'
#' @param rep_logical A logical vector aligned to \code{network$nodes}
#'   (\code{TRUE} = node on). May be empty.
#' @param csr_offsets,csr_nbrs CSR adjacency from \code{build_network()}.
#' @param z_vec Per-node z-scores aligned to the node order.
#' @param means,stds Score-context calibration vectors.
#'
#' @return A list with elements \code{rep} (the genome) and \code{scores} (the
#'   descending component scores).
.ga_make_individual <- function(rep_logical, csr_offsets, csr_nbrs, z_vec, means, stds) {
  scores <- component_scores_sorted(rep_logical, csr_offsets, csr_nbrs, z_vec, means, stds)
  list(rep = rep_logical, scores = scores)
}

#' Fitness of an individual
#'
#' The score of its highest scoring subnetwork, or \code{0} when it has none.
#'
#' @param ind A GA individual from \code{.ga_make_individual()}.
#'
#' @return A numeric score.
.ga_score <- function(ind) {
  if (length(ind$scores) == 0L) 0 else ind$scores[1L]
}

#' Compare two individuals rank-by-rank
#'
#' Compares the score-sorted subnetwork scores position by position; the first
#' strict difference decides. If one score vector is a prefix of the other, the
#' individual with more subnetworks wins.
#'
#' @param a,b GA individuals.
#'
#' @return \code{1} if \code{a} is better, \code{-1} if \code{b} is better,
#'   \code{0} if they are equivalent.
.ga_compare <- function(a, b) {
  sa <- a$scores
  sb <- b$scores
  m <- min(length(sa), length(sb))
  k <- 1L
  while (k <= m) {
    if (sa[k] > sb[k]) {
      return(1L)
    }
    if (sa[k] < sb[k]) {
      return(-1L)
    }
    k <- k + 1L
  }
  if (length(sa) > length(sb)) {
    return(1L)
  }
  if (length(sb) > length(sa)) {
    return(-1L)
  }
  0L
}

#' Sort a population from best to worst
#'
#' Implements the same ordering as \code{.ga_compare()} in a vectorised way:
#' score vectors are laid out in a matrix padded with \code{-Inf}, then ordered
#' lexicographically (descending), with the number of subnetworks as a final
#' tie-break so that more subnetworks ranks higher.
#'
#' @param pop A list of GA individuals.
#'
#' @return The population reordered, best individual first.
.ga_sort_desc <- function(pop) {
  n <- length(pop)
  if (n <= 1L) {
    return(pop)
  }

  lens <- vapply(pop, function(x) length(x$scores), integer(1))
  maxlen <- max(lens)
  if (maxlen == 0L) {
    return(pop)
  }

  mat <- matrix(-Inf, nrow = n, ncol = maxlen)
  for (i in seq_len(n)) {
    s <- pop[[i]]$scores
    if (length(s) > 0L) mat[i, seq_along(s)] <- s
  }

  keys <- c(lapply(seq_len(maxlen), function(j) mat[, j]), list(lens))
  ord <- do.call(order, c(keys, list(decreasing = TRUE)))
  pop[ord]
}

#' Pick one index by rank-proportional (roulette) selection
#'
#' @param weights Numeric weights aligned to the (sorted) population.
#' @param total Sum of \code{weights}.
#'
#' @return An integer index into the population.
.ga_pick <- function(weights, total) {
  rand <- stats::runif(1) * total
  np <- length(weights)
  idx <- -1L
  rr <- 1L
  while (rr <= np && idx == -1L) {
    rand <- rand - weights[rr]
    if (rand <= 0) idx <- rr
    rr <- rr + 1L
  }
  if (idx == -1L) idx <- np
  idx
}

#' Crossover and mutation of two parent genomes
#'
#' With probability \code{ga_crossover_rate} the parents are recombined
#' (uniform crossover by default); otherwise the children are empty genomes.
#' Mutation flips each bit with probability \code{ga_mutation_rate}.
#'
#' @param p1,p2 Parent genomes (logical vectors).
#' @param params A list of run parameters.
#' @param crossover_type One of \code{"UNIFORM"}, \code{"SINGLEPOINT"} or
#'   \code{"MULTIPOINT"}.
#'
#' @return A list of two child genomes (logical vectors).
.ga_crossover_mutation <- function(p1, p2, params, crossover_type = "UNIFORM") {
  n <- length(p1)
  c1 <- logical(0)
  c2 <- logical(0)

  if (stats::runif(1) < params$ga_crossover_rate) {
    if (identical(crossover_type, "SINGLEPOINT")) {
      cp <- sample.int(n, 1L) - 1L # crossover point in [0, n-1]
      if (cp >= 1L) {
        idx <- seq_len(cp)
        c1 <- c(c1, p1[idx])
        c2 <- c(c2, p2[idx])
      }
      if (cp < n) {
        idx <- (cp + 1L):n
        c1 <- c(c1, p2[idx])
        c2 <- c(c2, p1[idx])
      }
    } else if (identical(crossover_type, "MULTIPOINT")) {
      c1 <- logical(n)
      c2 <- logical(n)
      flag <- 0L
      count <- 0L
      for (i in seq_len(n)) {
        if (flag == 0L) {
          c1[i] <- p1[i]
          c2[i] <- p2[i]
        } else {
          c1[i] <- p2[i]
          c2[i] <- p1[i]
        }
        count <- count + 1L
        if (count == 10L) {
          count <- 0L
          flag <- (flag + 1L) %% 2L
        }
      }
    } else { # UNIFORM
      from_p1 <- stats::runif(n) < 0.5
      c1 <- ifelse(from_p1, p1, p2)
      c2 <- ifelse(from_p1, p2, p1)
    }
  }

  if (params$ga_mutation_rate > 0) {
    if (length(c1) > 0L) {
      flip <- stats::runif(length(c1)) < params$ga_mutation_rate
      c1[flip] <- !c1[flip]
    }
    if (length(c2) > 0L) {
      flip <- stats::runif(length(c2)) < params$ga_mutation_rate
      c2[flip] <- !c2[flip]
    }
  }

  list(c1, c2)
}

#' Run the genetic-algorithm active subnetwork search
#'
#' The initial population is seeded randomly (each gene switched on with
#' probability \code{params$gene_init_prob}); when
#' \code{params$start_with_all_positives} is set, one individual containing all
#' positive-z-score nodes is added. The population then evolves by rank-based
#' selection, uniform crossover and optional mutation. Every ten generations
#' the worst 10\% of the population is replaced with fresh random individuals,
#' and the previous best individual is preserved. The search stops after
#' \code{params$ga_iterations} generations or
#' once the best individual is unchanged for 50 generations.
#'
#' @param network A network from \code{build_network()}.
#' @param sc A score context from \code{build_score_context()}.
#' @param params A list of run parameters.
#' @param verbose Logical; emit progress messages.
#'
#' @return A list of subnetwork objects from the best individual found.
.genetic_algorithm <- function(network, sc, params, verbose = FALSE) {
  nodes <- network$nodes
  N <- length(nodes)
  if (N == 0L) {
    return(list())
  }

  set.seed(params$seed)

  # Precomputed inputs for the C++ component scorer (built once, reused for
  # every individual evaluation).
  z_vec <- as.numeric(sc$z[nodes])
  means <- as.numeric(sc$means)
  stds <- as.numeric(sc$stds)
  offsets <- network$csr_offsets
  nbrs <- network$csr_nbrs

  pop_size <- params$ga_population_size
  n_replace <- as.integer(pop_size * 0.1)

  # --- initial population -------------------------------------------------
  population <- vector("list", pop_size)
  for (i in seq_len(pop_size)) {
    rep_i <- stats::runif(N) < params$gene_init_prob
    population[[i]] <- .ga_make_individual(rep_i, offsets, nbrs, z_vec, means, stds)
  }
  if (isTRUE(params$start_with_all_positives)) {
    rep_pos <- as.numeric(sc$z[nodes]) > 0
    population[[pop_size]] <- .ga_make_individual(rep_pos, offsets, nbrs, z_vec, means, stds)
  }
  population <- .ga_sort_desc(population)

  weights <- pop_size - (seq_len(pop_size) - 1L) # popSize, popSize-1, ..., 1
  total_weight <- sum(weights)

  add_random_count <- 0L
  iter <- 0L
  last_best <- population[[1L]]
  last_best_repeat <- 0L
  running <- TRUE

  while (running) {
    # --- breed a new population ----------------------------------------
    new_pop <- vector("list", 0L)
    while (length(new_pop) < pop_size) {
      i1 <- .ga_pick(weights, total_weight)
      i2 <- .ga_pick(weights, total_weight)
      kids <- .ga_crossover_mutation(
        population[[i1]]$rep,
        population[[i2]]$rep,
        params, "UNIFORM"
      )
      new_pop[[length(new_pop) + 1L]] <- .ga_make_individual(kids[[1L]], offsets, nbrs, z_vec, means, stds)
      new_pop[[length(new_pop) + 1L]] <- .ga_make_individual(kids[[2L]], offsets, nbrs, z_vec, means, stds)
    }
    if (length(new_pop) > pop_size) new_pop <- new_pop[seq_len(pop_size)]
    new_pop <- .ga_sort_desc(new_pop)

    # --- every 10 generations replace the worst 10% with random ones ---
    if (add_random_count == 10L) {
      if (n_replace >= 1L) {
        for (i in seq_len(n_replace)) {
          rep_r <- stats::runif(N) < params$gene_init_prob
          new_pop[[pop_size - i + 1L]] <- .ga_make_individual(rep_r, offsets, nbrs, z_vec, means, stds)
        }
      }
      add_random_count <- 0L
    }

    # --- elitism: never let the best score regress ---------------------
    if (.ga_compare(new_pop[[1L]], population[[1L]]) < 0L) {
      new_pop[[pop_size]] <- population[[1L]]
    }
    new_pop <- .ga_sort_desc(new_pop)
    population <- new_pop

    add_random_count <- add_random_count + 1L
    iter <- iter + 1L

    if (.ga_compare(population[[1L]], last_best) == 0L) {
      last_best_repeat <- last_best_repeat + 1L
    } else {
      last_best <- population[[1L]]
      last_best_repeat <- 0L
    }

    if (verbose) {
      message(
        "GA generation ", iter,
        " best score=", signif(.ga_score(population[[1L]]), 4)
      )
    }

    if (last_best_repeat >= 50L) {
      running <- FALSE
      if (verbose) message("The score did not improve in 50 generations")
    }
    if (iter >= params$ga_iterations) running <- FALSE
  }

  # Materialise the best individual's subnetworks once (the search itself only
  # ever needed the scores held on each individual).
  best_rep <- population[[1L]]$rep
  .sort_subnetworks_desc(.find_subnetworks(network, sc, nodes[which(best_rep)]))
}
