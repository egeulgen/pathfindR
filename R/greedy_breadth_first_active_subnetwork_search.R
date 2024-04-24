#' Calculate Background Score
#'
#' @param pin_graph \code{\link[igraph]{igraph}} object containing the protein-protein
#' interaction network
#' @param number_of_iterations the number of iterations for background score calculations
#' @param pin_scores_vec vector of score for each gene in the PIN
#'
#' @return a list containing 2 elements \describe{
#'   \item{sampling_score_means}{vector of sampling score means}
#'   \item{sampling_score_stds}{vector of sampling score standard deviations}
#' }
calculate_background_score <- function(pin_graph, number_of_iterations, pin_scores_vec) {

    sampling_score_sums <- rep(0, igraph::vcount(pin_graph))
    sampling_score_square_sums <- rep(0, igraph::vcount(pin_graph))

    node_list_for_sampling <- igraph::V(pin_graph)

    for (iter in seq_len(number_of_iterations)) {

        # progress message
        if (iter %% 50 == 0) {
            message(round(iter / number_of_iterations, 2) * 100, "% of total iterations in background score calculation finished")
        }

        # shuffle nodes
        node_list_for_sampling <- sample(node_list_for_sampling)

        z_sum <- 0
        number_of_nodes_in_subnetwork <- 0

        for (node in node_list_for_sampling) {

            z_sum <- z_sum + pin_scores_vec[[node]]  #Much faster than using V(pin_graph)$score[node]

            number_of_nodes_in_subnetwork <- number_of_nodes_in_subnetwork + 1

            score <- z_sum/sqrt(number_of_nodes_in_subnetwork)

            sampling_score_sums[number_of_nodes_in_subnetwork] <- sampling_score_sums[number_of_nodes_in_subnetwork] + score
            sampling_score_square_sums[number_of_nodes_in_subnetwork] <- sampling_score_square_sums[number_of_nodes_in_subnetwork] + score ^ 2
        }
    }

    # These are vector operations
    sampling_score_means <- sampling_score_sums / number_of_iterations

    sampling_score_stds <- sampling_score_square_sums/number_of_iterations - sampling_score_means ^ 2
    sampling_score_stds <- sqrt(sampling_score_stds + 1e-10)

    # Addition of small number has two purposes: 1.Prevents division by zero in calculate_component_score() 2.Corrects for very small negative numbers that might appear,
    # like -3.388132e-21

    ## Explanation of the operation above var = SUM((x-xmean)^2) / N var = SUM(x^2 - 2*xmean*x + xmean^2)/N var = SUM(x^2)/N - (2*xmean*SUM(x))/N + (N*xmean^2)/N var =
    ## SUM(x^2)/N - 2*xmean^2 + xmean^2 var = SUM(x^2)/N - xmean^2
    return(list(sampling_score_means = sampling_score_means, sampling_score_stds = sampling_score_stds))
}


#' Calculate Component Score
#'
#' @inheritParams calculate_background_score
#' @param component component (vector)
#' @param sampling_score_means vector of background sampling score means
#' @param sampling_score_stds vector of background sampling score standard deviations
#'
#' @return the score of the component
calculate_component_score <- function(pin_scores_vec, component, sampling_score_means, sampling_score_stds) {
    numberOfNodes <- length(component)

    if (numberOfNodes == 0)
        return(0)

    score <- sum(pin_scores_vec[component]) / sqrt(numberOfNodes)
    score <- (score - sampling_score_means[numberOfNodes]) / sampling_score_stds[numberOfNodes]
    return(score)
}


#' Greedy Breadth-First Active Subnetwork search method implementation in R
#'
#' @param seed_node the seed node (an \code{\link[igraph]{igraph}} vs object)
#' @inheritParams calculate_background_score
#' @param input_scores_df data frame containing scores per each gene in the PIN
#' @param background_sampling_result list containing vector of background sampling score means
#' and vector of background sampling score standard deviations
#' @param max_depth maximum depth in search
#' @param check_second_neighbors boolean to indicate whether to check second
#' neighbors when any direct neighbor does not improve the score
#'
#' @return an active module (vector)
greedy_breadth_first_active_subnetwork_search <- function(seed_node, pin_graph, input_scores_df, pin_scores_vec, background_sampling_result, max_depth, check_second_neighbors) {

    sampling_score_means <- background_sampling_result$sampling_score_means
    sampling_score_stds <- background_sampling_result$sampling_score_stds

    comp <- c()

    queue <- c(seed_node)
    checked_in_greedy <- c(seed_node)
    distances_from_seed <- c(0)

    will_be_checked_for_neighbors <- c()
    # When check_second_neighbors==TRUE, a node, the nodes which do not increase the score by themselves are added to the end of the queue and they are checked when
    # there are only this kind of nodes.  The reason is, an important neighbor might already be added by the help of another important node and we may not need this
    # unimportant node

    while (length(queue) > 0) {
        if (setequal(queue, will_be_checked_for_neighbors)) {
            # all nodes in the queue are nodes whose neighbors will be checked for addition
            node <- queue[1]  # get next node
            queue <- queue[-1]  # remove the node

            node_distance <- distances_from_seed[1]
            distances_from_seed <- distances_from_seed[-1]

            neighbor_names <- names(igraph::neighbors(pin_graph, node))
            neighbor_scores_df <- input_scores_df[input_scores_df$Gene %in% neighbor_names, ]
            neighbor_names <- neighbor_scores_df$Gene[order(neighbor_scores_df$Score, decreasing = TRUE)]

            current_score <- calculate_component_score(pin_graph = pin_graph, pin_scores_vec = pin_scores_vec, compo = comp, sampling_score_means = sampling_score_means, sampling_score_stds = sampling_score_stds)
            comp <- c(comp, node)
            will_be_checked_for_neighbors <- will_be_checked_for_neighbors[will_be_checked_for_neighbors != node]  #remove from postponed list

            neighbor_added <- FALSE
            for (neighbor_name in neighbor_names) {
                if (pin_scores_vec[neighbor_name] > 0) {
                  neighbor_node <- igraph::V(pin_graph)[neighbor_name]  #getting igraph.vs from gene name
                  if (!neighbor_node %in% checked_in_greedy) {
                    checked_in_greedy <- c(checked_in_greedy, neighbor_node)
                    new_score <- calculate_component_score(pin_graph = pin_graph, pin_scores_vec = pin_scores_vec, compo = c(comp, neighbor_node), sampling_score_means = sampling_score_means,
                      sampling_score_stds = sampling_score_stds)
                    if (new_score > current_score) {
                      queue <- c(queue, neighbor_node)
                      distances_from_seed <- c(distances_from_seed, node_distance + 1)
                      neighbor_added <- TRUE
                    }
                  }
                }
            }

            if (!neighbor_added)
                comp <- comp[-length(comp)]  # Removing node from comp

        } else {
            node <- queue[1]  # get next node
            queue <- queue[-1]  # remove the node

            node_distance <- distances_from_seed[1]
            distances_from_seed <- distances_from_seed[-1]

            if (node %in% will_be_checked_for_neighbors) {
                # Sending to the end of the queue, will be checked when there are only this kind of nodes
                queue <- c(queue, neighbor_node)
                distances_from_seed <- c(distances_from_seed, node_distance)
            } else {
                current_score <- calculate_component_score(pin_graph = pin_graph, pin_scores_vec = pin_scores_vec, compo = comp, sampling_score_means = sampling_score_means, sampling_score_stds = sampling_score_stds)
                new_score <- calculate_component_score(pin_graph = pin_graph, pin_scores_vec = pin_scores_vec, compo = c(comp, node), sampling_score_means = sampling_score_means, sampling_score_stds = sampling_score_stds)


                if (new_score > current_score) {

                  comp <- c(comp, node)

                  if (node_distance < max_depth) {
                    # Its distance is less than max_depth, which means we can go further and check its neighbors
                    neighbor_names <- names(igraph::neighbors(pin_graph, node))
                    neighbor_scores_df <- input_scores_df[input_scores_df$Gene %in% neighbor_names, ]
                    neighbor_names <- neighbor_scores_df$Gene[order(neighbor_scores_df$Score, decreasing = TRUE)]

                    for (neighbor_name in neighbor_names) {
                      if (pin_scores_vec[neighbor_name] > 0)
                        {
                          neighbor_node <- igraph::V(pin_graph)[neighbor_name]  # getting igraph.vs from gene name
                          if (!neighbor_node %in% checked_in_greedy) {
                            checked_in_greedy <- c(checked_in_greedy, neighbor_node)
                            queue <- c(queue, neighbor_node)
                            distances_from_seed <- c(distances_from_seed, node_distance + 1)
                          }
                        }  #else break
                    }
                  }

                } else {
                  if (check_second_neighbors) {
                    if (node_distance < max_depth) {
                      # if we will be able to add neighbors of this node Postponing
                      will_be_checked_for_neighbors <- c(will_be_checked_for_neighbors, node)
                      queue <- c(queue, node)
                      distances_from_seed <- c(distances_from_seed, node_distance)
                    }
                  }
                }
            }
        }
    }

    return(comp)
}


#' Preprocessing and data generation for active subnetwork search functions
#'
#' @inheritParams active_snw_search
#'
#' @return list of processed data to be used in active subnetwork search functions
active_snw_search_preprocessing <- function(input_for_search, pin_name_path) {
  # read PIN
  pin_df <- utils::read.delim(pin_path, header = FALSE)
  pin_df <- subset(pin_df, select = -2)

  pin <- igraph::graph_from_data_frame(pin_df, directed = FALSE, vertices = NULL)

  # read scores file
  scores_df <- input_for_search
  colnames(scores_df) <- c("Gene", "Score")

  # qnorm that is used for z-score conversion gives Inf for less than (1 - 5e-17)
  # we are applying a threshold at 1E-15 and converting any smaller value to 1e-15
  scores_df$Score <- vapply(scores_df$Score, function(x) ifelse(x < 1e-15, 1e-15, x), 0.01)

  # p-values are converted to z-scores
  scores_df$Score <- stats::qnorm(1 - scores_df$Score)

  # add rows with a score of 0 for missing genes (PIN genes not in input)
  scores_df_to_append <- data.frame(Gene = setdiff(igraph::V(pin)$name,
                                                   scores_df$Gene),
                                    Score = 0)
  scores_df <- rbind(scores_df, scores_df_to_append)

  # assign scores to the PIN graph
  idx <- match(igraph::V(pin)$name, scores_df$Gene)
  igraph::V(pin)$score <- scores_df$Score[idx]

  # create a vector of scores (for faster execution)
  scores_vec <- scores_df$Score[idx]
  names(scores_vec) <- scores_df$Gene[idx]

  return(list(pin = pin, scores_df = scores_df, scores_vec = scores_vec))
}



#' Greed Active Subnetwork Search wrapper for pathfindR
#'
#' @inheritParams active_snw_search
#' @param max_depth Sets max depth in greedy search, 0 for no limit (default = 1) #TODO Move to `active_snw_search` doc
#' @param check_second_neighbors in greedy search, whether to check second
#' neighbors when any direct neighbor does not improve the score (default = TRUE) #TODO Default TRUE?  #TODO Move to `active_snw_search` doc
#'
#' @return list containing active subnetworks
#' @export
greedy_active_subnetwork_search <- function(input_for_search, pin_name_path, max_depth, check_second_neighbors) {

  res <- active_snw_search_preprocessing(input_for_search = input_for_search, pin_name_path = pin_name_path)

  input_scores_df <- res$scores_df
  pin_graph <- res$pin
  pin_scores_vec <- res$scores_vec
  sig_genes <- input_for_search$GENE

  background_sampling_result <- calculate_background_score(
    pin_graph = pin_graph,
    number_of_iterations = 1000,
    scores_vec = pin_scores_vec
  )

  active_modules <- list()
  num_of_seeds_used <- 0
  ratio_print_threshold <- 10
  for (seed in sig_genes) {
    # progress
    num_of_seeds_used <- num_of_seeds_used + 1
    used_seed_ratio <- round(num_of_seeds_used / length(sig_genes), 2) * 100
    if (used_seed_ratio >= ratio_print_threshold) {
      message(used_seed_ratio, "% of seeds are checked")
      ratio_print_threshold <- ratio_print_threshold + 10
    }

    seed_node <- igraph::V(pin_graph)[seed] # getting igraph.vs from gene name
    comp <- greedy_breadth_first_active_subnetwork_search(seed_node = seed_node,
                                                          pin_graph = pin_graph,
                                                          input_scores_df = input_scores_df,
                                                          pin_scores_vec = pin_scores_vec,
                                                          background_sampling_result = background_sampling_result,
                                                          max_depth = max_depth,
                                                          check_second_neighbors = check_second_neighbors)
    # check if the same module exists
    same_exists <- FALSE
    for (active_module in active_modules){
      if(setequal(comp, active_module)){
        same_exists <- TRUE
        break
      }
    }

    if(!same_exists)
      active_modules[[length(active_modules) + 1]] <- names(comp)
  }
  return(active_modules)
}
