iden_hub <- function(
  multi_layered_network_file,
  output_file,
  file_type = c("csv", "tsv"),
  top_n_hubs = NULL
) {
  file_type <- match.arg(file_type)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_file)) {
    dir.create(output_file, recursive = TRUE)
  }

  # Extract base name from input file for output specificity
  input_file_base_name <- tools::file_path_sans_ext(basename(multi_layered_network_file))
  # Clean the base name for use in filenames (e.g., replace non-alphanumeric with underscore)
  cleaned_input_file_name <- gsub("[^A-Za-z0-9_]", "", input_file_base_name)

  # 1. Load the multi-layered network file
  required_cols <- c("Feature1", "Feature2") # Only these are strictly needed for graph structure
  if (!file.exists(multi_layered_network_file)) {
    stop("Network file not found: '", multi_layered_network_file, "'")
  }

  network_data <- read_input_file(multi_layered_network_file, file_type = file_type, stringsAsFactors = FALSE)
  if (!all(required_cols %in% colnames(network_data))) {
    stop("Network file must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # 2. Create an undirected graph
  g <- igraph::graph_from_data_frame(d = network_data[, required_cols], directed = FALSE)

  # 3. Find all maximal cliques
  cliques <- tryCatch(
    igraph::max_cliques(g),
    error = function(e) {
      stop(paste("Error finding maximal cliques: ", e$message, ". This can be computationally intensive for large graphs.", sep = ""))
    }
  )

  # 4. Calculate MCC scores for each node based on the provided formula
  # Initialize MCC scores to 0 for all nodes
  mcc_scores <- setNames(numeric(igraph::vcount(g)), igraph::V(g)$name)

  # Sum (clique_size - 1)! for all maximal cliques containing the node
  for (clique_nodes_indices in cliques) {
    clique_size <- length(clique_nodes_indices)
    clique_score <- factorial(clique_size - 1) # factorial(0) = 1, factorial(1) = 1
    for (node_index in clique_nodes_indices) {
      node_name <- igraph::V(g)$name[node_index]
      mcc_scores[node_name] <- mcc_scores[node_name] + clique_score
    }
  }

  # Pre-calculate degrees and local clustering coefficients
  node_degrees <- igraph::degree(g)
  # transitivity(type="local") returns NaN for nodes with degree < 2.
  node_clustering_coeffs <- igraph::transitivity(g, type = "local", vids = igraph::V(g))
  names(node_clustering_coeffs) <- igraph::V(g)$name

  # Apply the special case: If there is no edge between the neighbors of the node v, MCC(v) = degree(v)
  for (node_name in igraph::V(g)$name) {
    current_degree <- node_degrees[node_name]

    # Determine if "no edge between neighbors" condition is met
    is_no_edge_between_neighbors <- FALSE
    if (current_degree == 0) {
      # Isolated node: no neighbors, so no edges between them. MCC should be 0 (its degree).
      is_no_edge_between_neighbors <- TRUE
    } else if (current_degree == 1) {
      # Node with one neighbor: trivially, no edges between multiple neighbors. MCC should be 1 (its degree).
      is_no_edge_between_neighbors <- TRUE
    } else { # current_degree >= 2
      # For nodes with 2 or more neighbors, check if local clustering coefficient is 0
      if (!is.na(node_clustering_coeffs[node_name]) && node_clustering_coeffs[node_name] == 0) {
        is_no_edge_between_neighbors <- TRUE
      }
    }

    # If the special condition is met, override the MCC score with the node's degree
    if (is_no_edge_between_neighbors) {
      mcc_scores[node_name] <- current_degree
    }
  }

  # 5. Create a data frame of results and rank
  hub_results_df <- dplyr::arrange(
    data.frame(
      Node = names(mcc_scores),
      MCC_score = mcc_scores,
      stringsAsFactors = FALSE
    ),
    dplyr::desc(MCC_score) # Sort in descending order of MCC_score
  )

  # 6. Apply top_n_hubs filter if specified
  if (!is.null(top_n_hubs) && is.numeric(top_n_hubs) && top_n_hubs > 0) {
    if (top_n_hubs > nrow(hub_results_df)) {
      warning("Requested top_n_hubs (", top_n_hubs, ") is greater than total nodes (", nrow(hub_results_df), "). Returning all nodes.")
    } else {
      hub_results_df <- head(hub_results_df, n = top_n_hubs)
    }
  }

  # 7. Save the results
  output_filename <- paste0("hub_mcc_", cleaned_input_file_name,
                            if(!is.null(top_n_hubs)) paste0("_top", top_n_hubs) else "", ".csv")
  output_filepath <- file.path(output_file, output_filename)
  write.csv(hub_results_df, output_filepath, row.names = FALSE)

  return(invisible(NULL))
}
