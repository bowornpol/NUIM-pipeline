find_path <- function(
  multi_layered_network_file,
  source_node,
  target_node,
  output_file,
  file_type = c("csv", "tsv")
) {
  file_type <- match.arg(file_type)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_file)) {
    dir.create(output_file, recursive = TRUE)
  }

  # 1. Load network
  required_cols <- c("Feature1", "Feature2", "edge_score", "edge_type")
  if (!file.exists(multi_layered_network_file)) {
    stop("Network file not found: '", multi_layered_network_file, "'")
  }

  network_data <- read_input_file(multi_layered_network_file, file_type = file_type, stringsAsFactors = FALSE)
  if (!all(required_cols %in% colnames(network_data))) {
    stop("Network file must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # Ensure edge_score is numeric
  network_data$edge_score <- as.numeric(network_data$edge_score)

  # --- Convert edge_score to absolute value before weight transformation ---
  network_data$edge_score_abs <- abs(network_data$edge_score)

  # 2. Create graph directly with edge attributes
  g <- igraph::graph_from_data_frame(d = network_data, directed = FALSE)

  # Apply weight transformation using the *absolute* edge_score
  igraph::E(g)$weight <- sapply(igraph::E(g)$edge_score_abs, function(w) {
    if (is.na(w)) {
      Inf # Treat NA scores as infinite weight (unreachable)
    } else if (w == 0) {
      0.001
    } else if (w == 1) {
      1 / (w + 0.1)
    } else {
      1 / w
    }
  })

  # 3. Validate nodes
  if (!source_node %in% igraph::V(g)$name) stop("Source node not found in network: ", source_node)
  if (!target_node %in% igraph::V(g)$name) stop("Target node not found in network: ", target_node)

  # 4. Find shortest path
  result <- igraph::shortest_paths(g, from = source_node, to = target_node, weights = igraph::E(g)$weight, output = "both")

  path_vertices <- result$vpath[[1]]
  path_edges <- result$epath[[1]]

  if (length(path_vertices) > 1 && length(path_edges) > 0) {
    edge_df <- igraph::as_data_frame(g, what = "edges")[path_edges, ]
    path_df <- dplyr::select(edge_df,
                             Source = from,
                             Target = to,
                             original_edge_score = edge_score,
                             transformed_weight = weight,
                             edge_type = edge_type)

    # Sanitize node names for filename
    safe_from <- gsub("[^A-Za-z0-9_.-]", "_", source_node)
    safe_to <- gsub("[^A-Za-z0-9_.-]", "_", target_node)
    output_filepath <- file.path(output_file, paste0("path_", safe_from, "_to_", safe_to, ".csv"))

    write.csv(path_df, output_filepath, row.names = FALSE)
  } else {
    # No finite path found, do not write a file.
  }
  return(invisible(NULL))
}
