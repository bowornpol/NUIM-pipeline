node_prior <- function(
  multi_layered_network_file,
  output_file,
  file_type = c("csv", "tsv"),
  time_step_interval = 0.01,
  stabilization_threshold = 0.0001,
  stabilization_window_size = 10,
  filter_other_metabolite_edges
) {
  file_type <- match.arg(file_type)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_file)) {
    dir.create(output_file, recursive = TRUE)
  }

  # Extract base name from input file for output specificity
  input_file_base_name <- tools::file_path_sans_ext(basename(multi_layered_network_file))
  cleaned_input_file_name <- gsub("[^A-Za-z0-9_]", "", input_file_base_name)

  # 1. Load the multi-layered network file (FULL network data)
  if (!file.exists(multi_layered_network_file)) {
    stop("Input network file not found: '", multi_layered_network_file, "'. Cannot proceed.")
  }

  combined_data_full <- tryCatch(
    read_input_file(multi_layered_network_file, file_type = file_type, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error reading multi-layered network file '", multi_layered_network_file, "': ", e$message, sep = ""))
    }
  )

  # Validate required columns in combined_data_full
  required_cols_network <- c("Feature1", "Feature2", "edge_score", "edge_type")
  if (!all(required_cols_network %in% colnames(combined_data_full))) {
    stop(paste("Input network file '", multi_layered_network_file, "' must contain columns: ", paste(required_cols_network, collapse = ", "), ". Please ensure it's the output from Module 2, Step 4.", sep = ""))
  }

  # Ensure edge_score is numeric
  if (!is.numeric(combined_data_full$edge_score)) {
    stop("Column 'edge_score' in the input network file must be numeric.")
  }

  # 2. Identify all unique metabolite nodes from the FULL network (for looping and for filtering)
  all_metabolite_nodes <- unique(
    dplyr::pull(
      dplyr::filter(combined_data_full, edge_type == "Pathway-Metabolite"),
      Feature2
    )
  )
  all_metabolite_nodes <- as.character(all_metabolite_nodes)

  if (length(all_metabolite_nodes) == 0) {
    stop("No metabolite nodes found in the network (based on 'Pathway-Metabolite' edge_type). Cannot proceed diffusion.") # Changed Edge_Type
  }

  # Define the function for heat vector at time t
  H_vector_func <- function(t_val, laplacian_matrix, initial_heat_vector) {
    # Ensure dimensions match before multiplication
    if (nrow(laplacian_matrix) != length(initial_heat_vector)) {
      stop("Dimension mismatch: Laplacian matrix rows (", nrow(laplacian_matrix), ") and initial heat vector length (", length(initial_heat_vector), ") do not match.")
    }
    expm::expm(-laplacian_matrix * t_val) %*% initial_heat_vector
  }

  # Function to find stabilization time and return correlations
  find_stabilization_data <- function(current_L, current_H0, time_interval, threshold, window_size) {
    time_steps_local <- seq(0, 1, by = time_interval)

    heat_vectors_over_time <- sapply(time_steps_local, function(t) {
      H_vector_func(t, current_L, current_H0)
    })

    if (ncol(heat_vectors_over_time) < 2) {
      # If only one time step or no progression, return empty correlations and last time step
      return(list(stabilization_time = time_steps_local[length(time_steps_local)], correlations_df = data.frame(Time = numeric(), Correlation = numeric())))
    }

    correlations <- numeric(ncol(heat_vectors_over_time) - 1)

    for (i in seq_along(correlations)) {
      H_current <- heat_vectors_over_time[, i + 1]
      H_prev <- heat_vectors_over_time[, i]
      correlations[i] <- ifelse(sd(H_current) == 0 || sd(H_prev) == 0, 1, stats::cor(H_current, H_prev, method = "spearman", use = "complete.obs"))
    }

    correlation_df <- data.frame(Time = time_steps_local[-1], Correlation = correlations)

    stabilization_t_found <- time_steps_local[length(correlations) + 1]

    # Ensure window_size is not larger than available correlations
    if (window_size >= length(correlations)) {
      stabilization_t_found <- time_steps_local[length(correlations) + 1]
    } else {
      for (i in seq(length(correlations) - window_size + 1)) {
        diffs <- abs(diff(correlations[i:(i + window_size - 1)]))
        if (all(diffs < threshold)) {
          stabilization_t_found <- time_steps_local[i+1]
          break
        }
      }
    }

    return(list(stabilization_time = stabilization_t_found, correlations_df = correlation_df))
  }

  # 3. Loop through each metabolite node to perform heat diffusion
  for (seed_metabolite_id in all_metabolite_nodes) { # Loop through all identified metabolites as seeds
    # --- Determine network for current seed based on filter_other_metabolite_edges ---
    if (filter_other_metabolite_edges) {
      # Identify all other metabolite nodes to exclude (all metabolites MINUS the current seed)
      other_metabolite_nodes_to_exclude <- setdiff(all_metabolite_nodes, seed_metabolite_id)

      # Filter edges: keep only edges where *neither* Feature1 nor Feature2 is one of the 'other metabolite' nodes
      current_combined_data <- dplyr::filter(
        combined_data_full,
        !(Feature1 %in% other_metabolite_nodes_to_exclude) &
          !(Feature2 %in% other_metabolite_nodes_to_exclude)
      )

      # Check if the seed node itself is still present after filtering
      if (!seed_metabolite_id %in% unique(c(current_combined_data$Feature1, current_combined_data$Feature2))) {
        warning("Seed metabolite '", seed_metabolite_id, "' is not present in the graph after filtering. This implies it was only connected to other metabolites, which were removed. Skipping diffusion for this seed.")
        next
      }

    } else {
      # Use the full network data (original behavior)
      current_combined_data <- combined_data_full
    }

    # --- Convert edge_score to absolute value ---
    current_combined_data$edge_score <- abs(current_combined_data$edge_score)

    # --- Create graph and Laplacian matrix for the CURRENT network configuration ---

    # Check for valid edges/nodes to create a graph
    if (nrow(current_combined_data) == 0) {
      warning("No edges found for graph construction for seed '", seed_metabolite_id, "' after filtering. Skipping diffusion.")
      next
    }

    # Get all unique nodes from the current filtered data to ensure the graph includes them
    all_nodes_in_current_data <- unique(c(current_combined_data$Feature1, current_combined_data$Feature2))

    g_current <- igraph::graph_from_data_frame(d = current_combined_data[, c("Feature1", "Feature2")], directed = FALSE,
                                               vertices = all_nodes_in_current_data)

    # Assign weights (handling potential mismatches as before)
    graph_edges_for_weighting <- dplyr::ungroup(
      dplyr::mutate(
        dplyr::rowwise(igraph::as_data_frame(g_current, what = "edges")),
        Node_A = min(from, to), Node_B = max(from, to)
      )
    )

    current_combined_data_sorted <- dplyr::ungroup(
      dplyr::distinct(
        dplyr::select(
          dplyr::mutate(
            dplyr::rowwise(current_combined_data),
            Node_A = min(Feature1, Feature2), Node_B = max(Feature1, Feature2)
          ),
          Node_A, Node_B, edge_score
        )
      )
    )

    # Create a unique identifier for edges in both dataframes for robust matching
    edge_id_graph <- paste0(graph_edges_for_weighting$Node_A, "_", graph_edges_for_weighting$Node_B)
    edge_id_data <- paste0(current_combined_data_sorted$Node_A, "_", current_combined_data_sorted$Node_B)

    matched_weights <- current_combined_data_sorted$edge_score[match(edge_id_graph, edge_id_data)]

    igraph::E(g_current)$weight <- matched_weights

    if (any(is.na(igraph::E(g_current)$weight))) {
      warning("Some edges in the graph for seed '", seed_metabolite_id, "' could not be matched to an 'edge_score' in the input data. Assigning 0 weight to unmatched edges.")
      igraph::E(g_current)$weight[is.na(igraph::E(g_current)$weight)] <- 0
    }

    # If the graph has become too small (e.g., only 1 node, no edges, or disconnected)
    if (igraph::vcount(g_current) < 2 || igraph::ecount(g_current) == 0 || !igraph::is.connected(g_current)) {
      warning("Graph for seed '", seed_metabolite_id, "' is too sparse, disconnected, or lacks sufficient nodes/edges after filtering (nodes: ", igraph::vcount(g_current), ", edges: ", igraph::ecount(g_current), ", connected: ", igraph::is.connected(g_current), "). Skipping diffusion for this seed as Laplacian cannot be computed meaningfully.")
      next
    }

    L_current <- igraph::laplacian_matrix(g_current, weights = igraph::E(g_current)$weight)
    L_current <- as.matrix(L_current)

    # Initialize the heat vector H_0 for the current seed within the CURRENT graph's nodes
    H_0_current_graph <- numeric(igraph::vcount(g_current))
    names(H_0_current_graph) <- igraph::V(g_current)$name
    seed_index_current_graph <- which(igraph::V(g_current)$name == seed_metabolite_id)

    if (length(seed_index_current_graph) == 0) {
      warning("Seed metabolite '", seed_metabolite_id, "' not found in the *filtered* graph's node set. This should not happen if previous checks are correct. Skipping diffusion.")
      next
    }
    H_0_current_graph[seed_index_current_graph] <- 1.0

    # Find stabilization time and get correlations data
    stabilization_data_result <- find_stabilization_data(L_current, H_0_current_graph, time_step_interval, stabilization_threshold, stabilization_window_size)
    stabilization_t <- stabilization_data_result$stabilization_time
    correlation_df <- stabilization_data_result$correlations_df

    # Calculate final heat scores at stabilization time
    final_heat_scores <- H_vector_func(stabilization_t, L_current, H_0_current_graph)

    # Create output data frame for this metabolite
    output_df <- dplyr::arrange(
      data.frame(
        Node = igraph::V(g_current)$name, # Nodes from the current filtered graph
        Heat_score = round(final_heat_scores, 10),
        stringsAsFactors = FALSE
      ),
      dplyr::desc(Heat_score) # Sort by heat score (descending)
    )

    # Generate and save output files with cleaned input file name
    cleaned_seed_id <- gsub("[^A-Za-z0-9_]", "", seed_metabolite_id)

    # Save heat scores to CSV
    output_file_name_heat <- paste0("heat_scores_", cleaned_seed_id, "_", cleaned_input_file_name, ".csv")
    output_path_heat <- file.path(output_file, output_file_name_heat)
    write.csv(output_df, file = output_path_heat, row.names = FALSE)

    # Save correlation data to CSV
    output_file_name_correlation_csv <- paste0("spearman_correlations_", cleaned_seed_id, "_", cleaned_input_file_name, ".csv")
    output_path_correlation_csv <- file.path(output_file, output_file_name_correlation_csv)
    write.csv(correlation_df, file = output_path_correlation_csv, row.names = FALSE)

    # Generate and save correlation plot
    correlation_plot <- ggplot2::ggplot(correlation_df, ggplot2::aes(x = Time, y = Correlation)) +
      ggplot2::geom_line() +
      ggplot2::geom_vline(xintercept = stabilization_t, linetype = "dashed", color = "#a62140") +
      ggplot2::geom_text(
        ggplot2::aes(
          x = stabilization_t,
          y = max(Correlation, na.rm = TRUE) * 0.5,
          label = paste("t =", round(stabilization_t, 4))
        ),
        color = "#a62140", hjust = -0.1, vjust = 0.5, size = 5, fontface = "bold"
      ) +
      ggplot2::xlab("Time step") +
      ggplot2::ylab("Spearman correlation with previous time step") +
      ggplot2::scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 1)
      )

    output_file_name_plot <- paste0("correlation_plot_", cleaned_seed_id, "_", cleaned_input_file_name, ".jpg")
    output_path_plot <- file.path(output_file, output_file_name_plot)
    ggplot2::ggsave(output_path_plot, plot = correlation_plot, width = 8, height = 5, dpi = 600)
  }

  return(invisible(NULL))
}
