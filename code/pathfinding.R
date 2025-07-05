library(igraph)
library(dplyr)
library(stringr)

pathfinding <- function(
  multi_layered_network_file,
  source_node,
  target_node,
  output_directory
) {
  message("Starting shortest pathfinding using Dijkstra's algorithm.")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
    message("Created output directory: ", output_directory)
  } else {
    message("Output directory already exists: ", output_directory)
  }
  
  # 1. Load network
  message("\n1. Loading multi-layered network from: ", multi_layered_network_file)
  
  required_cols <- c("Feature1", "Feature2", "Edge_Score", "Edge_Type")
  if (!file.exists(multi_layered_network_file)) {
    stop("Network file not found: '", multi_layered_network_file, "'")
  }
  
  network_data <- read.csv(multi_layered_network_file, stringsAsFactors = FALSE)
  if (!all(required_cols %in% colnames(network_data))) {
    stop("Network file must contain columns: ", paste(required_cols, collapse = ", "))
  }
  message("  Successfully loaded and validated network file.")
  
  if (!is.numeric(network_data$Edge_Score)) {
    stop("Column 'Edge_Score' must be numeric.")
  }
  
  # Convert Edge_Score to absolute
  network_data$Edge_Score <- abs(network_data$Edge_Score)
  
  # 2. Create graph directly with edge attributes
  message("  Creating graph and assigning weights.")
  g <- graph_from_data_frame(d = network_data, directed = FALSE)
  
  # Apply weight transformation
  E(g)$weight <- sapply(E(g)$Edge_Score, function(w) {
    if (is.na(w)) {
      Inf
    } else if (w < 1) {
      1 / w
    } else if (w == 1) {
      1 / (w + 0.1)
    } else {
      w
    }
  })
  
  # 3. Validate nodes
  message("\n2. Validating source and target nodes against the network.")
  if (!source_node %in% V(g)$name) stop("Source node not found in network: ", source_node)
  if (!target_node %in% V(g)$name) stop("Target node not found in network: ", target_node)
  message("  Both source and target nodes found.")
  
  # 4. Find shortest path
  message("\n3. Finding shortest path from '", source_node, "' to '", target_node, "'...")
  result <- shortest_paths(g, from = source_node, to = target_node, weights = E(g)$weight, output = "both")
  
  path_vertices <- result$vpath[[1]]
  path_edges <- result$epath[[1]]
  
  if (length(path_vertices) > 1 && length(path_edges) > 0) {
    message("  Path found with ", length(path_edges), " steps.")
    
    edge_df <- igraph::as_data_frame(g, what = "edges")[path_edges, ]
    path_df <- edge_df %>%
      select(Source = from, Target = to, Original_Edge_Score = Edge_Score, Transformed_Weight = weight, Edge_Type)
    
    safe_from <- gsub("[^A-Za-z0-9_.-]", "_", source_node)
    safe_to <- gsub("[^A-Za-z0-9_.-]", "_", target_node)
    output_file <- file.path(output_directory, paste0("path_", safe_from, "_to_", safe_to, ".csv"))
    
    write.csv(path_df, output_file, row.names = FALSE)
    message("  Saved path to: ", output_file)
    
    invisible(list(path_df = path_df))
    
    message("Node prioritization complete.")
    return(invisible(NULL))
  } else {
    message("  No finite path found between source and target.")
    invisible(NULL)
  }
}