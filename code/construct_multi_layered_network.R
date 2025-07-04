library(dplyr)
library(tidyr)
library(stringr)

construct_multi_layered_network <- function(
  gsea_results_file,
  microbe_pathway_file,
  pathway_jaccard_file,
  pathway_metabolite_file,
  output_directory 
) {
  message("Starting multi-layered network construction.")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
    message("Created output directory: ", output_directory)
  } else {
    message("Output directory already exists: ", output_directory)
  }
  
  # Determine gsea_suffix and target groups for processing 
  gsea_suffix <- NULL 
  gsea_target_group_from_filename <- NULL 
  
  if (!is.null(gsea_results_file)) {
    gsea_basename <- tools::file_path_sans_ext(basename(gsea_results_file))
    
    match_result <- stringr::str_match(gsea_basename, "^gsea_results_([^_]+)_vs_([^_]+)$")
    
    if (!is.na(match_result[1,1])) { 
      gsea_source_group <- match_result[1,2] 
      gsea_target_group_from_filename <- match_result[1,3] 
      
      gsea_suffix <- gsea_target_group_from_filename 
      message("Derived GSEA filename suffix: '", gsea_suffix, "' (from GSEA target group)")
      
    } else {
      message("  Warning: GSEA results file name '", basename(gsea_results_file), "' did not strictly match the 'gsea_results_*_vs_*.csv' pattern for precise group filtering. No GSEA-specific group filtering will be applied.")
    }
  } else {
    message("No GSEA results file provided, defaulting to processing all groups from metadata.")
  }
  
  # 1. Load GSEA results to identify all pathways to be included
  message("\n1. Identifying pathways from GSEA results for filtering...")
  if (!file.exists(gsea_results_file)) {
    stop("GSEA results file not found: '", gsea_results_file, "'. Cannot proceed.")
  }
  
  gsea_df <- tryCatch(
    read.csv(gsea_results_file, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error reading GSEA results file '", gsea_results_file, "': ", e$message, sep = ""))
    }
  )
  
  if (!"ID" %in% colnames(gsea_df)) {
    stop("GSEA results file '", gsea_results_file, "' must contain an 'ID' column for pathways.")
  }
  
  pathways_to_integrate_set <- unique(as.character(gsea_df$ID))
  
  if (length(pathways_to_integrate_set) == 0) {
    stop("No pathways identified from GSEA results '", gsea_results_file, "'. Cannot construct multi-layered network.")
  }
  message("    Identified ", length(pathways_to_integrate_set), " unique pathways from GSEA results for integration.")
  
  all_network_edges <- list()
  
  # 2. Load and process Microbe-Pathway Network edges
  message("\n2. Processing Microbe-Pathway network edges...")
  if (!file.exists(microbe_pathway_file)) {
    message("    Microbe-Pathway network file not found: '", microbe_pathway_file, "'. Skipping this layer.")
  } else {
    mp_df <- tryCatch(
      read.csv(microbe_pathway_file, stringsAsFactors = FALSE),
      error = function(e) {
        warning(paste("Error reading Microbe-Pathway file '", basename(microbe_pathway_file), "': ", e$message, ". Skipping this layer.", sep = ""))
        return(NULL)
      }
    )
    
    if (!is.null(mp_df) && all(c("TaxonID", "FunctionID", "relative_contribution") %in% colnames(mp_df))) {
      # Filter by pathways identified from GSEA files (FunctionID must be a GSEA pathway)
      mp_filtered <- mp_df %>%
        filter(FunctionID %in% pathways_to_integrate_set) %>%
        select(
          Feature1 = TaxonID,
          Feature2 = FunctionID,
          Edge_Score = relative_contribution
        ) %>%
        mutate(
          Edge_Type = "Microbe-Pathway"
        )
      if (nrow(mp_filtered) > 0) {
        all_network_edges[[length(all_network_edges) + 1]] <- mp_filtered
        message("    Added ", nrow(mp_filtered), " microbe-pathway edges from ", basename(microbe_pathway_file), ".")
      } else {
        message("    No relevant microbe-pathway edges found in ", basename(microbe_pathway_file), " after filtering by GSEA pathways.")
      }
    } else {
      warning("    Microbe-Pathway file '", basename(microbe_pathway_file), "' missing required columns ('TaxonID', 'FunctionID', 'relative_contribution') or is empty. Skipping this layer.")
    }
  }
  
  # 3. Load and process Pathway-Pathway Network edges
  message("\n3. Processing Pathway-Pathway network edges...")
  if (!file.exists(pathway_jaccard_file)) {
    message("    Pathway-Pathway network file not found: '", pathway_jaccard_file, "'. Skipping this layer.")
  } else {
    pp_df <- tryCatch(
      read.csv(pathway_jaccard_file, stringsAsFactors = FALSE),
      error = function(e) {
        warning(paste("Error reading Pathway-Pathway file '", basename(pathway_jaccard_file), "': ", e$message, ". Skipping this layer.", sep = ""))
        return(NULL)
      }
    )
    
    if (!is.null(pp_df) && all(c("FunctionID_1", "FunctionID_2", "jaccard_index") %in% colnames(pp_df))) {
      # Filter by pathways identified from GSEA files (both FunctionID_1 AND FunctionID_2 must be GSEA pathways)
      pp_filtered <- pp_df %>%
        filter(FunctionID_1 %in% pathways_to_integrate_set & FunctionID_2 %in% pathways_to_integrate_set) %>%
        select(
          Feature1 = FunctionID_1,
          Feature2 = FunctionID_2,
          Edge_Score = jaccard_index
        ) %>%
        mutate(
          Edge_Type = "Pathway-Pathway"
        )
      if (nrow(pp_filtered) > 0) {
        all_network_edges[[length(all_network_edges) + 1]] <- pp_filtered
        message("    Added ", nrow(pp_filtered), " pathway-pathway edges from ", basename(pathway_jaccard_file), ".")
      } else {
        message("    No relevant pathway-pathway edges found in ", basename(pathway_jaccard_file), " after filtering by GSEA pathways.")
      }
    } else {
      warning("    Pathway-Pathway file '", basename(pathway_jaccard_file), "' missing required columns ('FunctionID_1', 'FunctionID_2', 'jaccard_index') or is empty. Skipping this layer.")
    }
  }
  
  # 4. Load and process Pathway-Metabolite Network edges
  message("\n4. Processing Pathway-Metabolite network edges...")
  if (!file.exists(pathway_metabolite_file)) {
    message("    Pathway-Metabolite network file not found: '", pathway_metabolite_file, "'. Skipping this layer.")
  } else {
    pm_df <- tryCatch(
      read.csv(pathway_metabolite_file, stringsAsFactors = FALSE),
      error = function(e) {
        warning(paste("Error reading Pathway-Metabolite file '", basename(pathway_metabolite_file), "': ", e$message, ". Skipping this layer.", sep = ""))
        return(NULL)
      }
    )
    
    if (!is.null(pm_df) && all(c("FunctionID", "MetaboliteID", "Correlation") %in% colnames(pm_df))) {
      # Filter by pathways identified from GSEA files (FunctionID must be a GSEA pathway)
      pm_filtered <- pm_df %>%
        filter(FunctionID %in% pathways_to_integrate_set) %>%
        select(
          Feature1 = FunctionID,
          Feature2 = MetaboliteID,
          Edge_Score = Correlation
        ) %>%
        mutate(
          Edge_Type = "Pathway-Metabolite"
        )
      if (nrow(pm_filtered) > 0) {
        all_network_edges[[length(all_network_edges) + 1]] <- pm_filtered
        message("    Added ", nrow(pm_filtered), " pathway-metabolite edges from ", basename(pathway_metabolite_file), ".")
      } else {
        message("    No relevant pathway-metabolite edges found in ", basename(pathway_metabolite_file), " after filtering by GSEA pathways.")
      }
    } else {
      warning("    Pathway-Metabolite file '", basename(pathway_metabolite_file), "' missing required columns ('FunctionID', 'MetaboliteID', 'Correlation') or is empty. Skipping this layer.")
    }
  }
  
  # 5. Combine all network edges into a single data frame
  message("\n5. Combining all network edges...")
  if (length(all_network_edges) == 0) {
    stop("No network edges were collected from any layer after filtering by GSEA pathways. Please check input files and GSEA results.")
  }
  
  # Define the required output columns
  final_cols <- c("Feature1", "Feature2", "Edge_Score", "Edge_Type")
  
  final_network_df <- bind_rows(all_network_edges) %>%
    select(all_of(final_cols)) # Ensure only the specified columns are kept and in order
  
  message("    Total integrated edges: ", nrow(final_network_df))
  
  # 6. Save the final integrated network with a dynamic filename
  dynamic_output_filename <- if (!is.null(gsea_suffix)) {
    paste0("multi_layered_network_", gsea_suffix, ".csv")
  } else {
    # Fallback if GSEA suffix could not be determined
    paste0("multi_layered_network_overall.csv") 
  }
  
  final_output_filepath <- file.path(output_directory, dynamic_output_filename)
  message("\n6. Saving the final integrated multi-layered network to: ", final_output_filepath)
  write.csv(final_network_df, final_output_filepath, row.names = FALSE)
  message("Multi-layered network construction complete.")
  return(invisible(NULL))
}