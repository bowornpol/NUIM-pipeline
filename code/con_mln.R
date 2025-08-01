con_mln <- function(
  gsea_results_file,
  microbe_pathway_file,
  pathway_pathway_file,
  pathway_metabolite_file,
  output_file,
  file_type = c("csv", "tsv")
) {
  file_type <- match.arg(file_type)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_file)) {
    dir.create(output_file, recursive = TRUE)
  }

  # Initialize variables to hold GSEA derived names for filename and filtering
  gsea_suffix <- NULL
  gsea_comparison_group <- NULL
  gsea_pathways_to_filter <- NULL

  # 1. Load GSEA results and extract relevant pathways and filename components
  if (!is.null(gsea_results_file)) {
    if (file.exists(gsea_results_file)) {
      gsea_basename <- tools::file_path_sans_ext(basename(gsea_results_file))

      # Extract source and target groups, and the full comparison group
      match_result <- stringr::str_match(gsea_basename, "^gsea_results_([^_]+)_vs_([^_]+).*$")

      if (!is.na(match_result[1,1])) {
        gsea_source_group <- match_result[1,2]
        gsea_target_group <- match_result[1,3]
        gsea_suffix <- gsea_target_group # Assign target group to gsea_suffix
        gsea_comparison_group <- paste0(gsea_source_group, "_vs_", gsea_target_group) # Assign full comparison group
      }

      gsea_results_df <- tryCatch(
        read_input_file(gsea_results_file, file_type = file_type, header = TRUE, stringsAsFactors = FALSE),
        error = function(e) {
          stop(paste("Error loading GSEA results file '", gsea_results_file, "': ", e$message, sep = ""))
        }
      )

      if ("ID" %in% colnames(gsea_results_df)) {
        # Filter for significant pathways based on 'p.adjust' and 'NES'
        # Using a relaxed cutoff for initial pathway filtering for integration
        # (specific cutoffs for each layer are handled by their respective functions)
        gsea_pathways_to_filter <- unique(gsea_results_df$ID)
      } else {
        warning("GSEA results file must contain an 'ID' column with pathway identifiers. No GSEA pathway filtering will be applied.")
      }
    } else {
      warning("GSEA results file '", gsea_results_file, "' not found. No GSEA pathway filtering will be applied to layers.")
    }
  } else {
    # No GSEA results file provided. No GSEA pathway filtering will be applied to layers.
  }

  all_network_edges <- list()

  # 2. Process Microbe-Pathway Network
  if (!is.null(microbe_pathway_file) && file.exists(microbe_pathway_file)) {
    mp_df <- tryCatch(
      read_input_file(microbe_pathway_file, file_type = file_type, header = TRUE, stringsAsFactors = FALSE),
      error = function(e) {
        warning(paste("Error loading Microbe-Pathway file '", microbe_pathway_file, "': ", e$message, ". Skipping this layer.", sep = ""))
        return(NULL)
      }
    )
    if (!is.null(mp_df)) {
      # Ensure expected columns are present
      required_mp_cols <- c("TaxonID", "FunctionID", "relative_contribution")
      if (!all(required_mp_cols %in% colnames(mp_df))) {
        warning("Microbe-Pathway file missing required columns (TaxonID, FunctionID, relative_contribution). Skipping this layer.")
        mp_df <- NULL
      } else {
        # Filter by GSEA pathways if available
        if (!is.null(gsea_pathways_to_filter)) {
          mp_df <- dplyr::filter(mp_df, FunctionID %in% gsea_pathways_to_filter)
        }

        if (nrow(mp_df) > 0) {
          # Standardize column names
          mp_df_standardized <- dplyr::select(mp_df,
                                              Feature1 = TaxonID,
                                              Feature2 = FunctionID,
                                              edge_score = relative_contribution
          )
          mp_df_standardized$edge_type <- "Microbe-Pathway"
          all_network_edges[[length(all_network_edges) + 1]] <- mp_df_standardized
        }
      }
    }
  }

  # 3. Process Pathway-Pathway Network (Jaccard Index)
  if (!is.null(pathway_pathway_file) && file.exists(pathway_pathway_file)) {
    pj_df <- tryCatch(
      read_input_file(pathway_pathway_file, file_type = file_type, header = TRUE, stringsAsFactors = FALSE),
      error = function(e) {
        warning(paste("Error loading Pathway-Pathway Jaccard file '", pathway_pathway_file, "': ", e$message, ". Skipping this layer.", sep = ""))
        return(NULL)
      }
    )
    if (!is.null(pj_df)) {
      required_pj_cols <- c("FunctionID_1", "FunctionID_2", "jaccard_index")
      if (!all(required_pj_cols %in% colnames(pj_df))) {
        warning("Pathway-Pathway Jaccard file missing required columns (FunctionID_1, FunctionID_2, jaccard_index). Skipping this layer.")
        pj_df <- NULL
      } else {
        # Filter by GSEA pathways if available (both ends of the edge)
        if (!is.null(gsea_pathways_to_filter)) {
          pj_df <- dplyr::filter(pj_df, FunctionID_1 %in% gsea_pathways_to_filter & FunctionID_2 %in% gsea_pathways_to_filter)
        }

        if (nrow(pj_df) > 0) {
          # Standardize column names
          pj_df_standardized <- dplyr::select(pj_df,
                                              Feature1 = FunctionID_1,
                                              Feature2 = FunctionID_2,
                                              edge_score = jaccard_index
          )
          pj_df_standardized$edge_type <- "Pathway-Pathway"
          all_network_edges[[length(all_network_edges) + 1]] <- pj_df_standardized
        }
      }
    }
  }

  # 4. Process Pathway-Metabolite Network
  if (!is.null(pathway_metabolite_file) && file.exists(pathway_metabolite_file)) {
    pm_df <- tryCatch(
      read_input_file(pathway_metabolite_file, file_type = file_type, header = TRUE, stringsAsFactors = FALSE),
      error = function(e) {
        warning(paste("Error loading Pathway-Metabolite file '", pathway_metabolite_file, "': ", e$message, ". Skipping this layer.", sep = ""))
        return(NULL)
      }
    )
    if (!is.null(pm_df)) {
      required_pm_cols <- c("FunctionID", "MetaboliteID", "correlation")
      if (!all(required_pm_cols %in% colnames(pm_df))) {
        warning("Pathway-Metabolite file missing required columns (FunctionID, MetaboliteID, correlation). Skipping this layer.")
        pm_df <- NULL
      } else {
        # Filter by GSEA pathways if available
        if (!is.null(gsea_pathways_to_filter)) {
          pm_df <- dplyr::filter(pm_df, FunctionID %in% gsea_pathways_to_filter)
        }

        if (nrow(pm_df) > 0) {
          # Standardize column names
          pm_df_standardized <- dplyr::select(pm_df,
                                              Feature1 = FunctionID,
                                              Feature2 = MetaboliteID,
                                              edge_score = correlation
          )
          pm_df_standardized$edge_type <- "Pathway-Metabolite"
          all_network_edges[[length(all_network_edges) + 1]] <- pm_df_standardized
        }
      }
    }
  }

  # 5. Combine all network edges into a single data frame
  if (length(all_network_edges) == 0) {
    stop("No network edges were collected from any layer after filtering by GSEA pathways. Please check input files and GSEA results.")
  }

  # Define the required output columns
  final_cols <- c("Feature1", "Feature2", "edge_score", "edge_type") # Changed column names

  final_network_df <- dplyr::select(dplyr::bind_rows(all_network_edges), dplyr::all_of(final_cols)) # Ensure only the specified columns are kept and in order

  # 6. Save the final integrated network with a dynamic filename
  dynamic_output_filename <- if (!is.null(gsea_suffix) && !is.null(gsea_comparison_group)) {
    # Use the target group and the full comparison group from GSEA
    paste0("multi_layered_network_", gsea_suffix, "_from_gsea_", gsea_comparison_group, ".csv")
  } else {
    # Fallback if GSEA suffix or comparison group could not be determined
    paste0("multi_layered_network_overall.csv")
  }

  final_output_filepath <- file.path(output_file, dynamic_output_filename)
  write.csv(final_network_df, final_output_filepath, row.names = FALSE)

  return(final_output_filepath)
}
