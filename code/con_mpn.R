con_mpn <- function(
  contrib_file,
  metadata_file,
  taxonomy_file,
  output_file,
  file_type = c("csv", "tsv"),
  filtering = c("unfiltered", "mean", "median", "top10%", "top25%", "top50%", "top75%")
) {
  file_type <- match.arg(file_type)
  filtering <- match.arg(filtering)

  if (!dir.exists(output_file)) {
    dir.create(output_file, recursive = TRUE)
  }

  contrib <- tryCatch(
    read_input_file(contrib_file, file_type = file_type, stringsAsFactors = FALSE),
    error = function(e) stop(paste("Error loading contribution file: ", e$message))
  )
  metadata <- tryCatch(
    read_input_file(metadata_file, file_type = file_type, stringsAsFactors = FALSE),
    error = function(e) stop(paste("Error loading metadata file: ", e$message))
  )
  taxonomy <- tryCatch(
    read_input_file(taxonomy_file, file_type = file_type, stringsAsFactors = FALSE),
    error = function(e) stop(paste("Error loading taxonomy file: ", e$message))
  )

  if (!"class" %in% colnames(metadata)) {
    metadata$class <- "all"
  } else {
    metadata$class <- as.character(metadata$class)
  }

  contrib$SampleID <- as.character(contrib$SampleID)
  metadata$SampleID <- as.character(metadata$SampleID)
  contrib$FeatureID <- as.character(contrib$FeatureID)
  taxonomy$FeatureID <- as.character(taxonomy$FeatureID)
  contrib$FunctionID <- as.character(contrib$FunctionID)
  taxonomy$TaxonID <- as.character(taxonomy$TaxonID)

  merged <- merge(contrib, metadata, by = "SampleID", all.x = TRUE)
  if (nrow(merged) == 0) stop("Merging contribution and metadata resulted in an empty data frame.")

  cols_to_keep_from_taxonomy <- c("FeatureID", "TaxonID")
  if (!all(cols_to_keep_from_taxonomy %in% colnames(taxonomy))) {
    stop("Missing expected columns in taxonomy file: ", paste(setdiff(cols_to_keep_from_taxonomy, colnames(taxonomy)), collapse = ", "))
  }
  taxonomy_for_merge <- dplyr::select(taxonomy, dplyr::all_of(cols_to_keep_from_taxonomy))

  merged <- merge(merged, taxonomy_for_merge, by = "FeatureID", all.x = TRUE)
  if (!"taxon_function_abun" %in% colnames(merged)) {
    stop("Column 'taxon_function_abun' not found after merging.")
  }
  merged$taxon_function_abun <- as.numeric(merged$taxon_function_abun)

  unique_classes_clean <- unique(merged$class)
  unique_classes_clean <- unique_classes_clean[!is.na(unique_classes_clean)]

  if (length(unique_classes_clean) == 0) {
    warning("No valid (non-NA) classes found. Exiting.")
    return(character(0))
  }

  microbe_pathway_output_paths <- c()

  for (current_class in unique_classes_clean) {
    merged_class <- dplyr::filter(merged, class == current_class)
    if (nrow(merged_class) == 0) next

    taxon_function_total_class <- aggregate(
      taxon_function_abun ~ FunctionID + TaxonID,
      data = merged_class, sum, na.rm = TRUE
    )
    if (nrow(taxon_function_total_class) == 0) next

    function_total_class <- aggregate(
      taxon_function_abun ~ FunctionID,
      data = taxon_function_total_class, sum, na.rm = TRUE
    )
    colnames(function_total_class)[2] <- "total_abundance_all_taxa"
    taxon_function_total_class <- merge(taxon_function_total_class, function_total_class, by = "FunctionID")
    taxon_function_total_class$relative_contribution <- with(taxon_function_total_class,
                                                             ifelse(total_abundance_all_taxa == 0, 0, taxon_function_abun / total_abundance_all_taxa)
    )

    if (filtering != "unfiltered") {
      if (filtering %in% c("mean", "median")) {
        FUN_used <- if (filtering == "mean") mean else median
        threshold_df <- aggregate(relative_contribution ~ FunctionID,
                                  data = taxon_function_total_class,
                                  FUN = FUN_used, na.rm = TRUE)
        colnames(threshold_df)[2] <- "threshold"
        taxon_function_total_class <- merge(taxon_function_total_class, threshold_df, by = "FunctionID")
        taxon_function_total_class <- subset(taxon_function_total_class, relative_contribution >= threshold | is.na(threshold))
      } else {
        percent_map <- c("top10%" = 0.10, "top25%" = 0.25, "top50%" = 0.50, "top75%" = 0.75)
        taxon_function_total_class <- dplyr::ungroup(
          dplyr::filter(
            dplyr::mutate(
              dplyr::arrange(
                dplyr::group_by(taxon_function_total_class, FunctionID),
                dplyr::desc(relative_contribution)
              ),
              rank = dplyr::row_number(),
              n_taxa = dplyr::n(),
              cutoff = pmax(ceiling(percent_map[filtering] * dplyr::n()), 1)
            ),
            rank <= cutoff
          )
        )
        taxon_function_total_class <- dplyr::select(taxon_function_total_class, -rank, -n_taxa, -cutoff)
      }
      if (nrow(taxon_function_total_class) == 0) {
        next
      }
    }

    taxon_function_total_class <- taxon_function_total_class[
      order(taxon_function_total_class$FunctionID, -taxon_function_total_class$relative_contribution),
    ]
    file_suffix <- gsub("%", "", filtering)
    full_output_path <- file.path(output_file, paste0("microbe_pathway_network_", current_class, "_", file_suffix, ".csv"))
    write.csv(taxon_function_total_class, full_output_path, row.names = FALSE)
    microbe_pathway_output_paths <- c(microbe_pathway_output_paths, full_output_path)
  }

  return(microbe_pathway_output_paths)
}
