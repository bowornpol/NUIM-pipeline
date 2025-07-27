con_ppn <- function(
  abundance_file,
  metadata_file,
  map_file,
  output_file,
  file_type = c("csv", "tsv"),
  pvalueCutoff,
  pAdjustMethod = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none"),
  rank_by = c("signed_log_pvalue", "log2FoldChange")
) {
  # Validate inputs
  file_type <- match.arg(file_type)
  pAdjustMethod <- match.arg(pAdjustMethod)
  rank_by <- match.arg(rank_by)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_file)) {
    dir.create(output_file, recursive = TRUE)
  }

  # 1. Load abundance data from provided file
  gene_abundance <- tryCatch(
    read_input_file(abundance_file, file_type = file_type, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error loading abundance file '", abundance_file, "': ", e$message, sep = ""))
    }
  )

  # 2. Load sample metadata
  metadata <- tryCatch(
    read_input_file(metadata_file, file_type = file_type, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error loading metadata file '", metadata_file, "': ", e$message, sep = ""))
    }
  )

  # Ensure the correct SampleID and class columns are present
  if (!"SampleID" %in% colnames(metadata)) {
    stop("The 'sample_metadata.csv' file must contain a 'SampleID' column.")
  }
  if (!"class" %in% colnames(metadata)) {
    stop("The 'sample_metadata.csv' file must contain a 'class' column for group definition.")
  }

  # 3. Filter metadata samples present in abundance data
  sample_ids <- colnames(gene_abundance)
  metadata <- dplyr::filter(metadata, SampleID %in% sample_ids)

  # 4. Set 'condition' factor from 'class' column
  metadata$condition <- as.factor(metadata$class)
  rownames(metadata) <- metadata$SampleID

  # 5. Round abundance counts for DESeq2 compatibility and ensure sample order
  gene_abundance_rounded <- round(gene_abundance)
  gene_abundance_rounded <- gene_abundance_rounded[, rownames(metadata), drop = FALSE] # drop=FALSE to handle single sample case

  # 6. Create DESeq2 dataset and run differential expression analysis
  dds <- tryCatch(
    DESeq2::DESeqDataSetFromMatrix(countData = gene_abundance_rounded,
                                   colData = metadata,
                                   design = ~ condition),
    error = function(e) {
      stop(paste("Error creating DESeq2 dataset: ", e$message, ". Check abundance data (must be integers) and metadata consistency.", sep = ""))
    }
  )

  dds <- tryCatch(
    DESeq2::DESeq(dds),
    error = function(e) {
      stop(paste("Error running DESeq2 analysis: ", e$message, ". This might happen if groups have zero variance, or too few samples.", sep = ""))
    }
  )

  # 7. Get all pairwise condition comparisons
  conditions <- levels(metadata$condition)
  if (length(conditions) < 2) {
    stop("Less than two unique conditions found in metadata 'class' column. Cannot perform pairwise comparisons.")
  }
  comparisons <- combn(conditions, 2, simplify = FALSE)

  # 8. Load pathway-to-gene mapping and reshape into TERM2GENE format
  map_raw <- tryCatch(
    read_input_file(map_file, file_type = file_type, header = FALSE, fill = TRUE, stringsAsFactors = FALSE, skip = 1),
    error = function(e) {
      stop(paste("Error loading map file '", map_file, "': ", e$message, ". Ensure it's a valid CSV/TSV.", sep = ""))
    }
  )

  TERM2GENE <- dplyr::distinct(
    dplyr::filter(
      dplyr::select(
        tidyr::gather(map_raw, key = "temp_col", value = "gene", -V1),
        term = V1, gene
      ),
      gene != ""
    )
  )

  if (nrow(TERM2GENE) == 0) {
    stop("No valid pathway-gene mappings found after processing '", map_file, "'. Check file format.")
  }

  gsea_results_list <- list()
  gsea_output_paths <- c() # Initialize vector to store GSEA paths
  jaccard_output_paths <- c() # Initialize vector to store Jaccard paths

  # 9. Loop over each pairwise comparison to run GSEA
  for (i in seq_along(comparisons)) {
    comp <- comparisons[[i]]
    cond1 <- comp[1]
    cond2 <- comp[2]
    comparison_name <- paste0(cond1, "_vs_", cond2)

    # Get DESeq2 results for contrast cond2 vs cond1
    res <- DESeq2::results(dds, contrast = c("condition", cond2, cond1))

    # Prepare ranked gene list
    ranked_df <- as.data.frame(res[, c("log2FoldChange", "pvalue")])
    # Filter out NA values for relevant columns before ranking
    ranked_df <- ranked_df[!is.na(ranked_df$log2FoldChange) & !is.na(ranked_df$pvalue), ]

    if (nrow(ranked_df) == 0) {
      warning("No valid log2FoldChange or pvalue for comparison ", comparison_name, ". Skipping GSEA.")
      next
    }

    # Apply chosen ranking method
    if (rank_by == "signed_log_pvalue") {
      # Handle cases where pvalue might be 0, leading to -log10(0) = Inf
      min_pvalue_for_log <- min(ranked_df$pvalue[ranked_df$pvalue > 0], na.rm = TRUE) / 2
      ranked_df$pvalue[ranked_df$pvalue == 0] <- min_pvalue_for_log
      ranked_df$rank <- sign(ranked_df$log2FoldChange) * -log10(ranked_df$pvalue)
    } else if (rank_by == "log2FoldChange") {
      ranked_df$rank <- ranked_df$log2FoldChange
    }

    ranked_df <- ranked_df[order(ranked_df$rank, decreasing = TRUE), ]
    geneList <- setNames(ranked_df$rank, rownames(ranked_df))

    # Run GSEA using clusterProfiler
    gsea_res <- tryCatch(
      clusterProfiler::GSEA(geneList = geneList,
                            TERM2GENE = TERM2GENE,
                            pvalueCutoff = pvalueCutoff,
                            pAdjustMethod = pAdjustMethod,
                            seed = TRUE,
                            verbose = FALSE),
      error = function(e) {
        warning(paste("GSEA failed for comparison ", comparison_name, ": ", e$message, ". Skipping.", sep = ""))
        return(NULL)
      }
    )

    if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) {
      next
    }

    # Save GSEA results dataframe
    gsea_df <- as.data.frame(gsea_res)
    key <- comparison_name
    gsea_results_list[[key]] <- gsea_df

    # --- Construct full output path for GSEA results with parameters ---

    # Format parameters for filename
    pvalueCutoff_fname <- as.character(pvalueCutoff)
    pAdjustMethod_fname <- pAdjustMethod
    rank_by_fname <- rank_by

    gsea_output_filename <- paste0(
      "gsea_results_", key, "_",
      pvalueCutoff_fname, "_",
      pAdjustMethod_fname, "_",
      rank_by_fname,
      ".csv"
    )

    gsea_output_path <- file.path(output_file, gsea_output_filename)
    write.csv(gsea_df, gsea_output_path, row.names = FALSE)
    gsea_output_paths <- c(gsea_output_paths, gsea_output_path) # Store path
  }

  # 10. Compute Jaccard indices between pathways within each comparison's GSEA results
  jaccard_results_list <- list()

  if (length(gsea_results_list) == 0) {
    # No GSEA results to compute Jaccard indices. Skipping.
  } else {
    for (key in names(gsea_results_list)) {
      gsea_df <- gsea_results_list[[key]]

      # Filter for pathways with core enrichment genes (i.e., not empty or NA)
      gsea_df_filtered <- dplyr::filter(gsea_df, !is.na(core_enrichment) & core_enrichment != "")
      gene_sets <- strsplit(as.character(gsea_df_filtered$core_enrichment), "/")
      gene_sets <- lapply(gene_sets, function(x) unique(na.omit(x)))

      n <- length(gene_sets)
      res_list <- list()

      if (n < 2) { # Ensure there are at least two pathways to compare
        next
      }

      for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
          genes_i <- gene_sets[[i]]
          genes_j <- gene_sets[[j]]

          # Skip if either gene set is empty
          if (length(genes_i) == 0 || length(genes_j) == 0) {
            next
          }

          intersection <- length(intersect(genes_i, genes_j))
          union <- length(union(genes_i, genes_j))
          jaccard <- ifelse(union == 0, 0, intersection / union)

          if (jaccard > 0) {
            res_list[[length(res_list) + 1]] <- data.frame(
              FunctionID_1 = gsea_df_filtered$ID[i],
              FunctionID_2 = gsea_df_filtered$ID[j],
              jaccard_index = jaccard,
              comparison = key,
              stringsAsFactors = FALSE
            )
          }
        }
      }

      if (length(res_list) > 0) {
        jaccard_df <- do.call(rbind, res_list)
        # Construct full output path for Jaccard results using the new format
        jaccard_output_path <- file.path(output_file, paste0("pathway_pathway_network_from_gsea_", key, ".csv"))
        write.csv(jaccard_df, jaccard_output_path, row.names = FALSE)
        jaccard_output_paths <- c(jaccard_output_paths, jaccard_output_path)
      } else {
        # No Jaccard indices > 0 found for comparison. Skipping saving file.
      }
    }
  }

  return(list(gsea_paths = gsea_output_paths, jaccard_paths = jaccard_output_paths))
}
