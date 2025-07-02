# NUIM: Network-based utility for integrating microbiome and metabolome data

We developed NUIM, a modular, network-based framework for integrating microbiome and metabolome data systematically. NUIM consists of three modules: (1) data preparation and processing, (2) network construction, and (3) network analysis.

<p align="center">
  <img src="figures/NUIM_overview.png" width="500"/>
</p>

## Module 1: Data preparation and processing

This module defines the procedures required to prepare and process the input data for downstream network construction.

- Input data includes `microbial sequencing reads in FASTQ format` and `metabolite concentration table`.  
- Microbiome data processing involves the use of QIIME2 to generate a feature table and representative sequences. These outputs are subsequently processed using PICRUSt2 for functional prediction, yielding gene abundance, pathway abundance, and pathway contribution data.  
- Although metabolome data processing may vary depending on user preference and experimental design, NUIM assumes that metabolite concentrations have been appropriately processed by standard practice. For example, users may employ established platforms such as [Metabox](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giae005/7629842) or [MetaboAnalyst](https://academic.oup.com/nar/article/52/W1/W398/7642060) to perform metabolomics data processing.

### <ins>QIIME2 workflow</ins>

This section provides a general QIIME2 workflow for processing paired-end 16S rRNA sequencing data. The goal is to generate a feature table and representative sequences for PICRUSt2.

#### **Required inputs**

| File            | Description                         |
|-----------------|-----------------------------------|
| `FASTQ`         | Raw paired-end sequencing reads. |
| `manifest.tsv`  | Mapping of sample IDs to FASTQ files. |  

<p align="center">
  <img src="figures/QIIME2_overview.png" width="500"/>
</p>

```bash
# Activate QIIME2 environment
conda activate qiime2

# STEP 1: Import paired-end FASTQ files using a manifest file
# The manifest CSV should map sample IDs to file paths of forward and reverse reads.

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.csv \
  --input-format PairedEndFastqManifestPhred33 \
  --output-path demux_reads.qza

# STEP 2: Trim primers/adapters
# Replace the primer sequences below with the ones used in your sequencing protocol.

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux_reads.qza \
  --p-front-f <forward_primer_sequence> \
  --p-front-r <reverse_primer_sequence> \
  --o-trimmed-sequences trimmed_reads.qza

# STEP 3: Summarize read quality
# Review this visualization to determine appropriate truncation lengths for the next step.

qiime demux summarize \
  --i-data trimmed_reads.qza \
  --o-visualization quality_summary.qzv

# STEP 4: Denoise reads using DADA2
# Use the quality summary from Step 3 to choose truncation lengths and other parameters.

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed_reads.qza \
  --p-trunc-len-f <truncation_length_forward> \
  --p-trunc-len-r <truncation_length_reverse> \
  --o-table feature_table.qza \
  --o-representative-sequences rep_seqs.qza \
  --o-denoising-stats denoising_stats.qza

# STEP 5: Filter low-abundance features
# Adjust minimum sample count and read frequency as appropriate for your dataset.

qiime feature-table filter-features \
  --i-table feature_table.qza \
  --p-min-samples <minimum_number_of_samples> \
  --p-min-frequency <minimum_total_frequency> \
  --o-filtered-table filtered_table.qza

# STEP 6: Rarefy the table to a uniform sequencing depth
# Use rarefaction curves to help select the sampling depth.

qiime feature-table rarefy \
  --i-table filtered_table.qza \
  --p-sampling-depth <depth_to_rarefy> \
  --o-rarefied-table normalized_table.qza

# STEP 7: Taxonomic classification
# Use a pre-trained classifier appropriate for your 16S region (e.g., SILVA, Greengenes).

qiime feature-classifier classify-sklearn \
  --i-classifier <pretrained_classifier.qza> \
  --i-reads filtered_rep_seqs.qza \
  --o-classification taxonomy.qza

# STEP 8: Export data for PICRUSt2
# Export the final feature table and representative sequences for functional prediction.

qiime tools export \
  --input-path normalized_table.qza \
  --output-path feature-table

qiime tools export \
  --input-path rep_seqs.qza \
  --output-path rep_seqs

qiime tools export \
  --input-path taxonomy.qza \
  --output-path taxonomy
```

#### **QIIME2 outputs**

| File                | Description                      |
|---------------------|---------------------------------|
| `feature-table.biom` | Feature (ASV) count table. |
| `rep_seqs.fasta`    | Representative sequences. |
| `taxonomy.tsv`      | Taxonomic assignment of ASVs. | 

### <ins>PICRUSt2 workflow</ins>

PICRUSt2 predicts functional profiles from 16S rRNA data. This step uses a feature table and representative sequences from QIIME2.

#### **Required inputs**

| File                | Description                     |
|---------------------|---------------------------------|
| `feature-table.biom` | Feature table exported from QIIME2. |
| `rep_seqs.fasta`    | Representative sequences exported from QIIME2. |
| `pathway_gene_map.tsv` | Maps pathway IDs to their associated gene/KO IDs. **First column:** pathway ID; **other columns:** gene IDs. |

<p align="center">
  <img src="figures/PICRUSt2_overview.png" width="700"/>
</p>

```bash
# Activate PICRUSt2 environment
conda activate picrust2

# Run PICRUSt2 pipeline (no NSTI filtering)
picrust2_pipeline.py \
  -s rep_seqs.fasta \
  -i feature-table.biom \
  -o picrust2_out_nsti_default

# Optional: Add --max_nsti to filter low-confidence ASVs
# e.g., --max_nsti 0.5

# Infer pathway abundance and pathway contribution
# Users can choose to use either KEGG or MetaCyc pathways

# Generate pathway abundance data (sample × pathway)
pathway_pipeline.py \
  -i picrust2_out/<function_output_folder>/pred_metagenome_unstrat.tsv.gz \
  -o pathways_abundance \
  --no_regroup \
  --map <pathway_gene_map.tsv>

# Generate pathway contribution data (ASV × pathway)
pathway_pipeline.py \
  -i picrust2_out/<function_output_folder>/pred_metagenome_unstrat.tsv.gz \
  -o pathways_contrib \
  --per_sequence_contrib \
  --per_sequence_abun picrust2_out/<function_output_folder>/seqtab_norm.tsv.gz \
  --per_sequence_function <predicted_functions.tsv.gz> \
  --no_regroup \
  --map <pathway_gene_map.tsv>

# Gene abundance data (KEGG Orthologs) output by PICRUSt2:
# Located in: KO_metagenome_out/pred_metagenome_unstrat.tsv.gz
# This file contains predicted gene family (KO) abundances per sample.
```

#### **PICRUSt2 outputs**

| File                      | Description                            |
|---------------------------|--------------------------------------|
| `pred_metagenome_unstrat.tsv` | Predicted gene (KO) abundance. |
| `path_abun_unstrat.tsv`   | Predicted pathway abundance. |
| `path_abun_contrib.tsv`   | Predicted pathway contribution. |

## Module 2: Network construction

This module constructs a tripartite network linking microbial taxa, functional pathways, and metabolites. Follow the steps below to construct each network layer using the provided R functions:

### <ins>Step 1: Microbe–pathway network construction</ins>

The microbe–pathway network is constructed from pathway contribution data, with edges representing the relative contribution of each microbe to specific pathways.

#### **Required inputs**

| File | Description |
|---|---|
| `path_abun_contrib.csv` | Pathway contribution data from PICRUSt2. **Required columns**: `SampleID`, `FeatureID`, `FunctionID`, `taxon_function_abun`. |
| `sample_metadata.csv` | Sample metadata. **Required for group-specific analysis**; must contain `SampleID` and `class` columns. If not provided or `class` column is missing, data will be processed as one 'overall' group. |
| `taxonomy.csv` | Taxonomy annotations mapping `FeatureID` to taxonomic name (`TaxonID`) from QIIME2. **Required columns**: `FeatureID`, `TaxonID`. |

```r
library(dplyr)

construct_microbe_pathway_network <- function(
  contrib_file,
  metadata_file,
  taxonomy_file,
  output_file, # Now directly the output directory path
  filtering = c("unfiltered", "mean", "median", "top10%", "top25%", "top50%", "top75%")
) {
  filtering <- match.arg(filtering)
  message("Starting microbe-pathway network construction with filtering: ", filtering)
  
  # --- Create output directory if it doesn't exist ---
  # Now 'output_file' is treated directly as the directory path
  if (!dir.exists(output_file)) {
    dir.create(output_file, recursive = TRUE)
    message("Created output directory: ", output_file)
  } else {
    message("Output directory already exists: ", output_file)
  }
  # ----------------------------------------------------
  
  # Load data with robust error handling
  message("Loading data...")
  contrib <- tryCatch(
    read.csv(contrib_file, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error loading contribution file '", contrib_file, "': ", e$message, sep = ""))
    }
  )
  metadata <- tryCatch(
    read.csv(metadata_file, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error loading metadata file '", metadata_file, "': ", e$message, sep = ""))
    }
  )
  taxonomy <- tryCatch(
    read.csv(taxonomy_file, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error loading taxonomy file '", taxonomy_file, "': ", e$message, sep = ""))
    }
  )
  
  # --- Debugging: Check initial data frames dimensions and column names ---
  message("Dim contrib: ", paste(dim(contrib), collapse = "x"), " | Cols: ", paste(colnames(contrib), collapse = ", "))
  message("Dim metadata: ", paste(dim(metadata), collapse = "x"), " | Cols: ", paste(colnames(metadata), collapse = ", "))
  message("Dim taxonomy: ", paste(dim(taxonomy), collapse = "x"), " | Cols: ", paste(colnames(taxonomy), collapse = ", "))
  # -----------------------------------------------------------------------
  
  # Ensure 'class' column exists in metadata
  if (!"class" %in% colnames(metadata)) {
    message("'class' column not found in metadata. Creating a default 'all' class.")
    metadata$class <- "all"
  } else {
    message("'class' column found in metadata.")
    metadata$class <- as.character(metadata$class)
  }
  
  # Ensure ID columns are consistent types (character) for merging
  if ("SampleID" %in% colnames(contrib)) contrib$SampleID <- as.character(contrib$SampleID)
  if ("SampleID" %in% colnames(metadata)) metadata$SampleID <- as.character(metadata$SampleID)
  if ("FeatureID" %in% colnames(contrib)) contrib$FeatureID <- as.character(contrib$FeatureID)
  if ("FeatureID" %in% colnames(taxonomy)) taxonomy$FeatureID <- as.character(taxonomy$FeatureID)
  if ("FunctionID" %in% colnames(contrib)) contrib$FunctionID <- as.character(contrib$FunctionID)
  if ("TaxonID" %in% colnames(taxonomy)) taxonomy$TaxonID <- as.character(taxonomy$TaxonID)
  
  # Merge contribution and metadata
  message("Merging contribution and metadata...")
  merged <- merge(contrib, metadata, by = "SampleID", all.x = TRUE)
  message("  Dim after merging contrib and metadata: ", paste(dim(merged), collapse = "x"))
  if (nrow(merged) == 0) {
    stop("Merging contribution and metadata resulted in an empty data frame. Check 'SampleID' column consistency and data in both files.")
  }
  
  # Select only necessary columns from taxonomy (FeatureID, TaxonID) to avoid conflicts
  cols_to_keep_from_taxonomy <- c("FeatureID", "TaxonID")
  message("Pre-merge check with taxonomy:")
  message("  Columns in 'merged' before taxonomy merge: ", paste(colnames(merged), collapse = ", "))
  message("  Columns in 'taxonomy' before merge: ", paste(colnames(taxonomy), collapse = ", "))
  if (!all(cols_to_keep_from_taxonomy %in% colnames(taxonomy))) {
    missing_tax_cols <- setdiff(cols_to_keep_from_taxonomy, colnames(taxonomy))
    stop(paste("Error: Missing expected columns in taxonomy file for selection: ", paste(missing_tax_cols, collapse = ", ")))
  }
  taxonomy_for_merge <- taxonomy %>% select(all_of(cols_to_keep_from_taxonomy))
  message("  Columns selected from 'taxonomy' for merge: ", paste(colnames(taxonomy_for_merge), collapse = ", "))
  
  # Merge taxonomy information
  message("Merging with taxonomy data...")
  merged <- merge(merged, taxonomy_for_merge, by = "FeatureID", all.x = TRUE)
  message("  Dim after merging taxonomy: ", paste(dim(merged), collapse = "x"))
  if (nrow(merged) == 0) {
    stop("Merging taxonomy resulted in an empty data frame. Check 'FeatureID' column consistency in merged data and taxonomy file.")
  }
  if (!"taxon_function_abun" %in% colnames(merged)) {
    stop("Column 'taxon_function_abun' not found after merging. Please check your 'contrib_file' for this column.")
  }
  merged$taxon_function_abun <- as.numeric(merged$taxon_function_abun)
  
  # Clean unique classes (remove NAs)
  unique_classes_clean <- unique(merged$class)
  unique_classes_clean <- unique_classes_clean[!is.na(unique_classes_clean)]
  
  message("Unique classes identified for processing: ",
          if (length(unique_classes_clean) > 0) paste(sort(unique_classes_clean), collapse = ", ") else "None (or all are NA)")
  
  if (length(unique_classes_clean) == 0) {
    warning("No valid (non-NA) unique classes found in merged data. Output files will not be generated. Please check your 'metadata.csv' 'class' column and 'SampleID' matching.")
    return(invisible(NULL))
  }
  
  # Loop over each valid class
  for (current_class in unique_classes_clean) {
    message("\nProcessing class: '", current_class, "'")
    
    merged_class <- merged %>% filter(class == current_class)
    message("  Dim merged_class for '", current_class, "': ", paste(dim(merged_class), collapse = "x"))
    if (nrow(merged_class) == 0) {
      message("  No data for class '", current_class, "'. Skipping this class.")
      next
    }
    
    # Aggregate taxon-function abundance
    message("  Aggregating taxon-function abundance for '", current_class, "'...")
    if (!all(c("FunctionID", "TaxonID", "taxon_function_abun") %in% colnames(merged_class))) {
      warning("  Missing 'FunctionID', 'TaxonID', or 'taxon_function_abun' for class '", current_class, "'. Skipping aggregation.")
      next
    }
    taxon_function_total_class <- aggregate(
      taxon_function_abun ~ FunctionID + TaxonID,
      data = merged_class, sum, na.rm = TRUE
    )
    message("  Dim aggregated taxon-function data: ", paste(dim(taxon_function_total_class), collapse = "x"))
    if (nrow(taxon_function_total_class) == 0) {
      message("  No aggregated data (taxon-function abundance) for class '", current_class, "'. Skipping this class.")
      next
    }
    
    # Calculate total abundance per function and relative contribution
    message("  Calculating relative contributions for '", current_class, "'...")
    function_total_class <- aggregate(
      taxon_function_abun ~ FunctionID,
      data = taxon_function_total_class, sum, na.rm = TRUE
    )
    colnames(function_total_class)[2] <- "total_abundance_all_taxa"
    taxon_function_total_class <- merge(taxon_function_total_class, function_total_class, by = "FunctionID")
    taxon_function_total_class$relative_contribution <- with(taxon_function_total_class,
                                                             ifelse(total_abundance_all_taxa == 0, 0, taxon_function_abun / total_abundance_all_taxa)
    )
    message("  Dim after relative contribution calculation: ", paste(dim(taxon_function_total_class), collapse = "x"))
    
    # Apply filtering
    if (filtering != "unfiltered") {
      message("  Applying filtering: ", filtering, " for class '", current_class, "'...")
      if (filtering %in% c("mean", "median")) {
        threshold_df <- aggregate(relative_contribution ~ FunctionID, data = taxon_function_total_class,
                                  FUN = ifelse(filtering == "mean", mean, median), na.rm = TRUE)
        colnames(threshold_df)[2] <- "threshold"
        taxon_function_total_class <- merge(taxon_function_total_class, threshold_df, by = "FunctionID")
        taxon_function_total_class <- subset(taxon_function_total_class, relative_contribution >= threshold | is.na(threshold))
      } else {
        percent_map <- c("top10%" = 0.10, "top25%" = 0.25, "top50%" = 0.50, "top75%" = 0.75)
        top_percent <- percent_map[filtering]
        taxon_function_total_class <- taxon_function_total_class %>%
          group_by(FunctionID) %>% arrange(desc(relative_contribution)) %>%
          mutate(rank = row_number(), n_taxa = n(), cutoff = pmax(ceiling(top_percent * n_taxa), 1)) %>%
          filter(rank <= cutoff) %>% ungroup() %>% select(-rank, -n_taxa, -cutoff)
      }
      message("  Dim after filtering: ", paste(dim(taxon_function_total_class), collapse = "x"))
      if (nrow(taxon_function_total_class) == 0) {
        message("  No data remaining after filtering (", filtering, ") for class '", current_class, "'. Skipping output for this class.")
        next
      }
    }
    
    # Sort and write output
    message("  Saving results for class '", current_class, "'...")
    taxon_function_total_class <- taxon_function_total_class[order(taxon_function_total_class$FunctionID, -taxon_function_total_class$relative_contribution), ]
    file_suffix <- gsub("%", "", filtering)
    # Construct full_output_path directly using the provided output_file as the directory
    full_output_path <- file.path(output_file, paste0("microbe_pathway_network_", current_class, "_", file_suffix, ".csv"))
    write.csv(taxon_function_total_class, full_output_path, row.names = FALSE)
    message("  Saved results for class '", current_class, "' to: ", full_output_path)
  }
  message("Processing complete.")
  return(invisible(NULL))
}

# Example usage:
construct_microbe_pathway_network(
  contrib_file = "path_abun_contrib.csv",      
  metadata_file = "sample_metadata.csv",      
  taxonomy_file = "taxonomy.csv",      
  output_file = "microbe_pathway_network_results", # Output directory for results
  filtering = "median" # User can choose from "unfiltered", "mean", "median", "top10%", "top25%", "top50%", "top75%"
)
```

#### **Example output**

The output `microbe_pathway_network.csv` is a network table showing which microbes (taxa) contribute to which pathways, along with their contribution values.

| FunctionID | TaxonID | total_abundance | total_abundance_all_taxa | relative_contribution | median_contribution |
|------------|---------|------------------|----------------------------|------------------------|----------------------|
| ko00365    | g__Bilophila | 44.5 | 44.5 | 1.0 | 1.0 |
| ko00571    | g__Bifidobacterium | 1073.7 | 1076.3 | 0.998 | 0.5 |
| ko00720    | g__Blautia | 14084.3 | 16493.0 | 0.854 | 0.0055 |

Each row represents a weighted edge linking a microbial taxon (`TaxonID`) to a functional pathway (`FunctionID`) with the strength of the edge defined by the `relative_contribution`.
- `total_abundance`: The absolute contribution of a given taxon to a pathway summed across all samples.  
- `total_abundance_all_taxa`: The total combined contribution of all taxa to the same pathway, used as a baseline for normalization.
- `median_contribution` (and other thresholds like mean or top%) are shown in the table to indicate the filtering cutoff used for each pathway. This helps explain which taxa passed the filtering based on their relative contribution.

### <ins>Step 2: Pathway–pathway network construction</ins>

The pathway–pathway network is constructed using pathways identified as significant through Gene Set Enrichment Analysis (GSEA). Edges between pathways are defined based on shared genes, and Jaccard indices represent edge weights.

#### **Required inputs**

| File | Description |
|---|---|
| `pred_metagenome_unstrat.csv` | Gene abundance data from PICRUSt2. |
| `sample_metadata.csv` | Sample metadata with group or condition information. **Required columns**: `SampleID`, `class`. |
| `pathway_gene_map.csv` | Maps pathway IDs to their associated gene/KO IDs. **First column:** pathway ID; **other columns:** gene IDs. |

```r
library(DESeq2)
library(clusterProfiler)
library(dplyr)
library(tidyr)
set.seed(123) # Set seed for reproducibility

construct_pathway_pathway_network <- function(
  abundance_file,
  metadata_file,
  map_file,
  output_file,
  pvalueCutoff,
  pAdjustMethod = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none")
) {
  # Validate pAdjustMethod input
  pAdjustMethod <- match.arg(pAdjustMethod)
  message("Starting pathway-pathway network construction.")
  message("Using p-value cutoff: ", pvalueCutoff, " and p-adjustment method: ", pAdjustMethod, ".")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_file)) {
    dir.create(output_file, recursive = TRUE)
    message("Created output directory: ", output_file)
  } else {
    message("Output directory already exists: ", output_file)
  }
  
  # 1. Load abundance data from provided file
  message("1. Loading abundance data from: ", abundance_file)
  gene_abundance <- tryCatch(
    read.csv(abundance_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error loading abundance file '", abundance_file, "': ", e$message, sep = ""))
    }
  )
  message("   Loaded abundance data. Dimensions: ", paste(dim(gene_abundance), collapse = "x"))
  
  # 2. Load sample metadata
  message("2. Loading sample metadata from: ", metadata_file)
  metadata <- tryCatch(
    read.csv(metadata_file, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error loading metadata file '", metadata_file, "': ", e$message, sep = ""))
    }
  )
  message("   Loaded metadata. Dimensions: ", paste(dim(metadata), collapse = "x"))
  
  # Ensure the correct SampleID and class columns are present
  if (!"SampleID" %in% colnames(metadata)) {
    stop("The 'sample_metadata.csv' file must contain a 'SampleID' column.")
  }
  if (!"class" %in% colnames(metadata)) {
    stop("The 'sample_metadata.csv' file must contain a 'class' column for group definition.")
  }
  
  # 3. Filter metadata samples present in abundance data
  message("3. Filtering metadata to match samples in abundance data...")
  sample_ids <- colnames(gene_abundance)
  initial_metadata_rows <- nrow(metadata)
  metadata <- metadata %>% filter(SampleID %in% sample_ids)
  if (nrow(metadata) == 0) {
    stop("No matching samples found between abundance data and metadata after filtering.")
  }
  message("   Filtered metadata. Kept ", nrow(metadata), " out of ", initial_metadata_rows, " samples.")
  
  # 4. Set 'condition' factor from 'class' column
  message("4. Setting 'condition' factor from 'class' column and aligning data...")
  metadata$condition <- as.factor(metadata$class)
  rownames(metadata) <- metadata$SampleID # Use SampleID for row names
  
  # 5. Round abundance counts for DESeq2 compatibility and ensure sample order
  message("5. Rounding abundance counts and aligning sample order for DESeq2...")
  gene_abundance_rounded <- round(gene_abundance)
  gene_abundance_rounded <- gene_abundance_rounded[, rownames(metadata), drop = FALSE] # drop=FALSE to handle single sample case
  message("   Abundance data ready for DESeq2. Final dimensions: ", paste(dim(gene_abundance_rounded), collapse = "x"))
  
  # 6. Create DESeq2 dataset and run differential expression analysis
  message("6. Creating DESeq2 dataset and running differential expression analysis...")
  dds <- tryCatch(
    DESeqDataSetFromMatrix(countData = gene_abundance_rounded,
                           colData = metadata,
                           design = ~ condition),
    error = function(e) {
      stop(paste("Error creating DESeq2 dataset: ", e$message, ". Check abundance data (must be integers) and metadata consistency.", sep = ""))
    }
  )
  
  dds <- tryCatch(
    DESeq(dds),
    error = function(e) {
      stop(paste("Error running DESeq2 analysis: ", e$message, ". This might happen if groups have zero variance, or too few samples.", sep = ""))
    }
  )
  message("   DESeq2 analysis complete.")
  
  # 7. Get all pairwise condition comparisons
  conditions <- levels(metadata$condition)
  if (length(conditions) < 2) {
    stop("Less than two unique conditions found in metadata 'class' column. Cannot perform pairwise comparisons.")
  }
  comparisons <- combn(conditions, 2, simplify = FALSE)
  message("7. Identified ", length(comparisons), " pairwise comparisons: ",
          paste(sapply(comparisons, function(x) paste(x, collapse = " vs ")), collapse = ", "))
  
  # 8. Load pathway-to-gene mapping and reshape into TERM2GENE format
  message("8. Loading pathway-to-gene mapping from: ", map_file)
  map_raw <- tryCatch(
    read.csv(map_file, header = FALSE, fill = TRUE, stringsAsFactors = FALSE, skip = 1),
    error = function(e) {
      stop(paste("Error loading map file '", map_file, "': ", e$message, ". Ensure it's a valid CSV.", sep = ""))
    }
  )
  
  TERM2GENE <- map_raw %>%
    gather(key = "temp_col", value = "gene", -V1) %>% # V1 becomes 'term', others become 'gene'
    select(term = V1, gene) %>%
    filter(gene != "") %>% # Remove empty gene entries
    distinct() # Ensure unique pathway-gene pairs
  
  if (nrow(TERM2GENE) == 0) {
    stop("No valid pathway-gene mappings found after processing '", map_file, "'. Check file format.")
  }
  message("   Processed ", nrow(TERM2GENE), " unique pathway-gene mappings.")
  
  gsea_results_list <- list()
  
  # 9. Loop over each pairwise comparison to run GSEA
  message("9. Running GSEA for each pairwise comparison...")
  for (i in seq_along(comparisons)) {
    comp <- comparisons[[i]]
    cond1 <- comp[1]
    cond2 <- comp[2]
    comparison_name <- paste0(cond1, "_vs_", cond2)
    message("   Processing comparison (", i, "/", length(comparisons), "): ", cond2, " vs ", cond1)
    
    # Get DESeq2 results for contrast cond2 vs cond1
    res <- results(dds, contrast = c("condition", cond2, cond1))
    
    # Prepare ranked gene list: sign(log2FoldChange) * -log10(pvalue)
    ranked_df <- as.data.frame(res[, c("log2FoldChange", "pvalue")])
    ranked_df <- ranked_df[!is.na(ranked_df$log2FoldChange) & !is.na(ranked_df$pvalue), ]
    
    if (nrow(ranked_df) == 0) {
      warning("No valid log2FoldChange or pvalue for comparison ", comparison_name, ". Skipping GSEA.")
      next
    }
    
    # Handle cases where pvalue might be 0, leading to -log10(0) = Inf
    # A common approach is to set a minimum p-value
    min_pvalue_for_log <- min(ranked_df$pvalue[ranked_df$pvalue > 0], na.rm = TRUE) / 2
    ranked_df$pvalue[ranked_df$pvalue == 0] <- min_pvalue_for_log
    
    ranked_df$rank <- sign(ranked_df$log2FoldChange) * -log10(ranked_df$pvalue)
    ranked_df <- ranked_df[order(ranked_df$rank, decreasing = TRUE), ]
    geneList <- setNames(ranked_df$rank, rownames(ranked_df))
    
    # Run GSEA using clusterProfiler
    gsea_res <- tryCatch(
      GSEA(geneList = geneList,
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
      message("   No significant GSEA results found for ", comparison_name, ".")
      next
    }
    
    # Save GSEA results dataframe
    gsea_df <- as.data.frame(gsea_res)
    key <- comparison_name
    gsea_results_list[[key]] <- gsea_df
    
    # Construct full output path for GSEA results
    gsea_output_path <- file.path(output_file, paste0("gsea_results_", key, ".csv"))
    write.csv(gsea_df, gsea_output_path, row.names = FALSE)
    message("   Saved GSEA results for ", comparison_name, " to: ", gsea_output_path)
  }
  
  # 10. Compute Jaccard indices between pathways within each comparison's GSEA results
  message("10. Computing Jaccard indices for overlapping pathways...")
  jaccard_results_list <- list()
  
  if (length(gsea_results_list) == 0) {
    message("   No GSEA results to compute Jaccard indices. Skipping.")
  } else {
    for (key in names(gsea_results_list)) {
      message("   Calculating Jaccard index for comparison: ", key)
      gsea_df <- gsea_results_list[[key]]
      
      # Filter for pathways with core enrichment genes (i.e., not empty or NA)
      gsea_df_filtered <- gsea_df %>% filter(!is.na(core_enrichment) & core_enrichment != "")
      gene_sets <- strsplit(as.character(gsea_df_filtered$core_enrichment), "/")
      gene_sets <- lapply(gene_sets, function(x) unique(na.omit(x))) # Ensure unique genes per set and remove NAs
      
      n <- length(gene_sets)
      res_list <- list()
      
      if (n < 2) { # Ensure there are at least two pathways to compare
        message(paste("   Less than two significant pathways with core enrichment for comparison '", key, "'. Skipping Jaccard index calculation for this comparison.", sep = ""))
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
              pathway_1 = gsea_df_filtered$ID[i],
              pathway_2 = gsea_df_filtered$ID[j],
              jaccard_index = jaccard,
              comparison = key,
              stringsAsFactors = FALSE
            )
          }
        }
      }
      
      if (length(res_list) > 0) {
        jaccard_df <- do.call(rbind, res_list)
        # Construct full output path for Jaccard results
        jaccard_output_path <- file.path(output_file, paste0("pathway_jaccard_", key, ".csv"))
        write.csv(jaccard_df, jaccard_output_path, row.names = FALSE)
        message("   Saved Jaccard index results for ", key, " to: ", jaccard_output_path)
      } else {
        message("   No Jaccard indices > 0 found for comparison '", key, "'. Skipping saving file.")
      }
    }
  }
  
  message("Pathway-pathway network construction complete.")
}

# Example usage:
construct_pathway_pathway_network(
   abundance_file = "pred_metagenome_unstrat.csv", 
   metadata_file = "sample_metadata.csv",          
   map_file = "pathway_gene_map.csv",              
   output_file = "pathway_pathway_network_results", # Output directory for results       
   pvalueCutoff = 0.05, # User MUST specify this value, e.g., 0.05
   pAdjustMethod = "BH" # User can choose from "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
)
```
#### **Example output**

The function creates an output directory (e.g., `pathway_pathway_network_results`) containing `.csv` files for each pairwise group comparison (e.g., `G1_vs_G2`).

For each comparison, two types of files are generated:

1.  **`gsea_results_[group1]_vs_[group2].csv`**: Contains detailed Gene Set Enrichment Analysis (GSEA) results for pathways. Each row describes an enriched pathway, including its ID, description, enrichment score (NES), and adjusted p-value (`p.adjust`). The `core_enrichment` column lists the key genes driving the enrichment.

2.  **`pathway_jaccard_[group1]_vs_[group2].csv`**: Quantifies the similarity between *significant* pathways based on shared "core enriched" genes using the Jaccard index. This table defines the edges of the pathway-pathway network.

**Example table: `pathway_jaccard_G1_vs_G2.csv`**

| pathway_1 | pathway_2 | jaccard_index | comparison |
|:----------|:----------|:--------------|:-----------|
| ko00500   | ko00230   | 0.02173913    | G1_vs_G2   |
| ko00500   | ko00030   | 0.035714286   | G1_vs_G2   |
| ko00500   | ko00052   | 0.083333333   | G1_vs_G2   |
| ko00550   | ko00470   | 0.064516129   | G1_vs_G2   |

Each row represents a connection between two pathways (`pathway_1`, `pathway_2`). The `jaccard_index` (0-1) indicates the degree of shared genes between them; a higher value means more overlap and a stronger functional relationship. 

### <ins>Pathway–metabolite network construction</ins>

The pathway–metabolite network is constructed by calculating pairwise correlation (e.g., Spearman or Pearson) between pathway abundance and metabolite concentrations.  

### <ins>Multi-layered network</ins>

## Module 3: Network Analysis

This module provides three network analyses designed to identify context-specific associations:

- The hub identification uses the Maximal Clique Centrality (MCC) algorithm to identify key microbial pathways.  
- The pathfinding uses the Dijkstra's algorithm to identify the shortest path between the selected source and target nodes.  
- The node prioritization uses the Laplacian Heat Diffusion (LHD) algorithm to identify microbe-associated metabolites.
