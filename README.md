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
| `pathway_gene_map.tsv` | Maps pathway IDs to their associated gene/KO IDs. **First column:** pathway IDs; **other columns:** gene IDs. |

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

### <ins>STEP 1: Microbe–pathway network construction</ins>

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
  output_file, 
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
  message("Microbe-pathway network construction complete.")
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

### <ins>STEP 2: Pathway–pathway network construction</ins>

The pathway–pathway network is constructed using pathways identified as significant through Gene Set Enrichment Analysis (GSEA). Edges between pathways are defined based on shared genes, and Jaccard indices represent edge weights.

#### **Required inputs**

| File | Description |
|---|---|
| `pred_metagenome_unstrat.csv` | Gene abundance data from PICRUSt2. |
| `sample_metadata.csv` | Sample metadata with group or condition information. **Required columns**: `SampleID`, `class`. |
| `pathway_gene_map.csv` | Maps pathway IDs to their associated gene/KO IDs. **First column:** pathway IDs; **other columns:** gene IDs. |

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
  pAdjustMethod = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none"),
  rank_by = c("signed_log_pvalue", "log2FoldChange")
) {
  # Validate inputs
  pAdjustMethod <- match.arg(pAdjustMethod)
  rank_by <- match.arg(rank_by) # Validate new parameter
  
  message("Starting pathway-pathway network construction.")
  message("Using p-value cutoff: ", pvalueCutoff, " and p-adjustment method: ", pAdjustMethod, ".")
  message("Genes will be ranked by: ", rank_by, ".") # Status message for ranking method
  
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
      message("   Ranking by signed -log10(p-value) for ", comparison_name, ".")
    } else if (rank_by == "log2FoldChange") {
      ranked_df$rank <- ranked_df$log2FoldChange
      message("   Ranking by log2FoldChange for ", comparison_name, ".")
    }
    
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
              FunctionID_1 = gsea_df_filtered$ID[i], # Renamed from pathway_1
              FunctionID_2 = gsea_df_filtered$ID[j], # Renamed from pathway_2
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
   pAdjustMethod = "BH", # User can choose from "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
   rank_by = "signed_log_pvalue" # User can choose from "signed_log_pvalue" or "log2FoldChange"
)
```
#### **Example output**

The function creates an output directory (e.g., `pathway_pathway_network_results`) containing `.csv` files for each pairwise group comparison (e.g., `G1_vs_G2`).

For each comparison, two types of files are generated:

1.  **`gsea_results_[group1]_vs_[group2].csv`**: Contains detailed Gene Set Enrichment Analysis (GSEA) results for pathways. Each row describes an enriched pathway, including its ID, description, enrichment score (NES), and adjusted p-value (`p.adjust`). The `core_enrichment` column lists the key genes driving the enrichment.

2.  **`pathway_jaccard_[group1]_vs_[group2].csv`**: Quantifies the similarity between *significant* pathways based on shared "core enriched" genes using the Jaccard index. This table defines the edges of the pathway-pathway network.

**Example table: `pathway_jaccard_G1_vs_G2.csv`**

| FunctionID_1 | FunctionID_2 | jaccard_index | comparison |
|:----------|:----------|:--------------|:-----------|
| ko00500   | ko00230   | 0.02173913    | G1_vs_G2   |
| ko00500   | ko00030   | 0.035714286   | G1_vs_G2   |
| ko00500   | ko00052   | 0.083333333   | G1_vs_G2   |
| ko00550   | ko00470   | 0.064516129   | G1_vs_G2   |

Each row represents a connection between two pathways (`pathway_1`, `pathway_2`). The `jaccard_index` (0-1) indicates the degree of shared genes between them; a higher value means more overlap and a stronger functional relationship. 

### <ins>Pathway–metabolite network construction</ins>

The pathway–metabolite network is constructed by calculating pairwise correlation (e.g., Spearman or Pearson) between pathway abundance and metabolite concentrations.  

#### **Required inputs**

| File                        | Description                                                                                                                              |
|:----------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------|
| `path_abun_unstrat.csv`     | Pathway abundance data. First column: `SampleID`; other columns: `Function IDs`. Values will be converted to relative abundance within the function. |
| `metabolite_concentration.csv` | Metabolite concentration data. First column: `SampleID`; other columns: metabolite names.                                                |
| `sample_metadata.csv`       | (Optional) Sample metadata with group or condition information. Required columns: `SampleID`, `class`. If not provided or class column is missing, correlations will be performed on the overall dataset. |

```r
library(dplyr)
library(tidyr)

construct_pathway_metabolite_network <- function(
  pathway_abundance_file,
  metabolite_concentration_file,
  output_file,
  metadata_file,
  correlation_method = c("spearman", "pearson"),
  filter_by = c("none", "p_value", "q_value"),
  corr_cutoff,
  p_value_cutoff,
  q_value_cutoff,
  q_adjust_method = c("bonferroni", "fdr")
) {
  # Validate inputs
  correlation_method <- match.arg(correlation_method)
  filter_by <- match.arg(filter_by)
  q_adjust_method <- match.arg(q_adjust_method)

  message("Starting pathway-metabolite network construction.")
  message("Correlation method: ", correlation_method)
  message("Filtering results by: ", filter_by)

  if (filter_by == "p_value" && is.null(p_value_cutoff)) {
    stop("Error: 'p_value_cutoff' must be specified if 'filter_by' is 'p_value'.")
  }
  if (filter_by == "q_value" && is.null(q_value_cutoff)) {
    stop("Error: 'q_value_cutoff' must be specified if 'filter_by' is 'q_value'.")
  }

  message("Absolute correlation coefficient cutoff: ", corr_cutoff)
  if (filter_by == "p_value") message("P-value cutoff: ", p_value_cutoff)
  if (filter_by == "q_value") message("Q-value cutoff: ", q_value_cutoff, " (using ", q_adjust_method, " correction)")

  # Create output directory if it doesn't exist
  if (!dir.exists(output_file)) {
    dir.create(output_file, recursive = TRUE)
    message("Created output directory: ", output_file)
  } else {
    message("Output directory already exists: ", output_file)
  }

  # 1. Load data
  message("1. Loading pathway abundance and metabolite concentration data...")
  pathway_abun <- tryCatch(
    read.csv(pathway_abundance_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error loading pathway abundance file '", pathway_abundance_file, "': ", e$message, sep = ""))
    }
  )
  metabolite_conc <- tryCatch(
    read.csv(metabolite_concentration_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error loading metabolite concentration file '", metabolite_concentration_file, "': ", e$message, sep = ""))
    }
  )
  message("   Loaded pathway abundance. Dimensions: ", paste(dim(pathway_abun), collapse = "x"))
  message("   Loaded metabolite concentration. Dimensions: ", paste(dim(metabolite_conc), collapse = "x"))

  # 2. Load metadata if provided
  metadata <- NULL
  if (!is.null(metadata_file)) {
    message("2. Loading sample metadata from: ", metadata_file)
    metadata <- tryCatch(
      read.csv(metadata_file, header = TRUE, stringsAsFactors = FALSE),
      error = function(e) {
        warning(paste("Warning: Error loading metadata file '", metadata_file, "': ", e$message, ". Proceeding without group-specific analysis.", sep = ""))
        return(NULL)
      }
    )
    if (!is.null(metadata)) {
      message("   Loaded metadata. Dimensions: ", paste(dim(metadata), collapse = "x"))
      if (!"SampleID" %in% colnames(metadata)) {
        warning("Metadata file missing 'SampleID' column. Proceeding without group-specific analysis.")
        metadata <- NULL
      } else if (!"class" %in% colnames(metadata)) {
        warning("Metadata file missing 'class' column. Proceeding without group-specific analysis.")
        metadata <- NULL
      } else {
        metadata$SampleID <- as.character(metadata$SampleID)
        metadata$class <- as.factor(metadata$class)
        rownames(metadata) <- metadata$SampleID
      }
    }
  } else {
    message("2. No metadata file provided. Performing correlation on overall dataset.")
  }

  # 3. Data preparation: Ensure SampleIDs align and convert pathway abundance to relative values
  message("3. Preparing data for correlation...")

  # Merge data to get common samples and align
  common_samples <- intersect(rownames(pathway_abun), rownames(metabolite_conc))
  if (!is.null(metadata)) {
    common_samples <- intersect(common_samples, rownames(metadata))
  }

  if (length(common_samples) == 0) {
    stop("No common samples found between pathway abundance, metabolite concentration, and metadata (if provided).")
  }
  message("   Number of common samples: ", length(common_samples))

  pathway_abun_filtered <- pathway_abun[common_samples, , drop = FALSE]
  metabolite_conc_filtered <- metabolite_conc[common_samples, , drop = FALSE]
  if (!is.null(metadata)) {
    metadata_filtered <- metadata[common_samples, , drop = FALSE]
  } else {
    metadata_filtered <- data.frame(SampleID = common_samples, class = "overall", row.names = common_samples)
  }

  # Convert pathway abundance to relative values (row-wise normalization)
  # Sum of all pathway abundances for each sample
  pathway_abun_sum_per_sample <- rowSums(pathway_abun_filtered, na.rm = TRUE)
  # Replace 0 sums with 1 to avoid NaN/Inf for division, but keep values 0 in this case
  pathway_abun_relative <- pathway_abun_filtered / ifelse(pathway_abun_sum_per_sample == 0, 1, pathway_abun_sum_per_sample)
  # For rows that originally summed to zero, they should remain zero after division by 1.
  # Ensure if original value was NA, it remains NA.
  pathway_abun_relative[is.na(pathway_abun_filtered)] <- NA
  message("   Pathway abundance converted to relative values (sum per sample).")

  # Ensure all data are numeric
  pathway_abun_relative[] <- lapply(pathway_abun_relative, as.numeric)
  metabolite_conc_filtered[] <- lapply(metabolite_conc_filtered, as.numeric)

  # Prepare groups for iteration
  groups_to_process <- unique(metadata_filtered$class)
  all_correlation_results <- list()

  # 4. Loop over each group to perform correlations
  message("4. Performing correlations for each group...")
  for (current_group in groups_to_process) {
    message("   Processing group: '", current_group, "'")

    samples_in_group <- rownames(metadata_filtered[metadata_filtered$class == current_group, , drop = FALSE])

    path_data_group <- pathway_abun_relative[samples_in_group, , drop = FALSE]
    met_data_group <- metabolite_conc_filtered[samples_in_group, , drop = FALSE]

    # Check for sufficient samples for correlation
    # Filter out columns that are all NA or have zero variance within the group
    path_data_group_clean <- path_data_group[, colSums(is.na(path_data_group)) != nrow(path_data_group) & apply(path_data_group, 2, var, na.rm=TRUE) != 0, drop = FALSE]
    met_data_group_clean <- met_data_group[, colSums(is.na(met_data_group)) != nrow(met_data_group) & apply(met_data_group, 2, var, na.rm=TRUE) != 0, drop = FALSE]

    if (ncol(path_data_group_clean) == 0 || ncol(met_data_group_clean) == 0 || nrow(path_data_group_clean) < 3) {
      warning("   Not enough valid (non-constant, non-NA) data or samples (less than 3) in group '", current_group, "' to perform meaningful correlations. Skipping.")
      next
    }

    # Calculate correlation matrix using stats::cor
    corr_matrix <- cor(path_data_group_clean, met_data_group_clean, method = correlation_method, use = "pairwise.complete.obs")

    # Initialize p-value matrix
    p_matrix <- matrix(NA, nrow = nrow(corr_matrix), ncol = ncol(corr_matrix),
                       dimnames = dimnames(corr_matrix))

    # Calculate p-values for each correlation coefficient
    # Using t-distribution approximation for both Pearson and Spearman
    for (i in 1:nrow(corr_matrix)) { # Loop through pathways
      for (j in 1:ncol(corr_matrix)) { # Loop through metabolites
        r_val <- corr_matrix[i, j]
        if (!is.na(r_val)) {
          # Get actual number of complete observations for this specific pair
          valid_pairs_count <- sum(complete.cases(path_data_group_clean[, rownames(corr_matrix)[i]],
                                                  met_data_group_clean[, colnames(corr_matrix)[j]]))

          if (valid_pairs_count >= 3) {
            t_statistic <- r_val * sqrt((valid_pairs_count - 2) / (1 - r_val^2))
            p_matrix[i, j] <- 2 * stats::pt(abs(t_statistic), df = valid_pairs_count - 2, lower.tail = FALSE)
          }
        }
      }
    }


    # Reshape results into long format
    corr_long <- as.data.frame(corr_matrix) %>%
      tibble::rownames_to_column(var = "FunctionID") %>% # Changed to FunctionID
      pivot_longer(cols = -FunctionID, names_to = "MetaboliteID", values_to = "Correlation")

    p_long <- as.data.frame(p_matrix) %>%
      tibble::rownames_to_column(var = "FunctionID") %>% # Changed to FunctionID
      pivot_longer(cols = -FunctionID, names_to = "MetaboliteID", values_to = "P_value")

    # Combine correlation and p-value data
    combined_results <- left_join(corr_long, p_long, by = c("FunctionID", "MetaboliteID")) # Changed to FunctionID
    
    # Remove rows where correlation or p-value is NA (due to insufficient valid pairs)
    combined_results <- combined_results %>% filter(!is.na(Correlation) & !is.na(P_value))

    # 5. Apply filtering based on user choice
    message("   Applying filters for group '", current_group, "'...")

    # Filter by absolute correlation coefficient first
    combined_results_filtered <- combined_results %>%
      filter(abs(Correlation) >= corr_cutoff)

    if (filter_by == "p_value") {
      combined_results_filtered <- combined_results_filtered %>%
        filter(P_value <= p_value_cutoff)
      message("   Filtered by p-value <= ", p_value_cutoff)
    } else if (filter_by == "q_value") {
      # Calculate Q-values
      combined_results_filtered$Q_value <- stats::p.adjust(combined_results_filtered$P_value, method = q_adjust_method)
      combined_results_filtered <- combined_results_filtered %>%
        filter(Q_value <= q_value_cutoff)
      message("   Filtered by q-value <= ", q_value_cutoff, " (", q_adjust_method, " correction)")
    }

    if (nrow(combined_results_filtered) == 0) {
      message("   No significant correlations found after filtering for group '", current_group, "'. Skipping output for this group.")
      next
    }

    combined_results_filtered$Group <- current_group
    all_correlation_results[[current_group]] <- combined_results_filtered

    # Save results for the current group
    output_filename <- paste0("pathway_metabolite_network_", current_group, ".csv")
    output_filepath <- file.path(output_file, output_filename)
    write.csv(combined_results_filtered, output_filepath, row.names = FALSE)
    message("   Saved results for group '", current_group, "' to: ", output_filepath)
  }
  message("Pathway-metabolite network construction complete.")
}

# Example usage:
construct_pathway_metabolite_network(
  pathway_abundance_file = "path_abun_unstrat.csv", 
  metabolite_concentration_file = "metabolite_concentration.csv", 
  output_file = "pathway_metabolite_network_results", # Output directory for results
  metadata_file = "sample_metadata.csv", # Optional, set to NULL if no groups
  correlation_method = "pearson", # Choose "spearman" or "pearson"
  filter_by = "none", # Choose "none", "p_value", or "q_value"
  corr_cutoff = 0.5, # Absolute correlation coefficient cutoff (e.g., 0.5)
  p_value_cutoff = NULL, # Set if filter_by = "p_value"
  q_value_cutoff = NULL, # Set if filter_by = "q_value"
  q_adjust_method = "fdr" # Choose "bonferroni" or "fdr" if filter_by = "q_value"
)
```
#### **Example output**

The function creates an output directory (e.g., `pathway_metabolite_network_results`) containing `.csv` files for each group analyzed (e.g., `pathway_metabolite_network_G1.csv`, `pathway_metabolite_network_G2.csv`, or `pathway_metabolite_network_overall.csv` if no groups are defined).

Each output file is a table representing the pathway-metabolite network edges for that specific group.

**Example table: `pathway_metabolite_network_G1.csv`**

| FunctionID | MetaboliteID | Correlation | P_value | Q_value | Group |
|:-----------|:-------------|:------------|:--------|:--------|:------|
| ko00010    | M_glucose    | 0.82        | 0.001   | 0.005   | G1  |
| ko00020    | M_amino_acid | -0.75       | 0.003   | 0.008   | G1  |
| ko00300    | M_lipid      | 0.68        | 0.01    | 0.02    | G1  |
| ko00400    | M_vitB       | -0.62       | 0.025   | 0.04    | G1  |

Each row represents a correlation (edge) between a pathway (`FunctionID`) and a metabolite (`MetaboliteID`).
- `Correlation`: The correlation coefficient (Spearman or Pearson) indicating the strength and direction of the relationship.
- `P_value`: The raw p-value for the correlation.
- `Q_value`: The adjusted p-value (q-value), if q-value filtering was selected.
- `Group`: The specific group for which the correlation was calculated.

### <ins>Multi-layered network</ins>

## Module 3: Network Analysis

This module provides three network analyses designed to identify context-specific associations:

- The hub identification uses the Maximal Clique Centrality (MCC) algorithm to identify key microbial pathways.  
- The pathfinding uses the Dijkstra's algorithm to identify the shortest path between the selected source and target nodes.  
- The node prioritization uses the Laplacian Heat Diffusion (LHD) algorithm to identify microbe-associated metabolites.
