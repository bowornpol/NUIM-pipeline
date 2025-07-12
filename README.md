# NUIM: Network-based utility for integrating microbiome and metabolome data

We developed NUIM, a modular, network-based framework for integrating microbiome and metabolome data systematically. NUIM consists of three modules: (1) data preparation and processing, (2) network construction, and (3) network analysis. It provides a wide variety of network analyses to identify context-specific associations. These include hub identification, pathfinding, and node prioritization.

<p align="center">
  <img src="figures/NUIM_overview.png" width="500"/>
</p>

#### **Required tools and R packages**

The NUIM pipeline depends on a set of external bioinformatics tools and an R-based environment for network-based analyses. The versions of these tools and packages may vary depending on the user. The following table details the requirements: 

| Category          | Tool / R package         |
| :---------------- | :----------------------- |
| **External tools**| `QIIME2` – for microbiome data processing. |
|                   | `PICRUSt2` – for functional prediction.     |
|                   | `Cytoscape` + `cytoHubba` plugin – for hub identification. | 
| **R packages**    | `dplyr` – for data manipulation.              |
|                   | `tidyr` – for data reshaping.                 |
|                   | `stats` – for statistical testing.            |
|                   | `tibble` – for tidy data frame support.       |
|                   | `stringr` – for string operations.            |
|                   | `tools` – for file and object manipulation.   |
|                   | `DESeq2` – for differential expression analysis. |
|                   | `clusterProfiler` – for GSEA analysis.                 | 
|                   | `igraph` – for handling network data and pathfinding analysis. |
|                   | `expm` – for matrix exponential in node prioritization analysis. | 
|                   | `ggplot2` – for data visualization.           |

## Module 1: Data preparation and processing

This module defines the procedures required to prepare and process the input data for downstream network construction.

- Input data includes `microbial sequencing reads in FASTQ format` and `metabolite concentration table`.  
- Microbiome data processing involves the use of QIIME2 to generate a feature table and representative sequences. These outputs are subsequently processed using PICRUSt2 for functional prediction, yielding gene abundance, pathway abundance, and pathway contribution data.  
- Although metabolome data processing may vary depending on user preference and experimental design, NUIM assumes that metabolite concentrations have been appropriately processed by standard practice. For example, users may employ established platforms such as [Metabox](https://metsysbio.com/tools_protocols/metabox-2-0) or [MetaboAnalyst](https://www.metaboanalyst.ca/MetaboAnalyst/docs/RTutorial.xhtml#2.2%20Data%20Processing%20and%20Statistical%20Analysis) to perform metabolomics data processing.

### <ins>QIIME2 workflow</ins>

This section provides a general QIIME2 workflow for processing paired-end 16S rRNA sequencing data. The goal is to generate a feature table and representative sequences for PICRUSt2.

#### **Required inputs**

| File            | Description                         |
|-----------------|-----------------------------------|
| `FASTQ`         | 16S rRNA gene sequencing reads. |
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
| `pathway_gene_map.tsv` | Maps pathway IDs to their associated gene/KO IDs. **First column**: pathway IDs; **Other columns**: gene/KO IDs. |

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
  -o picrust2_out

# Optional: Add --max_nsti to filter low-confidence ASVs
# e.g., --max_nsti 0.5

# Infer pathway abundance and pathway contribution
# Users can choose to use either KEGG or MetaCyc pathways

# Generate pathway abundance data (sample × pathway)
pathway_pipeline.py \
  -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz \
  -o pathway_abundance \
  --no_regroup \
  --map pathway_gene_map.tsv

# Generate pathway contribution data (ASV × pathway)
pathway_pipeline.py \
  -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz \
  -o pathway_contrib \
  --per_sequence_contrib \
  --per_sequence_abun KO_metagenome_out/seqtab_norm.tsv.gz \
  --per_sequence_function KO_predicted.tsv.gz \
  --no_regroup \
  --map pathway_gene_map.tsv

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

This module constructs a tripartite network linking microbial taxa, functional pathways, and metabolites, based on the processed data from **Module 1**. Follow the steps below to construct each network layer using the provided R functions:

### <ins>STEP 1: Microbe–pathway network construction</ins>

The microbe–pathway network is constructed from pathway contribution data, with edges representing the relative contribution of each microbe to specific pathways.

#### **Required inputs**

| File                    | Description                                                                                                                                           | Required columns                               |
| :---------------------- | :---------------------------------------------------------------------------------------------------------------------------------------------------- | :--------------------------------------------- |
| `path_abun_contrib.csv` | Pathway contribution data from PICRUSt2.                                                                                                              | `SampleID`, `FeatureID`, `FunctionID`, `taxon_function_abun` |
| `sample_metadata.csv`   | Sample metadata. Required for group-specific analysis. If not provided or `class` column is missing, data will be processed as one 'overall' group. | `SampleID`, `class`                            |
| `taxonomy.csv`          | Taxonomy annotations mapping `FeatureID` to taxonomic name (`TaxonID`) from QIIME2.                                                                               | `FeatureID`, `TaxonID`                         |

<details>
<summary>Click to show the full R function</summary>

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

  # Create output directory if it doesn't exist
  if (!dir.exists(output_file)) {
    dir.create(output_file, recursive = TRUE)
    message("Created output directory: ", output_file)
  } else {
    message("Output directory already exists: ", output_file)
  }
  
  message("Loading data...")
  contrib <- tryCatch(
    read.csv(contrib_file, stringsAsFactors = FALSE),
    error = function(e) stop(paste("Error loading contribution file: ", e$message))
  )
  metadata <- tryCatch(
    read.csv(metadata_file, stringsAsFactors = FALSE),
    error = function(e) stop(paste("Error loading metadata file: ", e$message))
  )
  taxonomy <- tryCatch(
    read.csv(taxonomy_file, stringsAsFactors = FALSE),
    error = function(e) stop(paste("Error loading taxonomy file: ", e$message))
  )
  
  message("Dim contrib: ", paste(dim(contrib), collapse = "x"), " | Cols: ", paste(colnames(contrib), collapse = ", "))
  message("Dim metadata: ", paste(dim(metadata), collapse = "x"), " | Cols: ", paste(colnames(metadata), collapse = ", "))
  message("Dim taxonomy: ", paste(dim(taxonomy), collapse = "x"), " | Cols: ", paste(colnames(taxonomy), collapse = ", "))
  
  if (!"class" %in% colnames(metadata)) {
    message("'class' column not found in metadata. Creating a default 'all' class.")
    metadata$class <- "all"
  } else {
    message("'class' column found in metadata.")
    metadata$class <- as.character(metadata$class)
  }
  
  contrib$SampleID <- as.character(contrib$SampleID)
  metadata$SampleID <- as.character(metadata$SampleID)
  contrib$FeatureID <- as.character(contrib$FeatureID)
  taxonomy$FeatureID <- as.character(taxonomy$FeatureID)
  contrib$FunctionID <- as.character(contrib$FunctionID)
  taxonomy$TaxonID <- as.character(taxonomy$TaxonID)
  
  message("Merging contribution and metadata...")
  merged <- merge(contrib, metadata, by = "SampleID", all.x = TRUE)
  message("  Dim after merging contrib and metadata: ", paste(dim(merged), collapse = "x"))
  if (nrow(merged) == 0) stop("Merging contribution and metadata resulted in an empty data frame.")
  
  cols_to_keep_from_taxonomy <- c("FeatureID", "TaxonID")
  if (!all(cols_to_keep_from_taxonomy %in% colnames(taxonomy))) {
    stop("Missing expected columns in taxonomy file: ", paste(setdiff(cols_to_keep_from_taxonomy, colnames(taxonomy)), collapse = ", "))
  }
  taxonomy_for_merge <- taxonomy %>% select(all_of(cols_to_keep_from_taxonomy))
  
  message("Merging with taxonomy data...")
  merged <- merge(merged, taxonomy_for_merge, by = "FeatureID", all.x = TRUE)
  message("  Dim after merging taxonomy: ", paste(dim(merged), collapse = "x"))
  if (!"taxon_function_abun" %in% colnames(merged)) {
    stop("Column 'taxon_function_abun' not found after merging.")
  }
  merged$taxon_function_abun <- as.numeric(merged$taxon_function_abun)
  
  unique_classes_clean <- unique(merged$class)
  unique_classes_clean <- unique_classes_clean[!is.na(unique_classes_clean)]
  message("Unique classes identified for processing: ",
          if (length(unique_classes_clean) > 0) paste(sort(unique_classes_clean), collapse = ", ") else "None")
  
  if (length(unique_classes_clean) == 0) {
    warning("No valid (non-NA) classes found. Exiting.")
    return(invisible(NULL))
  }
  
  for (current_class in unique_classes_clean) {
    message("\nProcessing class: '", current_class, "'")
    merged_class <- merged %>% filter(class == current_class)
    message("  Dim merged_class: ", paste(dim(merged_class), collapse = "x"))
    if (nrow(merged_class) == 0) next
    
    message("  Aggregating taxon-function abundance...")
    taxon_function_total_class <- aggregate(
      taxon_function_abun ~ FunctionID + TaxonID,
      data = merged_class, sum, na.rm = TRUE
    )
    message("  Dim aggregated data: ", paste(dim(taxon_function_total_class), collapse = "x"))
    if (nrow(taxon_function_total_class) == 0) next
    
    message("  Calculating relative contributions...")
    function_total_class <- aggregate(
      taxon_function_abun ~ FunctionID,
      data = taxon_function_total_class, sum, na.rm = TRUE
    )
    colnames(function_total_class)[2] <- "total_abundance_all_taxa"
    taxon_function_total_class <- merge(taxon_function_total_class, function_total_class, by = "FunctionID")
    taxon_function_total_class$relative_contribution <- with(taxon_function_total_class,
                                                             ifelse(total_abundance_all_taxa == 0, 0, taxon_function_abun / total_abundance_all_taxa)
    )
    message("  Dim after contribution calc: ", paste(dim(taxon_function_total_class), collapse = "x"))
    
    if (filtering != "unfiltered") {
      message("  Applying filtering: ", filtering)
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
        top_percent <- percent_map[filtering]
        taxon_function_total_class <- taxon_function_total_class %>%
          group_by(FunctionID) %>% arrange(desc(relative_contribution)) %>%
          mutate(rank = row_number(), n_taxa = n(), cutoff = pmax(ceiling(top_percent * n_taxa), 1)) %>%
          filter(rank <= cutoff) %>%
          ungroup() %>% select(-rank, -n_taxa, -cutoff)
      }
      message("  Dim after filtering: ", paste(dim(taxon_function_total_class), collapse = "x"))
      if (nrow(taxon_function_total_class) == 0) {
        message("  No data left after filtering. Skipping.")
        next
      }
    }
    
    message("  Saving results for class: ", current_class)
    taxon_function_total_class <- taxon_function_total_class[
      order(taxon_function_total_class$FunctionID, -taxon_function_total_class$relative_contribution),
    ]
    file_suffix <- gsub("%", "", filtering)
    full_output_path <- file.path(output_file, paste0("microbe_pathway_network_", current_class, "_", file_suffix, ".csv"))
    write.csv(taxon_function_total_class, full_output_path, row.names = FALSE)
    message("  File saved to: ", full_output_path)
  }
  
  message("Microbe-pathway network construction complete.")
  return(invisible(NULL))
}
```
</details> 

#### **Example usage**

```r
source("https://raw.githubusercontent.com/bowornpol/NUIM-pipeline/main/code/construct_microbe_pathway_network.R")

construct_microbe_pathway_network(
  contrib_file = "path_abun_contrib.csv",      
  metadata_file = "sample_metadata.csv",      
  taxonomy_file = "taxonomy.csv",      
  output_file = "microbe_pathway_network_results", # Output directory for results
  filtering = "median" # User can choose from "unfiltered", "mean", "median", "top10%", "top25%", "top50%", "top75%"
)
```

#### **Example output**

The function generates a directory (e.g., `microbe_pathway_network_results`), which contains output files specific to the group analyzed and the filtering method applied (e.g., `microbe_pathway_network_[class]_[filtering].csv`).

**Example table: `microbe_pathway_network_G2_median.csv`**

| FunctionID | TaxonID            | total_abundance | total_abundance_all_taxa | relative_contribution | median_contribution |
|:-----------|:-------------------|:----------------|:-------------------------|:----------------------|:--------------------|
| ko00365    | g__Bilophila       | 44.5            | 44.5                     | 1.0                   | 1.0                 |
| ko00571    | g__Bifidobacterium | 1073.7          | 1076.3                   | 0.97                 | 0.5                 |
| ko00720    | g__Blautia         | 14084.3         | 16493.0                  | 0.85                 | 0.05              |
| ...        | ...                | ...             | ...                      | ...                   | ...                 |

Each row represents a weighted edge linking a microbial taxon (`TaxonID`) to a functional pathway (`FunctionID`) with the strength of the edge defined by the `relative_contribution`.
- `total_abundance`: The absolute contribution of a given taxon to a pathway summed across all samples.  
- `total_abundance_all_taxa`: The total combined contribution of all taxa to the same pathway, used as a baseline for normalization.
- `median_contribution` (and other thresholds like mean or top%) are shown in the table to indicate the filtering cutoff used for each pathway. This helps explain which taxa passed the filtering based on their relative contribution.

### <ins>STEP 2: Pathway–pathway network construction</ins>

The pathway–pathway network is constructed using pathways identified as significant through Gene Set Enrichment Analysis (GSEA). Edges between pathways are defined based on shared genes, and Jaccard indices represent edge weights.

#### **Required inputs**

| File                        | Description                                                                 | Required Columns                                           |
| :-------------------------- | :-------------------------------------------------------------------------- | :--------------------------------------------------------- |
| `pred_metagenome_unstrat.csv` | Gene abundance data from PICRUSt2.                                          | `SampleID`, Gene/KO IDs (as columns)                       |
| `sample_metadata.csv`       | Sample metadata with group or condition information.                        | `SampleID`, `class`                                        |
| `pathway_gene_map.csv`      | Maps pathway IDs to their associated gene/KO IDs.                           | **First column**: pathway IDs; **Other columns**: gene/KO IDs |

<details>
<summary>Click to show the full R function</summary>

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
  message("    Loaded abundance data. Dimensions: ", paste(dim(gene_abundance), collapse = "x"))
  
  # 2. Load sample metadata
  message("2. Loading sample metadata from: ", metadata_file)
  metadata <- tryCatch(
    read.csv(metadata_file, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error loading metadata file '", metadata_file, "': ", e$message, sep = ""))
    }
  )
  message("    Loaded metadata. Dimensions: ", paste(dim(metadata), collapse = "x"))
  
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
  message("    Filtered metadata. Kept ", nrow(metadata), " out of ", initial_metadata_rows, " samples.")
  
  # 4. Set 'condition' factor from 'class' column
  message("4. Setting 'condition' factor from 'class' column and aligning data...")
  metadata$condition <- as.factor(metadata$class)
  rownames(metadata) <- metadata$SampleID # Use SampleID for row names
  
  # 5. Round abundance counts for DESeq2 compatibility and ensure sample order
  message("5. Rounding abundance counts and aligning sample order for DESeq2...")
  gene_abundance_rounded <- round(gene_abundance)
  gene_abundance_rounded <- gene_abundance_rounded[, rownames(metadata), drop = FALSE] # drop=FALSE to handle single sample case
  message("    Abundance data ready for DESeq2. Final dimensions: ", paste(dim(gene_abundance_rounded), collapse = "x"))
  
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
  message("    DESeq2 analysis complete.")
  
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
  message("    Processed ", nrow(TERM2GENE), " unique pathway-gene mappings.")
  
  gsea_results_list <- list()
  
  # 9. Loop over each pairwise comparison to run GSEA
  message("9. Running GSEA for each pairwise comparison...")
  for (i in seq_along(comparisons)) {
    comp <- comparisons[[i]]
    cond1 <- comp[1]
    cond2 <- comp[2]
    comparison_name <- paste0(cond1, "_vs_", cond2)
    message("    Processing comparison (", i, "/", length(comparisons), "): ", cond2, " vs ", cond1)
    
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
      message("    Ranking by signed -log10(p-value) for ", comparison_name, ".")
    } else if (rank_by == "log2FoldChange") {
      ranked_df$rank <- ranked_df$log2FoldChange
      message("    Ranking by log2FoldChange for ", comparison_name, ".")
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
      message("    No significant GSEA results found for ", comparison_name, ".")
      next
    }
    
    # Save GSEA results dataframe
    gsea_df <- as.data.frame(gsea_res)
    key <- comparison_name
    gsea_results_list[[key]] <- gsea_df
    
    # --- Construct full output path for GSEA results with parameters ---
    
    # Format parameters for filename
    pvalueCutoff_fname <- as.character(pvalueCutoff) # Convert numeric to string
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
    message("    Saved GSEA results for ", comparison_name, " to: ", gsea_output_path)
  }
  
  # 10. Compute Jaccard indices between pathways within each comparison's GSEA results
  message("10. Computing Jaccard indices for overlapping pathways...")
  jaccard_results_list <- list()
  
  if (length(gsea_results_list) == 0) {
    message("    No GSEA results to compute Jaccard indices. Skipping.")
  } else {
    for (key in names(gsea_results_list)) {
      message("    Calculating Jaccard index for comparison: ", key)
      gsea_df <- gsea_results_list[[key]]
      
      # Filter for pathways with core enrichment genes (i.e., not empty or NA)
      gsea_df_filtered <- gsea_df %>% filter(!is.na(core_enrichment) & core_enrichment != "")
      gene_sets <- strsplit(as.character(gsea_df_filtered$core_enrichment), "/")
      gene_sets <- lapply(gene_sets, function(x) unique(na.omit(x))) # Ensure unique genes per set and remove NAs
      
      n <- length(gene_sets)
      res_list <- list()
      
      if (n < 2) { # Ensure there are at least two pathways to compare
        message(paste("    Less than two significant pathways with core enrichment for comparison '", key, "'. Skipping Jaccard index calculation for this comparison.", sep = ""))
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
        # Construct full output path for Jaccard results (unchanged from previous version)
        jaccard_output_path <- file.path(output_file, paste0("pathway_jaccard_", key, ".csv"))
        write.csv(jaccard_df, jaccard_output_path, row.names = FALSE)
        message("    Saved Jaccard index results for ", key, " to: ", jaccard_output_path)
      } else {
        message("    No Jaccard indices > 0 found for comparison '", key, "'. Skipping saving file.")
      }
    }
  }
  
  message("Pathway-pathway network construction complete.")
  return(invisible(NULL))
}
```

</details>

#### **Example usage**

```r
source("https://raw.githubusercontent.com/bowornpol/NUIM-pipeline/main/code/construct_pathway_pathway_network.R")

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

1.  **`gsea_results_[class]_vs_[class]_[pvalueCutoff]_[pAdjustMethod]_[rank_by].csv`**: Contains detailed Gene Set Enrichment Analysis (GSEA) results for pathways. Each row describes an enriched pathway, including its ID, description, enrichment score (NES), and adjusted p-value (`p.adjust`). The `core_enrichment` column lists the key genes driving the enrichment.

2.  **`pathway_jaccard_[class]_vs_[class].csv`**: Quantifies the similarity between *significant* pathways based on shared "core enriched" genes using the Jaccard index. This table defines the edges of the pathway-pathway network.

**Example table: `pathway_jaccard_G1_vs_G2.csv`**

| FunctionID_1 | FunctionID_2 | jaccard_index | comparison |
|:----------|:----------|:--------------|:-----------|
| ko00500    | ko00230    | 0.021         | G1_vs_G2   |
| ko00500    | ko00030    | 0.035         | G1_vs_G2   |
| ko00500    | ko00052    | 0.083         | G1_vs_G2   |
| ko00550    | ko00470    | 0.064         | G1_vs_G2   |
| ...        | ...        | ...           | ...        |

Each row represents a connection between two pathways (`pathway_1`, `pathway_2`). 
- `jaccard_index`: A value from 0 to 1 indicating the degree of shared genes between the two pathways; a higher value means more overlap and a stronger functional relationship.
- `comparison`: Specifies the pairwise group comparison for which this Jaccard index was calculated. 

### <ins>STEP 3: Pathway–metabolite network construction</ins>

The pathway–metabolite network is constructed by calculating pairwise correlation (e.g., Spearman or Pearson) between pathway abundance and metabolite concentrations.  

#### **Required inputs**

| File                           | Description                                                                                                                                                        | Required columns                             |
| :-----------------------------| :------------------------------------------------------------------------------------------------------------------------------------------------------------------ | :-------------------------------------------- |
| `path_abun_unstrat.csv`       | Pathway abundance data. Values will be converted to relative abundance within the function.                                                                        | `SampleID`, `FunctionID`                      |
| `metabolite_concentration.csv`| Metabolite concentration data.                                                                                                                                     | `SampleID`, Metabolite names (as columns)     |
| `sample_metadata.csv`         | Sample metadata with group or condition information. If not provided or `class` column is missing, correlations will be performed on the overall dataset. | `SampleID`, `class`                           |
| `gsea_results_*.csv`          | GSEA results containing identified pathways.                                                                                                                       | `ID` (pathway ID)                             |

<details>
<summary>Click to show the full R function</summary>

```r
library(dplyr)
library(tidyr)
library(stats)
library(tibble)
library(stringr)
library(tools) 

construct_pathway_metabolite_network <- function(
  pathway_abundance_file,
  metabolite_concentration_file,
  gsea_results_file,
  metadata_file,
  output_file, # This is used as the output directory
  correlation_method = c("spearman", "pearson"),
  filter_by = c("none", "p_value", "q_value"), # Still used for filtering logic, just not filename
  corr_cutoff,
  p_value_cutoff, # Still used for filtering logic, just not filename
  q_value_cutoff, # Still used for filtering logic, just not filename
  q_adjust_method = c("bonferroni", "fdr") # Still used for filtering logic, just not filename
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
  
  # Determine gsea_suffix for filename and group processing priority
  gsea_suffix <- NULL
  gsea_target_group_from_filename <- NULL
  
  if (!is.null(gsea_results_file)) {
    # Check if the GSEA file exists before trying to read its name
    if (file.exists(gsea_results_file)) {
      gsea_basename <- tools::file_path_sans_ext(basename(gsea_results_file))
      
      # Account for additional parameters after the group names
      match_result <- stringr::str_match(gsea_basename, "^gsea_results_([^_]+)_vs_([^_]+).*$")
      
      if (!is.na(match_result[1,1])) {
        gsea_source_group <- match_result[1,2]
        gsea_target_group_from_filename <- match_result[1,3]
        
        gsea_suffix <- gsea_target_group_from_filename
        message("Derived GSEA filename suffix: '", gsea_suffix, "' (from GSEA target group)")
        
      } else {
        message("  Warning: GSEA results file name '", basename(gsea_results_file), "' did not strictly match the 'gsea_results_SOURCE_vs_TARGET*.csv' pattern. No GSEA-specific group filtering will be applied based on filename.")
      }
    } else {
      message("  Warning: GSEA results file '", gsea_results_file, "' not found. No GSEA-specific group filtering will be applied based on filename.")
    }
  } else {
    message("No GSEA results file provided. Processing will default to groups from metadata or overall if no metadata.")
  }
  
  # 1. Load data
  message("1. Loading pathway abundance, metabolite concentration, and GSEA results...")
  
  # Load Pathway Abundance (Expected: Pathways as Rows, Samples as Columns)
  pathway_abun_absolute <- tryCatch(
    read.csv(pathway_abundance_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error loading pathway abundance file '", pathway_abundance_file, "': ", e$message, sep = ""))
    }
  )
  
  # Load Metabolite Concentration (Expected: Samples as Rows, Metabolites as Columns)
  metabolite_conc_absolute <- tryCatch(
    read.csv(metabolite_concentration_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error loading metabolite concentration file '", metabolite_concentration_file, "': ", e$message, sep = ""))
    }
  )
  
  message("  Loaded pathway abundance. Dimensions: ", paste(dim(pathway_abun_absolute), collapse = "x"))
  message("  Loaded metabolite concentration. Dimensions: ", paste(dim(metabolite_conc_absolute), collapse = "x"))
  
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
      message("  Loaded metadata. Dimensions: ", paste(dim(metadata), collapse = "x"))
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
    message("2. No metadata file provided. Performing correlation on overall dataset (single group).")
  }
  
  # 3. Data preparation: Ensure SampleIDs align
  message("3. Preparing data for correlation...")
  
  # Get common samples across all datasets
  common_samples <- intersect(colnames(pathway_abun_absolute), rownames(metabolite_conc_absolute))
  if (!is.null(metadata)) {
    common_samples <- intersect(common_samples, rownames(metadata))
  }
  
  if (length(common_samples) == 0) {
    stop("No common samples found between pathway abundance, metabolite concentration, and metadata (if provided). Please check SampleIDs for consistency across all input files.")
  }
  message("  Number of common samples: ", length(common_samples))
  
  # Filter abundance data to common samples
  pathway_abun_filtered_by_samples <- pathway_abun_absolute[, common_samples, drop = FALSE]
  
  # metabolite_conc_filtered_by_samples will be Samples x Metabolites
  metabolite_conc_filtered_by_samples <- metabolite_conc_absolute[common_samples, , drop = FALSE]
  
  # Filter metadata to common samples OR create an "overall" group if no metadata
  if (!is.null(metadata)) {
    metadata_filtered <- metadata[common_samples, , drop = FALSE]
  } else {
    metadata_filtered <- data.frame(SampleID = common_samples, class = "overall", row.names = common_samples)
    metadata_filtered$class <- as.factor(metadata_filtered$class)
    message("  No metadata provided, processing all samples as a single 'overall' group.")
  }
  
  # Convert pathway abundance to relative values (SAMPLE-WISE normalization)
  message("  Converting pathway abundance to relative values (sample-wise) using dplyr syntax...")
  pathway_abun_temp_tibble <- pathway_abun_filtered_by_samples %>%
    tibble::rownames_to_column(var = "FunctionID")
  
  relative_abundance_tibble <- pathway_abun_temp_tibble %>%
    mutate(across(-FunctionID, ~ {
      col_sum <- sum(., na.rm = TRUE)
      if (col_sum == 0) {
        0
      } else {
        . / col_sum
      }
    }, .names = "{col}"))
  
  pathway_abun_relative <- relative_abundance_tibble %>%
    tibble::column_to_rownames(var = "FunctionID")
  
  pathway_abun_relative[is.na(pathway_abun_filtered_by_samples)] <- NA
  message("  Pathway abundance converted to relative values.")
  
  # Filter normalized pathways based on GSEA results
  pathway_abun_for_correlation <- pathway_abun_relative
  if (!is.null(gsea_results_file) && file.exists(gsea_results_file)) { # Check file existence again before loading
    message("  Filtering normalized pathways based on GSEA results...")
    gsea_results_df <- tryCatch(
      read.csv(gsea_results_file, header = TRUE, stringsAsFactors = FALSE),
      error = function(e) {
        stop(paste("Error loading GSEA results file '", gsea_results_file, "': ", e$message, sep = ""))
      }
    )
    
    pathway_id_col <- NULL
    if ("ID" %in% colnames(gsea_results_df)) {
      pathway_id_col <- "ID"
    } else {
      stop("GSEA results file must contain an 'ID' column with pathway identifiers.")
    }
    
    pathways_from_gsea <- unique(gsea_results_df[[pathway_id_col]])
    
    initial_normalized_pathway_count <- nrow(pathway_abun_relative)
    pathway_abun_for_correlation <- pathway_abun_relative[rownames(pathway_abun_relative) %in% pathways_from_gsea, , drop = FALSE]
    
    if (nrow(pathway_abun_for_correlation) == 0) {
      stop("No pathways from GSEA results were found in the normalized pathway abundance file. Please check pathway identifiers in both files.")
    }
    message("  Filtered normalized pathway abundance to ", nrow(pathway_abun_for_correlation), " pathways (from ", initial_normalized_pathway_count, " initial) based on GSEA results.")
  } else {
    message("  No GSEA results file provided or file not found. Proceeding with all normalized pathways for correlation.")
  }
  
  # Ensure all data are numeric for correlation (convert data frames to numeric)
  pathway_abun_for_correlation[] <- lapply(pathway_abun_for_correlation, as.numeric)
  metabolite_conc_filtered_by_samples[] <- lapply(metabolite_conc_filtered_by_samples, as.numeric)
  
  # Determine groups to process based on metadata AND GSEA filename priority
  all_metadata_groups <- unique(metadata_filtered$class)
  groups_to_process <- all_metadata_groups # Default: process all groups from metadata
  
  if (!is.null(gsea_target_group_from_filename)) { # If a target group was successfully extracted from GSEA filename
    if (gsea_target_group_from_filename %in% all_metadata_groups) {
      # If the GSEA target group exists in metadata, process ONLY that group
      groups_to_process <- gsea_target_group_from_filename
      message("  Restricting processing to group: '", gsea_target_group_from_filename, "' as implied by GSEA results filename.")
    } else {
      warning("  GSEA target group '", gsea_target_group_from_filename, "' (from filename) not found in metadata groups. Processing all groups from metadata.")
      # Fallback to all groups from metadata if the target group from filename isn't in metadata
      groups_to_process <- all_metadata_groups
    }
  } else {
    message("  Processing all groups from metadata (no specific GSEA target group identified from filename).")
  }
  
  if (length(groups_to_process) == 0) {
    stop("No valid groups to process after filtering based on GSEA filename and metadata. Please check group names consistency.")
  }
  message("Final groups to process: ", paste(groups_to_process, collapse = ", "))
  
  all_correlation_results <- list()
  
  # 4. Loop over each (now potentially filtered) group to perform correlations
  message("4. Performing correlations for each selected group...")
  for (current_group in groups_to_process) {
    message("  Processing group: '", current_group, "'")
    
    # Get samples belonging to the current group
    samples_in_group <- rownames(metadata_filtered[metadata_filtered$class == current_group, , drop = FALSE])
    
    if (length(samples_in_group) < 3) {
      warning("  Not enough samples (less than 3) in group '", current_group, "' to perform meaningful correlations. Skipping.")
      next
    }
    
    # Subset data for the current group
    path_data_group <- pathway_abun_for_correlation[, samples_in_group, drop = FALSE]
    met_data_group <- metabolite_conc_filtered_by_samples[samples_in_group, , drop = FALSE]
    
    # Transpose pathway data for correlation: cor.test expects vectors or Samples as Rows
    path_data_group <- t(path_data_group)
    
    # Get column names for pathways and metabolites
    pathway_names <- colnames(path_data_group)
    metabolite_names <- colnames(met_data_group)
    
    # Initialize empty list to store correlation results in long format
    results_list_for_group <- list()
    
    message("    Calculating correlations using cor.test for each pair...")
    # Loop through each pathway and each metabolite to get individual correlation tests
    for (path_name in pathway_names) {
      for (met_name in metabolite_names) {
        vec_path <- path_data_group[, path_name]
        vec_met <- met_data_group[, met_name]
        
        # Calculate number of complete observations for this specific pair
        valid_obs_count <- sum(complete.cases(vec_path, vec_met))
        
        # Only perform cor.test if there are at least 3 valid pairs
        if (valid_obs_count >= 3) {
          cor_args <- list(x = vec_path, y = vec_met, method = correlation_method)
          if (correlation_method == "spearman") {
            cor_args$exact <- FALSE # For large N, exact p-value is computationally intensive for Spearman
          }
          
          cor_test_result <- tryCatch({
            do.call(stats::cor.test, cor_args)
          }, error = function(e) {
            return(NULL) # Return NULL if error occurs
          })
          
          if (!is.null(cor_test_result)) {
            results_list_for_group[[length(results_list_for_group) + 1]] <- data.frame(
              FunctionID = path_name,
              MetaboliteID = met_name,
              Correlation = cor_test_result$estimate,
              P_value = cor_test_result$p.value,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    
    if (length(results_list_for_group) == 0) {
      message("    No valid pathway-metabolite correlations found for group '", current_group, "'. Skipping output for this group.")
      next # Skip to next group if no results
    }
    
    # Combine all results for the current group into a single data frame
    combined_results <- do.call(rbind, results_list_for_group)
    
    # Calculate Q-values using the specified adjustment method for ALL p-values in this group
    combined_results$Q_value <- stats::p.adjust(combined_results$P_value, method = q_adjust_method)
    
    # 5. Apply filtering based on user choice (now including Q_value)
    message("  Applying filters for group '", current_group, "'...")
    
    # Filter by absolute correlation coefficient first
    combined_results_filtered <- combined_results %>%
      filter(abs(Correlation) >= corr_cutoff)
    
    if (filter_by == "p_value") {
      combined_results_filtered <- combined_results_filtered %>%
        filter(P_value <= p_value_cutoff)
      message("  Filtered by p-value <= ", p_value_cutoff)
    } else if (filter_by == "q_value") {
      combined_results_filtered <- combined_results_filtered %>%
        filter(Q_value <= q_value_cutoff)
      message("  Filtered by q-value <= ", q_value_cutoff, " (", q_adjust_method, " correction)")
    }
    
    if (nrow(combined_results_filtered) == 0) {
      message("  No significant correlations found after filtering for group '", current_group, "'. Skipping output for this group.")
      next
    }
    
    combined_results_filtered$Group <- current_group
    all_correlation_results[[current_group]] <- combined_results_filtered
    
    # --- Save results with specific filename based on parameters ---
    
    # Determine the group name for the filename
    # Prioritize gsea_suffix if it was successfully extracted and applies to this group
    group_name_for_file <- current_group
    if (!is.null(gsea_suffix) && current_group == gsea_suffix) {
      group_name_for_file <- gsea_suffix
    }
    # If current_group is "overall" (meaning no metadata was provided), and gsea_suffix exists,
    # then use gsea_suffix in the filename. Otherwise, "overall" is correct.
    if (current_group == "overall" && !is.null(gsea_suffix)) {
      group_name_for_file <- gsea_suffix
    }
    
    # Format cutoffs for filename (keeping decimals)
    # Use formatC to avoid scientific notation for very small numbers, use "f" for fixed
    corr_cutoff_fname <- formatC(corr_cutoff, format = "f", digits = 2) 
    
    # Construct the simplified output filename
    output_filename <- paste0(
      "pathway_metabolite_network_",
      group_name_for_file, "_",
      correlation_method, "_",
      corr_cutoff_fname,
      ".csv"
    )
    
    output_filepath <- file.path(output_file, output_filename)
    write.csv(combined_results_filtered, output_filepath, row.names = FALSE)
    message("  Saved results for group '", current_group, "' to: ", output_filepath)
  }
  
  message("Pathway-metabolite network construction complete.")
  return(invisible(NULL))
}
```

</details>

#### **Example usage**

```r
source("https://raw.githubusercontent.com/bowornpol/NUIM-pipeline/main/code/construct_pathway_metabolite_network.R")

construct_pathway_metabolite_network(
  pathway_abundance_file = "path_abun_unstrat.csv", 
  metabolite_concentration_file = "metabolite_concentration.csv",
  gsea_results_file = "pathway_pathway_network_results/gsea_results_G1_vs_G2_0.1_fdr_signed_log_pvalue.csv",
  metadata_file = "sample_metadata.csv", # Optional, set to NULL if no groups
  output_file = "pathway_metabolite_network_results", # Output directory for results
  correlation_method = "pearson", # User can choose from "spearman" or "pearson"
  filter_by = "none", # User can choose from "none", "p_value", or "q_value"
  corr_cutoff = 0.5, # Absolute correlation coefficient cutoff (e.g., 0.5)
  p_value_cutoff = NULL, # Set if filter_by = "p_value"
  q_value_cutoff = NULL, # Set if filter_by = "q_value"
  q_adjust_method = NULL # User can choose "bonferroni" or "fdr" if filter_by = "q_value"
)
```

#### **Example output**

The function creates an output directory (e.g., `pathway_metabolite_network_results`) containing `.csv` files for each group analyzed (e.g., `pathway_metabolite_network_[class]_[correlation_method]_[corr_cutoff].csv` or `pathway_metabolite_network_overall.csv` if no groups are defined).

**Example table: `pathway_metabolite_network_G2_pearson_corr_0.30.csv`**

| FunctionID | MetaboliteID | Correlation | P_value | Q_value | Group |
|:-----------|:-------------|:------------|:--------|:--------|:------|
| ko00010    | acetate       | 0.82        | 0.001   | 0.005   | G2    |
| ko00020    | butyrate      | -0.75       | 0.003   | 0.008   | G2    |
| ko00300    | propionate    | 0.68        | 0.015   | 0.023   | G2    |
| ko00400    | lactate       | -0.62       | 0.025   | 0.041   | G2    |
| ...        | ...           | ...         | ...     | ...     | ...   |

Each row represents a correlation (edge) between a pathway (`FunctionID`) and a metabolite (`MetaboliteID`).
- `Correlation`: The correlation coefficient (Spearman or Pearson) indicating the strength and direction of the relationship.
- `P_value`: The raw p-value for the correlation.
- `Q_value`: The adjusted p-value (q-value), if q-value filtering was selected.
- `Group`: The specific group for which the correlation was calculated.

### <ins>STEP 4: Multi-layered network construction</ins>

These networks are finally integrated through connected pathway nodes to construct a multi-layered network.

### **Required inputs**

| File | Description | Required columns |
| :------------- | :---------- | :--------------- |
| `microbe_pathway_network_*.csv` | Microbe-pathway network table. | `TaxonID`, `FunctionID`, `relative_contribution` |
| `pathway_jaccard_*.csv` | Pathway-pathway network table. | `FunctionID_1`, `FunctionID_2`, `jaccard_index` |
| `pathway_metabolite_network_*.csv` | Pathway-metabolite network table. | `FunctionID`, `MetaboliteID`, `Correlation` |
| `gsea_results_*.csv` | GSEA results containing identified pathways. | `ID` (pathway ID) |

<details>
<summary>Click to show the full R function</summary>

```r
library(dplyr)
library(tidyr)
library(stringr)
library(tools) 

construct_multi_layered_network <- function(
  gsea_results_file,
  microbe_pathway_file,
  pathway_jaccard_file,
  pathway_metabolite_file,
  output_directory # Renamed from output_file to clarify it's a directory
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
    # Check if the GSEA file actually exists before trying to parse its name
    if (file.exists(gsea_results_file)) {
      gsea_basename <- tools::file_path_sans_ext(basename(gsea_results_file))
      
      # Account for additional parameters after the group names
      match_result <- stringr::str_match(gsea_basename, "^gsea_results_([^_]+)_vs_([^_]+).*$")
      
      if (!is.na(match_result[1,1])) { # If the full pattern matches successfully
        gsea_source_group <- match_result[1,2] 
        gsea_target_group_from_filename <- match_result[1,3] 
        
        gsea_suffix <- gsea_target_group_from_filename # Use target group as suffix for filename
        message("Derived GSEA filename suffix: '", gsea_suffix, "' (from GSEA target group)")
        
      } else {
        message("  Warning: GSEA results file name '", basename(gsea_results_file), "' did not strictly match the 'gsea_results_SOURCE_vs_TARGET*.csv' pattern for precise group identification. Output filename will use 'overall' suffix.")
      }
    } else {
      message("  Warning: GSEA results file '", gsea_results_file, "' not found. Output filename will use 'overall' suffix.")
    }
  } else {
    message("No GSEA results file provided, defaulting to processing without GSEA-specific group identification.")
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
```

</details>

#### **Example usage**

```r
source("https://raw.githubusercontent.com/bowornpol/NUIM-pipeline/main/code/construct_multi_layered_network.R")

# Define the full path and filename for your input CSV file
my_microbe_pathway_file <- "microbe_pathway_network_results/microbe_pathway_network_G2_median.csv"
my_pathway_jaccard_file <- "pathway_pathway_network_results/pathway_jaccard_G1_vs_G2.csv"
my_pathway_metabolite_file <- "pathway_metabolite_network_results/pathway_metabolite_network_G2_pearson_corr_0.30.csv"
my_gsea_file <- "pathway_pathway_network_results/gsea_results_G1_vs_G2_0.1_fdr_signed_log_pvalue.csv"

# Define the full path and filename for your output CSV file
my_output_file <- "multi_layered_network_results"

# Call the function with your specific file paths 
construct_multi_layered_network(
  microbe_pathway_file = my_microbe_pathway_file,
  pathway_jaccard_file = my_pathway_jaccard_file,
  pathway_metabolite_file = my_pathway_metabolite_file,
  gsea_results_file = my_gsea_file,
  output_file = my_output_file
)
```

#### **Example output**

The `construct_multi_layered_network` function generates a single CSV file at the specified `output_file` path (e.g., `multi_layered_network_results/multi_layered_network_[class].csv`). This file integrates all specified network layers through connected GSEA-significant pathway nodes.

**Example table: `multi_layered_network_G2.csv`**

| Feature1      | Feature2      | Edge_Score      | Edge_Type           |
| :------------ | :------------ | :-------------- | :------------------ |
| g__Blautia  | ko00010       | 0.85            | Microbe-Pathway     |
| g__Bifidobacterium  | ko00020       | 0.90            | Microbe-Pathway     |
| ko00010       | ko00020       | 0.15            | Pathway-Pathway     |
| ko00020       | ko00052       | 0.21            | Pathway-Pathway     |
| ko00010       | butyrate     | 0.78            | Pathway-Metabolite  |
| ko00030       | acetate     | 0.55            | Pathway-Metabolite  |
| ...           | ...           | ...             | ...                 |

Each row represents a connection between two features (`Feature1`, `Feature2`).
- `Edge_Score`: Quantifies the strength or value of this connection between `Feature1` and `Feature2`.
- `Edge_Type`: Indicates the original network layer from which the connection originated.

## Module 3: Network Analysis

This module provides three network analyses designed to identify context-specific associations. It takes the multi-layered network generated from **Module 2** as input.

### <ins>(1) Hub identification</ins>

The hub identification uses the Maximal Clique Centrality (MCC) algorithm to identify key microbial pathways.  

#### **Required inputs**

| File | Description | Required columns |
| :------------- | :---------- | :--------------- |
| `multi_layered_network_*.csv` | Multi-layered network table. | `Feature1`, `Feature2`, `Edge_Score`, `Edge_Type` |

The hub identification will be perform on Cytoscape with cytoHubba plugin. User can follow these steps:

1. Go to **File → Import → Network from File**, and select your input file (e.g., `multi_layered_network_G2.csv`).
2. In the import settings, assign **Source Node** to the first column and **Target Node** to the second column.
3. Open the **cytoHubba** tab in the Control Panel. Under *Target Network*, select your imported network and click **Calculate** under the *Node's Score* section.
4. In the *Select nodes with Hubba nodes* section, specify the number of top nodes to rank using the **MCC** algorithm, then click **Submit**.
5. To export the result table, navigate to the **Unassigned Tables** tab, then click **Export Table to File**.

#### **Example output**

The example output below shows the result of applying the MCC algorithm using cytoHubba plugin in Cytoscape. The top-ranked hub pathways are highlighted in the network view and listed in a results table, along with various centrality scores used to assess their importance.

<p align="center">
  <img src="figures/cytohubba.png" width="850"/>
</p>

### <ins>(2) Pathfinding</ins>

The pathfinding uses the Dijkstra's algorithm to identify the shortest path between the selected source and target nodes.

#### **Required inputs**

| File | Description | Required columns |
| :------------- | :---------- | :--------------- |
| `multi_layered_network_*.csv` | Multi-layered network table. | `Feature1`, `Feature2`, `Edge_Score`, `Edge_Type` |

<details>
<summary>Click to show the full R function</summary>

```r
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
```

</details>

#### **Example usage**

```r
source("https://raw.githubusercontent.com/bowornpol/NUIM-pipeline/main/code/pathfinding.R")

# Define the full path and filename for your input CSV file
my_multi_layered_network_file <- "multi_layered_network_results/multi_layered_network_G2.csv" 

# Define the source node 
my_source_node <- "g__Megamonas"

# Define the target node 
my_target_node <- "acetate"

# Define the full path for your output directory
my_output_directory <- "pathfinding_results"

# Call the function with your specific file paths
pathfinding(
  multi_layered_network_file = my_multi_layered_network_file, 
  source_node = my_source_node,                            
  target_node = my_target_node,
  output_directory = my_output_directory
)
```

#### **Example output**

The `pathfinding` function generates a single CSV file named `path_[source]_to_[target].csv` for the shortest path found between the specified `source_node` and `target_node`. This file is saved in the specified `output_directory`.

**Example table: `path_g_Megamonas_to_acetate.csv`**

| Source     | Target    | Original_Edge_Score | Transformed_Weight | Edge_Type           |
|------------|-----------|----------------------|---------------------|----------------------|
| g_Megamonas  | ko00540   | 0.45                 | 2.21                | Microbe-Pathway      |
| ko00540    | acetate   | 0.59                 | 1.70                | Pathway-Metabolite   |

Each row in the ``path_<source>_to_<target>.csv` file represents an individual edge along the route from the source node to the target node.
- `Source`, `Target`: The two nodes connected by the edge.
- `Original_Edge_Score`: The input edge score from the multi-layered network.
- `Transformed_Weight`: The cost used for Dijkstra’s algorithm.
- `Edge_Type`: Indicates the original network layer from which the connection originated.

### <ins>(3) Node prioritization</ins>

The node prioritization uses the Laplacian Heat Diffusion (LHD) algorithm to identify microbe-associated metabolites.

#### **Required inputs**

| File | Description | Required columns |
| :------------- | :---------- | :--------------- |
| `multi_layered_network_*.csv` | Multi-layered network table. | `Feature1`, `Feature2`, `Edge_Score`, `Edge_Type` |

<details>
<summary>Click to show the full R function</summary>

```r
library(igraph)
library(expm)
library(dplyr)
library(ggplot2)

node_prioritization <- function(
  multi_layered_network_file,
  output_directory,
  time_step_interval = 0.01,
  stabilization_threshold = 0.0001,
  stabilization_window_size = 10,
  filter_other_metabolite_edges
) {
  message("Starting node prioritization using LHD algorithm.")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
    message("Created output directory: ", output_directory)
  } else {
    message("Output directory already exists: ", output_directory)
  }
  
  # Extract base name from input file for output specificity
  input_file_base_name <- tools::file_path_sans_ext(basename(multi_layered_network_file))
  # Clean the base name for use in filenames (e.g., replace non-alphanumeric with underscore)
  cleaned_input_file_name <- gsub("[^A-Za-z0-9_]", "", input_file_base_name)
  
  # 1. Load the multi-layered network file (FULL network data)
  message("\n1. Loading multi-layered network from: ", multi_layered_network_file)
  if (!file.exists(multi_layered_network_file)) {
    stop("Input network file not found: '", multi_layered_network_file, "'. Cannot proceed.")
  }
  
  combined_data_full <- tryCatch( # Renamed to combined_data_full to denote original data
    read.csv(multi_layered_network_file, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error reading multi-layered network file '", multi_layered_network_file, "': ", e$message, sep = ""))
    }
  )
  
  # Validate required columns in combined_data_full
  required_cols_network <- c("Feature1", "Feature2", "Edge_Score", "Edge_Type")
  if (!all(required_cols_network %in% colnames(combined_data_full))) {
    stop(paste("Input network file '", multi_layered_network_file, "' must contain columns: ", paste(required_cols_network, collapse = ", "), ". Please ensure it's the output from Module 2, Step 4.", sep = ""))
  }
  
  # Ensure Edge_Score is numeric
  if (!is.numeric(combined_data_full$Edge_Score)) {
    stop("Column 'Edge_Score' in the input network file must be numeric.")
  }
  
  # 2. Identify all unique metabolite nodes from the FULL network (for looping and for filtering)
  all_metabolite_nodes <- combined_data_full %>%
    filter(Edge_Type == "Pathway-Metabolite") %>%
    pull(Feature2) %>%
    unique() %>%
    as.character()
  
  if (length(all_metabolite_nodes) == 0) {
    stop("No metabolite nodes found in the network (based on 'Pathway-Metabolite' Edge_Type). Cannot proceed diffusion.")
  }
  message("Identified ", length(all_metabolite_nodes), " unique metabolite nodes for seeding and potential filtering.")
  
  # Define the function for heat vector at time t
  H_vector_func <- function(t_val, laplacian_matrix, initial_heat_vector) {
    # Ensure dimensions match before multiplication
    if (nrow(laplacian_matrix) != length(initial_heat_vector)) {
      stop("Dimension mismatch: Laplacian matrix rows (", nrow(laplacian_matrix), ") and initial heat vector length (", length(initial_heat_vector), ") do not match.")
    }
    expm(-laplacian_matrix * t_val) %*% initial_heat_vector
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
      correlations[i] <- ifelse(sd(H_current) == 0 || sd(H_prev) == 0, 1, cor(H_current, H_prev, method = "spearman", use = "complete.obs"))
    }
    
    correlation_df <- data.frame(Time = time_steps_local[-1], Correlation = correlations) # Time corresponds to the later time point of the pair
    
    stabilization_t_found <- time_steps_local[length(correlations) + 1] # Default to last time if no stabilization
    
    # Ensure window_size is not larger than available correlations
    if (window_size >= length(correlations)) {
      stabilization_t_found <- time_steps_local[length(correlations) + 1]
    } else {
      for (i in seq(length(correlations) - window_size + 1)) {
        diffs <- abs(diff(correlations[i:(i + window_size - 1)]))
        if (all(diffs < threshold)) {
          stabilization_t_found <- time_steps_local[i+1] # Time at the start of the stabilization *window* of the second vector in the pair.
          break
        }
      }
    }
    
    return(list(stabilization_time = stabilization_t_found, correlations_df = correlation_df))
  }
  
  # 3. Loop through each metabolite node to perform heat diffusion
  message("\nPerforming heat diffusion for each metabolite seed node...")
  for (seed_metabolite_id in all_metabolite_nodes) { # Loop through all identified metabolites as seeds
    message("  Processing seed metabolite: '", seed_metabolite_id, "'")
    
    # --- Determine network for current seed based on filter_other_metabolite_edges ---
    if (filter_other_metabolite_edges) {
      message("    Applying filtering: Excluding all nodes that are metabolites EXCEPT the current seed ('", seed_metabolite_id, "') and their connected edges.")
      
      # Identify all other metabolite nodes to exclude (all metabolites MINUS the current seed)
      other_metabolite_nodes_to_exclude <- setdiff(all_metabolite_nodes, seed_metabolite_id)
      
      # Filter edges: keep only edges where *neither* Feature1 nor Feature2 is one of the 'other metabolite' nodes
      current_combined_data <- combined_data_full %>%
        filter(
          !(Feature1 %in% other_metabolite_nodes_to_exclude) &
            !(Feature2 %in% other_metabolite_nodes_to_exclude)
        )
      
      # Check if the seed node itself is still present after filtering 
      if (!seed_metabolite_id %in% unique(c(current_combined_data$Feature1, current_combined_data$Feature2))) {
        warning("  Seed metabolite '", seed_metabolite_id, "' is not present in the graph after filtering. This implies it was only connected to other metabolites, which were removed. Skipping diffusion for this seed.")
        next
      }
      
    } else {
      # Use the full network data (original behavior)
      current_combined_data <- combined_data_full
      message("    Using the full network for diffusion.")
    }
    
    # --- Convert Edge_Score to absolute value ---
    current_combined_data$Edge_Score <- abs(current_combined_data$Edge_Score)
    message("    Edge_Scores converted to absolute values.")
    
    # --- Create graph and Laplacian matrix for the CURRENT network configuration ---
    # This block is now inside the loop, as the network structure changes per seed if filtered
    
    # Check for valid edges/nodes to create a graph
    if (nrow(current_combined_data) == 0) {
      warning("  No edges found for graph construction for seed '", seed_metabolite_id, "' after filtering. Skipping diffusion.")
      next
    }
    
    # Get all unique nodes from the current filtered data to ensure the graph includes them
    all_nodes_in_current_data <- unique(c(current_combined_data$Feature1, current_combined_data$Feature2))
    
    g_current <- graph_from_data_frame(d = current_combined_data[, c("Feature1", "Feature2")], directed = FALSE,
                                       vertices = all_nodes_in_current_data) # Explicitly define vertices
    
    # Assign weights (handling potential mismatches as before)
    graph_edges_for_weighting <- igraph::as_data_frame(g_current, what = "edges") %>%
      rowwise() %>%
      mutate(Node_A = min(from, to), Node_B = max(from, to)) %>%
      ungroup()
    
    current_combined_data_sorted <- current_combined_data %>%
      rowwise() %>%
      mutate(Node_A = min(Feature1, Feature2), Node_B = max(Feature1, Feature2)) %>%
      ungroup() %>%
      select(Node_A, Node_B, Edge_Score) %>%
      distinct()
    
    weighted_edges_df <- graph_edges_for_weighting %>%
      left_join(current_combined_data_sorted, by = c("Node_A", "Node_B")) %>%
      select(from, to, Edge_Score)
    
    E(g_current)$weight <- weighted_edges_df$Edge_Score[match(paste0(graph_edges_for_weighting$from, graph_edges_for_weighting$to),
                                                              paste0(weighted_edges_df$from, weighted_edges_df$to))]
    
    if (any(is.na(E(g_current)$weight))) {
      warning("  Some edges in the graph for seed '", seed_metabolite_id, "' could not be matched to an 'Edge_Score' in the input data. Assigning 0 weight to unmatched edges.")
      E(g_current)$weight[is.na(E(g_current)$weight)] <- 0
    }
    
    # If the graph has become too small (e.g., only 1 node, no edges, or disconnected)
    if (vcount(g_current) < 2 || ecount(g_current) == 0 || !is.connected(g_current)) {
      warning("  Graph for seed '", seed_metabolite_id, "' is too sparse, disconnected, or lacks sufficient nodes/edges after filtering (nodes: ", vcount(g_current), ", edges: ", ecount(g_current), ", connected: ", is.connected(g_current), "). Skipping diffusion for this seed as Laplacian cannot be computed meaningfully.")
      next
    }
    
    L_current <- laplacian_matrix(g_current, weights = E(g_current)$weight)
    L_current <- as.matrix(L_current)
    
    # Initialize the heat vector H_0 for the current seed within the CURRENT graph's nodes
    H_0_current_graph <- numeric(vcount(g_current))
    names(H_0_current_graph) <- V(g_current)$name
    seed_index_current_graph <- which(V(g_current)$name == seed_metabolite_id)
    
    if (length(seed_index_current_graph) == 0) {
      warning("  Seed metabolite '", seed_metabolite_id, "' not found in the *filtered* graph's node set. This should not happen if previous checks are correct. Skipping diffusion.")
      next
    }
    H_0_current_graph[seed_index_current_graph] <- 1.0
    
    # Find stabilization time and get correlations data
    stabilization_data_result <- find_stabilization_data(L_current, H_0_current_graph, time_step_interval, stabilization_threshold, stabilization_window_size)
    stabilization_t <- stabilization_data_result$stabilization_time
    correlation_df <- stabilization_data_result$correlations_df
    
    message("    Stabilization time for '", seed_metabolite_id, "': t = ", round(stabilization_t, 4))
    
    # Calculate final heat scores at stabilization time
    final_heat_scores <- H_vector_func(stabilization_t, L_current, H_0_current_graph)
    
    # Create output data frame for this metabolite
    output_df <- data.frame(
      Node = V(g_current)$name, # Nodes from the current filtered graph
      Heat_Score = round(final_heat_scores, 10),
      stringsAsFactors = FALSE
    ) %>%
      arrange(desc(Heat_Score)) # Sort by heat score (descending)
    
    # Generate and save output files with cleaned input file name
    cleaned_seed_id <- gsub("[^A-Za-z0-9_]", "", seed_metabolite_id) # Clean seed ID for filename
    
    # Save heat scores to CSV
    output_file_name_heat <- paste0("heat_scores_", cleaned_seed_id, "_", cleaned_input_file_name, ".csv")
    output_path_heat <- file.path(output_directory, output_file_name_heat)
    write.csv(output_df, file = output_path_heat, row.names = FALSE)
    message("    Saved heat scores for '", seed_metabolite_id, "' to: ", output_path_heat)
    
    # Save correlation data to CSV
    output_file_name_correlation_csv <- paste0("spearman_correlations_", cleaned_seed_id, "_", cleaned_input_file_name, ".csv")
    output_path_correlation_csv <- file.path(output_directory, output_file_name_correlation_csv)
    write.csv(correlation_df, file = output_path_correlation_csv, row.names = FALSE)
    message("    Saved Spearman correlation data for '", seed_metabolite_id, "' to: ", output_path_correlation_csv)
    
    # Generate and save correlation plot
    correlation_plot <- ggplot(correlation_df, aes(x = Time, y = Correlation)) + 
      geom_line() +
      geom_vline(xintercept = stabilization_t, linetype = "dashed", color = "#a62140") + 
      geom_text(
        aes(
          x = stabilization_t, 
          y = max(Correlation, na.rm = TRUE) * 0.5, # Position text dynamically
          label = paste("t =", round(stabilization_t, 4))
        ),
        color = "#a62140", hjust = -0.1, vjust = 0.5, size = 5, fontface = "bold"
      ) +
      xlab("Time step") + 
      ylab("Spearman correlation with previous time step") + 
      scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + 
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
      )
    
    output_file_name_plot <- paste0("correlation_plot_", cleaned_seed_id, "_", cleaned_input_file_name, ".jpg")
    output_path_plot <- file.path(output_directory, output_file_name_plot)
    ggsave(output_path_plot, plot = correlation_plot, width = 8, height = 5, dpi = 600)
    message("    Saved correlation plot for '", seed_metabolite_id, "' to: ", output_path_plot)
  }
  
  message("Node prioritization complete.")
  return(invisible(NULL))
}
```

</details>

#### **Example usage**

```r
source("https://raw.githubusercontent.com/bowornpol/NUIM-pipeline/main/code/node_prioritization.R")

# Define the full path and filename for your input CSV file
my_multi_layered_network_file <- "multi_layered_network_results/multi_layered_network_G2.csv"

# Define the full path for your output directory
my_output_directory <- "node_prioritization_results"

# Call the function with your specific file paths
node_prioritization(
  multi_layered_network_file = my_multi_layered_network_file,
  output_directory = my_output_directory,
  # filter_other_metabolite_edges: This parameter controls the scope of the diffusion network
  # for each metabolite seed.
  #   - Set to TRUE: For each metabolite seed, the diffusion will occur ONLY within
  #     the sub-network that does *not* contain any other metabolite nodes (except the seed itself).
  #     This focuses the diffusion on the pathways and microbes directly connected to the seed,
  #     without "leaking" heat to other metabolites.
  #   - Set to FALSE: The diffusion will occur across the ENTIRE multi-layered network,
  #     allowing heat to spread to all connected nodes, including other metabolite nodes.
  #     This provides a broader view of the seed's influence across the full network.
  filter_other_metabolite_edges = TRUE 
)
```

#### **Example output**

The `node_prioritization` function generates multiple output files for each diffusion set in the specified `output_directory`:

1.  **`heat_scores_[metabolite]_[network].csv`**: A CSV file containing the final heat scores for all nodes that participated in the diffusion from that specific metabolite seed, sorted by `Heat_Score` in descending order.
2.  **`spearman_correlations_[metabolite]_[network]*.csv`**: A CSV file detailing the Spearman correlations between heat vectors at consecutive time steps, used for stabilization assessment.
3.  **`correlation_plot_[metabolite]_[network]*.jpg`**: A JPEG image visualizing the Spearman correlations over time, with the identified stabilization time marked.

**Example table: `heat_scores_acetate_multi_layered_network_G2.csv`**

| Node | Heat_Score |
| :------------- | :---------- |
| acetate | 0.480 |
| ko00010 | 0.123 |
| g__Blautia | 0.106 |
| ko00400 | 0.088 |
| g__Bifidobacterium | 0.023 |
| ... | ... |

Each row in the `heat_scores_*.csv` file represents a node in the diffusion network and its calculated heat score.
- `Node`: The identifier for a node in the network (e.g., microbe, pathway, metabolite).
- `Heat_Score`: The calculated heat score for that node at the diffusion stabilization time, indicating its importance relative to the seed node.
- 

Here is a comparison of NUIM, MIMOSA2, MMIP, MintTea, and MiMeNet in a plain text table format:

```
| Feature                     | NUIM                                                                                                                              | MIMOSA2                                                                                                      | MMIP                                                                                                         | MintTea                                                                                                              | MiMeNet                                                                                                       |
|-----------------------------|-----------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------|
| **Primary Goal** | To explore and visualize potential mechanisms linking microbes and metabolites by constructing an integrated, multi-layered network. | To predict metabolite levels from microbial community data and identify the key microbial taxa responsible for these metabolic shifts. | To provide a web-based platform for comparing microbial communities and predicting their metabolic potential from 16S rRNA data. | To identify robust, co-varying sets of multi-omic features (modules) that are strongly associated with a specific host phenotype, such as disease. | To construct a bipartite network between microbial taxa and metabolites, identifying significant associations. |
| **Input Data** | 16S rRNA sequencing data and a corresponding metabolomics dataset.                                                                | Microbial composition data (e.g., 16S or metagenomics) and, for validation, a corresponding metabolomics dataset. | Primarily designed for 16S rRNA sequencing data, but can integrate metabolomics data if available.             | At least two omics datasets (e.g., metagenomics, metabolomics) and a corresponding categorical phenotype for each sample. | Microbiome (e.g., 16S rRNA or metagenomics) and metabolomics data.                                            |
| **Core Integration Method** | Builds a tripartite network (Microbe-Pathway-Metabolite) using functional prediction, correlation, and gene set enrichment.         | Utilizes a knowledge-based approach, creating a community metabolic model from reference genomes to simulate and predict metabolic outputs. | Integrates functional prediction (e.g., PICRUSt2), statistical comparisons, and machine learning.              | Employs sparse Generalized Canonical Correlation Analysis (sGCCA) and consensus clustering.                      | Often uses correlation-based approaches (e.g., Spearman, Pearson) to establish links, followed by network analysis. |
| **Primary Output** | A multi-layered network visualization, identified hub pathways, shortest paths, and prioritized nodes.                              | A list of predicted metabolite abundances and the specific microbes predicted to be responsible for their production or consumption. | Visualizations and statistical results comparing microbial composition and predicted metabolic profiles.       | A list of multi-omic modules, where each module contains features highly correlated and associated with the phenotype. | A bipartite microbe-metabolite network, lists of significant associations, and network statistics.             |
| **Key Strength** | Provides a structured, hierarchical framework that visualizes connections from microbes to pathways to metabolites, aiding in mechanistic interpretation. | Offers strong, mechanistically plausible predictions based on known metabolic pathways and genome annotations.     | Highly accessible and user-friendly web-based platform that simplifies the process of integrated analysis.     | Statistically robust method for identifying clinically relevant, co-varying sets of multi-omic features with high predictive power. | Relatively straightforward to implement and interpret for direct pairwise microbe-metabolite relationships. Provides a clear visual representation. |
| **Limitations** | Relies on the accuracy of functional prediction and correlation, which does not confirm causation. Complexity may require more user expertise. | Performance is highly dependent on the completeness and accuracy of underlying reference genome and pathway databases. | Predictions are based on inferred function, not direct measurement, which can be a significant simplification.      | The output modules are statistically derived and require biological validation to confirm underlying mechanisms. | Primarily focuses on pairwise associations, which may miss indirect or pathway-mediated interactions. Can be prone to false positives if not rigorously filtered. |
```
