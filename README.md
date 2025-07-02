# NUIM: Network-based utility for integrating microbiome and metabolome data

We developed NUIM, a modular, network-based framework for integrating microbiome and metabolome data systematically. NUIM consists of three modules: (1) data preparation and processing, (2) network construction, and (3) network analysis.

<p align="center">
  <img src="figures/NUIM_overview.png" width="500"/>
</p>

## Module 1: Data preparation and processing

This module defines the procedures required to prepare and process the input data for downstream network construction and analysis.

- Input data includes `microbial sequencing reads in FASTQ format` and `metabolite concentration table`.  
- Microbiome data processing involves the use of QIIME2 to generate a feature table and representative sequences. These outputs are subsequently processed using PICRUSt2 for functional prediction, yielding gene abundance, pathway abundance, and pathway contribution data.  
- Although metabolome data processing may vary depending on user preference and experimental design, NUIM assumes that metabolite concentrations have been appropriately processed by standard practice. For example, users may employ established platforms such as [Metabox](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giae005/7629842) or [MetaboAnalyst](https://academic.oup.com/nar/article/52/W1/W398/7642060) to perform normalization, transformation, and quality control of metabolomics data.

### <ins>QIIME2 workflow</ins>

This section provides a general QIIME2 workflow for processing paired-end 16S rRNA sequencing data. The goal is to generate a feature table and representative sequences for PICRUSt2.

#### **Required Inputs**

| File            | Description                         |
|-----------------|-----------------------------------|
| `FASTQ`         | Raw paired-end sequencing reads   |
| `manifest.tsv`  | Mapping of sample IDs to FASTQ files |  

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

#### **QIIME2 Outputs**

| File                | Description                      |
|---------------------|---------------------------------|
| `feature-table.biom` | Feature (ASV) count table        |
| `rep_seqs.fasta`    | Representative sequences         |
| `taxonomy.tsv`      | Taxonomic assignment of ASVs | 

### <ins>PICRUSt2 workflow</ins>

PICRUSt2 predicts functional profiles from 16S rRNA data. This step uses a feature table and representative sequences from QIIME2.

#### **Required Inputs**

| File                | Description                     |
|---------------------|---------------------------------|
| `feature-table.biom` | Feature table exported from QIIME2 (`normalized_table.qza`) |
| `rep_seqs.fasta`    | Representative sequences exported from QIIME2 (`rep_seqs.qza`) |

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
  --map <pathway_mapping_file.tsv>

# Generate pathway contribution data (ASV × pathway)
pathway_pipeline.py \
  -i picrust2_out/<function_output_folder>/pred_metagenome_unstrat.tsv.gz \
  -o pathways_contrib \
  --per_sequence_contrib \
  --per_sequence_abun picrust2_out/<function_output_folder>/seqtab_norm.tsv.gz \
  --per_sequence_function <predicted_functions.tsv.gz> \
  --no_regroup \
  --map <pathway_mapping_file.tsv>

# Gene abundance data (KEGG Orthologs) output by PICRUSt2:
# Located in: KO_metagenome_out/pred_metagenome_unstrat.tsv.gz
# This file contains predicted gene family (KO) abundances per sample.
```

#### **PICRUSt2 Outputs**

| File                      | Description                            |
|---------------------------|--------------------------------------|
| `pred_metagenome_unstrat.tsv` | Predicted gene (KO) abundance        |
| `path_abun_unstrat.tsv`   | Predicted pathway abundance |
| `path_abun_contrib.tsv`   | Predicted pathway contribution |

## Module 2: Network Construction

This module constructs a tripartite network linking microbial taxa, metabolic pathways, and metabolites. Below, we provide R functions for constructing each network layer.

### <ins>Microbe–pathway network construction</ins>

The microbe–pathway network is constructed from pathway contribution data, with edges representing the relative contribution of each microbe to specific pathways.

#### **Required Inputs**

| File                  | Description                                                                 |
|-----------------------|-----------------------------------------------------------------------------|
| `path_abun_contrib.csv` | Pathway contribution data fdrom PICRUSt2  
**Required columns**: `SampleID`, `FeatureID`, `FunctionID`, `taxon_function_abun`. |
| `sample_metadata.csv`   | Sample metadata. Used when there are multiple experimental groups.  
**Required column**: `SampleID`. |
| `taxa_name.csv`         | Taxonomy annotations mapping `FeatureID` to taxonomic name (`TaxonID`).  
**Required columns**: `FeatureID`, `TaxonID`. |

> `sample_metadata.csv` is optional if group-specific processing is not needed.  
> `taxa_name.csv` is required to assign taxonomic labels to features.


```r
library(dplyr)

construct_microbe_pathway_network <- function(
  contrib_file,
  metadata_file,
  taxonomy_file,
  output_file = "microbe_pathway_network.csv",
  filtering = c("unfiltered", "mean", "median", "top10%", "top25%", "top50%", "top75%")
) {
  # Match filtering argument to allowed options
  filtering <- match.arg(filtering)

  # Read input CSV files
  contrib <- read.csv(contrib_file)        # Columns: SampleID, FeatureID, FunctionID, taxon_function_abun, etc.
  metadata <- read.csv(metadata_file)      # Columns: SampleID, group info (optional)
  taxonomy <- read.csv(taxonomy_file)      # Columns: FeatureID, TaxonID, taxonomy details

  # Merge contribution data with metadata by SampleID
  merged <- merge(contrib, metadata, by = "SampleID")
  # Merge the above result with taxonomy data by FeatureID (microbial feature identifier)
  merged <- merge(merged, taxonomy, by = "FeatureID", all.x = TRUE)

  # Summarize total abundance of each function contributed by each taxon
  taxon_function_total <- aggregate(
    taxon_function_abun ~ FunctionID + TaxonID,
    data = merged,
    FUN = sum,
    na.rm = TRUE
  )

  # Calculate total abundance of each function across all taxa (denominator for relative contribution)
  function_total <- aggregate(
    taxon_function_abun ~ FunctionID,
    data = taxon_function_total,
    FUN = sum,
    na.rm = TRUE
  )
  colnames(function_total)[2] <- "total_abundance_all_taxa"

  # Merge total abundance with taxon-function abundance to compute relative contribution
  taxon_function_total <- merge(taxon_function_total, function_total, by = "FunctionID")
  taxon_function_total$relative_contribution <- taxon_function_total$taxon_function_abun / taxon_function_total$total_abundance_all_taxa

  # Filter data based on user-selected threshold method
  if (filtering != "unfiltered") {
    if (filtering %in% c("mean", "median")) {
      # Calculate mean or median relative contribution per function
      threshold_df <- aggregate(
        relative_contribution ~ FunctionID,
        data = taxon_function_total,
        FUN = ifelse(filtering == "mean", mean, median),
        na.rm = TRUE
      )
      colnames(threshold_df)[2] <- "threshold"
      # Filter to keep only taxa with relative contribution >= threshold
      taxon_function_total <- merge(taxon_function_total, threshold_df, by = "FunctionID")
      taxon_function_total <- subset(taxon_function_total, relative_contribution >= threshold | is.na(threshold))
    } else if (grepl("^top", filtering)) {
      # For topX% filters: keep taxa contributing cumulatively up to the specified top percentage
      top_percent <- as.numeric(sub("top(\\d+)%", "\\1", filtering)) / 100
      taxon_function_total <- do.call(rbind, lapply(split(taxon_function_total, taxon_function_total$FunctionID), function(df) {
        # Sort by relative contribution descending
        df <- df[order(-df$relative_contribution), ]
        # Calculate cumulative sum
        df$cum_sum <- cumsum(df$relative_contribution)
        # Keep rows with cumulative sum less than or equal to top_percent
        df[df$cum_sum <= top_percent, ]
      }))
    }
  }

  # Order final results by FunctionID and descending relative contribution
  taxon_function_total <- taxon_function_total[order(taxon_function_total$FunctionID, -taxon_function_total$relative_contribution), ]

  # Save filtered and ordered data to CSV
  write.csv(taxon_function_total, output_file, row.names = FALSE)
}

# Example usage:
construct_microbe_pathway_network(
  contrib_file = "path_abun_contrib.csv",     # Pathway contribution data file
  metadata_file = "sample_metadata.csv",      # Sample metadata (optional: group info)
  taxonomy_file = "taxa_name.csv",         # Taxonomy info per microbial feature
  output_file = "microbe_pathway_network.csv",# Output file path
  filtering = "median"                                 # Filtering threshold: options are "unfiltered", "mean", "median", "top10%", "top25%", "top50%", "top75%"
)
```

#### **Outputs**

| File                        | Description                                                                 |
|-----------------------------|-----------------------------------------------------------------------------|
| `microbe_pathway_network.csv` | A microbe-pathway network table with relative contribution values. |

### <ins>Pathway–pathway network construction</ins>

The pathway–pathway network is constructed using pathways identified as significant through Gene Set Enrichment Analysis (GSEA). Edges between pathways are defined based on shared genes, and Jaccard indices represent edge weights.  

### <ins>Pathway–metabolite network construction</ins>

The pathway–metabolite network is constructed by calculating pairwise correlation (e.g., Spearman or Pearson) between pathway abundance and metabolite concentrations.  


### <ins>Multi-layered network</ins>

## Module 3: Network Analysis

This module provides three network analyses designed to identify context-specific associations:

- The hub identification uses the Maximal Clique Centrality (MCC) algorithm to identify key microbial pathways.  
- The pathfinding uses the Dijkstra's algorithm to identify the shortest path between the selected source and target nodes.  
- The node prioritization uses the Laplacian Heat Diffusion (LHD) algorithm to identify microbe-associated metabolites.
