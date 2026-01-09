# Phylogeographic Analysis of Chikungunya Virus: A Comprehensive Genomic Surveillance Study

## Overview

This repository contains a comprehensive phylogeographic analysis of Chikungunya virus (CHIKV) based on genomic sequence data collected between 2017 and 2025. The study focuses on three major CHIKV genotypes (I-WestAfrica, II-ECSA, and III-Asian) and employs advanced phylogenetic methods to investigate global transmission patterns, evolutionary dynamics, and geographic spread.

## Abstract

Chikungunya virus (CHIKV) is a mosquito-borne alphavirus that has caused large-scale epidemics across tropical and subtropical regions worldwide. Despite its global health significance, the fine-scale transmission patterns and phylogeographic dynamics of CHIKV remain incompletely understood. This study presents a comprehensive genomic surveillance analysis of CHIKV using 1,500+ complete genome sequences collected from multiple continents between 2017 and 2025. We employed Bayesian phylogenetic methods (BEAST) to reconstruct time-calibrated phylogenies for three major CHIKV genotypes and analyzed transmission patterns using identical sequence pairs to characterize fine-scale geographic spread. Our analysis reveals distinct phylogeographic patterns across genotypes, with evidence of frequent inter-continental transmission events, particularly involving South America, Asia, and Africa. We applied relative risk (RR) analysis based on identical sequence pairs to identify transmission hotspots and quantify connectivity between geographic regions at multiple spatial scales (continent, country, and sub-national levels). Time-calibrated phylogenies indicate ongoing evolution and diversification of all three genotypes, with recent expansion events in Asia and South America. These findings provide critical insights into CHIKV transmission dynamics and inform public health surveillance and control strategies.

## Introduction

Chikungunya virus (CHIKV) is an arthropod-borne alphavirus (family *Togaviridae*) that causes chikungunya fever, a debilitating disease characterized by acute onset of fever and severe joint pain. First identified in Tanzania in 1952, CHIKV has since spread globally, with major outbreaks reported in Africa, Asia, the Indian Ocean islands, and the Americas. The virus is primarily transmitted by *Aedes aegypti* and *Aedes albopictus* mosquitoes, both of which are widely distributed across tropical and subtropical regions.

### Genomic Diversity and Lineages

CHIKV exhibits significant genetic diversity and can be classified into three main genotypes based on phylogenetic analysis of the envelope gene or complete genome sequences:

1. **Genotype I (West African)**: Originally restricted to West Africa, this genotype has been associated with sporadic cases and limited outbreaks in the region.

2. **Genotype II (East/Central/South African, ECSA)**: This is one of the most widespread genotypes, responsible for major epidemics in the Indian Ocean region, India, Southeast Asia, and parts of Europe. A key variant within this genotype, the E1-A226V mutation, enhances transmission by *Aedes albopictus*, contributing to its global spread.

3. **Genotype III (Asian)**: Originating in Asia, this genotype has been responsible for massive outbreaks in the Americas since its introduction in 2013. It has caused millions of cases across Central and South America, with ongoing transmission and local adaptation.

### Research Objectives

Despite decades of CHIKV research, several critical questions remain unanswered:

- What are the fine-scale transmission patterns between geographic regions at multiple spatial scales?
- How do transmission networks differ among the three major genotypes?
- What are the temporal dynamics of CHIKV spread and evolution?
- Which geographic regions serve as transmission hubs or bridges between continents?

To address these questions, we conducted a comprehensive phylogeographic analysis using genomic surveillance data. Our approach combines:

1. **Bayesian phylogenetic inference** using BEAST to reconstruct time-calibrated phylogenies
2. **Phylogeographic analysis** to track the geographic spread of viral lineages
3. **Transmission network analysis** using identical sequence pairs to identify fine-scale transmission patterns
4. **Relative risk (RR) analysis** to quantify connectivity between geographic regions

### Data Sources and Methods

We collected complete genome sequences from public databases (GISAID and GenBank) spanning the period from 2017 to 2025. After quality control and filtering (minimum length: 10,500 nucleotides), we compiled a dataset of 1,500+ sequences with comprehensive metadata including collection date, geographic location (continent, country, city), host type, and genotype classification.

**Phylogenetic Analysis Pipeline:**

1. **Sequence Alignment**: Multiple sequence alignment using MAFFT
2. **Phylogenetic Inference**: IQ-TREE for initial tree reconstruction, followed by BEAST analysis for time-calibrated phylogenies
3. **MCC Tree Extraction**: Maximum Clade Credibility (MCC) trees extracted from posterior tree distributions
4. **Time Calibration**: TreeTime for molecular clock calibration and divergence time estimation
5. **Phylogeographic Reconstruction**: Discrete trait analysis in BEAST to infer ancestral geographic states

**Transmission Analysis:**

1. **Identical Sequence Identification**: Identification of pairs of identical sequences within the dataset
2. **Geographic Annotation**: Annotation of sequence pairs with geographic metadata
3. **Relative Risk Calculation**: Computation of RR values to quantify transmission connectivity:
   ```
   RR = (observed pairs between regions) / (expected pairs based on random distribution)
   ```
4. **Uncertainty Quantification**: Bootstrap resampling to estimate confidence intervals for RR values

**Visualization:**

1. **Phylogenetic Trees**: Branch-colored trees showing geographic (country-level) distribution
2. **Geographic Maps**: Visualization of sequence distribution and transmission patterns
3. **Heatmaps**: Transmission rate matrices showing connectivity between regions

## Repository Structure

```
Chikungunya/
├── chikv_ALL/              # Combined analysis across all genotypes
│   ├── results/            # Analysis results (RR matrices, pair counts, etc.)
│   ├── scripts/            # Utility scripts for RR analysis
│   └── sub_trees/          # Sub-tree extraction and analysis
│
├── chikv_I_WestAfrica/     # Genotype I analysis
│   ├── West_African_aln_2MCC.tree        # MCC tree
│   ├── West_African_with_coordinates.tsv # Metadata with geographic coordinates
│   ├── plot_phylogenetic_tree.R          # Tree visualization script
│   └── results/            # Genotype-specific results
│
├── chikv_II_ECSA/          # Genotype II analysis
│   ├── ECSA_MCC.tree
│   ├── ECSA_with_coordinates.tsv
│   ├── plot_phylogenetic_tree.R
│   └── results/
│
├── chikv_III_Asian/        # Genotype III analysis
│   ├── Asian_alnMCC1.tree
│   ├── Asian_with_coordinates.tsv
│   ├── plot_phylogenetic_tree.R
│   └── results/
│
├── Genebank/               # Raw sequence data from GenBank
├── metadata.tsv            # Combined metadata file
└── README.md               # This file
```

## Key Features

### 1. Phylogenetic Tree Visualization

Scripts are available to generate publication-quality phylogenetic trees with:
- **Branch coloring** by geographic location (country-level)
- **Tip point coloring** by country
- **No tip labels** for cleaner visualization
- **Legend positioning** on the right side
- **Consistent styling** across all three genotypes

**Usage:**
```bash
cd chikv_III_Asian/
module load R/4.1.2
Rscript plot_phylogenetic_tree.R
```

This generates:
- `phylogenetic_tree.png` (high-resolution PNG)
- `phylogenetic_tree.pdf` (vector format PDF)

### 2. Transmission Analysis (RR Analysis)

The repository includes comprehensive scripts for analyzing transmission patterns using identical sequence pairs:

- **Pair identification**: `scripts/count_pairs.R`
- **RR calculation**: `scripts/RR_from_df_n_pairs.R`
- **Uncertainty estimation**: `scripts/uncertainty_RR.R`
- **Metadata integration**: `scripts/append_metadata_field.R`

**Outputs include:**
- Number of identical sequence pairs by geographic region
- Relative risk (RR) matrices for transmission connectivity
- Bootstrap confidence intervals for RR estimates
- Results at multiple spatial scales (continent, country, region)

### 3. Time-Calibrated Phylogenies

BEAST analyses were performed for each genotype to:
- Estimate divergence times
- Reconstruct ancestral states
- Infer phylogeographic history
- Quantify evolutionary rates

MCC trees are provided for each genotype, containing:
- Branch lengths in time units
- Posterior probabilities for nodes
- Ancestral geographic state reconstructions
- Height annotations

## Requirements

### Software Dependencies

- **R** (≥ 4.1.2) with packages:
  - `ape` - Phylogenetic analysis
  - `ggtree` - Phylogenetic tree visualization
  - `dplyr` - Data manipulation
  - `ggplot2` - Plotting
  - `RColorBrewer` - Color schemes
  - `tidytree` - Tree data manipulation

- **BEAST2** - Bayesian phylogenetic analysis
- **TreeAnnotator** - MCC tree extraction
- **IQ-TREE** - Maximum likelihood phylogeny
- **MAFFT** - Multiple sequence alignment
- **TreeTime** - Molecular clock calibration

### Installing R Packages

If `ggtree` is not available, the scripts will attempt to install it via BiocManager:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree")
```

Alternatively, install from GitHub:
```r
remotes::install_github("YuLab-SMU/ggtree")
```

## Usage Examples

### Generate Phylogenetic Tree

```bash
# For Genotype III (Asian)
cd /scr/u/dongw21/Chikungunya/chikv_III_Asian
module load R/4.1.2
Rscript plot_phylogenetic_tree.R
```

### Run Transmission Analysis

```r
# In R
source("chikv_ALL/scripts/utils_RR.R")
source("chikv_ALL/scripts/count_pairs.R")
source("chikv_ALL/scripts/RR_from_df_n_pairs.R")

# Example: Calculate RR by country
df_n_pairs <- count_identical_pairs(sequences, metadata, "Country")
df_RR <- get_df_RR(df_n_pairs, "Country")
```

## Results Summary

### Geographic Distribution

- **Genotype I (West Africa)**: Primarily restricted to West Africa and South America
- **Genotype II (ECSA)**: Widespread across Africa, Asia, and Europe
- **Genotype III (Asian)**: Dominant in Asia and the Americas, with recent expansion to new regions

### Transmission Patterns

- High connectivity observed between South America and Asia for Genotype III
- Frequent inter-continental transmission events identified
- Transmission hotspots identified in major urban centers
- Country-level RR analysis reveals asymmetric transmission flows

### Temporal Dynamics

- Continuous evolution and diversification across all genotypes
- Recent expansion events (2020-2025) in multiple regions
- Evidence of local adaptation and selection pressures

## Data Sources

### Sequence Data
- **GISAID**: EpiCoV database (acknowledgments provided in metadata)
- **GenBank**: NCBI Virus database

### Metadata Fields
- Accession ID
- Collection date
- Geographic location (Continent, Country, City, Coordinates)
- Host type (Human, Mosquito, etc.)
- Genotype classification

## Citation

If you use this code or data, please cite:

1. GISAID contributors (for sequences from GISAID)
2. Original GenBank submissions (for GenBank sequences)
3. This repository (if applicable)

## Acknowledgments

- Contributors to GISAID database
- GenBank submitters
- BEAST and TreeTime developers
- R package maintainers (ape, ggtree, etc.)

## Contact

For questions or issues, please contact the repository maintainer.

## License

[Specify license if applicable]

---

**Last Updated**: January 2025  
**Data Coverage**: 2017-2025  
**Total Sequences Analyzed**: 1,500+  
**Genotypes Analyzed**: I-WestAfrica, II-ECSA, III-Asian
