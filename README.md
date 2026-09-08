# 16S rRNA Amplicon Sequencing (V3–V4) Analysis Pipeline

[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![DADA2](https://img.shields.io/badge/Bioconductor-DADA2-2D728F.svg)](https://bioconductor.org/packages/release/bioc/html/dada2.html)
[![Phyloseq](https://img.shields.io/badge/Bioconductor-Phyloseq-2D728F.svg)](https://bioconductor.org/packages/release/bioc/html/phyloseq.html)
[![DESeq2](https://img.shields.io/badge/Bioconductor-DESeq2-2D728F.svg)](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
[![Vegan](https://img.shields.io/badge/CRAN-vegan-green.svg)](https://cran.r-project.org/package=vegan)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An end-to-end, publication-grade bioinformatic and statistical workflow for high-throughput **16S rRNA gene amplicon sequencing (hypervariable regions V3–V4)** in R. 

This pipeline covers everything from raw Illumina paired-end demultiplexed reads to high-resolution Amplicon Sequence Variant (ASV) inference, taxonomic assignment against SILVA, microbial community ecology (alpha & beta diversity, PERMANOVA), relative taxonomic profiling, and negative binomial differential abundance analysis (DESeq2).

---

## Table of Contents

- [Overview & Experimental Context](#overview--experimental-context)
- [Workflow Architecture](#workflow-architecture)
- [Pipeline Stages & Technical Details](#pipeline-stages--technical-details)
  - [1. Quality Control, Trimming, and Filtering](#1-quality-control-trimming-and-filtering)
  - [2. Error Modeling & Sample Inference (DADA2)](#2-error-modeling--sample-inference-dada2)
  - [3. Pair Merging, Length Validation, and Chimera Removal](#3-pair-merging-length-validation-and-chimera-removal)
  - [4. Taxonomic Assignment (SILVA v138.2)](#4-taxonomic-assignment-silva-v1382)
  - [5. Phyloseq Integration & Contaminant Filtering](#5-phyloseq-integration--contaminant-filtering)
  - [6. Alpha Diversity Analysis](#6-alpha-diversity-analysis)
  - [7. Beta Diversity & Multivariate Community Ecology](#7-beta-diversity--multivariate-community-ecology)
  - [8. Taxonomic Composition & Abundance Export](#8-taxonomic-composition--abundance-export)
  - [9. Differential Abundance Testing (DESeq2)](#9-differential-abundance-testing-deseq2)
- [Repository Structure](#repository-structure)
- [Installation & Prerequisites](#installation--prerequisites)
- [Input Data Requirements](#input-data-requirements)
- [Usage Guide](#usage-guide)
- [Expected Outputs & Deliverables](#expected-outputs--deliverables)
- [Author & Acknowledgments](#author--acknowledgments)
- [License](#license)

---

## Overview & Experimental Context

Amplicon sequencing of the 16S ribosomal RNA hypervariable regions **V3–V4** is a gold standard for profiling bacterial communities in complex environments. Unlike traditional 97% OTU clustering methods that obscure fine-scale variation, this workflow employs **DADA2** to resolve exact **Amplicon Sequence Variants (ASVs)** down to single-nucleotide differences.

### Experimental Study Design (CFI Study)
The pipeline is demonstrated on a gut microbiome dataset investigating the physiological and microbiome-modulating effects of a dietary intervention (**CFI** / Caffeine intervention) in mouse models under standard and high-fat diet conditions:
- **Control**: Baseline standard diet controls.
- **Control-CFI**: Healthy controls receiving dietary intervention (identifying baseline intervention effects).
- **HFD**: High-Fat Diet-induced dysbiosis model.
- **HFD-CFI**: High-Fat Diet mice receiving dietary intervention (evaluating microbiome rescue / therapeutic effects).
- **Negative & Positive Controls**: Sequencing/extraction controls filtered out prior to biological evaluations.

---

## Workflow Architecture

```mermaid
flowchart TD
    A[Raw Paired-End FASTQ<br/>R1 / R2 V3-V4 Reads] --> B[Quality Profiling<br/>plotQualityProfile]
    B --> C[Filtering & Trimming<br/>filterAndTrim: 240bp / 230bp, maxEE=2]
    C --> D[Parametric Error Learning<br/>learnErrors ML Model]
    D --> E[Sample Inference & Denoising<br/>dada algorithm]
    E --> F[Paired-End Merging<br/>mergePairs: 12bp min overlap]
    F --> G[ASV Table Construction<br/>makeSequenceTable]
    G --> H[V3-V4 Length Filtering<br/>399 bp to 430 bp]
    H --> I[Bimera / Chimera Removal<br/>removeBimeraDenovo: consensus]
    I --> J[Taxonomic Classification<br/>SILVA v138.2 to Genus & Species]
    J --> K[Phyloseq Object Assembly<br/>seqtab + taxa + metadata]
    K --> L[Contaminant Filtering<br/>Remove Mitochondria, Chloroplast, Controls]
    
    L --> M[Alpha Diversity<br/>Observed, Chao1, Shannon, Simpson<br/>Kruskal-Wallis & Wilcoxon FDR]
    L --> N[Beta Diversity<br/>Bray-Curtis PCoA + 95% Ellipses<br/>Global & Pairwise PERMANOVA]
    L --> O[Taxonomic Composition<br/>Relative Abundance Stacked Bars<br/>Phylum to Species Wide Matrix Export]
    L --> P[Differential Abundance<br/>DESeq2 poscounts + Local GLM<br/>HFD-CFI vs HFD | Control vs CFI]
```

---

## Pipeline Stages & Technical Details

### 1. Quality Control, Trimming, and Filtering
- Inspects per-base quality profiles for forward and reverse read pools (`plotQualityProfile`).
- Executes `filterAndTrim` with optimized cutoffs for V3–V4 amplicon sequencing:
  - Forward truncation: `truncLen = 240` bp.
  - Reverse truncation: `truncLen = 230` bp (accounting for lower reverse read quality).
  - Maximum expected errors: `maxEE = c(2, 2)`.
  - Quality score threshold: `truncQ = 2`.
  - PhiX spike-in removal: `rm.phix = TRUE`.

### 2. Error Modeling & Sample Inference (DADA2)
- Trains parametric error models on filtered reads (`learnErrors`) to differentiate true biological sequence variants from Illumina sequencing errors.
- Infers true sample sequences (`dada`) with exact single-nucleotide resolution.

### 3. Pair Merging, Length Validation, and Chimera Removal
- Merges forward and reverse denoised reads (`mergePairs`).
- Constructs the full Amplicon Sequence Variant (ASV) abundance table.
- Applies strict biological length filtering (`399:430 bp`), removing non-specific amplification products.
- Identifies and discards PCR chimeras de novo (`removeBimeraDenovo`, `method = "consensus"`).
- Tracks read retention metrics across each step of the pipeline (`track`).

### 4. Taxonomic Assignment (SILVA v138.2)
- Assigns bacterial taxonomy from Kingdom down to Genus using the SILVA v138.2 reference training set (`assignTaxonomy`, RDP naive Bayesian classifier).
- Assigns exact species matches using the SILVA species assignment database (`addSpecies`).

### 5. Phyloseq Integration & Contaminant Filtering
- Combines the ASV count matrix, sample metadata (`CFI metadata.csv`), and taxonomy table into a unified `phyloseq` object.
- Filters out non-bacterial or organellar DNA:
  - Excludes Family **Mitochondria**.
  - Excludes Order **Chloroplast**.
  - Removes unassigned phyla.
- Segregates experimental animals from sequencing control samples (`ps.mice`).

### 6. Alpha Diversity Analysis
- Calculates richness and evenness estimators:
  - **Observed ASVs** (direct richness)
  - **Chao1** (abundance-based richness estimator)
  - **Shannon Index** (richness & evenness)
  - **Simpson Index** (dominance & evenness)
- Visualizes distributions via publication-ready boxplots with `ggplot2`.
- Performs rigorous non-parametric hypothesis testing:
  - **Kruskal-Wallis test** across experimental groups.
  - **Pairwise Wilcoxon rank-sum tests** with Benjamini-Hochberg (BH/FDR) $p$-value adjustment.

### 7. Beta Diversity & Multivariate Community Ecology
- **Quality noise pruning**: Prunes rare ASVs (< 10 reads total across experiment) and low-prevalence artifacts (present in < 2 samples).
- **Ordination**: Computes Bray-Curtis dissimilarity distance matrices and performs Principal Coordinates Analysis (PCoA).
- Generates 2D PCoA ordination plots overlaid with 95% confidence multivariate normal ellipses (`stat_ellipse(type = "t")`).
- **Hypothesis Testing**:
  - Global **PERMANOVA** (`vegan::adonis2`, 999 permutations) testing community structure divergence across treatment groups.
  - **Pairwise PERMANOVA** comparisons:
    - `Control vs Control-CFI` ($p = 0.032$): Identifies the direct microbial modulation triggered by CFI in healthy hosts.
    - `HFD vs HFD-CFI` ($p = 0.040$): Evaluates the therapeutic capacity of CFI to rescue or remodel HFD-induced dysbiosis.

### 8. Taxonomic Composition & Abundance Export
- Computes relative abundance transformations via Total Sum Scaling (TSS).
- Agglomerates taxa at each canonical rank (`tax_glom`): **Phylum, Class, Order, Family, Genus, Species**.
- Generates faceted stacked bar plots (`ggplot2` + `facet_grid(~Group)`).
- Implements an automated custom export function (`export_wide_taxonomy`) generating tidy wide-format abundance matrices with sample metadata for all taxonomic ranks.

### 9. Differential Abundance Testing (DESeq2)
- Converts raw counts to DESeq2 format (`phyloseq_to_deseq2`).
- Implements `poscounts` size factor estimation, specifically formulated to handle sparsity and zero-inflation in microbiome count data.
- Fits generalized linear models (GLMs) using local dispersion estimation (`fitType = "local"`).
- Tests key biological contrasts:
  1. `HFD-CFI vs HFD`: Discovers specific bacterial biomarkers shifted or restored by the intervention.
  2. `Control-CFI vs Control`: Quantifies changes elicited in normal physiological status.
  3. `HFD vs Control-CFI`: Evaluates distance between disease state and healthy intervention state.
- Annotates statistically significant taxa ($\text{padj} < 0.05$) with full SILVA taxonomy and orders by $\log_2(\text{Fold Change})$.

---

## Repository Structure

```text
16S-rRNA-Seq-v3-v4-analysis/
├── 16S rRNA V3-V4 R codes/
│   └── CFI 16s data analysis.R      # Primary R analysis script (DADA2 to DESeq2)
├── .gitignore                       # Ignored cache, FASTQ, and large binaries
├── LICENSE                          # MIT License
└── README.md                        # Complete project documentation & workflow guide
```

---

## Installation & Prerequisites

### R Environment
Recommended: **R >= 4.0.0** (tested on R 4.2+ / 4.3+).

### Required R Packages

```r
# 1. Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 2. Install core Bioconductor packages
BiocManager::install(c("dada2", "phyloseq", "DESeq2"))

# 3. Install CRAN packages
install.packages(c("ggplot2", "dplyr", "tidyr", "vegan"))
```

### Reference Databases
Download the SILVA reference database files and place them in your analysis working directory:
- [silva_nr99_v138.2_toGenus_trainset.fa.gz](https://zenodo.org/records/10700595)
- [silva_v138.2_assignSpecies.fa.gz](https://zenodo.org/records/10700595)

---

## Input Data Requirements

1. **Demultiplexed Paired-End FASTQ Files**:
   - Forward reads: `<sample_id>_R1_001.fastq.gz`
   - Reverse reads: `<sample_id>_R2_001.fastq.gz`
2. **Sample Metadata File (`CFI metadata.csv`)**:
   - A CSV spreadsheet whose row names correspond exactly to the sample identifiers extracted from FASTQ filenames.
   - Recommended columns:
     - `Sample`: Unique identifier.
     - `Name`: Formatted display name for visualizations.
     - `Group`: Treatment group assignment (`Control`, `Control-CFI`, `HFD`, `HFD-CFI`, `Negative`, `Positive`).

---

## Usage Guide

1. Clone or download this repository:
   ```bash
   git clone https://github.com/taiwo65-codes/16S-rRNA-Seq-v3-v4-analysis.git
   cd 16S-rRNA-Seq-v3-v4-analysis
   ```

2. Open `16S rRNA V3-V4 R codes/CFI 16s data analysis.R` in RStudio or an R session.

3. Update the working path at the top of the script:
   ```r
   path <- "path/to/your/fastq_and_metadata_directory"
   ```

4. Run the script interactively or sequentially by section to inspect diagnostic quality plots and verify intermediate outputs.

---

## Expected Outputs & Deliverables

| Category | File Name | Description |
| :--- | :--- | :--- |
| **DADA2 Intermediates** | `seqtab_dada2_raw.rds` | Cleaned ASV count abundance matrix |
| **Taxonomy Table** | `taxa_dada2_raw.rds` | SILVA v138.2 taxonomic classifications |
| **Alpha Diversity** | `alpha_diversities_complete.csv` | Observed, Chao1, Shannon, Simpson indices per sample |
| **Taxonomic Profiles** | `phylum_rel_abundance_wide.csv` | Wide-format relative abundance across samples at Phylum level |
| | `*_rel_abundance_wide.txt` | Wide-format abundance tables for Class, Order, Family, Genus, Species |
| **Differential Abundance** | `DESeq2_HFD_vs_HFD-CFI.csv` | Differentially abundant taxa in HFD with CFI rescue ($\text{padj} < 0.05$) |
| | `DESeq2_Control_vs_Control-CFI.csv` | Differentially abundant taxa in healthy controls with CFI |
| | `DESeq2_Control-CFI_vs_HFD.csv` | Differentially abundant taxa between control-treated and HFD |

---

## Author & Acknowledgments

- **Author**: **Taiwo Bankole**
- **GitHub**: [@taiwo65-codes](https://github.com/taiwo65-codes)
- **Repository**: [16S-rRNA-Seq-v3-v4-analysis](https://github.com/taiwo65-codes/16S-rRNA-Seq-v3-v4-analysis)

---

## License

This project is licensed under the [MIT License](LICENSE) - feel free to adapt and use for academic and research workflows.
