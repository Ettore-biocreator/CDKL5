# CDKL5 Transcriptomics & Multi-Omics Analysis Pipeline

This repository contains a collection of general-purpose scripts used to perform RNA-seq preprocessing and downstream transcriptomic/multi-omics analyses in the context of CDKL5 deficiency disorder research.

> **Note:** These scripts represent the core analytical workflows and are intended as a general framework. They do not constitute the entirety of the analyses performed in the study. Parameters, file paths, and sample-specific settings will need to be adapted to your own data and environment.

---

## Data Availability

Raw and processed data will be deposited on [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/) upon publication of the associated article. This README will be updated with the accession number once available.

---

## Repository Structure

```
├── Preprocessing/
│   ├── build_STAR_index.sh           # Build STAR genome index from FASTA + GTF
│   └── preprocessing_fastq.sh        # Quality trimming (BBDuk), STAR alignment, and read counting (HTSeq)
│
└── R/
    ├── DESeq2_analysis.R             # Differential gene expression analysis (DESeq2)
    ├── DESeq2_exploration.R          # QC and data exploration: PCA, heatmaps, volcano plots, GSEA
    └── cross_omics_intersection.R    # Cross-species, multi-omics gene set intersection and Venn diagram
```

---

## Requirements

### Preprocessing (shell scripts)
- [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) (part of BBTools)
- [STAR](https://github.com/alexdobin/STAR)
- [SAMtools](http://www.htslib.org/)
- [HTSeq](https://htseq.readthedocs.io/)
- Conda (for environment management)

### R scripts
- R ≥ 4.0
- `DESeq2`, `ggplot2`, `pheatmap`, `RColorBrewer`
- `dplyr`, `tidyr`, `tibble`, `purrr`, `tools`
- `clusterProfiler`, `enrichplot`, `org.Hs.eg.db` / `org.Mm.eg.db`
- `VennDiagram`, `grid`, `futile.logger`
- `openxlsx`

---

## Usage

### 1. Build STAR Genome Index
Edit the path variables at the top of `build_STAR_index.sh` to point to your genome FASTA and annotation GTF, then run:
```bash
bash Preprocessing/build_STAR_index.sh
```

### 2. Preprocessing (Trimming → Alignment → Counting)
Edit the path and parameter variables at the top of `preprocessing_fastq.sh`, then run:
```bash
bash Preprocessing/preprocessing_fastq.sh
```
This will produce trimmed FASTQ files, STAR-aligned BAM files, and HTSeq read count tables.

### 3. Differential Expression Analysis
Open `R/DESeq2_analysis.R` and set the file paths and parameters in the configuration section at the top of the script. The script will output a tab-separated DEG results table.

### 4. Exploratory Analysis
`R/DESeq2_exploration.R` provides modular code blocks for:
- PCA plots
- Sample distance heatmaps
- MA plots
- Volcano plots (with optional gene highlighting)
- Expression heatmaps
- GSEA (Gene Ontology, Biological Process)

### 5. Cross-Omics Intersection
`R/cross_omics_intersection.R` integrates transcriptomics and proteomics data across multiple species (mouse, zebrafish, human cell lines), maps genes to human orthologs via the [Alliance of Genome Resources](https://www.alliancegenome.org/) orthology table, and generates a triple Venn diagram of overlapping gene sets.

---

## Notes

- All scripts are designed to be run after setting the variables in the dedicated configuration section at the top of each file. No hardcoded paths are present in the body of the scripts.
- SLURM job directives (`#SBATCH`) are included in the shell scripts for HPC cluster submission but can be ignored for local execution.
- The GSEA section in `DESeq2_exploration.R` assumes human gene symbols (`org.Hs.eg.db`). Change the `gsea_orgdb` parameter accordingly for other species.

---

## Citation

> Manuscript in preparation. Citation and GEO accession number will be added upon publication.

---

## Contact

For questions or issues regarding the scripts, please open an issue in this repository.
