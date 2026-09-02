# 🧪 NextLAMP Viral Benchmark & Validation Suite

This directory contains the isolated benchmark environment to validate **NextLAMP** performance, specificity, and sensitivity using viral genomes (**Zika virus - ZIKV** and **Dengue virus - DENV**).

## 📁 Directory Structure
- `data/`: Downloaded FASTA genomes, target/background lists, and Bowtie 2 indices.
- `configs/`: YAML configurations optimized for viral genome whole-genome LAMP design.
- `results/`: Publication packages, primer TSV tables, JSON bundles, and validation reports.
- `logs/`: Detailed execution logs.
- `scripts/`: Benchmark execution scripts.

## 🎯 Target & Background Setup

### Dataset 1: Zika Virus (ZIKV)
- **Target Taxon:** *Zika virus* (TaxID: `64320`)
- **Background Taxa:**
  - *Homo sapiens* (GRCh38 host)
  - *Aedes aegypti* (Vector)
  - *Dengue virus* (Flavivirus cross-reactivity check, TaxID: `12637`)
  - *Yellow fever virus* (TaxID: `11089`)

### Dataset 2: Dengue Virus Type-Specific (DENV-2)
- **Target Taxon:** *Dengue virus type 2* (TaxID: `11060`)
- **Background Taxa:**
  - *Dengue virus type 1, 3, 4*
  - *Homo sapiens*
