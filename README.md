# 🧬 NextLAMP: Whole-Genome LAMP Primer Design System

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Conda Environment](https://img.shields.io/badge/conda-environment.yml-green.svg)](environment.yml)
[![FAIR Compliant](https://img.shields.io/badge/FAIR-reproducible-orange.svg)](#-fair-principles--reproducibility)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

> **NextLAMP** is a modern, high-performance, FAIR-compliant software platform for automated **Loop-Mediated Isothermal Amplification (LAMP)** primer design at whole-genome scale.

---

## ⚡ Quick Setup (3 Simple Steps)

Open your terminal and run:

```bash
# 1. Clone the repository
git clone https://github.com/user/nextlamp.git
cd nextlamp

# 2. Create and activate the Conda environment (installs Python, Bowtie 2, and dependencies automatically)
conda env create -f environment.yml
conda activate nextlamp

# 3. Install NextLAMP in editable mode
pip install -e .
```

Done! The `nextlamp` and `nextlamp-prep` commands are now ready to use in your terminal.

---

## 🎯 Quick Start Guide (Choose Your Workflow)

### 🌟 Option A: NCBI Automated Download + Design (In 1 Command!)

Don't have local FASTA files downloaded? Specify taxonomic group names or TaxIDs and NextLAMP will fetch genomes from NCBI, format databases, build Bowtie 2 indices, and design primers automatically:

```bash
nextlamp-prep \
    --target-taxa "Babesia canis" \
    --common-taxa "Babesia" \
    --background-taxa "Apicomplexa" "Canis lupus familiaris" \
    --out-dir my_dataset \
    --run-nextlamp
```

---

### 📝 Option B: Using YAML Configuration Files (Recommended)

Ideal for reproducible, publishable scientific research:

```bash
# 1. Generate an annotated YAML configuration template
nextlamp --generate-config my_design.yaml

# 2. (Optional) Edit my_design.yaml if you wish to adjust temperatures or distance constraints

# 3. Run NextLAMP using the YAML file
nextlamp --config my_design.yaml
```

---

### 💻 Option C: Python API (Only 5 Lines of Code)

```python
from nextlamp import NextLampPipeline, export_results

# 1. Initialize the pipeline
pipeline = NextLampPipeline(
    target_fasta="data/target_babesia_canis.fa",
    index_prefix="data/db_completo_idx",
    targets_list_file="data/targets_list.txt",
    background_list_file="data/background_list.txt"
)

# 2. Run primer design (with 4 threads)
results, params, stats = pipeline.run(threads=4)

# 3. Export complete FAIR reproducibility bundle (JSON, TSV for Excel, TXT)
export_results(results, params, stats, out_json="results.json")
```

---

## 📂 Generated Output Artifacts

Each NextLAMP run produces **3 ready-to-use output files**:

| Output File | Purpose & Usage |
| :--- | :--- |
| **`results_primers.tsv`** | **For Laboratory Synthesis & Ordering / Excel:** Clean tabular spreadsheet with primer names, 5'→3' sequences, lengths, Tm, and GC%. |
| **`results_report.txt`** | **For Fast Human Interpretation:** Formatted summary report ranking top primer sets by thermal balance (`tm_balance`). |
| **`results.json`** | **For Scientific Reproducibility:** Immutable JSON bundle recording SHA-256 target hash, execution metadata, parameters, and selection funnels. |

---

## 🛡️ FAIR Principles & Reproducibility

NextLAMP generates a complete **FAIR Reproducibility Bundle** for every run:
- **Findable (F):** Digital provenance via target genome SHA-256 cryptographic hash (`target_fasta_sha256`) and ISO 8601 timestamp.
- **Accessible (A):** Three standard output artifacts (JSON metadata, TSV synthesis table, TXT report).
- **Interoperable (I):** Standard bioinformatic formats (FASTA, SAM, 5'→3' orientation, IUPAC codes).
- **Reusable (R):** Immutable record of all 17 thermodynamic, geometric, and distance parameters.

---

## 🛠️ CLI Reference Guide

### 1. Data Preparation Utility (`nextlamp-prep`)

| Argument | Example | Description |
| :--- | :--- | :--- |
| `--target-taxa` | `"Babesia canis"` | Primary target taxonomic group(s) |
| `--common-taxa` | `"Babesia"` | Common target taxonomic group(s) to conserve |
| `--background-taxa` | `"Apicomplexa" "Canis lupus familiaris"` | Background/host taxonomic group(s) to exclude |
| `--out-dir` | `my_dataset` | Directory to store genomes and Bowtie 2 index |
| `--run-nextlamp` | *(flag)* | Automatically trigger NextLAMP pipeline after data prep |

### 2. Primer Design Utility (`nextlamp`)

| Argument | Example | Description |
| :--- | :--- | :--- |
| `--config` | `config.yaml` | Run design pipeline using YAML config file |
| `--generate-config` | `config.yaml` | Generate annotated YAML configuration template |
| `--target-fasta` | `target.fa` | Path to target genome FASTA file |
| `--index-prefix` | `db_idx` | Path prefix of Bowtie 2 index |
| `--threads` | `4` | Number of CPU threads for parallel processing |
| `--out` | `results.json` | Output JSON file name |

---

## 🖥️ HPC Cluster Execution (SLURM)

To submit batch jobs to high-performance computing clusters:

```bash
sbatch scripts/nextlamp/run_nextlamp_slurm.sh
```

---

## ❓ Frequently Asked Questions (Troubleshooting)

- **Error: `bowtie2 binary not found`**
  - *Fix:* Make sure to activate the Conda environment (`conda activate nextlamp`). `bowtie2` is installed automatically inside the environment.
- **How to adjust Tm or GC thresholds?**
  - *Fix:* In your `config.yaml` file, adjust `min_gc`, `max_gc`, `min_tm`, and `max_tm`, or pass CLI arguments `--min-gc 30.0 --max-gc 70.0`.
