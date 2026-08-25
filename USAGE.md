# 📖 NextLAMP General User Guide & Application Manual

> **NextLAMP** is a universal, open-source bioinformatic software system designed to automate whole-genome **Loop-Mediated Isothermal Amplification (LAMP)** primer design for any target organism (viruses, bacteria, fungi, parasites, or host genes).

---

## 📌 Table of Contents
1. [Introduction & General Applications](#1-introduction--general-applications)
2. [Installation Guide](#2-installation-guide)
3. [Universal Workflows](#3-universal-workflows)
   - [Workflow 1: Automated NCBI Fetch & Design](#workflow-1-automated-ncbi-fetch--design)
   - [Workflow 2: Custom Local FASTA Files](#workflow-2-custom-local-fasta-files)
   - [Workflow 3: Reproducible YAML Configuration](#workflow-3-reproducible-yaml-configuration)
4. [Python API Integration](#4-python-api-integration)
5. [Understanding Output Results](#5-understanding-output-results)
6. [Parameters & Customization Reference](#6-parameters--customization-reference)

---

## 1. Introduction & General Applications

Isothermal LAMP assays require 6 distinct primer regions ($F3$, $F2$, $F1c$, $B1c$, $B2$, $B3$) with strict spatial distance constraints and precise melting temperatures ($T_m$). NextLAMP automates candidate generation, secondary structure checks (hairpins & self-dimers), specificity alignment against background genomes, and thermal balance optimization.

### Universal Use Cases:
- 🦠 **Virology:** Rapid diagnostic design for viral genomes (e.g., Dengue virus, Zika, Influenza, SARS-CoV-2).
- 🧫 **Bacteriology:** Species-specific or strain-specific assay design (e.g., *Escherichia coli*, *Staphylococcus aureus*, *Salmonella*).
- 🔬 **Parasitology & Mycology:** AT-rich or complex eukaryotic genomes (e.g., *Plasmodium*, *Babesia*, *Candida*).
- 🧬 **Genetics & Oncology:** Point mutation or gene-specific amplification.

---

## 2. Installation Guide

NextLAMP is packaged as a Conda environment containing all bioinformatic dependencies (`bowtie2`, `python`, `biopython`, `ncbi-datasets-cli`).

```bash
# Clone repository
git clone https://github.com/user/nextlamp.git
cd nextlamp

# Create Conda environment
conda env create -f environment.yml
conda activate nextlamp

# Install NextLAMP executable binaries
pip install -e .
```

Verify installation:
```bash
nextlamp --help
nextlamp-prep --help
```

---

## 3. Universal Workflows

### Workflow 1: Automated NCBI Fetch & Design

If you do not have local FASTA files, provide the NCBI organism names or TaxIDs. NextLAMP fetches genomes from NCBI, formats the databases, builds alignment indices, and executes primer design in one step:

```bash
nextlamp-prep \
    --target-taxa "Dengue virus 1" \
    --common-taxa "Dengue virus" \
    --background-taxa "Zika virus" "Homo sapiens" \
    --out-dir dengue_assay_dataset \
    --threads 4 \
    --run-nextlamp
```

---

### Workflow 2: Custom Local FASTA Files

When working with your own custom genomic data:

```bash
nextlamp \
    --target-fasta /path/to/my_target.fasta \
    --index-prefix /path/to/bowtie2_index_prefix \
    --targets-list /path/to/target_headers.txt \
    --background-list /path/to/background_headers.txt \
    --out my_custom_results.json \
    --threads 4
```

---

### Workflow 3: Reproducible YAML Configuration

For publishable and team-wide reproducible research, use declarative YAML files:

1. Generate a template:
```bash
nextlamp --generate-config my_experiment.yaml
```

2. Edit `my_experiment.yaml`:
```yaml
target_fasta: "data/my_organism.fa"
index_prefix: "data/genome_db_idx"
targets_list: "data/targets.txt"
background_list: "data/background.txt"
out: "experiment_results.json"
threads: 8

min_gc: 35.0
max_gc: 65.0
min_tm: 56.0
max_tm: 66.0
check_dimers: true
```

3. Execute:
```bash
nextlamp --config my_experiment.yaml
```

---

## 4. Python API Integration

Integrate NextLAMP directly into Python data science pipelines or Jupyter Notebooks:

```python
from nextlamp import prepare_nextlamp_dataset, NextLampPipeline, export_results

# Step 1 (Optional): Download dataset from NCBI
dataset = prepare_nextlamp_dataset(
    target_taxa=["Escherichia coli"],
    background_taxa=["Salmonella enterica", "Homo sapiens"],
    output_dir="ecoli_dataset"
)

# Step 2: Initialize & Run Pipeline
pipeline = NextLampPipeline(
    target_fasta=dataset["target_fasta"],
    index_prefix=dataset["index_prefix"],
    targets_list_file=dataset["targets_list"],
    background_list_file=dataset["background_list"]
)

results, params, stats = pipeline.run(min_gc=35.0, max_gc=65.0, threads=4)

# Step 3: Save Reproducibility Bundle
export_results(results, params, stats, out_json="ecoli_results.json")
```

---

## 5. Understanding Output Results

NextLAMP automatically outputs 3 formatted files:

1. **`results_primers.tsv` (Laboratory Ordering Table)**
   - Tab-separated spreadsheet containing primer names ($F3, F2, F1c, B1c, B2, B3$), 5'→3' sequences, positions, lengths, $T_m$ (°C), and GC%. Ready to copy-paste for oligo synthesis.

2. **`results_report.txt` (Human Summary Report)**
   - Human-readable summary detailing selection funnel statistics and ranking primer sets by thermal balance score ($tm\_balance = |T_m(F2) - T_m(B2)| + |T_m(F3) - T_m(B3)|$).

3. **`results.json` (Scientific Reproducibility Bundle)**
   - Machine-readable JSON recording cryptographic target FASTA SHA-256 hash, execution timestamp, complete parameter dictionary, and all candidate evaluation metrics.

---

## 6. Parameters & Customization Reference

| Parameter Key | CLI Flag | Default | Description |
| :--- | :--- | :---: | :--- |
| `target_fasta` | `--target-fasta` | *Required* | Path to primary target FASTA file |
| `index_prefix` | `--index-prefix` | *Required* | Bowtie 2 index path prefix |
| `targets_list` | `--targets-list` | *Required* | File with sequence headers of target genomes |
| `background_list` | `--background-list` | *Required* | File with sequence headers of background genomes |
| `min_gc` | `--min-gc` | `30.0` | Minimum primer GC content (%) |
| `max_gc` | `--max-gc` | `70.0` | Maximum primer GC content (%) |
| `min_tm` | `--min-tm` | `55.0` | Minimum primer melting temperature ($T_m$, °C) |
| `max_tm` | `--max-tm` | `68.0` | Maximum primer melting temperature ($T_m$, °C) |
| `dist_inner_min` | `--dist-inner-min` | `120` | Minimum core amplicon size ($F2$ to $B2$, bp) |
| `dist_inner_max` | `--dist-inner-max` | `180` | Maximum core amplicon size ($F2$ to $B2$, bp) |
| `threads` | `--threads` | `4` | Number of CPU threads for parallel alignment |
