# 🔬 Scientific Methodology: Genome-Scale LAMP Primer Design with NextLAMP

This document details the **experimental and bioinformatic methodology** for designing Loop-Mediated Isothermal Amplification (**LAMP**) primers and **Loop** primers (LoopF and LoopB) targeting the genus *_Babesia_* with multi-tiered specificity filtering against hosts, vectors, and related protozoa.

---

## 🎯 1. Scope & Assay Objectives

The LAMP assay aims for sensitive and specific detection of species within the genus *_Babesia_* (*B. canis*, *B. vogeli*, *B. gibsoni*, etc.), enabling real-time diagnostics in vertebrate host clinical samples as well as direct entomological/epidemiological surveillance in tick vectors.

---

## 🛡️ 2. Segmented Database Architecture & Early-Exit Filtering

To guarantee **zero cross-reactivity** in field Point-of-Care (PoC) assays or direct blood/vector extractions without processing giant monolithic files, NextLAMP's background exclusion database is structured into **segmented, independent database modules**:

### Pillar A: Vertebrate Hosts & Clinical Contaminants (Modular Indices)
Prevents non-specific binding with host genomic DNA and common kit/handling contaminants in blood and tissue samples:
- ***Canis lupus familiaris*** (`idx_dog` - Domestic Dog, RefSeq `GCF_011100685.1`)
- ***Felis catus*** (`idx_cat` - Domestic Cat)
- ***Homo sapiens*** (`idx_human` - Human Reference Genome GRCh38.p14)

### Pillar B: Tick Vectors (Vector Index - Ixodidae)
Prevents cross-reactivity in tick homogenates collected in the field or directly from animal skin:
- ***Rhipicephalus sanguineus***, ***Dermacentor reticulatus***, ***Rhipicephalus microplus***, ***Ixodes scapularis*** (`idx_ticks`)

### Pillar C: Phylogenetically Related Protozoa (Apicomplexa Index)
Prevents false positives due to evolutionary conservation with other hemoparasites and coccidia:
- *Theileria* spp., *Plasmodium* spp., *Toxoplasma gondii*, *Neospora caninum*, *Cryptosporidium* spp., *Eimeria* spp. (`idx_apicomplexa`)

### ⚡ Sequential Early-Exit Filtering
NextLAMP evaluates these indices sequentially using **early-exit short-circuiting**: if a candidate primer collides with the dog, cat, or human host genome in an early index, it is immediately discarded from the pipeline. This saves computational time and maintains a minimal RAM footprint (<500 MB).

---

## ⚙️ 3. Thermodynamic Parameters & Spatial Distance Constraints

The **NextLAMP** algorithm evaluates 8 oligonucleotides per candidate set (F3, F2, F1c, LoopF, B1c, LoopB, B2, B3) under the following constraints:

| Parameter / Oligonucleotide | Accepted Range / Threshold | Biological Function |
| :--- | :---: | :--- |
| **Melting Temperature (\(T_m\))** | \(55.0^\circ\text{C} - 68.0^\circ\text{C}\) | Optimal isothermal reaction range for *Bst* Polymerase |
| **GC Content (\(GC\%\))** | \(30.0\% - 70.0\%\) | Duplex stability |
| **Thermodynamic Balance (\(T_{m\text{ balance}}\))** | Minimize \(|T_m(F2)-T_m(B2)| + |T_m(F3)-T_m(B3)|\) | Symmetric and efficient amplicon synthesis |
| **Distance F3 → F2** | \(0 - 20\text{ bp}\) | Cleavage & displacement region for outer primer |
| **Distance F2 → F1c** | \(40 - 60\text{ bp}\) | Dumbbell loop structure formation |
| **Core Amplicon Size (F2 → B2)** | \(120 - 180\text{ bp}\) | Core isothermal amplification region |
| **Distance B1c → B2** | \(40 - 60\text{ bp}\) | Right dumbbell loop structure formation |
| **Distance B2 → B3** | \(0 - 20\text{ bp}\) | Right outer strand displacement |
| **Loop Primers (LoopF / LoopB)** | Regions between F1-F2 and B1-B2 | Reaction acceleration (amplification in < 20 min) |

---

## 🎯 4. Pan-Genome Target Conservation & Biological 3' Anchor Specificity

Instead of synthesizing a simple consensus sequence (which introduces IUPAC degenerate bases and reduces isothermal reaction sensitivity), NextLAMP implements **pan-genome target validation with biological 3' anchor specificity**:

1. **1-bp High-Resolution Candidate Generation:** Raw candidate primers are generated with base-by-base (`step = 1bp`) resolution from a high-quality reference genome (`target_fasta`) across both strands.
2. **Pan-Coverage Validation:** During Bowtie 2 alignment against target accessions (`targets_list.txt`), candidates are tested across target assemblies with a configurable coverage fraction (`min_target_coverage: 0.60`).
3. **Biological 3' End Anchor Rule:** To accommodate natural genomic divergence across strains without sacrificing amplification efficiency, NextLAMP permits up to **1 mismatch** in internal/5' regions, while strictly enforcing **0 mismatches in the 5-bp 3' end anchor region** required for *Bst* DNA polymerase extension.
4. **Biological Advantage:** Guarantees that designed primers are 100% conserved across all known pathogen isolates without requiring degenerate primers.

---

## ⚡ 5. Alignment Mechanism & FAIR Provenance

1. **Candidate Scanning:** Sliding-window evaluation across the reference genome, identifying viable kmers meeting thermodynamic constraints.
2. **Multi-Index Specificity Filtering:** Sequential early-exit Bowtie 2 alignment against segmented host, vector, and background indices, discarding non-specific candidates.
3. **Locus Deduplication:** Removal of spatial overlap redundancies to ensure broad coverage of unique genomic loci.
4. **FAIR Export Bundle:** Immutable export including JSON (SHA-256 target hash, parameter audit log), TSV (ready for laboratory ordering), and TXT (summary report).
