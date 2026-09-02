import json
import os
import re

DEFAULT_CONFIG = {
    "target_fasta": "data/target_babesia_canis.fa",
    "index_prefix": "data/db_completo_idx",
    "targets_list": "data/targets_list.txt",
    "background_list": "data/background_list.txt",
    "out": "nextlamp_success.json",
    "threads": 4,
    "max_sets": 10,
    "min_gc": 30.0,
    "max_gc": 70.0,
    "min_tm": 55.0,
    "max_tm": 68.0,
    "dist_f3_f2_min": 0,
    "dist_f3_f2_max": 20,
    "dist_f2_f1c_min": 40,
    "dist_f2_f1c_max": 60,
    "dist_inner_min": 120,
    "dist_inner_max": 180,
    "dist_b1c_b2_min": 40,
    "dist_b1c_b2_max": 60,
    "dist_b2_b3_min": 0,
    "dist_b2_b3_max": 20,
    "dist_f1c_b1c_max": 85,
    "min_tm_diff": 3.0,
    "check_dimers": True,
    "include_loops": True,
    "min_target_coverage": 1.0,
    "min_targets_count": None,
    "data_prep": {
        "target_taxa": ["Babesia canis"],
        "common_taxa": ["Babesia"],
        "background_taxa": ["Apicomplexa", "Canis lupus familiaris"],
        "max_genomes_per_taxon": 20,
        "out_dir": "dataset_babesia"
    }
}

def load_config(config_path: str) -> dict:
    """
    Loads configuration from a YAML or JSON file.
    """
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    ext = os.path.splitext(config_path)[1].lower()

    if ext in ('.json', '.js'):
        with open(config_path, 'r', encoding='utf-8') as f:
            return json.load(f)

    # Lightweight YAML / Key-Value Parser for YAML files without pyyaml dependency
    config_dict = {}
    current_section = None
    with open(config_path, 'r', encoding='utf-8') as f:
        for line in f:
            line_str = line.strip()
            if not line_str or line_str.startswith('#'):
                continue

            # Key-value mapping
            if ':' in line_str:
                parts = line_str.split(':', 1)
                key = parts[0].strip()
                val = parts[1].split('#')[0].strip()

                if not val:
                    current_section = key
                    config_dict[current_section] = {}
                    continue

                # Type casting
                if val.lower() in ('true', 'yes'):
                    typed_val = True
                elif val.lower() in ('false', 'no'):
                    typed_val = False
                elif val.startswith('[') and val.endswith(']'):
                    items = val[1:-1].split(',')
                    typed_val = [it.strip().strip('"').strip("'") for it in items if it.strip()]
                elif re.match(r'^-?\d+$', val):
                    typed_val = int(val)
                elif re.match(r'^-?\d+\.\d+$', val):
                    typed_val = float(val)
                else:
                    typed_val = val.strip('"').strip("'")

                if current_section and line.startswith(('  ', '\t')):
                    config_dict[current_section][key] = typed_val
                else:
                    current_section = None
                    config_dict[key] = typed_val
            elif line_str.startswith('- ') and current_section:
                item = line_str[2:].strip().strip('"').strip("'")
                if not isinstance(config_dict.get(current_section), list):
                    config_dict[current_section] = []
                config_dict[current_section].append(item)

    return config_dict

def generate_default_config_yaml(out_path: str):
    """
    Generates a clean, annotated YAML configuration file template for NextLAMP.
    """
    content = """# ========================================================================
#   NextLAMP Configuration Template (YAML)
#   Annotated parameters for reproducible Whole-Genome LAMP primer design
# ========================================================================

# --- Inputs & Outputs ---
target_fasta: "data/target_babesia_canis.fa"
# Specify a single index prefix or a list of segmented index prefixes for sequential early-exit filtering:
# index_prefix:
#   - "data/indexes/idx_dog"
#   - "data/indexes/idx_cat"
#   - "data/indexes/idx_ticks"
#   - "data/indexes/idx_apicomplexa"
index_prefix: "data/db_completo_idx"
targets_list: "data/targets_list.txt"
background_list: "data/background_list.txt"
out: "nextlamp_success.json"
threads: 4
max_sets: 10

# --- Thermodynamic Thresholds ---
min_gc: 30.0         # Minimum GC percentage for primers (%)
max_gc: 70.0         # Maximum GC percentage for primers (%)
min_tm: 55.0         # Minimum melting temperature (°C)
max_tm: 68.0         # Maximum melting temperature (°C)
min_tm_diff: 3.0     # Minimum Tm difference between F1c and B1c (°C)
check_dimers: true   # Enable self-dimer and hairpin filtering

# --- Spatial & Amplicon Distance Constraints (bp) ---
dist_f3_f2_min: 0    # F3 to F2 min distance
dist_f3_f2_max: 20   # F3 to F2 max distance
dist_f2_f1c_min: 40  # F2 to F1c min distance
dist_f2_f1c_max: 60  # F2 to F1c max distance
dist_inner_min: 120  # Core amplicon min distance (F2 to B2)
dist_inner_max: 180  # Core amplicon max distance (F2 to B2)
dist_b1c_b2_min: 40  # B1c to B2 min distance
dist_b1c_b2_max: 60  # B1c to B2 max distance
dist_b2_b3_min: 0    # B2 to B3 min distance
dist_b2_b3_max: 20   # B2 to B3 max distance
dist_f1c_b1c_max: 85 # Max distance between F1c and B1c

# --- Automated NCBI Data Prep (Optional Initial Step) ---
data_prep:
  target_taxa: ["Babesia canis"]
  common_taxa: ["Babesia"]
  background_taxa: ["Apicomplexa", "Canis lupus familiaris"]
  max_genomes_per_taxon: 20
  out_dir: "dataset_babesia"
"""
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"[OK] Generated default configuration template: {out_path}")

def generate_default_prep_config_yaml(out_path: str):
    """
    Generates an annotated YAML configuration file template specifically for Data Prep.
    """
    content = """# ========================================================================
#   NextLAMP Data Preparation Configuration (YAML)
# ========================================================================

out_dir: "dataset_babesia"
threads: 4
max_genomes_per_taxon: 20
run_nextlamp: true

# --- Option A: Specify Taxonomic Group Names or TaxIDs ---
target_taxa:
  - "Babesia canis"

common_taxa:
  - "Babesia"

background_taxa:
  - "Apicomplexa"
  - "Canis lupus familiaris"

# --- Option B: Specify Local Accession List Files (CSV / TSV / TXT) ---
# targets_list: "data/targets_list.txt"
# common_list: "data/common_list.txt"
# background_list: "data/background_list.txt"
"""
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"[OK] Generated Data Prep configuration template: {out_path}")
