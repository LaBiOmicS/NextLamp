import os
import zipfile
import glob
import subprocess
import json

# Define paths
base_dir = "/home/fabiano.menegidio/workdir/Omics/genomics/babesia"
data_dir = os.path.join(base_dir, "data")
os.makedirs(data_dir, exist_ok=True)

datasets_bin = "/home/fabiano.menegidio/miniforge3/envs/ncbi_env/bin/datasets"
apicomplexa_extract = os.path.join(data_dir, "apicomplexa_raw")
dog_extract = os.path.join(data_dir, "dog_raw")

# 1. Load assembly report to identify all Babesia accessions
report_path = os.path.join(apicomplexa_extract, "ncbi_dataset", "data", "assembly_data_report.jsonl")
babesia_accessions = set()
if os.path.exists(report_path):
    print("Loading assembly data report to identify Babesia genomes...")
    with open(report_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            data = json.loads(line)
            org = data.get("organism", {})
            org_name = org.get("organismName", "")
            if "Babesia" in org_name:
                acc = data.get("accession")
                if acc:
                    babesia_accessions.add(acc)
    print(f"Identified {len(babesia_accessions)} Babesia genomes/accessions.")
else:
    print(f"Error: {report_path} not found. Please ensure Apicomplexa genomes are downloaded and extracted.")
    exit(1)

import re

def normalize_acc(acc_str):
    if acc_str.startswith("GCA_"):
        return "GCF_" + acc_str[4:]
    return acc_str

# Find all FNA files and deduplicate GCA/GCF pairs (prefer RefSeq GCF_)
all_apicomplexa_fnas = glob.glob(os.path.join(apicomplexa_extract, "ncbi_dataset", "data", "GCA_*", "*.fna")) + \
                       glob.glob(os.path.join(apicomplexa_extract, "ncbi_dataset", "data", "GCF_*", "*.fna"))

fna_by_num_id = {}
for fna in all_apicomplexa_fnas:
    match = re.search(r'G[CFL]_[0-9]{9}\.[0-9]+', fna)
    acc = match.group(0) if match else os.path.basename(fna)
    num_id = re.sub(r'^(GCA_|GCF_)', '', acc)
    if num_id not in fna_by_num_id or acc.startswith("GCF_"):
        fna_by_num_id[num_id] = fna

apicomplexa_fnas = sorted(list(fna_by_num_id.values()))
dog_fnas = glob.glob(os.path.join(dog_extract, "ncbi_dataset", "data", "GCF_*", "*.fna"))

target_fasta_out = os.path.join(data_dir, "target_babesia_canis.fa")  # Template sequence for candidate generation (B. canis)
db_completo_out = os.path.join(data_dir, "db_completo.fa")

targets_list_file = os.path.join(data_dir, "targets_list.txt")
background_list_file = os.path.join(data_dir, "background_list.txt")

target_headers = []
background_headers = []

print("Processing and sorting genomes into Babesia (targets) and non-Babesia (background)...")

def read_fasta_as_pseudogenome(fna_path, spacer_len=100):
    contigs = []
    current_seq = []
    with open(fna_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current_seq:
                    contigs.append("".join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line.strip())
        if current_seq:
            contigs.append("".join(current_seq))
    return ("N" * spacer_len).join(contigs)

# Open files to write target list, background list, and full database
with open(db_completo_out, "w") as db_out, open(target_fasta_out, "w") as template_out:
    # Process Apicomplexa
    for fna in apicomplexa_fnas:
        # Check if this fna belongs to Babesia and capture accession
        matched_acc = None
        for acc in babesia_accessions:
            if acc in fna:
                matched_acc = acc
                break
        
        if not matched_acc:
            # Extract assembly accession from directory name
            parts = fna.split(os.sep)
            for part in parts:
                if part.startswith("GCA_") or part.startswith("GCF_"):
                    matched_acc = part
                    break
            if not matched_acc:
                matched_acc = os.path.basename(fna)

        is_babesia = matched_acc in babesia_accessions
        is_template = "GCA_045269395.1" in fna
        norm_acc = normalize_acc(matched_acc)
        
        pseudo_seq = read_fasta_as_pseudogenome(fna, spacer_len=100)
        
        if is_babesia:
            target_headers.append((norm_acc, norm_acc))
            db_out.write(f">{norm_acc}\n{pseudo_seq}\n")
            if is_template:
                template_out.write(f">{norm_acc}\n{pseudo_seq}\n")
        else:
            background_headers.append((norm_acc, norm_acc))
            db_out.write(f">{norm_acc}\n{pseudo_seq}\n")
            if is_template:
                template_out.write(f">{norm_acc}\n{pseudo_seq}\n")
                        
    # Process Dog (Dog is always background)
    for fna in dog_fnas:
        matched_acc = "Dog_GCF_011100685.1"
        pseudo_seq = read_fasta_as_pseudogenome(fna, spacer_len=100)
        background_headers.append((matched_acc, matched_acc))
        db_out.write(f">{matched_acc}\n{pseudo_seq}\n")

# Write lists
with open(targets_list_file, "w") as f:
    for h, acc in target_headers:
        f.write(f"{h}\t{acc}\n")

with open(background_list_file, "w") as f:
    for h, acc in background_headers:
        f.write(f"{h}\t{acc}\n")

unique_target_genomes = len(set(acc for _, acc in target_headers))
print(f"Total target sequences: {len(target_headers)} pseudo-genomes spanning {unique_target_genomes} unique Babesia assemblies, written to {targets_list_file}")
print(f"Total background sequences: {len(background_headers)} pseudo-genomes written to {background_list_file}")
print("All tasks completed successfully!")
