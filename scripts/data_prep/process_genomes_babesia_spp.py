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

# Find all FNA files
apicomplexa_fnas = glob.glob(os.path.join(apicomplexa_extract, "ncbi_dataset", "data", "GCA_*", "*.fna")) + \
                   glob.glob(os.path.join(apicomplexa_extract, "ncbi_dataset", "data", "GCF_*", "*.fna"))
dog_fnas = glob.glob(os.path.join(dog_extract, "ncbi_dataset", "data", "GCF_*", "*.fna"))

target_fasta_out = os.path.join(data_dir, "target_babesia_canis.fa")  # Template sequence for candidate generation (B. canis)
db_completo_out = os.path.join(data_dir, "db_completo.fa")

targets_list_file = os.path.join(data_dir, "targets_list.txt")
background_list_file = os.path.join(data_dir, "background_list.txt")

target_headers = []
background_headers = []

print("Processing and sorting genomes into Babesia (targets) and non-Babesia (background)...")

# Open files to write target list, background list, and full database
with open(db_completo_out, "w") as db_out, open(target_fasta_out, "w") as template_out:
    # Process Apicomplexa
    for fna in apicomplexa_fnas:
        # Check if this fna belongs to Babesia
        is_babesia = False
        for acc in babesia_accessions:
            if acc in fna:
                is_babesia = True
                break
        
        # Check if this is the template Babesia canis (GCA_045269395.1)
        is_template = "GCA_045269395.1" in fna
        
        with open(fna, "r") as f:
            for line in f:
                if line.startswith(">"):
                    header = line[1:].strip().split()[0]
                    if is_babesia:
                        target_headers.append(header)
                    else:
                        background_headers.append(header)
                    db_out.write(line)
                    if is_template:
                        template_out.write(line)
                else:
                    db_out.write(line)
                    if is_template:
                        template_out.write(line)
                        
    # Process Dog (Dog is always background)
    for fna in dog_fnas:
        with open(fna, "r") as f:
            for line in f:
                if line.startswith(">"):
                    header = line[1:].strip().split()[0]
                    background_headers.append(header)
                    db_out.write(line)
                else:
                    db_out.write(line)

# Write lists
with open(targets_list_file, "w") as f:
    for h in target_headers:
        f.write(h + "\n")

with open(background_list_file, "w") as f:
    for h in background_headers:
        f.write(h + "\n")

print(f"Total target sequences (all Babesia): {len(target_headers)} written to {targets_list_file}")
print(f"Total background sequences (non-Babesia + Dog): {len(background_headers)} written to {background_list_file}")
print("All tasks completed successfully!")
