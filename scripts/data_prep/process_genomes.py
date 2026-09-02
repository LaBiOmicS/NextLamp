import os
import zipfile
import glob
import subprocess

# Define paths
base_dir = "/home/fabiano.menegidio/workdir/Omics/genomics/babesia"
data_dir = os.path.join(base_dir, "data")
os.makedirs(data_dir, exist_ok=True)

datasets_bin = "/home/fabiano.menegidio/miniforge3/envs/ncbi_env/bin/datasets"
bowtie_build_bin = os.path.join(base_dir, "GLAPD", "bowtie", "bowtie-build")

apicomplexa_zip = os.path.join(data_dir, "apicomplexa.zip")
dog_zip = os.path.join(data_dir, "dog.zip")
human_zip = os.path.join(data_dir, "human.zip")

apicomplexa_extract = os.path.join(data_dir, "apicomplexa_raw")
dog_extract = os.path.join(data_dir, "dog_raw")
human_extract = os.path.join(data_dir, "human_raw")
ticks_extract = os.path.join(data_dir, "ticks_raw")

# 1. Download Apicomplexa (dehydrated if zip exists)
if not os.path.exists(apicomplexa_extract) and not os.path.exists(apicomplexa_zip):
    print("Downloading all Apicomplexa genomes (dehydrated)...")
    subprocess.run([
        datasets_bin, "download", "genome", "taxon", "Apicomplexa",
        "--exclude-atypical", "--dehydrated", "--filename", apicomplexa_zip
    ], check=True)

# 2. Download Dog (dehydrated if zip exists)
if not os.path.exists(dog_extract) and not os.path.exists(dog_zip):
    print("Downloading Dog genome (dehydrated)...")
    subprocess.run([
        datasets_bin, "download", "genome", "accession", "GCF_011100685.1",
        "--dehydrated", "--filename", dog_zip
    ], check=True)

# 2b. Download Human (dehydrated)
if not os.path.exists(human_extract) and not os.path.exists(human_zip):
    print("Downloading Human genome (GRCh38.p14 dehydrated)...")
    subprocess.run([
        datasets_bin, "download", "genome", "accession", "GCF_000001405.40",
        "--dehydrated", "--filename", human_zip
    ], check=True)

# 3. Extract & Rehydrate Apicomplexa
if os.path.exists(apicomplexa_zip):
    print("Extracting Apicomplexa...")
    with zipfile.ZipFile(apicomplexa_zip, 'r') as zip_ref:
        zip_ref.extractall(apicomplexa_extract)
    os.remove(apicomplexa_zip)

if os.path.exists(apicomplexa_extract):
    print("Rehydrating Apicomplexa genomes...")
    subprocess.run([datasets_bin, "rehydrate", "--directory", apicomplexa_extract], check=True)

# 4. Extract & Rehydrate Dog
if os.path.exists(dog_zip):
    print("Extracting Dog genome...")
    with zipfile.ZipFile(dog_zip, 'r') as zip_ref:
        zip_ref.extractall(dog_extract)
    os.remove(dog_zip)

if os.path.exists(dog_extract):
    print("Rehydrating Dog genome...")
    subprocess.run([datasets_bin, "rehydrate", "--directory", dog_extract], check=True)

# 4b. Extract & Rehydrate Human
if os.path.exists(human_zip):
    print("Extracting Human genome...")
    with zipfile.ZipFile(human_zip, 'r') as zip_ref:
        zip_ref.extractall(human_extract)
    os.remove(human_zip)

if os.path.exists(human_extract):
    print("Rehydrating Human genome...")
    subprocess.run([datasets_bin, "rehydrate", "--directory", human_extract], check=True)

# 5. Process and sort genomes into segmented FASTA files
target_fasta_out = os.path.join(data_dir, "target_babesia_canis.fa")
dog_fasta_out = os.path.join(data_dir, "host_dog.fa")
cat_fasta_out = os.path.join(data_dir, "host_cat.fa")
human_fasta_out = os.path.join(data_dir, "host_human.fa")
ticks_fasta_out = os.path.join(data_dir, "vectors_ticks.fa")
apicomplexa_bg_fasta_out = os.path.join(data_dir, "bg_apicomplexa.fa")

targets_list_file = os.path.join(data_dir, "targets_list.txt")
background_list_file = os.path.join(data_dir, "background_list.txt")

import re

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

# Select primary representative tick reference assemblies (R. microplus and I. scapularis)
all_ticks_fnas = glob.glob(os.path.join(ticks_extract, "ncbi_dataset", "data", "GCA_*", "*.fna")) + \
                 glob.glob(os.path.join(ticks_extract, "ncbi_dataset", "data", "GCF_*", "*.fna"))

primary_tick_accessions = {"GCF_013339725.1", "GCF_016920785.2", "GCA_031841145.2", "GCF_000208615.1"}
ticks_fnas = [fna for fna in all_ticks_fnas if any(acc in fna for acc in primary_tick_accessions)]
if not ticks_fnas and all_ticks_fnas:
    ticks_fnas = all_ticks_fnas[:2]  # Fallback to top 2 if accessions differ

cat_extract = os.path.join(data_dir, "cat_raw")
cat_fnas = glob.glob(os.path.join(cat_extract, "ncbi_dataset", "data", "GCF_*", "*.fna")) + \
           glob.glob(os.path.join(cat_extract, "ncbi_dataset", "data", "GCA_*", "*.fna"))

target_headers = []
background_headers = []

print(f"Found {len(apicomplexa_fnas)} Apicomplexa genomic FASTA files.")
print(f"Found {len(dog_fnas)} Dog genomic FASTA files.")
print(f"Found {len(cat_fnas)} Cat genomic FASTA files.")
print(f"Found {len(ticks_fnas)} Ticks genomic FASTA files.")

# Load all Babesia accession IDs and filter 1 representative accession per non-Babesia Apicomplexa species
babesia_accessions = set()
apicomplexa_bg_selected_accs = set()
species_selected = {}

jsonl_path = os.path.join(apicomplexa_extract, "ncbi_dataset", "data", "assembly_data_report.jsonl")
if os.path.exists(jsonl_path):
    import json
    with open(jsonl_path) as f:
        for line in f:
            data = json.loads(line)
            acc = data.get("accession")
            org = data.get("organism", {}).get("organismName", "")
            level = data.get("assemblyInfo", {}).get("assemblyLevel", "")
            refseq_category = data.get("assemblyInfo", {}).get("refseqCategory", "")
            
            if "Babesia" in org and acc:
                babesia_accessions.add(acc)
            elif acc:
                # Priority: RefSeq representative/reference > RefSeq (GCF) > Assembly Level
                is_gcf = 1 if acc.startswith("GCF_") else 0
                is_refseq_cat = 2 if refseq_category in ["reference genome", "representative genome"] else 0
                level_score = 3 if level == "Complete Genome" else (2 if level == "Chromosome" else 1)
                score = (is_refseq_cat, is_gcf, level_score)
                
                if org not in species_selected or score > species_selected[org]["score"]:
                    species_selected[org] = {"acc": acc, "score": score}

apicomplexa_bg_selected_accs = {info["acc"] for info in species_selected.values()}

print(f"Identified {len(babesia_accessions)} Babesia genome accessions as targets.")
print(f"Deduplicated Apicomplexa background to {len(apicomplexa_bg_selected_accs)} representative species accessions.")

# Filter apicomplexa_fnas to include all Babesia targets + only selected representative background genomes
filtered_apicomplexa_fnas = []
for fna in apicomplexa_fnas:
    is_target = any(acc in fna for acc in babesia_accessions)
    is_selected_bg = any(acc in fna for acc in apicomplexa_bg_selected_accs)
    if is_target or is_selected_bg:
        filtered_apicomplexa_fnas.append(fna)

print(f"Filtering reduced Apicomplexa FASTA files from {len(apicomplexa_fnas)} to {len(filtered_apicomplexa_fnas)} (42 targets + {len(apicomplexa_bg_selected_accs)} background).")

def normalize_acc(acc_str):
    if acc_str.startswith("GCA_"):
        return "GCF_" + acc_str[4:]
    return acc_str

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

# Write segmented files
with open(target_fasta_out, "w") as target_out, \
     open(apicomplexa_bg_fasta_out, "w") as bg_api_out, \
     open(dog_fasta_out, "w") as dog_out, \
     open(cat_fasta_out, "w") as cat_out, \
     open(human_fasta_out, "w") as human_out, \
     open(ticks_fasta_out, "w") as ticks_out:
    
    # Process Apicomplexa (Targets vs Background)
    for fna in filtered_apicomplexa_fnas:
        matched_acc = None
        for acc in babesia_accessions:
            if acc in fna:
                matched_acc = acc
                break
        if not matched_acc:
            parts = fna.split(os.sep)
            for part in parts:
                if part.startswith("GCA_") or part.startswith("GCF_"):
                    matched_acc = part
                    break
            if not matched_acc:
                matched_acc = os.path.basename(fna)

        is_target = matched_acc in babesia_accessions
        norm_acc = normalize_acc(matched_acc)
        pseudo_seq = read_fasta_as_pseudogenome(fna, spacer_len=100)

        if is_target:
            target_headers.append((norm_acc, norm_acc))
            target_out.write(f">{norm_acc}\n{pseudo_seq}\n")
        else:
            background_headers.append((norm_acc, norm_acc))
            bg_api_out.write(f">{norm_acc}\n{pseudo_seq}\n")
                        
    # Process Dog
    for fna in dog_fnas:
        matched_acc = "Dog_GCF_011100685.1"
        pseudo_seq = read_fasta_as_pseudogenome(fna, spacer_len=100)
        background_headers.append((matched_acc, matched_acc))
        dog_out.write(f">{matched_acc}\n{pseudo_seq}\n")

    # Process Cat
    for fna in cat_fnas:
        matched_acc = "Cat_GCF_018350175.1"
        pseudo_seq = read_fasta_as_pseudogenome(fna, spacer_len=100)
        background_headers.append((matched_acc, matched_acc))
        cat_out.write(f">{matched_acc}\n{pseudo_seq}\n")

    # Process Human
    human_fnas = glob.glob(os.path.join(human_extract, "ncbi_dataset", "data", "GCF_*", "*.fna")) + \
                 glob.glob(os.path.join(human_extract, "ncbi_dataset", "data", "GCA_*", "*.fna"))
    print(f"Found {len(human_fnas)} Human genomic FASTA files.")
    for fna in human_fnas:
        matched_acc = "Human_GRCh38.p14"
        pseudo_seq = read_fasta_as_pseudogenome(fna, spacer_len=100)
        background_headers.append((matched_acc, matched_acc))
        human_out.write(f">{matched_acc}\n{pseudo_seq}\n")

    # Process Ticks
    for fna in ticks_fnas:
        matched_acc = None
        parts = fna.split(os.sep)
        for part in parts:
            if part.startswith("GCA_") or part.startswith("GCF_"):
                matched_acc = part
                break
        if not matched_acc:
            matched_acc = os.path.basename(fna)
        norm_acc = normalize_acc(matched_acc)
        pseudo_seq = read_fasta_as_pseudogenome(fna, spacer_len=100)
        background_headers.append((norm_acc, norm_acc))
        ticks_out.write(f">{norm_acc}\n{pseudo_seq}\n")

# Write sequence lists
with open(targets_list_file, "w") as f:
    for h, acc in target_headers:
        f.write(f"{h}\t{acc}\n")

with open(background_list_file, "w") as f:
    for h, acc in background_headers:
        f.write(f"{h}\t{acc}\n")

unique_target_genomes = len(set(acc for _, acc in target_headers))
print(f"Target sequences written: {len(target_headers)} pseudo-genomes spanning {unique_target_genomes} unique Babesia assemblies to {targets_list_file}")
print(f"Background sequences written: {len(background_headers)} pseudo-genomes to {background_list_file}")

# 6. Build Bowtie index (Commented out to run on SLURM later)
# index_prefix = os.path.join(data_dir, "db_completo_idx")
# print("Building Bowtie index...")
# subprocess.run([
#     bowtie_build_bin, db_completo_out, index_prefix
# ], check=True)

print("All tasks completed successfully!")
