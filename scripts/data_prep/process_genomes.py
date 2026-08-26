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

# B. canis genome identifier
apicomplexa_fnas = glob.glob(os.path.join(apicomplexa_extract, "ncbi_dataset", "data", "GCA_*", "*.fna")) + \
                   glob.glob(os.path.join(apicomplexa_extract, "ncbi_dataset", "data", "GCF_*", "*.fna"))

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

# Write segmented files
with open(target_fasta_out, "w") as target_out, \
     open(apicomplexa_bg_fasta_out, "w") as bg_api_out, \
     open(dog_fasta_out, "w") as dog_out, \
     open(cat_fasta_out, "w") as cat_out, \
     open(human_fasta_out, "w") as human_out, \
     open(ticks_fasta_out, "w") as ticks_out:
    
    # Process Apicomplexa (Targets vs Background)
    for fna in filtered_apicomplexa_fnas:
        is_target = any(acc in fna for acc in babesia_accessions)
        with open(fna, "r") as f:
            for line in f:
                if line.startswith(">"):
                    header = line[1:].strip().split()[0]
                    if is_target:
                        target_headers.append(header)
                        target_out.write(line)
                    else:
                        background_headers.append(header)
                        bg_api_out.write(line)
                else:
                    if is_target:
                        target_out.write(line)
                    else:
                        bg_api_out.write(line)
                        
    # Process Dog
    for fna in dog_fnas:
        with open(fna, "r") as f:
            for line in f:
                if line.startswith(">"):
                    header = line[1:].strip().split()[0]
                    background_headers.append(header)
                    dog_out.write(line)
                else:
                    dog_out.write(line)

    # Process Cat
    for fna in cat_fnas:
        with open(fna, "r") as f:
            for line in f:
                if line.startswith(">"):
                    header = line[1:].strip().split()[0]
                    background_headers.append(header)
                    cat_out.write(line)
                else:
                    cat_out.write(line)

    # Process Human
    human_fnas = glob.glob(os.path.join(human_extract, "ncbi_dataset", "data", "GCF_*", "*.fna")) + \
                 glob.glob(os.path.join(human_extract, "ncbi_dataset", "data", "GCA_*", "*.fna"))
    print(f"Found {len(human_fnas)} Human genomic FASTA files.")
    for fna in human_fnas:
        with open(fna, "r") as f:
            for line in f:
                if line.startswith(">"):
                    header = line[1:].strip().split()[0]
                    background_headers.append(header)
                    human_out.write(line)
                else:
                    human_out.write(line)

    # Process Ticks
    for fna in ticks_fnas:
        with open(fna, "r") as f:
            for line in f:
                if line.startswith(">"):
                    header = line[1:].strip().split()[0]
                    background_headers.append(header)
                    ticks_out.write(line)
                else:
                    ticks_out.write(line)

# Write sequence lists
with open(targets_list_file, "w") as f:
    for h in target_headers:
        f.write(h + "\n")

with open(background_list_file, "w") as f:
    for h in background_headers:
        f.write(h + "\n")

print(f"Target sequences written: {len(target_headers)} to {targets_list_file}")
print(f"Background sequences written: {len(background_headers)} to {background_list_file}")

# 6. Build Bowtie index (Commented out to run on SLURM later)
# index_prefix = os.path.join(data_dir, "db_completo_idx")
# print("Building Bowtie index...")
# subprocess.run([
#     bowtie_build_bin, db_completo_out, index_prefix
# ], check=True)

print("All tasks completed successfully!")
