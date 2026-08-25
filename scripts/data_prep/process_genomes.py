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

apicomplexa_extract = os.path.join(data_dir, "apicomplexa_raw")
dog_extract = os.path.join(data_dir, "dog_raw")

# 1. Download Apicomplexa (dehydrated)
if not os.path.exists(apicomplexa_zip):
    print("Downloading all Apicomplexa genomes (dehydrated)...")
    subprocess.run([
        datasets_bin, "download", "genome", "taxon", "Apicomplexa",
        "--exclude-atypical", "--dehydrated", "--filename", apicomplexa_zip
    ], check=True)
else:
    print("Apicomplexa zip already exists.")

# 2. Download Dog (dehydrated)
if not os.path.exists(dog_zip):
    print("Downloading Dog genome (dehydrated)...")
    subprocess.run([
        datasets_bin, "download", "genome", "accession", "GCF_011100685.1",
        "--dehydrated", "--filename", dog_zip
    ], check=True)
else:
    print("Dog zip already exists.")

# 3. Extract Apicomplexa
if not os.path.exists(apicomplexa_extract):
    print("Extracting Apicomplexa...")
    with zipfile.ZipFile(apicomplexa_zip, 'r') as zip_ref:
        zip_ref.extractall(apicomplexa_extract)
else:
    print("Apicomplexa already extracted.")

# 3b. Rehydrate Apicomplexa
print("Rehydrating Apicomplexa genomes (downloading individual FASTAs)...")
subprocess.run([
    datasets_bin, "rehydrate", "--directory", apicomplexa_extract
], check=True)

# 4. Extract Dog
if not os.path.exists(dog_extract):
    print("Extracting Dog genome...")
    with zipfile.ZipFile(dog_zip, 'r') as zip_ref:
        zip_ref.extractall(dog_extract)
else:
    print("Dog already extracted.")

# 4b. Rehydrate Dog
print("Rehydrating Dog genome...")
subprocess.run([
    datasets_bin, "rehydrate", "--directory", dog_extract
], check=True)

# 5. Process and sort genomes
target_fasta_out = os.path.join(data_dir, "target_babesia_canis.fa")
background_fasta_out = os.path.join(data_dir, "background_db.fa")
db_completo_out = os.path.join(data_dir, "db_completo.fa")

targets_list_file = os.path.join(data_dir, "targets_list.txt")
background_list_file = os.path.join(data_dir, "background_list.txt")

# B. canis genome identifier (BrossiM_PMB)
# Let's search all .fna files in apicomplexa
apicomplexa_fnas = glob.glob(os.path.join(apicomplexa_extract, "ncbi_dataset", "data", "GCA_*", "*.fna")) + \
                   glob.glob(os.path.join(apicomplexa_extract, "ncbi_dataset", "data", "GCF_*", "*.fna"))

dog_fnas = glob.glob(os.path.join(dog_extract, "ncbi_dataset", "data", "GCF_*", "*.fna"))

target_headers = []
background_headers = []

print(f"Found {len(apicomplexa_fnas)} Apicomplexa genomic FASTA files.")
print(f"Found {len(dog_fnas)} Dog genomic FASTA files.")

# Write merged files
with open(target_fasta_out, "w") as target_out, \
     open(background_fasta_out, "w") as bg_out:
    
    # Process Apicomplexa
    for fna in apicomplexa_fnas:
        # Check if this is Babesia canis (GCA_045269395.1)
        is_target = "GCA_045269395.1" in fna
        
        with open(fna, "r") as f:
            for line in f:
                if line.startswith(">"):
                    header = line[1:].strip().split()[0]
                    if is_target:
                        target_headers.append(header)
                        target_out.write(line)
                    else:
                        background_headers.append(header)
                        bg_out.write(line)
                else:
                    if is_target:
                        target_out.write(line)
                    else:
                        bg_out.write(line)
                        
    # Process Dog
    for fna in dog_fnas:
        with open(fna, "r") as f:
            for line in f:
                if line.startswith(">"):
                    header = line[1:].strip().split()[0]
                    background_headers.append(header)
                    bg_out.write(line)
                else:
                    bg_out.write(line)

# Create db_completo.fa (targets + background)
with open(db_completo_out, "w") as db_out:
    for fpath in [target_fasta_out, background_fasta_out]:
        with open(fpath, "r") as f:
            for line in f:
                db_out.write(line)

# Write lists
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
