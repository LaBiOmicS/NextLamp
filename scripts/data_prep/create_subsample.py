import os
from Bio import SeqIO

base_dir = "/home/fabiano.menegidio/workdir/Omics/genomics/babesia"
data_dir = os.path.join(base_dir, "data")
sub_dir = os.path.join(data_dir, "subsample")
os.makedirs(sub_dir, exist_ok=True)

# 1. Subsample target (Babesia canis) - first 50,000 bp
target_in = os.path.join(data_dir, "target_babesia_canis.fa")
target_out = os.path.join(sub_dir, "target.fa")

target_headers = []
with open(target_out, "w") as out_f:
    for record in SeqIO.parse(target_in, "fasta"):
        sub_record = record[:50000]
        SeqIO.write(sub_record, out_f, "fasta")
        target_headers.append(record.id)
        break  # Just take the first chromosome/contig

# 2. Subsample background - first 200,000 bp of the background_db.fa
bg_in = os.path.join(data_dir, "background_db.fa")
bg_out = os.path.join(sub_dir, "background.fa")

bg_headers = []
with open(bg_out, "w") as out_f:
    count = 0
    for record in SeqIO.parse(bg_in, "fasta"):
        sub_record = record[:100000]
        SeqIO.write(sub_record, out_f, "fasta")
        bg_headers.append(record.id)
        count += 1
        if count >= 2:  # Just take first two contigs
            break

# 3. Create db_completo.fa for subsample (target + background)
db_out = os.path.join(sub_dir, "db_completo.fa")
with open(db_out, "w") as out_f:
    for fpath in [target_out, bg_out]:
        with open(fpath, "r") as in_f:
            out_f.write(in_f.read())

# 4. Write targets and background lists
with open(os.path.join(sub_dir, "targets_list.txt"), "w") as f:
    for h in target_headers:
        f.write(h + "\n")

with open(os.path.join(sub_dir, "background_list.txt"), "w") as f:
    for h in bg_headers:
        f.write(h + "\n")

print("Subsample created successfully under data/subsample/")
print(f"Targets: {target_headers}")
print(f"Backgrounds: {bg_headers}")
