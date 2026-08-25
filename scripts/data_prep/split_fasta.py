import os
from Bio import SeqIO

base_dir = "/home/fabiano.menegidio/workdir/Omics/genomics/babesia"
data_dir = os.path.join(base_dir, "data")
input_fasta = os.path.join(data_dir, "db_completo.fa")

# Target size per part: 3.0 GB
MAX_PART_SIZE = 3.0 * 1024 * 1024 * 1024

print("Splitting db_completo.fa into chunks of <= 3.0 GB...")

part_idx = 1
current_size = 0
current_records = []

def write_part(records, index):
    out_path = os.path.join(data_dir, f"db_part{index}.fa")
    print(f"Writing {len(records)} records to {out_path}...")
    with open(out_path, "w") as f:
        for rec in records:
            f.write(f">{rec.description}\n{str(rec.seq)}\n")
    return out_path

for record in SeqIO.parse(input_fasta, "fasta"):
    # Estimate size in bytes
    record_size = len(record.description) + len(record.seq) + 10  # header + sequence + newlines
    if current_size + record_size > MAX_PART_SIZE and current_records:
        write_part(current_records, part_idx)
        part_idx += 1
        current_records = [record]
        current_size = record_size
    else:
        current_records.append(record)
        current_size += record_size

if current_records:
    write_part(current_records, part_idx)

print(f"Successfully split into {part_idx} parts.")
with open(os.path.join(data_dir, "parts_count.txt"), "w") as f:
    f.write(str(part_idx) + "\n")
