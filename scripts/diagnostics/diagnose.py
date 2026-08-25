import json
import os
from nextlamp.pipeline import NextLampPipeline
from nextlamp.candidates import generate_candidates
from nextlamp.alignment import filter_by_specificity

target_fasta = "GLAPD/example/example.fa"
bowtie_path = "GLAPD/bowtie/bowtie"
index_prefix = "GLAPD/example/index"
targets_list = "GLAPD/example/target-list.txt"
background_list = "GLAPD/example/background-list.txt"

pipeline = NextLampPipeline(
    target_fasta=target_fasta,
    bowtie_path=bowtie_path,
    index_prefix=index_prefix,
    targets_list_file=targets_list,
    background_list_file=background_list
)

candidates = generate_candidates(pipeline.target_fasta)
filtered = filter_by_specificity(
    candidates_dict=candidates,
    bowtie_path=pipeline.bowtie_path,
    index_prefix=pipeline.index_prefix,
    targets_list=pipeline.targets,
    background_list=pipeline.backgrounds
)

f3_list = [c for c in filtered["F3_B3"] if c.strand == 1]
f2_list = [c for c in filtered["F2_B2"] if c.strand == 1]
f1c_list = [c for c in filtered["F1c_B1c"] if c.strand == -1]
b1c_list = [c for c in filtered["F1c_B1c"] if c.strand == 1]
b2_list = [c for c in filtered["F2_B2"] if c.strand == -1]
b3_list = [c for c in filtered["F3_B3"] if c.strand == -1]

print(f"Sizes - f3: {len(f3_list)}, f2: {len(f2_list)}, f1c: {len(f1c_list)}, b1c: {len(b1c_list)}, b2: {len(b2_list)}, b3: {len(b3_list)}")

f3_f2_ok = 0
f2_f1c_ok = 0
f1c_tm_ok = 0
f2_b2_ok = 0
b2_tm_ok = 0
b1c_b2_ok = 0
b3_b2_ok = 0
b1c_tm_ok = 0
f1c_b1c_dist_ok = 0

for f2 in f2_list:
    valid_f3s = [f3 for f3 in f3_list if 0 <= (f2.start - f3.end) <= 20]
    if valid_f3s:
        f3_f2_ok += 1
        
    valid_f1cs = [f1c for f1c in f1c_list if 40 <= (f1c.start - f2.start - 1) <= 60]
    if valid_f1cs:
        f2_f1c_ok += 1
        
    for f3 in valid_f3s:
        for f1c in valid_f1cs:
            if f1c.tm - f3.tm >= 3.0 and f1c.tm - f2.tm >= 3.0:
                f1c_tm_ok += 1
                
                valid_b2s = [b2 for b2 in b2_list if 120 <= (b2.end - f2.start - 2) <= 180]
                if valid_b2s:
                    f2_b2_ok += 1
                    
                for b2 in valid_b2s:
                    if f1c.tm - b2.tm >= 3.0:
                        b2_tm_ok += 1
                        
                        valid_b1cs = [b1c for b1c in b1c_list if 40 <= (b2.end - b1c.end - 1) <= 60]
                        if valid_b1cs:
                            b1c_b2_ok += 1
                            
                        valid_b3s = [b3 for b3 in b3_list if 0 <= (b3.start - b2.end) <= 20]
                        if valid_b3s:
                            b3_b2_ok += 1
                            
                        for b1c in valid_b1cs:
                            if b1c.tm - f3.tm >= 3.0 and b1c.tm - f2.tm >= 3.0 and b1c.tm - b2.tm >= 3.0:
                                b1c_tm_ok += 1
                                
                                if b1c.start - f1c.start <= 85 and b1c.start >= f1c.end:
                                    f1c_b1c_dist_ok += 1

print(f"f3_f2_ok: {f3_f2_ok}")
print(f"f2_f1c_ok: {f2_f1c_ok}")
print(f"f1c_tm_ok: {f1c_tm_ok}")
print(f"f2_b2_ok: {f2_b2_ok}")
print(f"b2_tm_ok: {b2_tm_ok}")
print(f"b1c_b2_ok: {b1c_b2_ok}")
print(f"b3_b2_ok: {b3_b2_ok}")
print(f"b1c_tm_ok: {b1c_tm_ok}")
print(f"f1c_b1c_dist_ok: {f1c_b1c_dist_ok}")
