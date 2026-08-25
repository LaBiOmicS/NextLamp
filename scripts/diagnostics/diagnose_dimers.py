import json
import os
from nextlamp.pipeline import NextLampPipeline
from nextlamp.candidates import generate_candidates
from nextlamp.alignment import filter_by_specificity
from nextlamp.combination import check_heterodimer

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

spatial_passed = 0
min_max_matches = 100

for f2 in f2_list:
    valid_f3s = [f3 for f3 in f3_list if 0 <= (f2.start - f3.end) <= 20]
    valid_f1cs = [f1c for f1c in f1c_list if 40 <= (f1c.start - f2.start - 1) <= 60]
    
    for f3 in valid_f3s:
        for f1c in valid_f1cs:
            if f1c.tm - f3.tm < 3.0 or f1c.tm - f2.tm < 3.0:
                continue
                
            valid_b2s = [b2 for b2 in b2_list if 120 <= (b2.end - f2.start - 2) <= 180]
            
            for b2 in valid_b2s:
                if f1c.tm - b2.tm < 3.0:
                    continue
                    
                valid_b1cs = [b1c for b1c in b1c_list if 40 <= (b2.end - b1c.end - 1) <= 60]
                valid_b3s = [b3 for b3 in b3_list if 0 <= (b3.start - b2.end) <= 20]
                
                for b1c in valid_b1cs:
                    if b1c.tm - f3.tm < 3.0 or b1c.tm - f2.tm < 3.0 or b1c.tm - b2.tm < 3.0:
                        continue
                        
                    if b1c.start - f1c.start > 85 or b1c.start < f1c.end:
                        continue
                        
                    for b3 in valid_b3s:
                        if f1c.tm - b3.tm < 3.0 or b1c.tm - b3.tm < 3.0:
                            continue
                            
                        spatial_passed += 1
                        
                        # Let's find the maximum dimer matches for this set
                        primer_seqs = [f3.seq, f2.seq, f1c.seq, b1c.seq, b2.seq, b3.seq]
                        max_match_in_set = 0
                        for i in range(6):
                            for j in range(i+1, 6):
                                # Check contiguous matches
                                seq1 = primer_seqs[i]
                                from Bio.Seq import Seq
                                seq2_rev_comp = str(Seq(primer_seqs[j]).reverse_complement())
                                len1 = len(seq1)
                                len2 = len(seq2_rev_comp)
                                for shift in range(-len1 + 1, len2):
                                    matches = 0
                                    for k in range(len1):
                                        idx2 = k + shift
                                        if 0 <= idx2 < len2:
                                            if seq1[k] == seq2_rev_comp[idx2]:
                                                matches += 1
                                                if matches > max_match_in_set:
                                                    max_match_in_set = matches
                                            else:
                                                matches = 0
                        if max_match_in_set < min_max_matches:
                            min_max_matches = max_match_in_set
                            
                        if spatial_passed >= 100:
                            break

print(f"Spatial passed sets checked: {spatial_passed}")
print(f"Minimum max contiguous matches in any set: {min_max_matches}")
