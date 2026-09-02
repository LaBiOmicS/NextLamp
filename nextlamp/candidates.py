import os
from Bio import SeqIO
from .thermo import calculate_tm, calculate_gc, has_hairpin, has_self_dimer

class CandidatePrimer:
    __slots__ = ('seq', 'start', 'end', 'strand', 'tm', 'gc', 'matched_targets')
    
    def __init__(self, seq: str, start: int, end: int, strand: int, tm: float, gc: float, matched_targets: list[str] = None):
        self.seq = seq
        self.start = start
        self.end = end
        self.strand = strand  # 1 for forward, -1 for reverse complement
        self.tm = tm
        self.gc = gc
        self.matched_targets = matched_targets if matched_targets is not None else []

    def to_dict(self):
        d = {
            "seq": self.seq,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "tm": self.tm,
            "gc": self.gc
        }
        if self.matched_targets:
            d["matched_targets"] = self.matched_targets
            d["matched_targets_count"] = len(self.matched_targets)
        return d

def find_candidates_in_sequence(seq_str: str, 
                                min_len: int = 18, 
                                max_len: int = 28,
                                min_tm: float = 55.0, 
                                max_tm: float = 68.0,
                                min_gc: float = 30.0, 
                                max_gc: float = 70.0) -> list[CandidatePrimer]:
    """
    Scans a single sequence string on both forward and reverse strands to find candidate primers.
    """
    candidates = []
    seq_str_upper = seq_str.upper()
    seq_len = len(seq_str_upper)

    _comp_table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    rev_comp_seq = seq_str_upper.translate(_comp_table)[::-1]

    # Fast single pass for forward and reverse strands
    for strand, s_str in [(1, seq_str_upper), (-1, rev_comp_seq)]:
        for i in range(0, seq_len - min_len + 1):  # Step by 1bp for complete candidate coverage
            for length in range(min_len, max_len + 1): # Step by 1bp for all primer lengths
                if i + length > seq_len:
                    break
                sub_seq = s_str[i : i + length]
                if 'N' in sub_seq:
                    break
                
                gc = calculate_gc(sub_seq)
                if not (min_gc <= gc <= max_gc):
                    continue

                tm = calculate_tm(sub_seq)
                if not (min_tm <= tm <= max_tm):
                    continue

                if has_hairpin(sub_seq) or has_self_dimer(sub_seq):
                    continue

                orig_start = i if strand == 1 else seq_len - (i + length)
                orig_end = i + length if strand == 1 else seq_len - i

                candidates.append(CandidatePrimer(
                    seq=sub_seq,
                    start=orig_start,
                    end=orig_end,
                    strand=strand,
                    tm=tm,
                    gc=gc
                ))

    return candidates
            
    return candidates

def generate_candidates(fasta_path: str,
                        min_gc: float = 30.0,
                        max_gc: float = 70.0) -> dict[str, list[CandidatePrimer]]:
    """
    Reads a FASTA file and generates candidate lists categorized by target regions.
    """
    all_candidates = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_str = str(record.seq)
        all_candidates.extend(find_candidates_in_sequence(seq_str, min_gc=min_gc, max_gc=max_gc))
        
    # Categorize candidates based on Tm and length
    categorized = {
        "F3_B3": [],   # Outer: shorter (15-25 bp, Tm 54-63°C)
        "F2_B2": [],   # Inner F2/B2: medium (15-25 bp, Tm 55-65°C)
        "F1c_B1c": [], # Inner F1c/B1c: longer (18-28 bp, Tm 58-68°C)
        "Loop": []     # Loop primers: medium (15-25 bp, Tm 54-65°C)
    }
    
    for cand in all_candidates:
        length = len(cand.seq)
        # Outer primers F3/B3: 15-25 bp, Tm 54-63
        if 15 <= length <= 25 and 54.0 <= cand.tm <= 63.0:
            categorized["F3_B3"].append(cand)
            
        # F2/B2 primers: 15-25 bp, Tm 55-65
        if 15 <= length <= 25 and 55.0 <= cand.tm <= 65.0:
            categorized["F2_B2"].append(cand)
            
        # F1c/B1c primers: 18-28 bp, Tm 58-68
        if 18 <= length <= 28 and 58.0 <= cand.tm <= 68.0:
            categorized["F1c_B1c"].append(cand)
            
        # Loop primers: 15-25 bp, Tm 54-65
        if 15 <= length <= 25 and 54.0 <= cand.tm <= 65.0:
            categorized["Loop"].append(cand)
            
    return categorized
