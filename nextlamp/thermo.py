import functools
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

def _translate(char: str) -> int:
    if char == 'A': return 0
    if char == 'T': return 1
    if char == 'C': return 2
    return 3

_doublets = {
    0: (7.9, 22.2),
    1: (7.2, 20.4),
    2: (8.4, 22.4),
    3: (7.8, 21.0),
    4: (7.2, 21.3),
    5: (7.9, 22.2),
    6: (8.2, 22.2),
    7: (8.5, 22.7),
    8: (8.5, 22.7),
    9: (7.8, 21.0),
    10: (8.0, 19.9),
    11: (10.6, 27.2),
    12: (8.2, 22.2),
    13: (8.4, 22.4),
    14: (9.8, 24.4),
    15: (8.0, 19.9),
}

@functools.lru_cache(maxsize=65536)
def calculate_tm(sequence: str) -> float:
    """
    Calculate the melting temperature (Tm) of a sequence using GLAPD's exact nearest-neighbor formula.
    """
    seq_upper = sequence.upper()
    length = len(seq_upper)
    total_h = 0.0
    total_s = 0.0
    for i in range(length - 1):
        pos = _translate(seq_upper[i]) * 4 + _translate(seq_upper[i+1])
        h, s = _doublets[pos]
        total_h += h
        total_s += s
        
    total_h = -total_h
    total_s = -total_s
    
    # 5' terminal correction
    if seq_upper[0] in ('A', 'T'):
        total_h += 2.3
        total_s += 4.1
    else:
        total_h += 0.1
        total_s -= 2.8
        
    # 3' terminal correction
    if seq_upper[-1] in ('A', 'T'):
        total_h += 2.3
        total_s += 4.1
    else:
        total_h += 0.1
        total_s -= 2.8
        
    # GLAPD custom salt/entropy correction
    return 1000.0 * total_h / (total_s - 0.51986 * (length - 1) - 36.70381) - 273.15

@functools.lru_cache(maxsize=65536)
def calculate_gc(sequence: str) -> float:
    """
    Calculate GC percentage of a sequence.
    """
    g_count = sequence.count('G') + sequence.count('g')
    c_count = sequence.count('C') + sequence.count('c')
    return (g_count + c_count) / len(sequence) * 100

_comp_table = str.maketrans("ATCGatcg", "TAGCtagc")

@functools.lru_cache(maxsize=65536)
def has_hairpin(sequence: str, min_stem: int = 4, min_loop: int = 3, max_loop: int = 12) -> bool:
    """
    Check if the sequence can form a hairpin structure with a complementary stem.
    Fast vectorized slice comparison.
    """
    seq_len = len(sequence)
    if seq_len < 2 * min_stem + min_loop:
        return False
    seq_upper = sequence.upper()
    comp_seq = seq_upper.translate(_comp_table)
    
    # Quick scan for complementary 4-bp stems
    max_idx = seq_len - 2 * min_stem - min_loop + 1
    for i in range(max_idx):
        rev_comp_stem5 = comp_seq[i : i + min_stem][::-1]
        max_j = min(i + min_stem + max_loop, seq_len - min_stem)
        for j in range(i + min_stem + min_loop, max_j + 1):
            if seq_upper[j : j + min_stem] == rev_comp_stem5:
                return True
    return False

@functools.lru_cache(maxsize=131072)
def has_self_dimer(sequence: str, max_contiguous_matches: int = 6) -> bool:
    """
    Check if the sequence can dimerize with itself (Watson-Crick pairing).
    Rejects primers with long self-complementary regions (>=6 bp) or 3'-end self-dimers (>=4 bp).
    """
    seq_upper = sequence.upper()
    rev_comp = seq_upper.translate(_comp_table)[::-1]
    seq_len = len(sequence)

    # 1. Contiguous match check across entire primer length
    for i in range(seq_len - max_contiguous_matches + 1):
        if seq_upper[i:i + max_contiguous_matches] in rev_comp:
            return True

    # 2. 3'-end dimer check (4bp match involving 3' tail)
    tail3 = seq_upper[-4:]
    if tail3 in rev_comp or rev_comp[-4:] in seq_upper:
        return True

    return False
