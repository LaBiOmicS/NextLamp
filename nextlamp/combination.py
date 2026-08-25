from .thermo import has_self_dimer

def check_heterodimer(seq1: str, seq2: str, max_contiguous_matches: int = 8) -> bool:
    """
    Check if two different primers can dimerize (cross-dimerization).
    Aligns seq1 against the reverse complement of seq2.
    """
    from Bio.Seq import Seq
    seq1_upper = seq1.upper()
    seq2_rev_comp = str(Seq(seq2.upper()).reverse_complement())
    
    len1 = len(seq1_upper)
    len2 = len(seq2_rev_comp)
    
    # Check for contiguous matches between seq1 and reverse complement of seq2
    for shift in range(-len1 + 1, len2):
        matches = 0
        for i in range(len1):
            j = i + shift
            if 0 <= j < len2:
                if seq1_upper[i] == seq2_rev_comp[j]:
                    matches += 1
                    if matches >= max_contiguous_matches:
                        return True
                else:
                    matches = 0
    return False

def check_set_dimers(primers: list[str]) -> bool:
    """
    Checks if any pair of primers in the set has dimerizing tendencies.
    """
    n = len(primers)
    for i in range(n):
        for j in range(i + 1, n):
            if check_heterodimer(primers[i], primers[j]):
                return True
    return False

def _classify_quality(tm_balance: float) -> str:
    """
    Classifies a LAMP primer set quality based on thermal balance (Tm imbalance).
    
    tm_balance = |Tm(F2) - Tm(B2)| + |Tm(F3) - Tm(B3)|
    
    Lower values mean the forward and backward primer pairs have more similar
    melting temperatures, leading to more uniform amplification.
    
    Biological rationale (Eiken Chemical / PrimerExplorer guidelines):
    - F3, B3, F2, B2 should each have Tm ~59-61°C (target 60°C).
    - Primers within the same functional pair (F3/B3 or F2/B2) should
      have Tm as close as possible; differences > 5°C between any two
      primers in the set may cause biased amplification or failure.
    - Since tm_balance sums TWO pair-wise |ΔTm| values, the thresholds
      below reflect the combined tolerance for both pairs:
    
    Thresholds:
        <= 2.0  ->  Excellent  (avg <=1°C per pair; within ideal Tm window)
        <= 5.0  ->  Good       (within the 5°C max single-pair tolerance)
        <= 8.0  ->  Acceptable (moderate imbalance; may reduce efficiency)
        >  8.0  ->  Poor       (exceeds recommended Tm tolerance)
    
    References:
        - Notomi et al. (2000) Nucleic Acids Res. 28(12):e63
        - Eiken Chemical Co. PrimerExplorer V5 Manual
        - Nagamine et al. (2002) Mol Cell Probes 16(3):223-9
    """
    if tm_balance <= 2.0:
        return "Excellent"
    elif tm_balance <= 5.0:
        return "Good"
    elif tm_balance <= 8.0:
        return "Acceptable"
    else:
        return "Poor"

def _add_ranks(results: list[dict]) -> list[dict]:
    """
    Adds a 1-based rank to each primer set (1 = best, lowest tm_balance).
    Results must already be sorted by tm_balance ascending.
    """
    for i, r in enumerate(results):
        r["rank"] = i + 1
    return results

def assemble_sets(filtered_candidates: dict, 
                  max_sets: int = 10,
                  min_dist_f3_f2: int = 0,
                  max_dist_f3_f2: int = 20,
                  min_dist_f2_f1c: int = 40,
                  max_dist_f2_f1c: int = 60,
                  min_dist_inner: int = 120,
                  max_dist_inner: int = 180,
                  min_dist_b1c_b2: int = 40,
                  max_dist_b1c_b2: int = 60,
                  min_dist_b2_b3: int = 0,
                  max_dist_b2_b3: int = 20,
                  max_dist_f1c_b1c: int = 85,
                  min_tm_diff_f1c_b1c: float = 3.0,
                  check_dimers: bool = True) -> list[dict]:
    """
    Assembles F3, F2, F1c, B1c, B2, B3 candidates into valid LAMP primer sets
    using spatial distance pruning and heterodimer checks.
    """
    f3_b3_list = filtered_candidates.get("F3_B3", [])
    f2_b2_list = filtered_candidates.get("F2_B2", [])
    f1c_b1c_list = filtered_candidates.get("F1c_B1c", [])
    loop_list = filtered_candidates.get("Loop", [])
    
    # Group candidates by strand for positive strand synthesis
    f3_list = [c for c in f3_b3_list if c.strand == 1]
    f2_list = [c for c in f2_b2_list if c.strand == 1]
    f1c_list = [c for c in f1c_b1c_list if c.strand == -1]
    
    b1c_list = [c for c in f1c_b1c_list if c.strand == 1]
    b2_list = [c for c in f2_b2_list if c.strand == -1]
    b3_list = [c for c in f3_b3_list if c.strand == -1]
    
    results = []
    
    # We loop over F2
    for f2 in f2_list:
        # 1. Find F3 (F3-F2 distance)
        valid_f3s = [f3 for f3 in f3_list if min_dist_f3_f2 <= (f2.start - f3.end) <= max_dist_f3_f2]
        if not valid_f3s:
            continue
            
        # 2. Find F1c (F2-F1c distance)
        valid_f1cs = [f1c for f1c in f1c_list if min_dist_f2_f1c <= (f1c.start - f2.start - 1) <= max_dist_f2_f1c]
        if not valid_f1cs:
            continue
            
        for f3 in valid_f3s:
            for f1c in valid_f1cs:
                # Tm Checks for F1c relative to F3 and F2
                if f1c.tm - f3.tm < min_tm_diff_f1c_b1c or f1c.tm - f2.tm < min_tm_diff_f1c_b1c:
                    continue
                    
                # 3. Find B2 (F2-B2 inner amplicon size)
                valid_b2s = [b2 for b2 in b2_list if min_dist_inner <= (b2.end - f2.start - 2) <= max_dist_inner]
                
                for b2 in valid_b2s:
                    # Tm Checks for F1c relative to B2
                    if f1c.tm - b2.tm < min_tm_diff_f1c_b1c:
                        continue
                        
                    # 4. Find B1c (B1c-B2 distance)
                    valid_b1cs = [b1c for b1c in b1c_list if min_dist_b1c_b2 <= (b2.end - b1c.end - 1) <= max_dist_b1c_b2]
                    if not valid_b1cs:
                        continue
                        
                    # 5. Find B3 (B2-B3 distance)
                    valid_b3s = [b3 for b3 in b3_list if min_dist_b2_b3 <= (b3.start - b2.end) <= max_dist_b2_b3]
                    
                    for b1c in valid_b1cs:
                        # Tm Checks for B1c relative to F3, F2, B2
                        if b1c.tm - f3.tm < min_tm_diff_f1c_b1c or b1c.tm - f2.tm < min_tm_diff_f1c_b1c or b1c.tm - b2.tm < min_tm_diff_f1c_b1c:
                            continue
                            
                        # 6. Find F1c-B1c distance limit
                        if b1c.start - f1c.start > max_dist_f1c_b1c or b1c.start < f1c.end:
                            continue
                            
                        for b3 in valid_b3s:
                            # Tm Checks for F1c and B1c relative to B3
                            if f1c.tm - b3.tm < min_tm_diff_f1c_b1c or b1c.tm - b3.tm < min_tm_diff_f1c_b1c:
                                continue
                                
                            # 7. Check for heterodimer formations in this candidate set
                            if check_dimers:
                                primer_seqs = [f3.seq, f2.seq, f1c.seq, b1c.seq, b2.seq, b3.seq]
                                if check_set_dimers(primer_seqs):
                                    continue
                                
                            # 8. Identify optional Loop primers (LoopF and LoopB)
                            loop_f_cand = None
                            loop_b_cand = None
                            for lcand in loop_list:
                                # LoopF (strand -1, between F2 and F1c)
                                if lcand.strand == -1 and f2.end <= lcand.start and lcand.end <= f1c.start:
                                    if loop_f_cand is None or abs(lcand.tm - f2.tm) < abs(loop_f_cand.tm - f2.tm):
                                        loop_f_cand = lcand
                                # LoopB (strand 1, between B1c and B2)
                                elif lcand.strand == 1 and b1c.end <= lcand.start and lcand.end <= b2.start:
                                    if loop_b_cand is None or abs(lcand.tm - b2.tm) < abs(loop_b_cand.tm - b2.tm):
                                        loop_b_cand = lcand

                            # If all tests pass, we have a valid LAMP set!
                            tm_balance = abs(f2.tm - b2.tm) + abs(f3.tm - b3.tm)
                            lamp_set = {
                                "F3": f3.to_dict(),
                                "F2": f2.to_dict(),
                                "F1c": f1c.to_dict(),
                                "B1c": b1c.to_dict(),
                                "B2": b2.to_dict(),
                                "B3": b3.to_dict(),
                                "tm_balance": round(tm_balance, 4),
                                "quality": _classify_quality(tm_balance)
                            }
                            if loop_f_cand:
                                lamp_set["LoopF"] = loop_f_cand.to_dict()
                            if loop_b_cand:
                                lamp_set["LoopB"] = loop_b_cand.to_dict()

                            # Locus Deduplication: avoid returning near-identical sets for the same F2-B2 locus
                            locus_key = (f2.start, b2.start)
                            existing_idx = None
                            for idx, existing in enumerate(results):
                                existing_key = (existing["F2"]["start"], existing["B2"]["start"])
                                if existing_key == locus_key:
                                    existing_idx = idx
                                    break

                            if existing_idx is not None:
                                # Replace if new set has a better tm_balance
                                if tm_balance < results[existing_idx]["tm_balance"]:
                                    results[existing_idx] = lamp_set
                            else:
                                results.append(lamp_set)
                            
                            results.sort(key=lambda x: x["tm_balance"])
                            if len(results) == max_sets:
                                return _add_ranks(results)

    return _add_ranks(results[:max_sets])
