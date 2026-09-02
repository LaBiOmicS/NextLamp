import functools
from bisect import bisect_left, bisect_right
from .thermo import has_self_dimer

_COMPLEMENT_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

@functools.lru_cache(maxsize=262144)
def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT_TABLE)[::-1]

@functools.lru_cache(maxsize=262144)
def check_heterodimer(seq1: str, seq2: str, max_contiguous_matches: int = 8) -> bool:
    """
    Check if two different primers can dimerize (cross-dimerization).
    Aligns seq1 against the reverse complement of seq2 using fast substring matching.
    """
    s1 = seq1.upper()
    s2_rc = _reverse_complement(seq2.upper())
    len1 = len(s1)
    
    # Fast check: any contiguous match of length >= max_contiguous_matches?
    for i in range(len1 - max_contiguous_matches + 1):
        if s1[i:i + max_contiguous_matches] in s2_rc:
            return True
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
    using spatial distance pruning (binary search) and heterodimer checks.
    """
    f3_b3_list = filtered_candidates.get("F3_B3", [])
    f2_b2_list = filtered_candidates.get("F2_B2", [])
    f1c_b1c_list = filtered_candidates.get("F1c_B1c", [])
    loop_list = filtered_candidates.get("Loop", [])
    
    # Group and sort candidates by spatial coordinates for binary search
    f3_list = sorted([c for c in f3_b3_list if c.strand == 1], key=lambda x: x.end)
    f3_ends = [c.end for c in f3_list]

    f2_list = sorted([c for c in f2_b2_list if c.strand == 1], key=lambda x: x.start)
    
    f1c_list = sorted([c for c in f1c_b1c_list if c.strand == -1], key=lambda x: x.start)
    f1c_starts = [c.start for c in f1c_list]
    
    b1c_list = sorted([c for c in f1c_b1c_list if c.strand == 1], key=lambda x: x.end)
    b1c_ends = [c.end for c in b1c_list]

    b2_list = sorted([c for c in f2_b2_list if c.strand == -1], key=lambda x: x.end)
    b2_ends = [c.end for c in b2_list]

    b3_list = sorted([c for c in f3_b3_list if c.strand == -1], key=lambda x: x.start)
    b3_starts = [c.start for c in b3_list]

    # Group and sort loop candidates for fast spatial bisect
    loop_f_list = sorted([c for c in loop_list if c.strand == -1], key=lambda x: x.start)
    loop_f_starts = [c.start for c in loop_f_list]

    loop_b_list = sorted([c for c in loop_list if c.strand == 1], key=lambda x: x.start)
    loop_b_starts = [c.start for c in loop_b_list]

    # Map for locus deduplication: (f2_start, b2_start) -> best lamp_set
    locus_results = {}

    for f2 in f2_list:
        # 1. Binary search for B2 (min_dist_inner <= b2.end - f2.start - 2 <= max_dist_inner)
        min_b2_end = f2.start + 2 + min_dist_inner
        max_b2_end = f2.start + 2 + max_dist_inner
        i_b2_min = bisect_left(b2_ends, min_b2_end)
        i_b2_max = bisect_right(b2_ends, max_b2_end)
        if i_b2_min >= i_b2_max:
            continue
        valid_b2s = b2_list[i_b2_min:i_b2_max]

        # 2. Binary search for F3 (min_dist_f3_f2 <= f2.start - f3.end <= max_dist_f3_f2)
        min_f3_end = f2.start - max_dist_f3_f2
        max_f3_end = f2.start - min_dist_f3_f2
        i_f3_min = bisect_left(f3_ends, min_f3_end)
        i_f3_max = bisect_right(f3_ends, max_f3_end)
        if i_f3_min >= i_f3_max:
            continue
        valid_f3s = f3_list[i_f3_min:i_f3_max]

        # 3. Binary search for F1c (min_dist_f2_f1c <= f1c.start - f2.start - 1 <= max_dist_f2_f1c)
        min_f1c_start = f2.start + 1 + min_dist_f2_f1c
        max_f1c_start = f2.start + 1 + max_dist_f2_f1c
        i_f1c_min = bisect_left(f1c_starts, min_f1c_start)
        i_f1c_max = bisect_right(f1c_starts, max_f1c_start)
        if i_f1c_min >= i_f1c_max:
            continue
        valid_f1cs = f1c_list[i_f1c_min:i_f1c_max]

        for b2 in valid_b2s:
            # 4. Binary search for B1c (min_dist_b1c_b2 <= b2.end - b1c.end - 1 <= max_dist_b1c_b2)
            min_b1c_end = b2.end - 1 - max_dist_b1c_b2
            max_b1c_end = b2.end - 1 - min_dist_b1c_b2
            i_b1c_min = bisect_left(b1c_ends, min_b1c_end)
            i_b1c_max = bisect_right(b1c_ends, max_b1c_end)
            if i_b1c_min >= i_b1c_max:
                continue
            valid_b1cs = b1c_list[i_b1c_min:i_b1c_max]

            # 5. Binary search for B3 (min_dist_b2_b3 <= b3.start - b2.end <= max_dist_b2_b3)
            min_b3_start = b2.end + min_dist_b2_b3
            max_b3_start = b2.end + max_dist_b2_b3
            i_b3_min = bisect_left(b3_starts, min_b3_start)
            i_b3_max = bisect_right(b3_starts, max_b3_start)
            if i_b3_min >= i_b3_max:
                continue
            valid_b3s = b3_list[i_b3_min:i_b3_max]

            # LoopB candidates for b1c-b2 region
            loop_b_cand = None
            lb_min = bisect_left(loop_b_starts, f2.start)
            lb_max = bisect_right(loop_b_starts, b2.start)
            for lcand in loop_b_list[lb_min:lb_max]:
                if lcand.end <= b2.start:
                    if loop_b_cand is None or abs(lcand.tm - b2.tm) < abs(loop_b_cand.tm - b2.tm):
                        loop_b_cand = lcand

            for f3 in valid_f3s:
                if check_dimers and (check_heterodimer(f3.seq, f2.seq) or check_heterodimer(f3.seq, b2.seq)):
                    continue

                for f1c in valid_f1cs:
                    if f1c.tm - f3.tm < min_tm_diff_f1c_b1c or f1c.tm - f2.tm < min_tm_diff_f1c_b1c or f1c.tm - b2.tm < min_tm_diff_f1c_b1c:
                        continue

                    if check_dimers and (check_heterodimer(f1c.seq, f2.seq) or check_heterodimer(f1c.seq, f3.seq) or check_heterodimer(f1c.seq, b2.seq)):
                        continue

                    # LoopF candidates for f2-f1c region
                    loop_f_cand = None
                    lf_min = bisect_left(loop_f_starts, f2.end)
                    lf_max = bisect_right(loop_f_starts, f1c.start)
                    for lcand in loop_f_list[lf_min:lf_max]:
                        if lcand.end <= f1c.start:
                            if loop_f_cand is None or abs(lcand.tm - f2.tm) < abs(loop_f_cand.tm - f2.tm):
                                loop_f_cand = lcand

                    for b1c in valid_b1cs:
                        if b1c.tm - f3.tm < min_tm_diff_f1c_b1c or b1c.tm - f2.tm < min_tm_diff_f1c_b1c or b1c.tm - b2.tm < min_tm_diff_f1c_b1c:
                            continue

                        if b1c.start - f1c.start > max_dist_f1c_b1c or b1c.start < f1c.end:
                            continue

                        if check_dimers and (check_heterodimer(b1c.seq, f2.seq) or check_heterodimer(b1c.seq, f3.seq) or 
                                            check_heterodimer(b1c.seq, f1c.seq) or check_heterodimer(b1c.seq, b2.seq)):
                            continue

                        for b3 in valid_b3s:
                            if f1c.tm - b3.tm < min_tm_diff_f1c_b1c or b1c.tm - b3.tm < min_tm_diff_f1c_b1c:
                                continue

                            if check_dimers and (check_heterodimer(b3.seq, f2.seq) or check_heterodimer(b3.seq, f3.seq) or 
                                                check_heterodimer(b3.seq, f1c.seq) or check_heterodimer(b3.seq, b2.seq) or 
                                                check_heterodimer(b3.seq, b1c.seq)):
                                continue

                            # Compute shared target coverage across all core primers in set
                            common_targets = None
                            for cand_obj in [f3, f2, f1c, b1c, b2, b3]:
                                m_targets = getattr(cand_obj, 'matched_targets', [])
                                if m_targets:
                                    if common_targets is None:
                                        common_targets = set(m_targets)
                                    else:
                                        common_targets &= set(m_targets)
                            
                            common_target_list = sorted(list(common_targets)) if common_targets else []

                            tm_balance = abs(f2.tm - b2.tm) + abs(f3.tm - b3.tm)
                            lamp_set = {
                                "F3": f3.to_dict(),
                                "F2": f2.to_dict(),
                                "F1c": f1c.to_dict(),
                                "B1c": b1c.to_dict(),
                                "B2": b2.to_dict(),
                                "B3": b3.to_dict(),
                                "tm_balance": round(tm_balance, 4),
                                "quality": _classify_quality(tm_balance),
                                "target_coverage_count": len(common_target_list),
                                "target_coverage_list": common_target_list
                            }
                            if loop_f_cand:
                                lamp_set["LoopF"] = loop_f_cand.to_dict()
                            if loop_b_cand:
                                lamp_set["LoopB"] = loop_b_cand.to_dict()

                            locus_key = (f2.start, b2.start)
                            if locus_key not in locus_results or tm_balance < locus_results[locus_key]["tm_balance"]:
                                locus_results[locus_key] = lamp_set
                                if len(locus_results) >= max_sets:
                                    sorted_results = sorted(locus_results.values(), key=lambda x: x["tm_balance"])
                                    return _add_ranks(sorted_results[:max_sets])

    sorted_results = sorted(locus_results.values(), key=lambda x: x["tm_balance"])
    return _add_ranks(sorted_results[:max_sets])

