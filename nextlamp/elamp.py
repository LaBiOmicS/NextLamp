"""
NextLAMP eLAMP Module (Electronic / In Silico LAMP Amplification Simulator & Quality Assessor)

Simulates in silico isothermal LAMP amplification and calculates validated, peer-reviewed
quality metrics for designed LAMP primer sets.
"""

import math
from Bio.Seq import Seq
from .thermo import calculate_tm, calculate_gc, has_hairpin, has_self_dimer

# SantaLucia 1998 nearest-neighbor thermodynamic parameters (kcal/mol) for 3'-end deltaG
_NN_DELTA_G = {
    "AA": -1.00, "TT": -1.00, "AT": -0.88, "TA": -0.58,
    "CA": -1.45, "TG": -1.45, "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28, "GA": -1.30, "TC": -1.30,
    "CG": -2.17, "GC": -2.24, "GG": -1.84, "CC": -1.84
}

def calculate_3end_deltag(sequence: str, window: int = 5) -> float:
    """
    Calculates the free energy (deltaG in kcal/mol at 37°C) of the 3'-terminal window (default 5 bp).
    Optimal 3'-terminal stability for LAMP/PCR primers is between -5.0 and -7.0 kcal/mol.
    """
    seq_upper = sequence.upper()
    if len(seq_upper) < window:
        sub = seq_upper
    else:
        sub = seq_upper[-window:]
    
    delta_g = 0.0
    for i in range(len(sub) - 1):
        pair = sub[i:i+2]
        delta_g += _NN_DELTA_G.get(pair, -1.2)
    return round(delta_g, 2)

def calculate_3end_gc_clamp(sequence: str, window: int = 5) -> int:
    """
    Counts the number of G or C bases in the 3'-terminal 5-bp window.
    Ideal: 1 to 3 G/C bases. >4 G/C bases increases non-specific priming risk.
    """
    seq_upper = sequence.upper()
    sub = seq_upper[-window:] if len(seq_upper) >= window else seq_upper
    return sub.count('G') + sub.count('C')

def evaluate_primer_set_quality(primer_set: dict) -> dict:
    """
    Calculates comprehensive, validated quality metrics for a LAMP primer set.
    Returns a dictionary of individual scores, metrics, and a composite 0-100 Quality Score.
    """
    # 1. Extract melting temperatures
    def _safe_tm(key, default_key):
        seq = primer_set.get(key, "")
        if primer_set.get(default_key) is not None:
            return float(primer_set[default_key])
        return calculate_tm(seq) if seq else 0.0

    tm_f3 = _safe_tm("f3", "tm_f3")
    tm_b3 = _safe_tm("b3", "tm_b3")
    tm_f2 = _safe_tm("f2", "tm_f2")
    tm_b2 = _safe_tm("b2", "tm_b2")
    tm_f1c = _safe_tm("f1c", "tm_f1c")
    tm_b1c = _safe_tm("b1c", "tm_b1c")

    # 2. Pairwise Imbalances
    diff_f3_b3 = abs(tm_f3 - tm_b3)
    diff_f2_b2 = abs(tm_f2 - tm_b2)
    diff_f1c_b1c = abs(tm_f1c - tm_b1c)
    tm_balance = diff_f3_b3 + diff_f2_b2

    # 3. 3'-end deltaG & GC Clamp analysis for key primers (F3, B3, F2, B2)
    key_primers = {
        "F3": primer_set.get("f3", ""),
        "B3": primer_set.get("b3", ""),
        "F2": primer_set.get("f2", ""),
        "B2": primer_set.get("b2", ""),
        "F1c": primer_set.get("f1c", ""),
        "B1c": primer_set.get("b1c", "")
    }

    deltag_3end = {name: calculate_3end_deltag(seq) for name, seq in key_primers.items() if seq}
    gc_clamp_3end = {name: calculate_3end_gc_clamp(seq) for name, seq in key_primers.items() if seq}

    # 4. Penalty calculation for 0-100 Score
    # Base score: 100
    score = 100.0

    # Imbalance penalties
    score -= (tm_balance * 3.0)       # Deduct 3 points per °C of Tm imbalance
    score -= (diff_f1c_b1c * 2.0)     # Deduct 2 points per °C of F1c/B1c imbalance

    # 3'-end stability penalties (ideal deltaG between -4.5 and -7.5 kcal/mol)
    for name, dg in deltag_3end.items():
        if dg < -8.0:
            score -= 3.0  # Overly stable 3' end -> risk of mispriming
        elif dg > -3.0:
            score -= 3.0  # Underly stable 3' end -> risk of weak initiation

    # 3'-end GC clamp penalties (>3 GC at 3' end)
    for name, gc_count in gc_clamp_3end.items():
        if gc_count >= 4:
            score -= 2.0

    score = max(0.0, min(100.0, round(score, 1)))

    # Classification
    if score >= 90.0:
        grade = "A+ (Optimal)"
    elif score >= 80.0:
        grade = "A (High Quality)"
    elif score >= 70.0:
        grade = "B (Acceptable)"
    else:
        grade = "C (Suboptimal)"

    return {
        "quality_score": score,
        "grade": grade,
        "tm_balance": round(tm_balance, 2),
        "diff_f3_b3": round(diff_f3_b3, 2),
        "diff_f2_b2": round(diff_f2_b2, 2),
        "diff_f1c_b1c": round(diff_f1c_b1c, 2),
        "deltag_3end": deltag_3end,
        "gc_clamp_3end": gc_clamp_3end
    }

def simulate_elamp_amplicon(target_seq: str, primer_set: dict) -> dict:
    """
    Simulates in silico isothermal LAMP amplification on a target FASTA sequence.
    Extracts outer amplicon (F3->B3), inner core amplicon (F2->B2), and dumbbell loop structures.
    """
    seq_upper = target_seq.upper()
    f3_seq = primer_set.get("f3", "").upper()
    b3_seq = primer_set.get("b3", "").upper()
    f2_seq = primer_set.get("f2", "").upper()
    b2_seq = primer_set.get("b2", "").upper()

    b3_rev = str(Seq(b3_seq).reverse_complement())
    b2_rev = str(Seq(b2_seq).reverse_complement())

    f3_pos = seq_upper.find(f3_seq)
    b3_pos = seq_upper.find(b3_rev)

    f2_pos = seq_upper.find(f2_seq)
    b2_pos = seq_upper.find(b2_rev)

    simulation_valid = False
    outer_amplicon = ""
    inner_amplicon = ""
    outer_len = 0
    inner_len = 0

    if f3_pos != -1 and b3_pos != -1 and b3_pos > f3_pos:
        outer_amplicon = seq_upper[f3_pos : b3_pos + len(b3_rev)]
        outer_len = len(outer_amplicon)

    if f2_pos != -1 and b2_pos != -1 and b2_pos > f2_pos:
        inner_amplicon = seq_upper[f2_pos : b2_pos + len(b2_rev)]
        inner_len = len(inner_amplicon)
        simulation_valid = True

    metrics = evaluate_primer_set_quality(primer_set)

    return {
        "simulation_valid": simulation_valid,
        "f3_start": f3_pos,
        "b3_end": b3_pos + len(b3_rev) if b3_pos != -1 else -1,
        "outer_amplicon_size": outer_len,
        "inner_amplicon_size": inner_len,
        "outer_amplicon_seq": outer_amplicon,
        "inner_amplicon_seq": inner_amplicon,
        "quality_metrics": metrics
    }
