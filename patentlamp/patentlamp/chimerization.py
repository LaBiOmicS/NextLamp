"""
Chimerization and Synthetic Linker Engineering Engine.
Modifies natural primers by introducing non-natural synthetic linkers, LNA/phosphorothioate backbone tags,
and non-nucleosidic chemical spacers (HEG, Spacer C3) to surpass LPI Art. 10 (natural sequence patentability exclusion).
Uses SantaLucia (1998) Nearest-Neighbor thermodynamics for full-length matrix alignment dG calculation.
Supports case-insensitive resolution of linkers, synonym mapping, and flexible dictionary key lookups.
"""

import math
from typing import Dict, List, Any, Tuple

SYNTHETIC_LINKERS = {
    "TAAA": {
        "type": "flexible_loop",
        "description": "4-nt synthetic flexible loop linker (5'-TAAA-3')",
        "is_non_nucleosidic": False,
        "st26_sequence": "TAAA",
        "human_readable": "TAAA"
    },
    "TTTT": {
        "type": "poly_t_spacer",
        "description": "Poly-T synthetic spacer tag (5'-TTTT-3')",
        "is_non_nucleosidic": False,
        "st26_sequence": "TTTT",
        "human_readable": "TTTT"
    },
    "TTTTAAAATTTT": {
        "type": "rigid_loop",
        "description": "12-nt rigid loop-forming synthetic linker",
        "is_non_nucleosidic": False,
        "st26_sequence": "TTTTAAAATTTT",
        "human_readable": "TTTTAAAATTTT"
    },
    "HEG": {
        "type": "non_nucleosidic",
        "description": "Hexaethylene glycol non-nucleosidic synthetic spacer (Spacer 18 / HEG)",
        "is_non_nucleosidic": True,
        "st26_sequence": "n",
        "human_readable": "-(HEG)-"
    },
    "SPACER_C3": {
        "type": "non_nucleosidic",
        "description": "Propane-1,3-diol non-nucleosidic synthetic spacer (Spacer C3)",
        "is_non_nucleosidic": True,
        "st26_sequence": "n",
        "human_readable": "-(C3)-"
    },
    "LNA_STABILIZED": {
        "type": "chemical_modification",
        "description": "Bases with Locked Nucleic Acid (LNA) modifications at 3' terminal",
        "is_non_nucleosidic": False,
        "st26_sequence": "TAAA",
        "human_readable": "TAAA[LNA]"
    }
}

LINKER_SYNONYMS = {
    "HEG": "HEG",
    "SPACER18": "HEG",
    "SPACER_18": "HEG",
    "C3": "SPACER_C3",
    "SPACER_C3": "SPACER_C3",
    "SPACERC3": "SPACER_C3",
    "TAAA": "TAAA",
    "TTTT": "TTTT",
    "TTTTAAAATTTT": "TTTTAAAATTTT",
    "LNA": "LNA_STABILIZED",
    "LNA_STABILIZED": "LNA_STABILIZED"
}

# SantaLucia (1998) Nearest-Neighbor parameters: dH (kcal/mol), dS (cal/K*mol)
NN_PARAMS = {
    "AA": (-7.6, -21.3), "TT": (-7.6, -21.3),
    "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9)
}

IUPAC_MAP = {
    'R': 'A', 'Y': 'C', 'S': 'G', 'W': 'A', 'K': 'G', 'M': 'A',
    'B': 'C', 'D': 'A', 'H': 'A', 'V': 'G', 'N': 'A'
}

def resolve_linker_info(linker_arg: str) -> Dict[str, Any]:
    """Resolves linker argument case-insensitively with synonym mapping and custom DNA support."""
    if not linker_arg:
        linker_arg = "TAAA"

    key_upper = linker_arg.strip().upper()
    canonical_key = LINKER_SYNONYMS.get(key_upper, key_upper)

    if canonical_key in SYNTHETIC_LINKERS:
        return SYNTHETIC_LINKERS[canonical_key]

    # Custom synthetic DNA linker
    return {
        "type": "custom_synthetic",
        "description": f"Custom synthetic oligonucleotide linker ({key_upper})",
        "is_non_nucleosidic": False,
        "st26_sequence": key_upper,
        "human_readable": key_upper
    }

def clean_iupac_sequence(sequence: str) -> str:
    """Replaces IUPAC degenerate base codes with unambiguous primary bases for thermodynamic modeling."""
    seq = sequence.upper()
    return "".join([IUPAC_MAP.get(b, b) for b in seq if b != 'N' and b != 'n'])

def calculate_gc(sequence: str) -> float:
    """Calculates GC percentage of a DNA sequence ignoring non-nucleosidic codes (n)."""
    if not sequence:
        return 0.0
    seq = clean_iupac_sequence(sequence)
    valid_bases = [b for b in seq if b in 'ATGC']
    if not valid_bases:
        return 0.0
    g_c = valid_bases.count('G') + valid_bases.count('C')
    return round((g_c / len(valid_bases)) * 100.0, 2)

def calculate_nearest_neighbor_dg(sequence: str, temp_c: float = 65.0) -> float:
    """
    Calculates free energy (dG in kcal/mol) using SantaLucia (1998) Nearest-Neighbor model at reaction temperature.
    """
    seq = clean_iupac_sequence(sequence)
    if len(seq) < 2:
        return 0.0

    temp_k = temp_c + 273.15
    dH_sum = 0.0
    dS_sum = 0.0

    if seq[0] in ['G', 'C']:
        dH_sum += 0.1
        dS_sum += -2.8
    else:
        dH_sum += 2.3
        dS_sum += 4.1

    if seq[-1] in ['G', 'C']:
        dH_sum += 0.1
        dS_sum += -2.8
    else:
        dH_sum += 2.3
        dS_sum += 4.1

    rev_comp = "".join([{'A':'T','T':'A','G':'C','C':'G'}.get(b,'N') for b in reversed(seq)])
    if seq == rev_comp:
        dS_sum += -1.4

    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        if dinuc in NN_PARAMS:
            dh, ds = NN_PARAMS[dinuc]
            dH_sum += dh
            dS_sum += ds

    dG = dH_sum - (temp_k * (dS_sum / 1000.0))
    return round(dG, 2)

def calculate_full_matrix_homodimer_dg(sequence: str, temp_c: float = 65.0) -> float:
    """
    Calculates global self-dimerization free energy (dG in kcal/mol) using full-length
    offset alignment scanning across all sliding overlap windows (SantaLucia NN).
    """
    seq = clean_iupac_sequence(sequence)
    if len(seq) < 4:
        return 0.0

    rev_comp = "".join([{'A':'T','T':'A','G':'C','C':'G'}.get(b,'N') for b in reversed(seq)])
    min_overlap = 4
    most_negative_dG = 0.0

    for offset in range(min_overlap, len(seq) + 1):
        sub_seq = seq[-offset:]
        sub_rc = rev_comp[:offset]
        
        matches = [sub_seq[k] == sub_rc[k] for k in range(len(sub_seq))]
        current_block = []
        for idx, match in enumerate(matches):
            if match:
                current_block.append(sub_seq[idx])
            else:
                if len(current_block) >= min_overlap:
                    block_str = "".join(current_block)
                    dg_val = calculate_nearest_neighbor_dg(block_str, temp_c=temp_c)
                    if dg_val < most_negative_dG:
                        most_negative_dG = dg_val
                current_block = []
        if len(current_block) >= min_overlap:
            block_str = "".join(current_block)
            dg_val = calculate_nearest_neighbor_dg(block_str, temp_c=temp_c)
            if dg_val < most_negative_dG:
                most_negative_dG = dg_val

    return round(most_negative_dG, 2)

def estimate_tm(sequence: str) -> float:
    """Estimates melting temperature (°C) using nearest-neighbor / GC approximation."""
    if not sequence:
        return 0.0
    seq = clean_iupac_sequence(sequence)
    if len(seq) < 14:
        return round((seq.count('A') + seq.count('T')) * 2 + (seq.count('G') + seq.count('C')) * 4, 2)
    return round(64.9 + 41.0 * (seq.count('G') + seq.count('C') - 16.4) / len(seq), 2)

def get_case_insensitive_val(d: Dict[str, Any], key: str, default: str = "") -> str:
    """Case-insensitive dictionary value lookup helper."""
    if not d:
        return default
    if key in d:
        return str(d[key])
    for k, v in d.items():
        if k.lower() == key.lower():
            return str(v)
    return default

def parse_and_split_fip_bip(primer_set: Dict[str, str]) -> Tuple[str, str, str, str]:
    """
    Intelligently extracts or infers F1c, F2, B1c, B2 components from a primer set dictionary.
    Supports case-insensitive key lookups and full fused FIP/BIP strings.
    """
    f1c = get_case_insensitive_val(primer_set, "F1c")
    f2 = get_case_insensitive_val(primer_set, "F2")
    b1c = get_case_insensitive_val(primer_set, "B1c")
    b2 = get_case_insensitive_val(primer_set, "B2")

    if not f1c or not f2:
        fip_full = get_case_insensitive_val(primer_set, "FIP")
        if fip_full:
            half = len(fip_full) // 2
            f1c = fip_full[:half]
            f2 = fip_full[half:]

    if not b1c or not b2:
        bip_full = get_case_insensitive_val(primer_set, "BIP")
        if bip_full:
            half = len(bip_full) // 2
            b1c = bip_full[:half]
            b2 = bip_full[half:]

    return f1c, f2, b1c, b2

class ChimerizationEngine:
    def __init__(self, default_linker: str = "TAAA"):
        self.default_linker = default_linker

    def engineer_chimeric_fip(self, f1c: str, f2: str, linker: str = None) -> Dict[str, Any]:
        """
        Engineers a synthetic non-natural FIP primer (F1c + Synthetic Linker + F2).
        Supports case-insensitive linker resolution, synonyms, and custom sequence tags.
        """
        linker_arg = linker if linker is not None else self.default_linker
        linker_info = resolve_linker_info(linker_arg)

        linker_st26 = linker_info["st26_sequence"]
        linker_human = linker_info["human_readable"]
        
        natural_fip = f1c + f2
        chimeric_fip_st26 = f1c + linker_st26 + f2
        chimeric_fip_human = f1c + linker_human + f2

        dg_nat = calculate_full_matrix_homodimer_dg(natural_fip)
        dg_chim = calculate_full_matrix_homodimer_dg(chimeric_fip_st26)

        tm_nat = estimate_tm(natural_fip)
        tm_chim = estimate_tm(chimeric_fip_st26)

        linker_start = len(f1c) + 1
        linker_end = len(f1c) + len(linker_st26)

        if dg_nat < 0.0:
            if dg_chim < 0.0:
                dimer_reduction_pct = max(0.0, min(100.0, ((abs(dg_nat) - abs(dg_chim)) / abs(dg_nat)) * 100.0))
            else:
                dimer_reduction_pct = 100.0
        else:
            dimer_reduction_pct = 0.0

        return {
            "natural_sequence": natural_fip,
            "synthetic_chimeric_sequence": chimeric_fip_st26,
            "synthetic_chimeric_human_readable": chimeric_fip_human,
            "f1c_length": len(f1c),
            "f2_length": len(f2),
            "linker_used": linker_arg,
            "linker_location": {"start": linker_start, "end": linker_end},
            "linker_info": linker_info,
            "natural_metrics": {"length": len(natural_fip), "tm": tm_nat, "gc": calculate_gc(natural_fip), "dG_dimer": dg_nat},
            "chimeric_metrics": {"length": len(chimeric_fip_st26), "tm": tm_chim, "gc": calculate_gc(chimeric_fip_st26), "dG_dimer": dg_chim},
            "dimer_reduction_percentage": round(dimer_reduction_pct, 2),
            "patentability_justification": f"Molecula quimerica sintetica nao-natural contendo ligante artificial especifico {linker_arg} (LPI Art. 10 superado)."
        }

    def process_primer_set(self, primer_set: Dict[str, str], linker: str = "TAAA") -> Dict[str, Any]:
        """
        Processes a full 6-primer LAMP set, converting natural FIP/BIP into patentable synthetic chimeras.
        Supports case-insensitive key lookups for F3, B3, LoopF, LoopB.
        """
        f1c, f2, b1c, b2 = parse_and_split_fip_bip(primer_set)

        fip_eng = self.engineer_chimeric_fip(f1c, f2, linker)
        bip_eng = self.engineer_chimeric_fip(b1c, b2, linker)

        chimeric_set = {
            "F3": get_case_insensitive_val(primer_set, "F3"),
            "B3": get_case_insensitive_val(primer_set, "B3"),
            "FIP_synthetic": fip_eng["synthetic_chimeric_sequence"],
            "BIP_synthetic": bip_eng["synthetic_chimeric_sequence"],
            "LoopF": get_case_insensitive_val(primer_set, "LoopF"),
            "LoopB": get_case_insensitive_val(primer_set, "LoopB"),
        }

        return {
            "original_set": primer_set,
            "synthetic_chimeric_set": chimeric_set,
            "fip_engineering": fip_eng,
            "bip_engineering": bip_eng
        }
