import hashlib
import json
import os
import time
from datetime import datetime

def compute_file_hash(filepath: str) -> str:
    """Computes SHA-256 hash of a file for reproducibility metadata."""
    if not os.path.isfile(filepath):
        return "N/A"
    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            sha256.update(chunk)
    return sha256.hexdigest()

from .elamp import evaluate_primer_set_quality, simulate_elamp_amplicon

def export_results(results: list[dict],
                   params: dict,
                   stats: dict,
                   out_json: str,
                   out_tsv: str = None,
                   out_txt: str = None):
    """
    Exports NextLAMP design results into JSON, TSV, and Text report formats.
    Enriches each primer set with in silico eLAMP quality evaluation metrics.
    """
    timestamp = datetime.now().isoformat()
    target_hash = compute_file_hash(params.get("target_fasta", ""))

    # Enrich results with eLAMP quality metrics
    target_seq = ""
    target_fasta_path = params.get("target_fasta")
    if target_fasta_path and os.path.isfile(target_fasta_path):
        from Bio import SeqIO
        try:
            record = next(SeqIO.parse(target_fasta_path, "fasta"))
            target_seq = str(record.seq)
        except Exception:
            pass

    for pset in results:
        if "elamp_metrics" not in pset:
            if target_seq:
                flat_set = {}
                for pname in ["F3", "F2", "F1c", "B1c", "B2", "B3"]:
                    if pname in pset:
                        flat_set[pname.lower()] = pset[pname]["seq"]
                        flat_set[f"tm_{pname.lower()}"] = pset[pname]["tm"]
                sim = simulate_elamp_amplicon(target_seq, flat_set)
                pset["elamp_metrics"] = sim["quality_metrics"]
                pset["elamp_amplicon"] = {
                    "outer_size": sim["outer_amplicon_size"],
                    "inner_size": sim["inner_amplicon_size"]
                }
            else:
                flat_set = {}
                for pname in ["F3", "F2", "F1c", "B1c", "B2", "B3"]:
                    if pname in pset:
                        flat_set[pname.lower()] = pset[pname]["seq"]
                        flat_set[f"tm_{pname.lower()}"] = pset[pname]["tm"]
                pset["elamp_metrics"] = evaluate_primer_set_quality(flat_set)

    metadata = {
        "tool": "NextLAMP",
        "version": "1.0.0",
        "timestamp": timestamp,
        "reproducibility": {
            "target_fasta": params.get("target_fasta"),
            "target_fasta_sha256": target_hash,
            "bowtie2_path": params.get("bowtie2_path"),
            "index_prefix": params.get("index_prefix"),
            "targets_count": len(params.get("targets", [])),
            "backgrounds_count": len(params.get("backgrounds", []))
        },
        "parameters": params,
        "statistics": stats,
        "designed_sets_count": len(results)
    }

    # 1. Export JSON (Complete Reproducibility Bundle)
    json_bundle = {
        "metadata": metadata,
        "primer_sets": results
    }
    with open(out_json, "w") as f:
        json.dump(json_bundle, f, indent=4)

    # 2. Export TSV (Lab ordering table)
    if not out_tsv:
        base, _ = os.path.splitext(out_json)
        out_tsv = f"{base}_primers.tsv"

    total_targets_count = len(params.get("targets", []))

    # 2. Export TSV (Lab ordering table with Target Specificity)
    if not out_tsv:
        base, _ = os.path.splitext(out_json)
        out_tsv = f"{base}_primers.tsv"

    with open(out_tsv, "w") as f:
        f.write("Rank\tQuality_Grade\tQuality_Score_100\tTm_Balance\tSet_Matched_Targets_Count\tSet_Matched_Targets\tPrimer_Name\tSequence_5to3\tStart_Pos\tEnd_Pos\tLength_bp\tTm_C\tGC_percent\tStrand\tPrimer_Matched_Targets_Count\tPrimer_Matched_Targets\n")
        for pset in results:
            rank = pset.get("rank", 0)
            metrics = pset.get("elamp_metrics", {})
            grade = metrics.get("grade", pset.get("quality", "N/A"))
            qscore = metrics.get("quality_score", 0.0)
            tm_bal = pset.get("tm_balance", 0.0)
            set_tgt_count = pset.get("target_coverage_count", total_targets_count)
            set_tgt_list = ",".join(pset.get("target_coverage_list", []))

            for pname in ["F3", "F2", "F1c", "LoopF", "B1c", "LoopB", "B2", "B3"]:
                if pname in pset:
                    p = pset[pname]
                    strand_str = "+" if p.get("strand", 1) == 1 else "-"
                    p_tgt_list = p.get("matched_targets", [])
                    p_tgt_count = p.get("matched_targets_count", len(p_tgt_list))
                    p_tgt_str = ",".join(p_tgt_list)
                    f.write(f"{rank}\t{grade}\t{qscore}\t{tm_bal}\t{set_tgt_count}\t{set_tgt_list}\t{pname}\t{p['seq']}\t{p['start']}\t{p['end']}\t{len(p['seq'])}\t{p['tm']}\t{p['gc']}\t{strand_str}\t{p_tgt_count}\t{p_tgt_str}\n")

    # 3. Export Formatted Summary Text Report (Human Interpretation)
    if not out_txt:
        base, _ = os.path.splitext(out_json)
        out_txt = f"{base}_report.txt"

    with open(out_txt, "w") as f:
        f.write("========================================================================\n")
        f.write("  NextLAMP: Whole-Genome LAMP Primer Design Report (with Target Specificity)\n")
        f.write("========================================================================\n\n")
        f.write(f"Execution Date:       {timestamp}\n")
        f.write(f"Target FASTA:         {params.get('target_fasta')}\n")
        f.write(f"Target SHA256:        {target_hash}\n")
        f.write(f"Target Genomes DB:    {total_targets_count}\n")
        f.write(f"Background Genomes DB:{len(params.get('backgrounds', []))}\n")
        f.write(f"Min Target Coverage:  {params.get('min_target_coverage', 1.0)*100:.1f}%\n")
        f.write(f"Threads Used:         {params.get('threads', 4)}\n\n")

        f.write("--- Design Parameters ---\n")
        f.write(f"  GC Content Range:   {params.get('min_gc')}% - {params.get('max_gc')}%\n")
        f.write(f"  Tm Range:           {params.get('min_tm')}°C - {params.get('max_tm')}°C\n")
        f.write(f"  F3-F2 Distance:     {params.get('min_dist_f3_f2')} - {params.get('max_dist_f3_f2')} bp\n")
        f.write(f"  F2-F1c Distance:    {params.get('min_dist_f2_f1c')} - {params.get('max_dist_f2_f1c')} bp\n")
        f.write(f"  F2-B2 Amplicon:     {params.get('min_dist_inner')} - {params.get('max_dist_inner')} bp\n")
        f.write(f"  B1c-B2 Distance:    {params.get('min_dist_b1c_b2')} - {params.get('max_dist_b1c_b2')} bp\n")
        f.write(f"  B2-B3 Distance:     {params.get('min_dist_b2_b3')} - {params.get('max_dist_b2_b3')} bp\n")
        f.write(f"  Dimers Filtered:    {params.get('check_dimers', True)}\n\n")

        f.write("--- Selection Funnel Statistics ---\n")
        if stats:
            f.write(f"  Raw Candidates Generated:     F3/B3={stats.get('raw_F3_B3',0)}, F2/B2={stats.get('raw_F2_B2',0)}, F1c/B1c={stats.get('raw_F1c_B1c',0)}\n")
            f.write(f"  Specific Candidates Passed:    F3/B3={stats.get('filt_F3_B3',0)}, F2/B2={stats.get('filt_F2_B2',0)}, F1c/B1c={stats.get('filt_F1c_B1c',0)}\n")
        f.write(f"  Final Designed Sets Output:    {len(results)}\n\n")

        f.write("========================================================================\n")
        f.write("  Designed LAMP Primer Sets (Ranked by eLAMP Quality Score & Target Coverage)\n")
        f.write("========================================================================\n")
        for pset in results:
            rank = pset.get("rank", 0)
            metrics = pset.get("elamp_metrics", {})
            grade = metrics.get("grade", pset.get("quality", "N/A"))
            qscore = metrics.get("quality_score", 0.0)
            tm_bal = pset.get("tm_balance", 0.0)
            set_tgt_count = pset.get("target_coverage_count", total_targets_count)
            set_tgt_list = pset.get("target_coverage_list", [])

            f.write(f"\n[SET #{rank}] Grade: {grade} | eLAMP Score: {qscore}/100 | Tm Balance: {tm_bal:.4f}\n")
            if total_targets_count > 0:
                pct = (set_tgt_count / total_targets_count) * 100
                f.write(f"  Target Specificity Coverage: {set_tgt_count}/{total_targets_count} targets ({pct:.1f}%)\n")
            if set_tgt_list:
                f.write(f"  Specific Target Genomes Matched: {', '.join(set_tgt_list)}\n")

            for pname in ["F3", "F2", "F1c", "LoopF", "B1c", "LoopB", "B2", "B3"]:
                if pname in pset:
                    p = pset[pname]
                    p_tgt_list = p.get("matched_targets", [])
                    p_tgt_str = f" [Matched: {len(p_tgt_list)} targets]" if p_tgt_list else ""
                    f.write(f"  {pname:5s}: 5'- {p['seq']} -3'  (pos: {p['start']}, Tm: {p['tm']:.1f}°C, GC: {p['gc']:.1f}%){p_tgt_str}\n")
        f.write("\n========================================================================\n")
