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

def export_results(results: list[dict],
                   params: dict,
                   stats: dict,
                   out_json: str,
                   out_tsv: str = None,
                   out_txt: str = None):
    """
    Exports NextLAMP design results into JSON, TSV, and Text report formats.
    Ensures complete scientific reproducibility and biological interpretation.
    """
    timestamp = datetime.now().isoformat()
    target_hash = compute_file_hash(params.get("target_fasta", ""))

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

    with open(out_tsv, "w") as f:
        f.write("Rank\tQuality\tTm_Balance\tPrimer_Name\tSequence_5to3\tStart_Pos\tEnd_Pos\tLength_bp\tTm_C\tGC_percent\tStrand\n")
        for pset in results:
            rank = pset.get("rank", 0)
            quality = pset.get("quality", "N/A")
            tm_bal = pset.get("tm_balance", 0.0)

            for pname in ["F3", "F2", "F1c", "B1c", "B2", "B3"]:
                if pname in pset:
                    p = pset[pname]
                    strand_str = "+" if p.get("strand", 1) == 1 else "-"
                    f.write(f"{rank}\t{quality}\t{tm_bal}\t{pname}\t{p['seq']}\t{p['start']}\t{p['end']}\t{len(p['seq'])}\t{p['tm']}\t{p['gc']}\t{strand_str}\n")

    # 3. Export Formatted Summary Text Report (Human Interpretation)
    if not out_txt:
        base, _ = os.path.splitext(out_json)
        out_txt = f"{base}_report.txt"

    with open(out_txt, "w") as f:
        f.write("========================================================================\n")
        f.write("  NextLAMP: Whole-Genome LAMP Primer Design Report\n")
        f.write("========================================================================\n\n")
        f.write(f"Execution Date:       {timestamp}\n")
        f.write(f"Target FASTA:         {params.get('target_fasta')}\n")
        f.write(f"Target SHA256:        {target_hash}\n")
        f.write(f"Target Genomes:       {len(params.get('targets', []))}\n")
        f.write(f"Background Genomes:   {len(params.get('backgrounds', []))}\n")
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
        f.write("  Designed LAMP Primer Sets (Ranked by Thermal Balance)\n")
        f.write("========================================================================\n")
        for pset in results:
            rank = pset.get("rank", 0)
            quality = pset.get("quality", "N/A")
            tm_bal = pset.get("tm_balance", 0.0)

            f.write(f"\n[SET #{rank}] Quality: {quality} | Tm Balance Score: {tm_bal:.4f}\n")
            for pname in ["F3", "F2", "F1c", "LoopF", "B1c", "LoopB", "B2", "B3"]:
                if pname in pset:
                    p = pset[pname]
                    f.write(f"  {pname:5s}: 5'- {p['seq']} -3'  (pos: {p['start']}, Tm: {p['tm']:.1f}°C, GC: {p['gc']:.1f}%)\n")
        f.write("\n========================================================================\n")
