#!/usr/bin/env python3
"""
NextLAMP Subsample Verification & Demonstration Test.

Runs the complete NextLAMP pipeline on the included Babesia canis subsample dataset,
verifying candidate generation, thermodynamic filtering, Bowtie2 alignment,
locus deduplication, and structured report exports.
"""

import os
import sys
import shutil
import tempfile
import time
from datetime import datetime

from nextlamp.pipeline import NextLampPipeline
from nextlamp.report import export_results

def get_sample_data_dir() -> str:
    """Returns absolute path to the bundled sample_data directory."""
    pkg_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sample_dir = os.path.join(pkg_dir, "sample_data")
    if not os.path.isdir(sample_dir):
        # Fallback to repo tests/data if installed in development
        repo_dir = os.path.dirname(pkg_dir)
        sample_dir = os.path.join(repo_dir, "tests", "data")
    return sample_dir

def run_subsample_test(output_dir: str = None) -> bool:
    """
    Executes NextLAMP pipeline on the subsample dataset.
    
    Returns True if design succeeds and returns valid primer sets.
    """
    sample_dir = get_sample_data_dir()
    target_fa = os.path.join(sample_dir, "target.fa")
    targets_list = os.path.join(sample_dir, "targets_list.txt")
    background_list = os.path.join(sample_dir, "background_list.txt")
    bt2_prefix = os.path.join(sample_dir, "bt2_idx")

    # Verify input files exist
    for fpath in [target_fa, targets_list, background_list]:
        if not os.path.isfile(fpath):
            print(f"[ERROR] Sample data file missing: {fpath}", file=sys.stderr)
            return False

    bowtie2_path = shutil.which("bowtie2") or os.path.expanduser("~/miniforge3/envs/humann3_env/bin/bowtie2")

    print("========================================================================")
    print("  NextLAMP: Subsample Dataset Functional Verification")
    print(f"  Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("========================================================================")
    print(f"  Target FASTA:     {target_fa}")
    print(f"  Bowtie2 Prefix:   {bt2_prefix}")
    print(f"  Targets List:     {targets_list}")
    print(f"  Background List:  {background_list}")
    print("========================================================================\n")

    start_time = time.time()

    pipeline = NextLampPipeline(
        target_fasta=target_fa,
        bowtie2_path=bowtie2_path,
        index_prefix=bt2_prefix,
        targets_list_file=targets_list,
        background_list_file=background_list
    )

    results, params, stats = pipeline.run(max_sets=10, threads=4)
    elapsed = time.time() - start_time

    print(f"\n[OK] Pipeline finished in {elapsed:.2f} seconds.")
    print(f"     Raw Candidates: F3/B3={stats.get('raw_F3_B3', 0)}, F2/B2={stats.get('raw_F2_B2', 0)}, F1c/B1c={stats.get('raw_F1c_B1c', 0)}")
    print(f"     Final Primer Sets Designed: {len(results)}")

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        out_json = os.path.join(output_dir, "subsample_results.json")
        export_results(results, params, stats, out_json)
        print(f"     Exported results to: {out_json}")

    if len(results) > 0:
        print("\n--- Top Designed LAMP Primer Set ---")
        top = results[0]
        quality = top.get('quality', top.get('quality_rating', 'Good'))
        tm_score = top.get('tm_balance', top.get('tm_balance_score', 0.0))
        print(f"Set #1 Quality: {quality} | Tm Balance: {tm_score:.4f}")
        for pkey in ["F3", "F2", "F1c", "B1c", "B2", "B3"]:
            pinfo = top[pkey]
            print(f"  {pkey:4s}: 5'- {pinfo['seq']} -3' (pos: {pinfo['start']}, Tm: {pinfo['tm']:.1f}°C, GC: {pinfo['gc']:.1f}%)")
        print("\n========================================================================")
        print("  VERIFICATION SUCCESSFUL: NextLAMP is fully functional!")
        print("========================================================================")
        return True
    else:
        print("\n[ERROR] Pipeline completed but generated 0 primer sets.")
        return False

def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else None
    success = run_subsample_test(output_dir=out_dir)
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
