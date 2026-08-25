#!/usr/bin/env python
"""
Local test for GLAPD and NextLAMP using the Babesia canis subsample dataset.

Dataset (tests/data/):
  - target.fa            (1 B. canis sequence: CM098081.1)
  - background.fa        (2 background sequences: CM002034.1, CM002035.1)
  - db_completo.fa       (target + background combined, 3 sequences)
  - targets_list.txt     / background_list.txt
  - bt1_idx.*            (Bowtie 1 index for GLAPD)
  - bt2_idx.*            (Bowtie 2 index for NextLAMP)

Usage:
  python tests/test_local.py
"""

import json
import os
import subprocess
import sys
import time

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTS_DIR = os.path.join(BASE_DIR, "tests")
DATA_DIR = os.path.join(TESTS_DIR, "data")
RESULTS_DIR = os.path.join(TESTS_DIR, "results")
GLAPD_DIR = os.path.join(BASE_DIR, "GLAPD")

# Binaries
BOWTIE1_PATH = os.path.join(GLAPD_DIR, "bowtie", "bowtie")
BOWTIE2_PATH = os.path.expanduser("~/miniforge3/envs/humann3_env/bin/bowtie2")
SINGLE_BIN = os.path.join(GLAPD_DIR, "Single")
LAMP_BIN = os.path.join(GLAPD_DIR, "LAMP")
PAR_PL = os.path.join(GLAPD_DIR, "par.pl")

# Data paths
TARGET_FA = os.path.join(DATA_DIR, "target.fa")
DB_FA = os.path.join(DATA_DIR, "db_completo.fa")
TARGETS_LIST = os.path.join(DATA_DIR, "targets_list.txt")
BACKGROUND_LIST = os.path.join(DATA_DIR, "background_list.txt")
BT1_INDEX = os.path.join(DATA_DIR, "bt1_idx")
BT2_INDEX = os.path.join(DATA_DIR, "bt2_idx")

# Ensure nextlamp module is importable
sys.path.insert(0, BASE_DIR)

os.makedirs(RESULTS_DIR, exist_ok=True)


def separator(title: str):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}\n")


# ======================================================================
#  TEST 1: GLAPD
# ======================================================================
def test_glapd():
    separator("TEST 1: GLAPD (Babesia subsample)")

    if not os.path.isfile(SINGLE_BIN) or not os.path.isfile(LAMP_BIN):
        print("[SKIP] GLAPD binaries not found (Single / LAMP).")
        print("       Run 'make' in GLAPD/ first.")
        return False

    out_name = "BabesiaTest"
    os.makedirs(os.path.join(GLAPD_DIR, "Inner"), exist_ok=True)
    os.makedirs(os.path.join(GLAPD_DIR, "Outer"), exist_ok=True)
    glapd_out = os.path.join(RESULTS_DIR, "glapd_success.txt")

    # Step 1: Single primer identification
    print("[GLAPD] Step 1: Identifying single primer regions...")
    t0 = time.time()
    r = subprocess.run(
        [SINGLE_BIN, "-in", TARGET_FA, "-out", out_name, "-check", "0"],
        cwd=GLAPD_DIR,
        capture_output=True, text=True
    )
    if r.returncode != 0:
        print(f"[FAIL] Single failed:\n{r.stderr}")
        return False
    print(f"[OK]   Single completed in {time.time()-t0:.1f}s")

    # Step 2: Commonality + Specificity with par.pl
    print("[GLAPD] Step 2: Running par.pl (commonality + specificity)...")
    t0 = time.time()
    r = subprocess.run(
        [
            "perl", PAR_PL,
            "--in", out_name,
            "--ref", TARGET_FA,
            "--dir", ".",
            "--bowtie", BOWTIE1_PATH,
            "--index", BT1_INDEX,
            "--common", TARGETS_LIST,
            "--specific", BACKGROUND_LIST,
        ],
        cwd=GLAPD_DIR,
        capture_output=True, text=True
    )
    if r.returncode != 0:
        print(f"[FAIL] par.pl failed:\n{r.stderr}")
        return False
    print(f"[OK]   par.pl completed in {time.time()-t0:.1f}s")

    # Step 3: LAMP primer set design (common + specific)
    print("[GLAPD] Step 3: Designing LAMP primer sets...")
    t0 = time.time()
    r = subprocess.run(
        [
            LAMP_BIN,
            "-in", out_name,
            "-ref", TARGET_FA,
            "-out", glapd_out,
            "-common", "-specific",
            "-check", "0",
        ],
        cwd=GLAPD_DIR,
        capture_output=True, text=True
    )
    if r.returncode != 0:
        print(f"[FAIL] LAMP failed:\n{r.stderr}")
        return False
    print(f"[OK]   LAMP completed in {time.time()-t0:.1f}s")

    # Show results
    if os.path.isfile(glapd_out) and os.path.getsize(glapd_out) > 0:
        with open(glapd_out) as f:
            content = f.read()
        n_sets = content.count("LAMP primers:")
        print(f"[OK]   GLAPD designed {n_sets} primer sets")
        print(f"       Output: {glapd_out}")
    else:
        print("[WARN] GLAPD output file is empty.")
        return False

    return True


# ======================================================================
#  TEST 2: NextLAMP
# ======================================================================
def test_nextlamp():
    separator("TEST 2: NextLAMP (Babesia subsample)")

    from nextlamp.pipeline import NextLampPipeline

    nextlamp_out = os.path.join(RESULTS_DIR, "nextlamp_success.json")

    # Check bowtie2 binary
    if not os.path.isfile(BOWTIE2_PATH):
        print(f"[SKIP] bowtie2 not found at: {BOWTIE2_PATH}")
        return False

    # Check bowtie2 index
    if not os.path.isfile(BT2_INDEX + ".1.bt2"):
        print(f"[SKIP] Bowtie 2 index not found: {BT2_INDEX}")
        return False

    print(f"[INFO] Target:     {TARGET_FA}")
    print(f"[INFO] Bowtie2:    {BOWTIE2_PATH}")
    print(f"[INFO] Index:      {BT2_INDEX}")

    pipeline = NextLampPipeline(
        target_fasta=TARGET_FA,
        bowtie2_path=BOWTIE2_PATH,
        index_prefix=BT2_INDEX,
        targets_list_file=TARGETS_LIST,
        background_list_file=BACKGROUND_LIST
    )

    from nextlamp.report import export_results

    t0 = time.time()
    results, params, stats = pipeline.run(max_sets=10)
    elapsed = time.time() - t0

    # Export reproducibility bundle (JSON, TSV, TXT report)
    export_results(results, params, stats, nextlamp_out)

    print(f"\n[OK]   NextLAMP completed in {elapsed:.1f}s")
    print(f"[OK]   Designed {len(results)} primer sets")
    print(f"       Output: {nextlamp_out}")

    # Display ranked results
    if results:
        print(f"\n{'─'*60}")
        print(f"  Ranked Results (best first, lower tm_balance = better)")
        print(f"{'─'*60}")
        for pset in results:
            rank = pset['rank']
            quality = pset['quality']
            tm_bal = pset['tm_balance']
            print(f"\n  #{rank} | tm_balance: {tm_bal:.4f} | quality: {quality}")
            print(f"    F3:  {pset['F3']['seq']}  (Tm: {pset['F3']['tm']:.1f}°C)")
            print(f"    F2:  {pset['F2']['seq']}  (Tm: {pset['F2']['tm']:.1f}°C)")
            print(f"    F1c: {pset['F1c']['seq']}  (Tm: {pset['F1c']['tm']:.1f}°C)")
            print(f"    B1c: {pset['B1c']['seq']}  (Tm: {pset['B1c']['tm']:.1f}°C)")
            print(f"    B2:  {pset['B2']['seq']}  (Tm: {pset['B2']['tm']:.1f}°C)")
            print(f"    B3:  {pset['B3']['seq']}  (Tm: {pset['B3']['tm']:.1f}°C)")
    else:
        print("[WARN] NextLAMP returned 0 primer sets.")
        return False

    return True


# ======================================================================
#  MAIN
# ======================================================================
if __name__ == "__main__":
    separator("NextLAMP + GLAPD — Local Test (Babesia subsample)")
    print(f"Base directory:  {BASE_DIR}")
    print(f"Test data:       {DATA_DIR}")
    print(f"Results output:  {RESULTS_DIR}")

    glapd_ok = test_glapd()
    nextlamp_ok = test_nextlamp()

    separator("SUMMARY")
    print(f"  GLAPD:    {'PASS ✓' if glapd_ok else 'FAIL / SKIPPED ✗'}")
    print(f"  NextLAMP: {'PASS ✓' if nextlamp_ok else 'FAIL / SKIPPED ✗'}")
    print()

    if not (glapd_ok and nextlamp_ok):
        sys.exit(1)
