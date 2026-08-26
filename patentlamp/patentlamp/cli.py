"""
PatentLAMP CLI: Command-Line Interface for Autonomous Patent Package Generation.
Supports multi-set processing, robust JSON error handling, UTF-8 output encoding, and complete patent artifact generation.
"""

import argparse
import json
import os
import sys

from patentlamp.chimerization import ChimerizationEngine
from patentlamp.inventive_step import InventiveStepProofEngine
from patentlamp.wipo_st26 import WIPOSequenceListingGenerator

def main():
    parser = argparse.ArgumentParser(description="PatentLAMP: Autonomous Patent-Ready Molecular Engineering & Intellectual Property Engine")
    parser.add_argument("--input-json", required=True, help="Path to NextLAMP results JSON file")
    parser.add_argument("--out-dir", required=True, help="Directory to save the patent package")
    parser.add_argument("--target-species", default="Babesia canis", help="Target organism name for patent text")
    parser.add_argument("--applicant", default="LaBiOmicS / Universidade", help="Applicant name for WIPO ST.26 listing")
    parser.add_argument("--linker", default="TAAA", help="Synthetic linker sequence to insert in FIP/BIP chimeras")

    args = parser.parse_args()

    if not os.path.exists(args.input_json):
        print(f"[ERROR] Input JSON file not found: {args.input_json}")
        sys.exit(1)

    os.makedirs(args.out_dir, exist_ok=True)

    print("==========================================================================")
    print("  PatentLAMP: Autonomous Patent-Ready Molecular Engineering Engine")
    print(f"  Processing NextLAMP Input: {args.input_json}")
    print(f"  Target Organism: {args.target_species}")
    print("==========================================================================")

    try:
        with open(args.input_json, "r", encoding="utf-8") as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        print(f"[ERROR] Failed to parse JSON input file '{args.input_json}': {e}")
        sys.exit(1)

    primer_sets = data.get("primer_sets", [])
    if not primer_sets:
        if isinstance(data, list):
            primer_sets = data
        else:
            print("[WARNING] No 'primer_sets' key found in JSON input. Checking root data...")

    if not primer_sets:
        print("[ERROR] No valid primer sets found to process.")
        sys.exit(1)

    print(f"--> Loaded {len(primer_sets)} candidate primer sets from NextLAMP.")

    # 1. Chimerization Engine across all primer sets
    print("\n--- Step 1: Executing Synthetic Chimerization Engine (LPI Art. 10 Compliance) ---")
    chim_engine = ChimerizationEngine(default_linker=args.linker)
    processed_sets = []
    for idx, pset in enumerate(primer_sets, start=1):
        proc = chim_engine.process_primer_set(pset, linker=args.linker)
        processed_sets.append(proc)
        fip_seq = proc['fip_engineering']['synthetic_chimeric_sequence']
        print(f"  [Set {idx}] Synthetic FIP Chimera: {fip_seq}")

    # 2. Inventive Step Engine across all sets
    print("\n--- Step 2: Evaluating Inventive Step & Technical Proofs (LPI Art. 13 Compliance) ---")
    proof_engine = InventiveStepProofEngine(target_species=args.target_species)
    all_evaluations = []
    for proc in processed_sets:
        ev = proof_engine.evaluate_comparative_advantage(proc)
        all_evaluations.append(ev)

    eval_json_path = os.path.join(args.out_dir, "patentlamp_inventive_step_proofs.json")
    with open(eval_json_path, "w", encoding="utf-8") as f:
        json.dump(all_evaluations, f, ensure_ascii=False, indent=4)
    print(f"  [OUTPUT] Multi-set comparative proof report saved to: {eval_json_path}")

    # 3. WIPO ST.26 & Patent Draft Package Generator
    print("\n--- Step 3: Generating WIPO ST.26 XML Listing & INPI Patent Specification ---")
    st26_gen = WIPOSequenceListingGenerator(applicant_name=args.applicant, title=f"Conjunto Oligonucleotidico Sintetico para Detecçao Isotermica de {args.target_species}")
    
    st26_xml_path = os.path.join(args.out_dir, "WIPO_ST26_Sequence_Listing.xml")
    st26_gen.generate_st26_xml(processed_sets, st26_xml_path)
    print(f"  [OUTPUT] WIPO ST.26 XML Sequence Listing: {st26_xml_path}")

    patent_md_path = os.path.join(args.out_dir, "INPI_Patent_Specification_Draft.md")
    st26_gen.generate_patent_draft_md(all_evaluations, patent_md_path)
    print(f"  [OUTPUT] INPI Patent Specification Draft: {patent_md_path}")

    print("\n==========================================================================")
    print("  PatentLAMP Package Completed Successfully!")
    print(f"  All patent artifacts saved in: {args.out_dir}")
    print("==========================================================================")

if __name__ == "__main__":
    main()
