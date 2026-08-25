#!/usr/bin/env python
import argparse
import os
import sys

# Ensure nextlamp module is in path
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)

from nextlamp.data_prep import prepare_nextlamp_dataset

def main():
    parser = argparse.ArgumentParser(
        description="NextLAMP Data Prep: Optional NCBI Genomes Downloader & Automated Bowtie 2 Indexer"
    )
    parser.add_argument("--config", default=None, help="Path to YAML or JSON configuration file")
    parser.add_argument("--generate-config", default=None, help="Path to generate an annotated Data Prep YAML configuration template")
    parser.add_argument("--targets-list", default=None, help="Text/CSV/TSV file with target genome accession IDs")
    parser.add_argument("--background-list", default=None, help="Text/CSV/TSV file with background genome accession IDs")
    parser.add_argument("--common-list", default=None, help="Optional file with common target accession IDs")
    parser.add_argument("--target-taxa", nargs="+", default=None, help="Taxonomic group names or TaxIDs for primary targets (e.g., 'Babesia canis')")
    parser.add_argument("--common-taxa", nargs="+", default=None, help="Taxonomic group names or TaxIDs for common targets (e.g., 'Babesia')")
    parser.add_argument("--background-taxa", nargs="+", default=None, help="Taxonomic group names or TaxIDs for background (e.g., 'Apicomplexa' 'Canis lupus familiaris')")
    parser.add_argument("--max-genomes-per-taxon", type=int, default=None, help="Maximum genomes to fetch per taxonomic group (default: 20)")
    parser.add_argument("--out-dir", default=None, help="Output directory to store downloaded genomes, formatted FASTAs, and Bowtie 2 index")
    parser.add_argument("--bowtie2-build-path", default=None, help="Path to bowtie2-build binary (auto-detected if omitted)")
    parser.add_argument("--threads", type=int, default=None, help="Number of CPU threads for Bowtie 2 indexing (default: 4)")
    parser.add_argument("--run-nextlamp", action="store_true", help="Automatically trigger NextLAMP pipeline after data preparation")

    args = parser.parse_args()

    if args.generate_config:
        from nextlamp.config import generate_default_prep_config_yaml
        generate_default_prep_config_yaml(args.generate_config)
        return

    # Configuration loading logic
    cfg = {
        "output_dir": "dataset_babesia",
        "threads": 4,
        "max_genomes_per_taxon": 20,
        "run_nextlamp": False
    }

    if args.config:
        from nextlamp.config import load_config
        print(f"[CONFIG] Loading Data Prep settings from YAML file: {args.config}")
        file_cfg = load_config(args.config)
        cfg.update(file_cfg)

    # CLI Overrides
    if args.out_dir: cfg["output_dir"] = args.out_dir
    if args.targets_list: cfg["target_list_file"] = args.targets_list
    if args.background_list: cfg["background_list_file"] = args.background_list
    if args.common_list: cfg["common_list_file"] = args.common_list
    if args.target_taxa: cfg["target_taxa"] = args.target_taxa
    if args.common_taxa: cfg["common_taxa"] = args.common_taxa
    if args.background_taxa: cfg["background_taxa"] = args.background_taxa
    if args.max_genomes_per_taxon is not None: cfg["max_genomes_per_taxon"] = args.max_genomes_per_taxon
    if args.threads is not None: cfg["threads"] = args.threads
    if args.run_nextlamp: cfg["run_nextlamp"] = True

    try:
        prep_results = prepare_nextlamp_dataset(
            output_dir=cfg["output_dir"],
            target_list_file=cfg.get("target_list_file"),
            background_list_file=cfg.get("background_list_file"),
            common_list_file=cfg.get("common_list_file"),
            target_taxa=cfg.get("target_taxa"),
            common_taxa=cfg.get("common_taxa"),
            background_taxa=cfg.get("background_taxa"),
            max_genomes_per_taxon=cfg.get("max_genomes_per_taxon", 20),
            bowtie2_build_path=cfg.get("bowtie2_build_path"),
            threads=cfg.get("threads", 4)
        )

        if args.run_nextlamp:
            print("\n--- Triggering NextLAMP Pipeline Automatically ---")
            from nextlamp.pipeline import NextLampPipeline
            from nextlamp.report import export_results

            pipeline = NextLampPipeline(
                target_fasta=prep_results["target_fasta"],
                index_prefix=prep_results["index_prefix"],
                targets_list_file=prep_results["targets_list"],
                background_list_file=prep_results["background_list"]
            )

            out_json = os.path.join(prep_results["output_dir"], "nextlamp_success.json")
            results, params, stats = pipeline.run(threads=args.threads)
            export_results(results, params, stats, out_json)
            print(f"\n[OK] NextLAMP Pipeline finished! Results stored in {prep_results['output_dir']}")

    except Exception as e:
        print(f"\n[ERROR] Data preparation failed: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
