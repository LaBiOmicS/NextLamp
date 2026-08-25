#!/usr/bin/env python
import argparse
import json
import sys
from nextlamp.pipeline import NextLampPipeline
from nextlamp.report import export_results

def main():
    parser = argparse.ArgumentParser(description="NextLAMP: A Modern, High-Performance Whole-Genome LAMP Primer Design Tool")
    parser.add_argument("--config", default=None, help="Path to YAML or JSON configuration file")
    parser.add_argument("--generate-config", default=None, help="Path to generate a clean, annotated YAML configuration template")
    parser.add_argument("--target-fasta", default=None, help="Path to target genome FASTA file")
    parser.add_argument("--bowtie2-path", default=None, help="Path to Bowtie 2 binary (auto-detected if omitted)")
    parser.add_argument("--index-prefix", default=None, help="Prefix of the Bowtie 2 index")
    parser.add_argument("--targets-list", default=None, help="File containing sequence headers of targets (one per line)")
    parser.add_argument("--background-list", default=None, help="File containing sequence headers of backgrounds (one per line)")
    parser.add_argument("--out", default=None, help="Path to save designed primer sets (JSON format)")
    parser.add_argument("--max-sets", type=int, default=None, help="Maximum number of primer sets to output")
    parser.add_argument("--min-gc", type=float, default=None, help="Minimum GC percentage for primers (default: 30.0)")
    parser.add_argument("--max-gc", type=float, default=None, help="Maximum GC percentage for primers (default: 70.0)")
    parser.add_argument("--min-tm", type=float, default=None, help="Minimum Tm for primers in °C (default: 55.0)")
    parser.add_argument("--max-tm", type=float, default=None, help="Maximum Tm for primers in °C (default: 68.0)")
    parser.add_argument("--dist-f3-f2-min", type=int, default=None, help="Minimum distance between F3 and F2 (default: 0)")
    parser.add_argument("--dist-f3-f2-max", type=int, default=None, help="Maximum distance between F3 and F2 (default: 20)")
    parser.add_argument("--dist-f2-f1c-min", type=int, default=None, help="Minimum distance between F2 and F1c (default: 40)")
    parser.add_argument("--dist-f2-f1c-max", type=int, default=None, help="Maximum distance between F2 and F1c (default: 60)")
    parser.add_argument("--dist-inner-min", type=int, default=None, help="Minimum inner amplicon size F2-B2 (default: 120)")
    parser.add_argument("--dist-inner-max", type=int, default=None, help="Maximum inner amplicon size F2-B2 (default: 180)")
    parser.add_argument("--dist-b1c-b2-min", type=int, default=None, help="Minimum distance between B1c and B2 (default: 40)")
    parser.add_argument("--dist-b1c-b2-max", type=int, default=None, help="Maximum distance between B1c and B2 (default: 60)")
    parser.add_argument("--dist-b2-b3-min", type=int, default=None, help="Minimum distance between B2 and B3 (default: 0)")
    parser.add_argument("--dist-b2-b3-max", type=int, default=None, help="Maximum distance between B2 and B3 (default: 20)")
    parser.add_argument("--dist-f1c-b1c-max", type=int, default=None, help="Maximum distance between F1c and B1c (default: 85)")
    parser.add_argument("--min-tm-diff", type=float, default=None, help="Minimum Tm difference for F1c/B1c over outer/inner primers (default: 3.0)")
    parser.add_argument("--no-check-dimers", action="store_true", help="Disable primer set heterodimer filtering")
    parser.add_argument("--include-loops", action="store_true", default=None, help="Enable searching and displaying LoopF and LoopB primers (default: True)")
    parser.add_argument("--no-loops", action="store_true", help="Disable Loop primer design")
    parser.add_argument("--threads", type=int, default=None, help="Number of CPU threads to use for parallel processing (default: 4)")
    parser.add_argument("--gpu", action="store_true", help="Enable NVIDIA CUDA GPU acceleration for candidate scanning and alignment")

    args = parser.parse_args()

    if args.generate_config:
        from nextlamp.config import generate_default_config_yaml
        generate_default_config_yaml(args.generate_config)
        return

    # Load configuration from file or use defaults
    from nextlamp.config import load_config, DEFAULT_CONFIG
    cfg = DEFAULT_CONFIG.copy()

    if args.config:
        print(f"[CONFIG] Loading parameters from configuration file: {args.config}")
        file_cfg = load_config(args.config)
        cfg.update(file_cfg)

    # CLI arguments override configuration file values if explicitly provided
    if args.target_fasta: cfg["target_fasta"] = args.target_fasta
    if args.bowtie2_path: cfg["bowtie2_path"] = args.bowtie2_path
    if args.index_prefix: cfg["index_prefix"] = args.index_prefix
    if args.targets_list: cfg["targets_list"] = args.targets_list
    if args.background_list: cfg["background_list"] = args.background_list
    if args.out: cfg["out"] = args.out
    if args.max_sets is not None: cfg["max_sets"] = args.max_sets
    if args.min_gc is not None: cfg["min_gc"] = args.min_gc
    if args.max_gc is not None: cfg["max_gc"] = args.max_gc
    if args.min_tm is not None: cfg["min_tm"] = args.min_tm
    if args.max_tm is not None: cfg["max_tm"] = args.max_tm
    if args.dist_f3_f2_min is not None: cfg["dist_f3_f2_min"] = args.dist_f3_f2_min
    if args.dist_f3_f2_max is not None: cfg["dist_f3_f2_max"] = args.dist_f3_f2_max
    if args.dist_f2_f1c_min is not None: cfg["dist_f2_f1c_min"] = args.dist_f2_f1c_min
    if args.dist_f2_f1c_max is not None: cfg["dist_f2_f1c_max"] = args.dist_f2_f1c_max
    if args.dist_inner_min is not None: cfg["dist_inner_min"] = args.dist_inner_min
    if args.dist_inner_max is not None: cfg["dist_inner_max"] = args.dist_inner_max
    if args.dist_b1c_b2_min is not None: cfg["dist_b1c_b2_min"] = args.dist_b1c_b2_min
    if args.dist_b1c_b2_max is not None: cfg["dist_b1c_b2_max"] = args.dist_b1c_b2_max
    if args.dist_b2_b3_min is not None: cfg["dist_b2_b3_min"] = args.dist_b2_b3_min
    if args.dist_b2_b3_max is not None: cfg["dist_b2_b3_max"] = args.dist_b2_b3_max
    if args.dist_f1c_b1c_max is not None: cfg["dist_f1c_b1c_max"] = args.dist_f1c_b1c_max
    if args.min_tm_diff is not None: cfg["min_tm_diff"] = args.min_tm_diff
    if args.no_check_dimers: cfg["check_dimers"] = False
    if args.include_loops: cfg["include_loops"] = True
    if args.no_loops: cfg["include_loops"] = False
    if args.threads is not None: cfg["threads"] = args.threads
    if args.gpu: cfg["gpu"] = True

    pipeline = NextLampPipeline(
        target_fasta=cfg["target_fasta"],
        bowtie2_path=cfg.get("bowtie2_path"),
        index_prefix=cfg["index_prefix"],
        targets_list_file=cfg["targets_list"],
        background_list_file=cfg["background_list"]
    )

    try:
        results, params, stats = pipeline.run(
            max_sets=cfg["max_sets"],
            min_gc=cfg["min_gc"],
            max_gc=cfg["max_gc"],
            min_tm=cfg["min_tm"],
            max_tm=cfg["max_tm"],
            min_dist_f3_f2=cfg["dist_f3_f2_min"],
            max_dist_f3_f2=cfg["dist_f3_f2_max"],
            min_dist_f2_f1c=cfg["dist_f2_f1c_min"],
            max_dist_f2_f1c=cfg["dist_f2_f1c_max"],
            min_dist_inner=cfg["dist_inner_min"],
            max_dist_inner=cfg["dist_inner_max"],
            min_dist_b1c_b2=cfg["dist_b1c_b2_min"],
            max_dist_b1c_b2=cfg["dist_b1c_b2_max"],
            min_dist_b2_b3=cfg["dist_b2_b3_min"],
            max_dist_b2_b3=cfg["dist_b2_b3_max"],
            max_dist_f1c_b1c=cfg["dist_f1c_b1c_max"],
            min_tm_diff_f1c_b1c=cfg["min_tm_diff"],
            check_dimers=cfg["check_dimers"],
            threads=cfg["threads"],
            use_gpu=cfg.get("gpu", False)
        )
        
        export_results(
            results=results,
            params=params,
            stats=stats,
            out_json=cfg["out"]
        )
        
        out_name = cfg["out"]
        base, _ = out_name.rsplit(".", 1) if "." in out_name else (out_name, "")
        print(f"\n[OK] Results exported successfully for reproducibility and interpretation:")
        print(f"     - JSON Bundle (Reproducibility & Metadata): {cfg['out']}")
        print(f"     - TSV Table (Lab & Excel Ordering):          {base}_primers.tsv")
        print(f"     - Text Summary Report (Human Interpretation): {base}_report.txt")
        
        # Display summary of designed sets (ordered best → worst)
        print(f"\n{'='*60}")
        print(f"  LAMP Primer Sets — Ranked by Thermal Balance (best first)")
        print(f"  tm_balance = |Tm(F2)−Tm(B2)| + |Tm(F3)−Tm(B3)|")
        print(f"  Lower value = more balanced Tm = better amplification")
        print(f"{'='*60}")
        for pset in results:
            rank = pset['rank']
            quality = pset['quality']
            tm_bal = pset['tm_balance']
            print(f"\n--- #{rank} | Tm Balance: {tm_bal:.4f} | Quality: {quality} ---")
            print(f"  F3:  {pset['F3']['seq']}  (pos: {pset['F3']['start']}, Tm: {pset['F3']['tm']:.1f}°C, GC: {pset['F3']['gc']:.1f}%)")
            print(f"  F2:  {pset['F2']['seq']}  (pos: {pset['F2']['start']}, Tm: {pset['F2']['tm']:.1f}°C, GC: {pset['F2']['gc']:.1f}%)")
            print(f"  F1c: {pset['F1c']['seq']}  (pos: {pset['F1c']['start']}, Tm: {pset['F1c']['tm']:.1f}°C, GC: {pset['F1c']['gc']:.1f}%)")
            print(f"  B1c: {pset['B1c']['seq']}  (pos: {pset['B1c']['start']}, Tm: {pset['B1c']['tm']:.1f}°C, GC: {pset['B1c']['gc']:.1f}%)")
            print(f"  B2:  {pset['B2']['seq']}  (pos: {pset['B2']['start']}, Tm: {pset['B2']['tm']:.1f}°C, GC: {pset['B2']['gc']:.1f}%)")
            print(f"  B3:  {pset['B3']['seq']}  (pos: {pset['B3']['start']}, Tm: {pset['B3']['tm']:.1f}°C, GC: {pset['B3']['gc']:.1f}%)")
            
    except Exception as e:
        print(f"Error running NextLAMP pipeline: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
