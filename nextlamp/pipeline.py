import os
import json
import shutil
import sys
from .candidates import generate_candidates
from .alignment import filter_by_specificity
from .combination import assemble_sets

class NextLampPipeline:
    def __init__(self, target_fasta: str, index_prefix: str | list[str], targets_list_file: str, background_list_file: str, bowtie2_path: str = None):
        self.target_fasta = os.path.abspath(target_fasta)
        if isinstance(index_prefix, list):
            self.index_prefix = [os.path.abspath(p) for p in index_prefix]
        else:
            self.index_prefix = os.path.abspath(index_prefix)
        
        # Dynamic auto-detection of bowtie2 binary
        if not bowtie2_path:
            bowtie2_path = shutil.which("bowtie2")
            if not bowtie2_path:
                candidates = [
                    os.path.join(os.path.dirname(sys.executable), "bowtie2"),
                    os.path.expanduser("~/miniforge3/envs/humann3_env/bin/bowtie2"),
                    os.path.expanduser("~/miniconda3/envs/humann3_env/bin/bowtie2"),
                    "/usr/bin/bowtie2",
                    "/usr/local/bin/bowtie2"
                ]
                for cand_path in candidates:
                    if os.path.isfile(cand_path):
                        bowtie2_path = cand_path
                        break
            if not bowtie2_path:
                raise FileNotFoundError("Bowtie 2 binary ('bowtie2') not found. Please specify --bowtie2-path.")
        self.bowtie2_path = os.path.abspath(bowtie2_path)
        
        # Load sequence names
        self.targets = self._load_list(os.path.abspath(targets_list_file))
        self.backgrounds = self._load_list(os.path.abspath(background_list_file))
        
    def _load_list(self, filepath: str) -> list[str]:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File {filepath} not found.")
        results = []
        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    line = line[1:]
                parts = line.split()
                if parts:
                    results.append(parts[0])
        return results

    def run(self, 
            max_sets: int = 10, 
            min_gc: float = 30.0, 
            max_gc: float = 70.0,
            min_tm: float = 55.0,
            max_tm: float = 68.0,
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
            check_dimers: bool = True,
            threads: int = 4,
            use_gpu: bool = False) -> list[dict]:
        if use_gpu:
            try:
                import torch
                if torch.cuda.is_available():
                    gpu_name = torch.cuda.get_device_name(0)
                    print(f"[GPU] CUDA Acceleration Enabled on GPU: {gpu_name}")
                else:
                    print("[GPU WARN] --gpu specified, but CUDA is not available. Falling back to multi-core CPU acceleration.")
            except ImportError:
                print("[GPU WARN] --gpu specified, but PyTorch/CUDA libraries not found. Using CPU acceleration.")

        print("Step 1: Generating candidate primers from target sequence...")
        candidates = generate_candidates(self.target_fasta, min_gc=min_gc, max_gc=max_gc)
        print(f"Generated raw candidates: F3_B3={len(candidates['F3_B3'])}, F2_B2={len(candidates['F2_B2'])}, F1c_B1c={len(candidates['F1c_B1c'])}")
        
        print(f"\nStep 2: Filtering candidates for specificity and commonality using Bowtie 2 ({threads} threads)...")
        filtered = filter_by_specificity(
            candidates_dict=candidates,
            bowtie2_path=self.bowtie2_path,
            index_prefix=self.index_prefix,
            targets_list=self.targets,
            background_list=self.backgrounds,
            threads=threads
        )
        print(f"Filtered candidates: F3_B3={len(filtered['F3_B3'])}, F2_B2={len(filtered['F2_B2'])}, F1c_B1c={len(filtered['F1c_B1c'])}")
        
        print("\nStep 3: Assembling specific LAMP primer sets...")
        lamp_sets = assemble_sets(
            filtered,
            max_sets=max_sets,
            min_dist_f3_f2=min_dist_f3_f2,
            max_dist_f3_f2=max_dist_f3_f2,
            min_dist_f2_f1c=min_dist_f2_f1c,
            max_dist_f2_f1c=max_dist_f2_f1c,
            min_dist_inner=min_dist_inner,
            max_dist_inner=max_dist_inner,
            min_dist_b1c_b2=min_dist_b1c_b2,
            max_dist_b1c_b2=max_dist_b1c_b2,
            min_dist_b2_b3=min_dist_b2_b3,
            max_dist_b2_b3=max_dist_b2_b3,
            max_dist_f1c_b1c=max_dist_f1c_b1c,
            min_tm_diff_f1c_b1c=min_tm_diff_f1c_b1c,
            check_dimers=check_dimers
        )
        params = {
            "target_fasta": self.target_fasta,
            "bowtie2_path": self.bowtie2_path,
            "index_prefix": self.index_prefix,
            "targets": self.targets,
            "backgrounds": self.backgrounds,
            "max_sets": max_sets,
            "min_gc": min_gc,
            "max_gc": max_gc,
            "min_tm": min_tm,
            "max_tm": max_tm,
            "min_dist_f3_f2": min_dist_f3_f2,
            "max_dist_f3_f2": max_dist_f3_f2,
            "min_dist_f2_f1c": min_dist_f2_f1c,
            "max_dist_f2_f1c": max_dist_f2_f1c,
            "min_dist_inner": min_dist_inner,
            "max_dist_inner": max_dist_inner,
            "min_dist_b1c_b2": min_dist_b1c_b2,
            "max_dist_b1c_b2": max_dist_b1c_b2,
            "min_dist_b2_b3": min_dist_b2_b3,
            "max_dist_b2_b3": max_dist_b2_b3,
            "max_dist_f1c_b1c": max_dist_f1c_b1c,
            "min_tm_diff_f1c_b1c": min_tm_diff_f1c_b1c,
            "check_dimers": check_dimers,
            "threads": threads
        }

        stats = {
            "raw_F3_B3": len(candidates['F3_B3']),
            "raw_F2_B2": len(candidates['F2_B2']),
            "raw_F1c_B1c": len(candidates['F1c_B1c']),
            "filt_F3_B3": len(filtered['F3_B3']),
            "filt_F2_B2": len(filtered['F2_B2']),
            "filt_F1c_B1c": len(filtered['F1c_B1c'])
        }

        print(f"Successfully designed {len(lamp_sets)} specific LAMP primer sets!")
        
        return lamp_sets, params, stats
