import os
import subprocess
import tempfile
from collections import defaultdict

def filter_by_specificity(candidates_dict: dict,
                          bowtie2_path: str,
                          index_prefix: str | list[str],
                          targets_list: list[str],
                          background_list: list[str],
                          max_target_mismatches: int = 0,
                          max_background_mismatches: int = 2,
                          threads: int = 4) -> dict:
    """
    Aligns candidate primers against one or multiple target & background databases using Bowtie 2.
    Supports sequential early-exit filtering over a list of index prefixes with minimal RAM footprint.
    """
    targets_set = set(targets_list)
    background_set = set(background_list)

    index_prefixes = [index_prefix] if isinstance(index_prefix, str) else index_prefix

    # 1. Deduplicate candidate sequences to avoid redundant Bowtie 2 queries
    unique_seqs = {}  # { seq: [ (category, candidate_obj), ... ] }
    for category, cand_list in candidates_dict.items():
        for cand in cand_list:
            if cand.seq not in unique_seqs:
                unique_seqs[cand.seq] = []
            unique_seqs[cand.seq].append((category, cand))

    seq_to_id = {seq: f"Q{i}" for i, seq in enumerate(unique_seqs.keys())}
    id_to_seq = {f"Q{i}": seq for i, seq in enumerate(unique_seqs.keys())}

    # 2. Write unique queries to temporary FASTA file
    temp_fasta = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
    for seq, qid in seq_to_id.items():
        temp_fasta.write(f">{qid}\n{seq}\n")
    temp_fasta.close()

    target_matches = defaultdict(set)
    background_hits = set()

    try:
        for idx in index_prefixes:
            # Check if all candidates were already eliminated
            active_queries = [qid for qid in id_to_seq.keys() if qid not in background_hits]
            if not active_queries:
                break

            cmd = [
                bowtie2_path,
                "-p", str(threads),
                "-f",
                "-x", idx,
                "-U", temp_fasta.name,
                "-k", "100",
                "--end-to-end",
                "--no-hd",
                "--no-sq",
                "--very-sensitive"
            ]

            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                text=True,
                bufsize=1048576  # 1MB buffer
            )

            for line in proc.stdout:
                if not line or line.startswith("@"):
                    continue
                parts = line.split("\t")
                if len(parts) < 11:
                    continue

                flag = int(parts[1])
                if flag & 4:  # Unmapped
                    continue

                query_id = parts[0]
                if query_id in background_hits:
                    continue  # Early exit: skip if query hit background in previous SAM line or index

                ref_name = parts[2]

                # Parse NM (edit distance/mismatches) tag
                mismatch_count = 0
                for part in parts[11:]:
                    if part.startswith("NM:i:"):
                        try:
                            mismatch_count = int(part.split(":")[-1])
                        except ValueError:
                            pass
                        break

                if ref_name in background_set and mismatch_count <= max_background_mismatches:
                    background_hits.add(query_id)
                    target_matches.pop(query_id, None)

                elif ref_name in targets_set and mismatch_count <= max_target_mismatches:
                    target_matches[query_id].add(ref_name)

            proc.stdout.close()
            proc.wait()

    finally:
        if os.path.exists(temp_fasta.name):
            os.remove(temp_fasta.name)

    filtered = {
        "F3_B3": [],
        "F2_B2": [],
        "F1c_B1c": [],
        "Loop": []
    }

    req_targets_count = len(targets_set)

    for qid, seq in id_to_seq.items():
        if qid in background_hits:
            continue
        if len(target_matches[qid]) >= req_targets_count:
            for category, cand in unique_seqs[seq]:
                filtered[category].append(cand)

    return filtered
