import os
import subprocess
import tempfile
from collections import defaultdict

import re

def normalize_assembly_id(acc: str) -> str:
    """Normalizes assembly accessions by mapping GenBank (GCA_) to RefSeq (GCF_) counterparts."""
    if acc.startswith("GCA_"):
        return "GCF_" + acc[4:]
    return acc

def filter_by_specificity(candidates_dict: dict,
                          bowtie2_path: str,
                          index_prefix: str | list[str],
                          targets_list: list[str],
                          background_list: list[str],
                          max_target_mismatches: int = 1,
                          max_background_mismatches: int = 2,
                          threads: int = 4,
                          min_target_coverage: float = 1.0,
                          min_targets_count: int = None) -> dict:
    """
    Aligns candidate primers against one or multiple target & background databases using Bowtie 2.
    Supports sequential early-exit filtering over a list of index prefixes with minimal RAM footprint.
    """
    # Build targets_map: { contig_header: group_or_assembly_id }
    targets_map = {}
    if isinstance(targets_list, dict):
        targets_map = {k: normalize_assembly_id(v) for k, v in targets_list.items()}
    elif isinstance(targets_list, (list, set)):
        for item in targets_list:
            if isinstance(item, (tuple, list)):
                targets_map[item[0]] = normalize_assembly_id(item[1])
            elif isinstance(item, str):
                parts = item.strip().split()
                if len(parts) >= 2:
                    targets_map[parts[0]] = normalize_assembly_id(parts[1])
                elif parts:
                    targets_map[parts[0]] = normalize_assembly_id(parts[0])

    background_set = set()
    if isinstance(background_list, (list, set, dict)):
        for item in background_list:
            if isinstance(item, (tuple, list)):
                background_set.add(item[0])
            elif isinstance(item, str):
                background_set.add(item.strip().split()[0])

    targets_set = set(targets_map.keys())
    total_target_groups = len(set(targets_map.values()))

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

                # Parse NM (mismatch count) and MD (mismatch positions) tags
                mismatch_count = 0
                md_tag = None
                for part in parts[11:]:
                    if part.startswith("NM:i:"):
                        try:
                            mismatch_count = int(part.split(":")[-1])
                        except ValueError:
                            pass
                    elif part.startswith("MD:Z:"):
                        md_tag = part[5:]

                if ref_name in background_set and mismatch_count <= max_background_mismatches:
                    background_hits.add(query_id)
                    target_matches.pop(query_id, None)

                elif ref_name in targets_set and mismatch_count <= max_target_mismatches:
                    # Enforce strict 3' end anchor rule: 0 mismatches allowed in the last 5 bp at 3' end
                    has_3prime_mismatch = False
                    if mismatch_count > 0 and md_tag:
                        # Extract query sequence length
                        qseq = parts[9]
                        qlen = len(qseq)
                        # Check if last 5 bp contains mismatch via MD tag
                        # MD tag contains numbers (matches) and letters (ref bases for mismatches)
                        # e.g. "15A4" -> 15 matches, 1 mismatch, 4 matches
                        md_matches = re.findall(r'(\d+)|([A-Za-z]|\^[A-Za-z]+)', md_tag)
                        pos = 0
                        for num_str, mismatch_str in md_matches:
                            if num_str:
                                pos += int(num_str)
                            elif mismatch_str:
                                if not mismatch_str.startswith("^"):
                                    if pos >= qlen - 5:  # Mismatch within last 5bp of 3' end
                                        has_3prime_mismatch = True
                                        break
                                    pos += len(mismatch_str)
                                else:
                                    pos += len(mismatch_str) - 1

                    if not has_3prime_mismatch:
                        target_matches[query_id].add(targets_map[ref_name])

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

    if min_targets_count is not None:
        req_targets_count = max(1, min(min_targets_count, total_target_groups))
    elif min_target_coverage is not None and 0.0 < min_target_coverage <= 1.0:
        req_targets_count = max(1, int(total_target_groups * min_target_coverage))
    else:
        req_targets_count = total_target_groups

    for qid, seq in id_to_seq.items():
        if qid in background_hits:
            continue
        matched_list = sorted(list(target_matches[qid]))
        if len(matched_list) >= req_targets_count:
            for category, cand in unique_seqs[seq]:
                cand.matched_targets = matched_list
                filtered[category].append(cand)

    return filtered
