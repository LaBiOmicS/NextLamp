import hashlib
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict

def parse_accession_list(filepath: str) -> list[str]:
    """
    Parses accession IDs from a text, TSV, or CSV file.
    Extracts valid NCBI Accession numbers (e.g., GCF_..., GCA_..., NC_..., CM..., etc.).
    """
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"Input file not found: {filepath}")

    accessions = []
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Extract first column if CSV/TSV
            parts = re.split(r'[,\t\s]+', line)
            acc = parts[0].strip()
            if acc.startswith(">"):
                acc = acc[1:].split()[0]
            if acc:
                accessions.append(acc)
    return list(dict.fromkeys(accessions))  # Deduplicate keeping order

def download_ncbi_sequence(accession: str, out_dir: str, retries: int = 3) -> str:
    """
    Downloads sequence in FASTA format from NCBI E-utilities (NucCore / Assembly)
    or NCBI Datasets REST API.
    """
    os.makedirs(out_dir, exist_ok=True)
    out_fasta = os.path.join(out_dir, f"{accession}.fa")

    if os.path.isfile(out_fasta) and os.path.getsize(out_fasta) > 100:
        print(f"  [CACHE] {accession} already downloaded: {out_fasta}")
        return out_fasta

    print(f"  [DOWNLOADING] {accession} from NCBI...")

    # Method 0: NCBI Datasets CLI (if installed in Conda environment via ncbi-datasets-cli)
    datasets_bin = shutil.which("datasets")
    if datasets_bin and accession.startswith(("GCF_", "GCA_")):
        zip_path = os.path.join(out_dir, f"{accession}_ds.zip")
        try:
            cmd = [datasets_bin, "download", "genome", "accession", accession, "--include", "genome", "--filename", zip_path]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode == 0 and os.path.isfile(zip_path):
                import zipfile
                with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                    for member in zip_ref.namelist():
                        if member.endswith('.fna') or member.endswith('.fa'):
                            with zip_ref.open(member) as src, open(out_fasta, 'wb') as dst:
                                shutil.copyfileobj(src, dst)
                            break
                if os.path.exists(zip_path):
                    os.remove(zip_path)
                if os.path.isfile(out_fasta) and os.path.getsize(out_fasta) > 100:
                    print(f"  [OK] Downloaded via Conda ncbi-datasets-cli {accession}")
                    return out_fasta
        except Exception:
            pass

    # Method 1: NCBI E-utilities efetch (for RefSeq/GenBank nucleotide accession IDs like NC_..., CM...)
    eutils_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    
    for attempt in range(retries):
        try:
            req = urllib.request.Request(eutils_url, headers={'User-Agent': 'NextLAMP/1.0'})
            with urllib.request.urlopen(req, timeout=30) as response, open(out_fasta, "w", encoding="utf-8") as f:
                content = response.read().decode('utf-8', errors='replace')
                if content.startswith(">"):
                    f.write(content)
                    print(f"  [OK] Downloaded {accession} ({len(content)} bytes)")
                    return out_fasta
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(2)
            else:
                print(f"  [WARN] E-utilities download failed for {accession}: {e}")

    # Method 2: Fallback NCBI Datasets API for Assembly Accessions (GCF_..., GCA_...)
    if accession.startswith(("GCF_", "GCA_")):
        datasets_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession}/download?include_annotation_type=GENOME_FASTA"
        try:
            zip_path = os.path.join(out_dir, f"{accession}.zip")
            urllib.request.urlretrieve(datasets_url, zip_path)
            import zipfile
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                for member in zip_ref.namelist():
                    if member.endswith('.fna') or member.endswith('.fa'):
                        with zip_ref.open(member) as src, open(out_fasta, 'wb') as dst:
                            shutil.copyfileobj(src, dst)
                        break
            if os.path.exists(zip_path):
                os.remove(zip_path)
            if os.path.isfile(out_fasta) and os.path.getsize(out_fasta) > 100:
                print(f"  [OK] Downloaded NCBI Assembly {accession}")
                return out_fasta
        except Exception as e:
            print(f"  [FAIL] Assembly download failed for {accession}: {e}")

    raise RuntimeError(f"Could not download genome for accession {accession} from NCBI.")

def resolve_taxonomic_query(taxa: list[str], max_genomes_per_taxon: int = 20) -> list[str]:
    """
    Queries NCBI Assembly DB by taxonomic group names or TaxIDs (e.g. 'Babesia canis', 'Apicomplexa', 'Canis lupus familiaris').
    Returns a list of assembly accessions (GCF_... / GCA_...).
    """
    import json
    accessions = []
    for taxon in taxa:
        taxon_str = taxon.strip()
        if not taxon_str:
            continue
        print(f"  [NCBI TAXONOMY SEARCH] Querying NCBI for taxon: '{taxon_str}'...")

        encoded_taxon = urllib.parse.quote(taxon_str)
        found_accs = []

        # Method 1: NCBI Datasets REST API
        datasets_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/{encoded_taxon}/dataset_report?page_size={max_genomes_per_taxon}"
        try:
            req = urllib.request.Request(datasets_url, headers={'User-Agent': 'NextLAMP/1.0'})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode('utf-8'))
                reports = data.get('reports', [])
                for r in reports:
                    acc = r.get('accession')
                    if acc:
                        found_accs.append(acc)
        except Exception:
            pass

        # Method 2: Fallback Entrez ESearch for NCBI Assembly
        if not found_accs:
            esearch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term={encoded_taxon}[Organism]+AND+%22latest%20refseq%22[filter]&retmode=json&retmax={max_genomes_per_taxon}"
            try:
                req = urllib.request.Request(esearch_url, headers={'User-Agent': 'NextLAMP/1.0'})
                with urllib.request.urlopen(req, timeout=30) as resp:
                    data = json.loads(resp.read().decode('utf-8'))
                    id_list = data.get('esearchresult', {}).get('idlist', [])
                    if id_list:
                        id_str = ",".join(id_list)
                        esummary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id={id_str}&retmode=json"
                        with urllib.request.urlopen(esummary_url, timeout=30) as s_resp:
                            s_data = json.loads(s_resp.read().decode('utf-8'))
                            result = s_data.get('result', {})
                            for uid in id_list:
                                item = result.get(uid, {})
                                acc = item.get('assemblyaccession')
                                if acc:
                                    found_accs.append(acc)
            except Exception as e:
                print(f"  [WARN] Entrez search failed for '{taxon_str}': {e}")

        if found_accs:
            print(f"  [OK] Found {len(found_accs)} NCBI Assemblies for '{taxon_str}': {found_accs[:3]}...")
            accessions.extend(found_accs)
        else:
            print(f"  [WARN] No NCBI assemblies found for taxon query '{taxon_str}'")

    return list(dict.fromkeys(accessions))

def prepare_nextlamp_dataset(output_dir: str,
                            target_list_file: str = None,
                            background_list_file: str = None,
                            common_list_file: str = None,
                            target_taxa: list[str] = None,
                            common_taxa: list[str] = None,
                            background_taxa: list[str] = None,
                            max_genomes_per_taxon: int = 20,
                            bowtie2_build_path: str = None,
                            threads: int = 4) -> dict:
    """
    Automated NCBI Data Preparation & Indexing pipeline for NextLAMP.
    Downloads genomes specified by accession files OR taxonomic group names (e.g., 'Babesia canis', 'Apicomplexa', 'Canis lupus familiaris'),
    generates formatted FASTA/list files, and builds Bowtie 2 index.
    """
    out_dir = os.path.abspath(output_dir)
    os.makedirs(out_dir, exist_ok=True)
    downloads_dir = os.path.join(out_dir, "downloads")

    print(f"\n{'='*70}")
    print(f"  NextLAMP Automated NCBI Data Preparation Pipeline")
    print(f"{'='*70}\n")
    print(f"Output Directory: {out_dir}")

    # 1. Gather Accessions from Files and Taxonomic Group Searches
    target_accs = parse_accession_list(target_list_file) if target_list_file else []
    background_accs = parse_accession_list(background_list_file) if background_list_file else []
    common_accs = parse_accession_list(common_list_file) if common_list_file else []

    if target_taxa:
        print("\n--- Resolving Target Taxonomic Groups ---")
        target_accs.extend(resolve_taxonomic_query(target_taxa, max_genomes_per_taxon))

    if common_taxa:
        print("\n--- Resolving Common Target Taxonomic Groups ---")
        common_accs.extend(resolve_taxonomic_query(common_taxa, max_genomes_per_taxon))

    if background_taxa:
        print("\n--- Resolving Background Taxonomic Groups ---")
        background_accs.extend(resolve_taxonomic_query(background_taxa, max_genomes_per_taxon))

    # Exclude target accessions from background accessions automatically
    target_set = set(target_accs) | set(common_accs)
    background_accs = [acc for acc in background_accs if acc not in target_set]

    target_accs = list(dict.fromkeys(target_accs))
    common_accs = list(dict.fromkeys(common_accs))
    background_accs = list(dict.fromkeys(background_accs))

    if not target_accs:
        raise ValueError("No target accessions or taxonomic groups provided. Please specify --targets-list or --target-taxa.")

    if not background_accs:
        raise ValueError("No background accessions or taxonomic groups provided. Please specify --background-list or --background-taxa.")

    print(f"[INFO] Parsed Accessions:")
    print(f"       - Primary Targets:  {len(target_accs)}")
    print(f"       - Common Targets:   {len(common_accs)}")
    print(f"       - Backgrounds:      {len(background_accs)}")

    # 2. Download target sequences
    target_files = []
    print("\n--- Downloading Target Genomes ---")
    for acc in target_accs:
        target_files.append(download_ncbi_sequence(acc, downloads_dir))

    common_files = []
    if common_accs:
        print("\n--- Downloading Common Target Genomes ---")
        for acc in common_accs:
            common_files.append(download_ncbi_sequence(acc, downloads_dir))

    # 3. Download background sequences
    background_files = []
    print("\n--- Downloading Background Genomes ---")
    for acc in background_accs:
        background_files.append(download_ncbi_sequence(acc, downloads_dir))

    # 4. Prepare combined FASTA files and headers list
    target_fasta = os.path.join(out_dir, "target.fa")
    db_completo_fasta = os.path.join(out_dir, "db_completo.fa")
    targets_list_txt = os.path.join(out_dir, "targets_list.txt")
    background_list_txt = os.path.join(out_dir, "background_list.txt")

    targets_headers = []
    background_headers = []

    # Write target.fa
    with open(target_fasta, "w", encoding="utf-8") as out_f:
        for fpath in target_files:
            with open(fpath, "r", encoding="utf-8") as in_f:
                for line in in_f:
                    if line.startswith(">"):
                        hdr = line[1:].strip().split()[0]
                        targets_headers.append(hdr)
                    out_f.write(line)

    # Write common targets into target headers list
    for fpath in common_files:
        with open(fpath, "r", encoding="utf-8") as in_f:
            for line in in_f:
                if line.startswith(">"):
                    hdr = line[1:].strip().split()[0]
                    targets_headers.append(hdr)

    # Write db_completo.fa (targets + common + backgrounds)
    with open(db_completo_fasta, "w", encoding="utf-8") as out_f:
        # Copy target fasta content
        with open(target_fasta, "r", encoding="utf-8") as in_f:
            out_f.write(in_f.read())
        for fpath in common_files:
            with open(fpath, "r", encoding="utf-8") as in_f:
                out_f.write(in_f.read())

        # Append background fastas
        for fpath in background_files:
            with open(fpath, "r", encoding="utf-8") as in_f:
                for line in in_f:
                    if line.startswith(">"):
                        hdr = line[1:].strip().split()[0]
                        background_headers.append(hdr)
                    out_f.write(line)

    # Write targets_list.txt and background_list.txt
    with open(targets_list_txt, "w", encoding="utf-8") as f:
        for hdr in dict.fromkeys(targets_headers):
            f.write(f"{hdr}\n")

    with open(background_list_txt, "w", encoding="utf-8") as f:
        for hdr in dict.fromkeys(background_headers):
            f.write(f"{hdr}\n")

    print("\n--- Data Formatting Complete ---")
    print(f"  [CREATED] Target FASTA:       {target_fasta}")
    print(f"  [CREATED] Complete DB FASTA:   {db_completo_fasta}")
    print(f"  [CREATED] Targets List:        {targets_list_txt} ({len(targets_headers)} sequences)")
    print(f"  [CREATED] Background List:     {background_list_txt} ({len(background_headers)} sequences)")

    # 5. Auto-detect & run bowtie2-build
    index_prefix = os.path.join(out_dir, "db_completo_idx")
    
    if not bowtie2_build_path:
        bowtie2_build_path = shutil.which("bowtie2-build")
        if not bowtie2_build_path:
            candidates = [
                os.path.join(os.path.dirname(sys.executable), "bowtie2-build"),
                os.path.expanduser("~/miniforge3/envs/humann3_env/bin/bowtie2-build"),
                "/usr/bin/bowtie2-build",
                "/usr/local/bin/bowtie2-build"
            ]
            for cand in candidates:
                if os.path.isfile(cand):
                    bowtie2_build_path = cand
                    break

    if bowtie2_build_path and os.path.isfile(bowtie2_build_path):
        print(f"\n--- Building Bowtie 2 Index ({threads} threads) ---")
        cmd = [bowtie2_build_path, "--threads", str(threads), db_completo_fasta, index_prefix]
        print(f"Executing: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        print(f"  [OK] Bowtie 2 Index created: {index_prefix}.*.bt2")
    else:
        print("\n[WARN] bowtie2-build binary not found. Skipping auto-indexing.")

    result_bundle = {
        "target_fasta": target_fasta,
        "db_completo_fasta": db_completo_fasta,
        "index_prefix": index_prefix,
        "targets_list": targets_list_txt,
        "background_list": background_list_txt,
        "output_dir": out_dir
    }

    print(f"\n{'='*70}")
    print(f"  NextLAMP Dataset Ready for Run!")
    print(f"{'='*70}\n")

    return result_bundle
