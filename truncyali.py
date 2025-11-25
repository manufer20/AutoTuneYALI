#!/usr/bin/env python3

"""
TruncYALI-style truncation cassette designer for Yarrowia lipolytica (W29).
UPDATED:
- Fully Automated: Uses Selenium Robot to check CHOPCHOP.
- Fallback: Local Doench 2014 if robot fails.
- Backbone: Restored to original specifications.
- This script designs truncation cassettes (e.g., deleting the first N bp of an ORF) rather than full-gene knockouts.
"""
from datetime import datetime
from dataclasses import dataclass
from typing import Dict, Tuple, Optional, List
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import argparse
import requests
import zipfile
import tempfile
import shutil
import math
import csv
import matplotlib.pyplot as plt
import time

# --- Selenium for Automation ---
try:
    from selenium import webdriver
    from selenium.webdriver.chrome.service import Service
    from selenium.webdriver.chrome.options import Options
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait
    from selenium.webdriver.support import expected_conditions as EC
    from webdriver_manager.chrome import ChromeDriverManager
    HAS_SELENIUM = True
except ImportError:
    HAS_SELENIUM = False
    print("Warning: 'selenium' not installed. CHOPCHOP automation will be skipped.")

# Optional: DNA Features Viewer
try:
    from dna_features_viewer import GraphicFeature, GraphicRecord
    HAS_DNA_VIEWER = True
except ImportError:
    HAS_DNA_VIEWER = False

# Optional: Pydna for Gel Simulation
try:
    from pydna.gel import gel
    from pydna.dseqrecord import Dseqrecord
    from pydna.ladders import GeneRuler_1kb
    HAS_PYDNA = True
except ImportError:
    HAS_PYDNA = False
    print("Warning: 'pydna' not installed. Gel simulation will be skipped.")


# =========================
# Constants
# =========================

GA_LEFT = "TTCGATTCCGGGTCGGCGCA"
GA_RIGHT = "ATCGCGTGCATTCGCGGCCG"

TRACR_SCAFFOLD = (
    "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGA"
    "AAAAGTGGCACCGAGTCGGTGGTGC"
)

TER_RPR1 = "TTTTTTTGTGTGGGTACGGT"


# =========================
# Data Structures
# =========================

@dataclass
class GeneInfo:
    gene_id: str
    chrom: str
    start: int  # 1-based, inclusive
    end: int    # 1-based, inclusive
    strand: str # '+' or '-'

@dataclass
class SgRNA:
    protospacer: str
    pam: str
    chrom: str
    cut_pos: int
    strand: str
    score: float = 0.0  # Added to track quality
    method: str = "Local" # Added to track source

@dataclass
class PrimerInfo:
    name: str
    seq: str
    genomic_start: int # 0-based coordinate on genome (for math)
    length: int

@dataclass
class KnockoutDesign:
    gene: GeneInfo
    sgrna: SgRNA
    upstream_arm: str
    downstream_arm: str
    fragment_seq: str
    primers: List[PrimerInfo]
    band_size_unedited: int
    band_size_ko: int
    search_window_seq: str = ""


# =========================
# Reference Files (Restored Original Logic)
# =========================

def ensure_reference_files(genome_path: str, gff_path: str) -> None:
    missing = []
    if not os.path.exists(genome_path): missing.append(genome_path)
    if not os.path.exists(gff_path): missing.append(gff_path)

    if not missing: return

    print("Missing reference file(s). Attempting automatic download...")
    url = (
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
        "GCF_001761485.1/download"
        "?include_annotation_type=GENOME_FASTA"
        "&include_annotation_type=GENOME_GFF"
        "&hydrated=FULLY_HYDRATED"
    )
    
    try:
        with tempfile.NamedTemporaryFile(suffix=".zip", delete=False) as tmp:
            tmp_path = tmp.name
            resp = requests.get(url, stream=True)
            resp.raise_for_status()
            for chunk in resp.iter_content(chunk_size=8192):
                if chunk: tmp.write(chunk)

        with zipfile.ZipFile(tmp_path, "r") as zf:
            fasta, gff = None, None
            for name in zf.namelist():
                if name.lower().endswith((".fna", ".fasta")): fasta = name
                if name.lower().endswith(".gff"): gff = name
            
            if not fasta or not gff: raise RuntimeError("Files not found in zip.")

            with zf.open(fasta) as src, open(genome_path, "wb") as dst:
                shutil.copyfileobj(src, dst)
            with zf.open(gff) as src, open(gff_path, "wb") as dst:
                shutil.copyfileobj(src, dst)
        os.remove(tmp_path)
    except Exception as e:
        raise FileNotFoundError(f"Download failed: {e}")


def load_genome_fasta(path: str) -> Dict[str, str]:
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(path, "fasta")}


def load_annotation_gff(gff_path: str) -> Dict[str, GeneInfo]:
    genes = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] not in ("gene", "CDS"): continue
            
            # Robust ID parsing from original script
            attrs = {k: v for k, v in [x.split("=", 1) for x in cols[8].split(";") if "=" in x]}
            gene_id = None
            for key in ("ID", "locus_tag", "gene", "Name"):
                val = attrs.get(key, "")
                if val.startswith("YALI1_"):
                    gene_id = val
                    break
            
            if not gene_id: continue

            start, end = int(cols[3]), int(cols[4])
            if gene_id in genes:
                g = genes[gene_id]
                genes[gene_id] = GeneInfo(gene_id, cols[0], min(g.start, start), max(g.end, end), cols[6])
            else:
                genes[gene_id] = GeneInfo(gene_id, cols[0], start, end, cols[6])
    return genes


# =========================
# Core Logic
# =========================

def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]

def gc_content(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0

# --- NEW: Doench 2014 Scoring ---
def calculate_doench_score(seq_30bp: str) -> float:
    seq = seq_30bp.upper()
    if len(seq) != 30: return 50.0 
    score = 0.59763615
    proto = seq[4:24]; gc = proto.count('G') + proto.count('C')
    score += -0.2026259 if gc < 10 else 0
    score += -0.1665835 if gc > 13 else 0
    weights = {('A',2):-0.04558,('C',2):-0.09094,('C',5):-0.05169,('T',5):-0.0597,('A',6):-0.06,('C',6):-0.02053,('G',6):0.0573,('A',10):-0.0072,('C',11):0.0295,('A',12):0.013,('C',13):-0.03784,('G',14):-0.01266,('A',15):0.0611,('C',15):-0.04518,('A',16):-0.0935,('C',16):0.0609,('G',16):-0.08655,('A',18):0.0534,('C',18):-0.05818,('G',18):-0.02908,('T',18):-0.05597,('G',19):0.0454,('G',20):-0.0759,('T',20):0.0649,('C',21):0.041,('T',21):-0.04217,('A',22):0.005,('C',22):-0.03666,('G',22):0.0601,('G',23):-0.0102,('C',24):0.0368,('T',24):-0.05458,('C',27):0.0854,('G',27):-0.06074,('A',28):0.0954,('C',28):0.0334}
    for i, n in enumerate(seq):
        if (n, i) in weights: score += weights[(n, i)]
    return int((1/(1+math.exp(-score*4)))*100)

# --- NEW: Selenium Robot ---
def get_chopchop_guide(sequence: str, gene_label: Optional[str] = None) -> Tuple[Optional[str], float]:
    """Query CHOPCHOP on the web using Selenium.

    - Sends the provided DNA sequence in FASTA mode.
    - Selects the Y. lipolytica W29 genome.
    - Waits for the results table.
    - Saves the full table as CSV (if gene_label is provided).
    - Saves a PNG screenshot of the results page (to visually inspect green bars).
    - Returns (20nt protospacer, CHOPCHOP score) for the top-ranked guide.
    """
    if not HAS_SELENIUM:
        return None, 0.0

    print("--> Launching CHOPCHOP Robot (Invisible)...")

    chrome_options = Options()
    chrome_options.add_argument("--headless")
    chrome_options.add_argument("--no-sandbox")
    chrome_options.add_argument("--disable-dev-shm-usage")
    chrome_options.add_argument("--log-level=3")

    service = Service(ChromeDriverManager().install())
    driver = webdriver.Chrome(service=service, options=chrome_options)

    try:
        # Open CHOPCHOP main page
        driver.get("https://chopchop.cbu.uib.no/")
        wait = WebDriverWait(driver, 20)

        # --- Paste sequence into FASTA input ---
        wait.until(EC.element_to_be_clickable((By.ID, "fasta"))).click()
        input_box = wait.until(EC.visibility_of_element_located((By.ID, "fastaInput")))
        input_box.clear()
        input_box.send_keys(sequence)

        # --- Select Y. lipolytica W29 genome ---
        # This uses the same genome ID we tested manually in the browser
        driver.execute_script("$('#genomeSelect').selectpicker('val', 'Ylipolytica_W29');")

        # --- Submit the job ---
        submit_btn = driver.find_element(By.ID, "searchRequest")
        driver.execute_script("arguments[0].click();", submit_btn)

        # --- Wait for the results table to appear ---
        rows = WebDriverWait(driver, 60).until(
            EC.presence_of_all_elements_located((By.CSS_SELECTOR, "tbody tr"))
        )

        # ============================
        # 1) Export full guide table as CSV
        # ============================
        try:
            header_cells = driver.find_elements(By.CSS_SELECTOR, "table thead tr th")
            if header_cells:
                headers = [h.text.strip() or f"col{i+1}" for i, h in enumerate(header_cells)]
            else:
                first_cols = rows[0].find_elements(By.TAG_NAME, "td")
                headers = [f"col{i+1}" for i in range(len(first_cols))]

            table_data = []
            for r in rows:
                cols2 = r.find_elements(By.TAG_NAME, "td")
                table_data.append([c.text.strip() for c in cols2])

            if gene_label:
                csv_name = f"{gene_label}_chopchop.csv"
                with open(csv_name, "w", newline="") as csvfile:
                    writer = csv.writer(csvfile)
                    writer.writerow(headers)
                    writer.writerows(table_data)
                print(f"    [Robot] Saved CHOPCHOP table as {csv_name}")
        except Exception as e:
            print(f"    [Robot] Could not export CHOPCHOP table: {e}")

        # ============================
        # 2) Save a PNG screenshot of the results page
        #    (top-of-page view with colored bars + rank-1 hit)
        # ============================
        try:
            # Make the window reasonably large for a nicer capture
            driver.set_window_size(1400, 1000)

            # Go to the very top of the page so the view matches what you see
            # when you open the results in a browser (summary plot + top rows)
            driver.execute_script("window.scrollTo(0, 0);")

            # Give the browser a moment to re-render
            time.sleep(1.0)

            png_name = f"{gene_label}_chopchop.png" if gene_label else "chopchop_results.png"

            # Capture the full browser window (including colored score bars
            # and the first-ranked guide row)
            driver.save_screenshot(png_name)
            print(f"    [Robot] Saved CHOPCHOP screenshot as {png_name}")
        except Exception as e:
            print(f"    [Robot] Could not save screenshot: {e}")

        # ============================
        # 3) Extract the top-ranked guide from row 1
        # ============================
        cols = rows[0].find_elements(By.TAG_NAME, "td")
        raw_seq = cols[1].text.strip().upper()  # usually 23 bp = 20-mer + NGG
        raw_score = cols[-1].text.strip()

        try:
            score_val = float(raw_score)
        except Exception:
            score_val = 0.0

        # CHOPCHOP usually reports 20-mer+NGG (23 nt). We only keep the 20-mer.
        if len(raw_seq) == 23:
            final_guide = raw_seq[:20]
        else:
            final_guide = raw_seq

        print(f"    [Robot] Found guide: {final_guide} (Score: {score_val})")
        return final_guide, score_val

    except Exception as e:
        print(f"    [Robot] Failed: {e}")
        return None, 0.0

    finally:
        driver.quit()

from typing import Tuple

def find_sgrna_in_orf(genome: Dict[str, str], gene: GeneInfo, window_size: int = 200) -> Tuple[Optional[SgRNA], str]:
    chrom_seq = genome[gene.chrom]

    # 1. PREPARE SEARCH WINDOW AROUND THE 5′ THIRD OF THE ORF
    window_half = window_size // 2
    orf_len = gene.end - gene.start + 1

    # Five-prime end depends on strand
    if gene.strand == '+':
        # 5′ end is at gene.start; move ~30% into the ORF
        five_prime_center = gene.start + int(orf_len * 0.3)
    else:
        # 5′ end is at gene.end on the reverse strand
        five_prime_center = gene.end - int(orf_len * 0.3)

    # Clamp to chromosome bounds (five_prime_center is 1-based)
    five_prime_center = max(1, min(len(chrom_seq), five_prime_center))

    # Convert to 0-based indices for slicing
    center_idx = five_prime_center - 1
    search_start = max(0, center_idx - window_half)
    search_end = min(len(chrom_seq), center_idx + window_half)
    search_seq = chrom_seq[search_start:search_end]
    
    # 2. TRY CHOPCHOP ROBOT
    guide_seq, robot_score = get_chopchop_guide(search_seq, gene.gene_id)
    
    if guide_seq:
        # Verify against local genome to get coordinates
        idx_fwd = chrom_seq.find(guide_seq)
        orf_start_0 = gene.start - 1
        orf_end_0 = gene.end  # one past the last index in 0-based slicing

        if idx_fwd != -1:
            # Enforce: inside ORF and at least 50 bp away from 5′ end on plus strand
            if (
                orf_start_0 <= idx_fwd and
                idx_fwd + 20 <= orf_end_0 and
                idx_fwd >= orf_start_0 + 50 and
                chrom_seq[idx_fwd + 20: idx_fwd + 23].endswith("GG")
            ):
                cut_pos = idx_fwd + 18  # SpCas9 cut 3 bp upstream of PAM
                return SgRNA(
                    guide_seq,
                    chrom_seq[idx_fwd + 20: idx_fwd + 23],
                    gene.chrom,
                    cut_pos,
                    "+",
                    robot_score,
                    "CHOPCHOP-Web",
                ), search_seq

        # Try reverse complement mapping
        rc = revcomp(guide_seq)
        idx_rev = chrom_seq.find(rc)
        if idx_rev != -1:
            # Enforce: inside ORF and at least 50 bp away from 5′ end on minus strand
            # ORF spans [orf_start_0, orf_end_0-1] in 0-based indices.
            # The 5′ end on the minus strand corresponds to the top 50 bp near gene.end.
            first50_threshold = gene.end - 50  # 0-based index where "safe zone" begins
            if (
                orf_start_0 <= idx_rev and
                idx_rev + 20 <= orf_end_0 and
                idx_rev + 20 <= first50_threshold and
                idx_rev >= 3 and
                chrom_seq[idx_rev - 3: idx_rev].startswith("CC")
            ):
                cut_pos = idx_rev + 4  # your previous minus-strand cut position
                return SgRNA(
                    guide_seq,
                    "NGG",
                    gene.chrom,
                    cut_pos,
                    "-",
                    robot_score,
                    "CHOPCHOP-Web",
                ), search_seq

        print("    [Robot] Guide violates ORF/5′-50bp constraints. Switching to local calculation...")

    # 3. FALLBACK: LOCAL SEARCH (Doench 2014)
    print("--> Running Local Doench 2014...")
    if gene.strand == '+':
        orf_seq = chrom_seq[gene.start - 1 : gene.end]
        abs_start = gene.start - 1
        strand = '+'
    else:
        orf_seq = revcomp(chrom_seq[gene.start - 1 : gene.end])
        abs_start = gene.start - 1
        strand = '-'

    if len(orf_seq) < 23:
        return None, search_seq
    
    # Center window on the 5′ third of the ORF (in ORF orientation)
    half = window_size // 2
    center = int(len(orf_seq) * 0.3)
    start_idx = max(0, center - half)
    end_idx = min(len(orf_seq), start_idx + window_size)
    search_region = orf_seq[start_idx:end_idx]
    
    candidates: List[Tuple[float, SgRNA]] = []
    for i in range(len(search_region) - 23):
        j = start_idx + i

        # Skip guides in the first 50 bp of the ORF to avoid cutting too close to the start
        if j < 50:
            continue

        protospacer = orf_seq[j:j+20]
        pam = orf_seq[j+20:j+23]
        
        if pam[1:] != "GG":
            continue
        if "TTTT" in protospacer:
            continue
        
        # Calculate Doench 2014 score (need 30 bp context)
        if strand == '+':
            ctx_start = abs_start + j - 4
            ctx_end = abs_start + j + 26
            ctx = chrom_seq[max(0, ctx_start):ctx_end]
        else:
            if j >= 4 and (j + 26) <= len(orf_seq):
                ctx = orf_seq[j-4:j+26]
            else:
                ctx = ""
            
        score = calculate_doench_score(str(ctx)) if len(ctx) == 30 else 50
        
        if gene.strand == '+':
            cut_pos = abs_start + j + 18
        else:
            cut_pos = gene.start + (len(orf_seq) - (j + 20 + 3)) + 3
        
        candidates.append(
            (score, SgRNA(protospacer, pam, gene.chrom, cut_pos, strand, score, "Local-Doench"))
        )

    if not candidates:
        return None, search_seq
    candidates.sort(key=lambda x: x[0], reverse=True)
    return candidates[0][1], search_seq

def extract_homology_arms(genome: Dict[str, str], gene: GeneInfo, arm_len: int) -> Tuple[str, str]:
    chrom_seq = genome[gene.chrom]
    if gene.strand == '+':
        up = chrom_seq[max(0, gene.start - 1 - arm_len) : gene.start - 1]
        down = chrom_seq[gene.end : gene.end + arm_len]
    else:
        up = revcomp(chrom_seq[gene.end : gene.end + arm_len])
        down = revcomp(chrom_seq[max(0, gene.start - 1 - arm_len) : gene.start - 1])
    return up, down


# =========================
# Primer Design + Math
# =========================

def find_best_primer_with_pos(sequence: str, direction: str) -> Tuple[str, int]:
    """
    Returns (primer_seq, start_index_in_sequence).
    """
    seq_to_scan = sequence.upper()
    
    # 1. Scan for ideal primer
    for i in range(0, len(seq_to_scan) - 18):
        for length in range(18, 25): 
            if i + length > len(seq_to_scan): break
            cand = seq_to_scan[i : i+length]
            gc = gc_content(cand)
            
            if 0.4 <= gc <= 0.6 and "TTTT" not in cand and "AAAA" not in cand:
                final_seq = cand if direction == 'fw' else revcomp(cand)
                return final_seq, i
    
    # 2. Fallback
    fallback = seq_to_scan[:20]
    final_seq = fallback if direction == 'fw' else revcomp(fallback)
    return final_seq, 0
def design_truncation(
    genome: Dict[str, str],
    gene: GeneInfo,
    trunc_len: int,
    arm_len: int = 162
) -> KnockoutDesign:
    """
    TruncYALI design:
      - Assume a 1 kb promoter upstream of the ATG (or downstream for '-' strand).
      - User chooses trunc_len (50–500 bp), which is the length of promoter
        we KEEP closest to the ATG.
      - We delete the distal part of the promoter and use:
          UP  arm:  162 bp outside the -1000 border
          DOWN arm: 162 bp starting at the truncation point X (-trunc_len)
      - For '-' strand, the same logic is applied in genomic coordinates,
        then arms are reverse-complemented for the donor.
    """
    chrom_seq = genome[gene.chrom]
    chrom_len = len(chrom_seq)

    if trunc_len < 50 or trunc_len > 500:
        raise ValueError("trunc_len should be between 50 and 500 bp.")

    # -----------------------------
    # Define promoter and truncation in genomic coordinates (0-based)
    # -----------------------------
    if gene.strand == '+':
        # ORF start (ATG) in 0-based coordinates
        tss0 = gene.start - 1  # position of the first coding base
        # Full 1 kb promoter upstream of ATG (clamped to contig start)
        prom_start0 = max(0, tss0 - 1000)   # corresponds to "-1000"
        prom_end0 = tss0                    # corresponds to "-1" (exclusive)
        prom_len = max(0, prom_end0 - prom_start0)

        # Clamp trunc_len so it never exceeds the usable promoter or the arm length
        max_trunc = max(1, min(prom_len, arm_len - 1))
        if trunc_len > max_trunc:
            print(
                f"[TruncYALI] Requested trunc_len {trunc_len} bp is too long for "
                f"{gene.gene_id} on '+' strand (promoter {prom_len} bp, arm_len {arm_len} bp). "
                f"Using trunc_len={max_trunc} bp instead."
            )
            trunc_len = max_trunc

        # Truncation point X = -trunc_len → genomic coordinate:
        # keep [-trunc_len .. -1], so X = tss0 - trunc_len
        trunc0 = tss0 - trunc_len  # this is the first base we KEEP

        # UP arm: -1000 - 162 → 162 bp immediately upstream of prom_start0
        up_start = max(0, prom_start0 - arm_len)
        up_end = prom_start0
        up_seq = chrom_seq[up_start:up_end]

        # DOWN arm: from -trunc_len (e.g. -50) across the ATG into the ORF by (arm_len - trunc_len)
        down_start = max(0, tss0 - trunc_len)
        down_end = min(chrom_len, tss0 + (arm_len - trunc_len))
        down_seq = chrom_seq[down_start:down_end]

        # Region actually deleted for PCR math: [prom_start0 .. trunc0)
        del_start = prom_start0
        del_end = trunc0

    else:
        # '-' strand: treat everything in an ORF-relative frame, then map to genomic.
        # r = 0 at the first coding base (ATG), promoter is r in [-1000 .. -1].
        tss0 = gene.end - 1  # 0-based index of first coding base on '-' strand

        def r_to_g(r: int) -> int:
            """
            Map ORF-relative coordinate r to genomic 0-based index.
            For the minus strand, going 'upstream' in ORF space (more negative r)
            means increasing genomic coordinate.
            """
            return tss0 - r

        # Approximate 1 kb promoter: r = -1000 .. -1.
        # Compute its span in genomic coordinates and clamp to contig bounds.
        prom_start_g_raw = r_to_g(-1000)
        prom_end_g_raw = r_to_g(-1) + 1  # exclusive

        # Normalise order and clamp to chromosome bounds
        prom_left = max(0, min(prom_start_g_raw, prom_end_g_raw))
        prom_right = min(chrom_len, max(prom_start_g_raw, prom_end_g_raw))
        prom_start_g, prom_end_g = prom_left, prom_right
        prom_len = max(0, prom_end_g - prom_start_g)

        # Clamp trunc_len so it never exceeds usable promoter or arm length
        max_trunc = max(1, min(prom_len, arm_len - 1))
        if trunc_len > max_trunc:
            print(
                f"[TruncYALI] Requested trunc_len {trunc_len} bp is too long for "
                f"{gene.gene_id} on '-' strand (promoter {prom_len} bp, arm_len {arm_len} bp). "
                f"Using trunc_len={max_trunc} bp instead."
            )
            trunc_len = max_trunc

        # Truncation boundary at r = -trunc_len (this is the first bp we KEEP in promoter).
        r_trunc = -trunc_len

        # UP arm: 162 bp just outside the -1000 border → r in [-1000-162 .. -1001]
        r_up_start = -1000 - arm_len
        r_up_end = -1000  # exclusive
        g_up_1 = r_to_g(r_up_start)
        g_up_2 = r_to_g(r_up_end - 1)
        up_start_gen = max(0, min(g_up_1, g_up_2))
        up_end_gen = min(chrom_len, max(g_up_1, g_up_2) + 1)
        up_seq_gen = chrom_seq[up_start_gen:up_end_gen]
        # Donor carries reverse complement of genomic segment
        up_seq = revcomp(up_seq_gen)

        # DOWN arm: 162 bp starting at r = -trunc_len, extending across the ATG
        # into the ORF (e.g. for trunc_len=50 → [-50 .. +112] in ORF coordinates).
        r_dn_start = r_trunc
        r_dn_end = r_trunc + arm_len  # exclusive
        g_dn_1 = r_to_g(r_dn_start)
        g_dn_2 = r_to_g(r_dn_end - 1)
        down_start_gen = max(0, min(g_dn_1, g_dn_2))
        down_end_gen = min(chrom_len, max(g_dn_1, g_dn_2) + 1)
        down_seq_gen = chrom_seq[down_start_gen:down_end_gen]
        down_seq = revcomp(down_seq_gen)

        # Deleted promoter segment corresponds to r in [-1000 .. -trunc_len-1].
        g_del_1 = r_to_g(-1000)
        g_del_2 = r_to_g(-trunc_len - 1)
        del_start = min(g_del_1, g_del_2)
        del_end = max(g_del_1, g_del_2) + 1

    # -----------------------------
    # Length of the deleted piece (for band size difference)
    # -----------------------------
    deleted_len = max(0, del_end - del_start)

    # -----------------------------
    # sgRNA selection in the DELETED PROMOTER (not the ORF)
    # -----------------------------
    # We want the cut to fall roughly in the first third of the DELETED promoter,
    # counted from the UP (5′) homology side:
    #   - For the '+' strand, the UP side is at the more upstream index (del_start,
    #     corresponding to ~-1000 relative to the ATG).
    #   - For the '-' strand, the UP side is at the more upstream index in gene
    #     space, which corresponds to the larger genomic index (del_end - 1).
    # We define a 200 bp window centered around that 5′-third position and send
    # ONLY that fragment to CHOPCHOP so that the cut is always close to the UP arm.
    prom_window = 200
    window_half = prom_window // 2

    if deleted_len <= 0:
        raise ValueError(f"Deleted promoter segment has non‑positive length for {gene.gene_id}.")

    if gene.strand == "+":
        # 5′/UP side is del_start (most upstream base of the deleted promoter).
        up_index = del_start
        target_center = up_index + deleted_len // 3
    else:
        # 5′/UP side is at the larger genomic index (del_end - 1).
        up_index = del_end - 1
        target_center = up_index - deleted_len // 3

    # Clamp target_center to lie strictly inside the deleted segment.
    target_center = max(del_start, min(del_end - 1, target_center))

    # Define the genomic window used for CHOPCHOP and local fallback.
    search_start = max(del_start, target_center - window_half)
    search_end = min(del_end, target_center + window_half)
    search_window = chrom_seq[search_start:search_end]

    # --- Primary: CHOPCHOP over the deleted promoter window ---
    sgrna = None
    guide_seq, robot_score = get_chopchop_guide(search_window, gene.gene_id + "_trunc")

    if guide_seq:
        # Map the guide back to the genome but RESTRICTED to the promoter window,
        # so we never pick a hit in the ORF.
        idx_fwd = chrom_seq.find(guide_seq, search_start, search_end)
        if idx_fwd != -1 and chrom_seq[idx_fwd + 20: idx_fwd + 23].endswith("GG"):
            cut_pos = idx_fwd + 18  # same convention as KOYALI
            sgrna = SgRNA(
                guide_seq,
                chrom_seq[idx_fwd + 20: idx_fwd + 23],
                gene.chrom,
                cut_pos,
                "+",
                robot_score,
                "CHOPCHOP-Web",
            )
        else:
            # Try reverse‑complement orientation inside the same window.
            rc = revcomp(guide_seq)
            idx_rev = chrom_seq.find(rc, search_start, search_end)
            if (
                idx_rev != -1
                and idx_rev >= 3
                and chrom_seq[idx_rev - 3: idx_rev].startswith("CC")
            ):
                cut_pos = idx_rev + 4  # use same minus‑strand convention as in KOYALI
                sgrna = SgRNA(
                    guide_seq,
                    "NGG",
                    gene.chrom,
                    cut_pos,
                    "-",
                    robot_score,
                    "CHOPCHOP-Web",
                )

    # --- Fallback: if CHOPCHOP fails or cannot be mapped, revert to ORF‑based search ---
    if sgrna is None:
        print("    [TruncYALI] CHOPCHOP promoter search failed; falling back to ORF‑based guide search.")
        sgrna, fallback_window = find_sgrna_in_orf(genome, gene)
        if not sgrna:
            raise ValueError(f"No suitable sgRNA found for {gene.gene_id}")
        # For the report, still keep the promoter window if it exists.
        if search_window:
            search_window = search_window
        else:
            search_window = fallback_window

    # -----------------------------
    # Primer design OUTSIDE the homology arms (as in KOYALI)
    # -----------------------------
    # We'll place primers flanking the entire edit region: from just outside
    # the UP arm to just outside the DOWN arm, analogous to KOYALI's "out-out" PCR.

    # For simplicity, define the bounding region of the HDR event on the genome.
    # On '+' this is [up_start .. down_end), on '-' we use the genomic coordinates
    # we actually took the arms from.
    if gene.strand == '+':
        hdr_left = up_start
        hdr_right = down_end
    else:
        # For '-' strand, arms were taken from up_start_gen/up_end_gen and
        # down_start_gen/down_end_gen.
        hdr_left = min(up_start_gen, down_start_gen)
        hdr_right = max(up_end_gen, down_end_gen)

    # Define left and right search windows for primers (outside the HDR region)
    chrom_seq = genome[gene.chrom]
    search_window_len = 200
    buffer = 10

    # Left side
    p_left_end = max(0, hdr_left - buffer)
    p_left_start = max(0, p_left_end - search_window_len)
    left_search_seq = chrom_seq[p_left_start:p_left_end]

    # Right side
    p_right_start = min(chrom_len, hdr_right + buffer)
    p_right_end = min(chrom_len, p_right_start + search_window_len)
    right_search_seq = chrom_seq[p_right_start:p_right_end]

    fw_seq, fw_rel_pos = find_best_primer_with_pos(left_search_seq, 'fw')
    rv_seq, rv_rel_pos = find_best_primer_with_pos(right_search_seq, 'rv')

    fw_genomic_start = p_left_start + fw_rel_pos
    rv_genomic_end = p_right_start + rv_rel_pos + len(rv_seq)

    band_unedited = rv_genomic_end - fw_genomic_start
    band_trunc = band_unedited - deleted_len

    # -----------------------------
    # Build synthetic fragment:
    # GA_LEFT – sgRNA – scaffold – terminator – UP – DOWN – GA_RIGHT
    # -----------------------------
    fragment_seq = (
        GA_LEFT
        + sgrna.protospacer.upper()
        + TRACR_SCAFFOLD
        + TER_RPR1
        + up_seq
        + down_seq
        + GA_RIGHT
    )

    primers = [
        PrimerInfo("TRUNC_VERIFY_FWD", fw_seq, fw_genomic_start, len(fw_seq)),
        PrimerInfo("TRUNC_VERIFY_REV", rv_seq, rv_genomic_end, len(rv_seq)),
    ]

    return KnockoutDesign(
        gene=gene,
        sgrna=sgrna,
        upstream_arm=up_seq,
        downstream_arm=down_seq,
        fragment_seq=fragment_seq,
        primers=primers,
        band_size_unedited=band_unedited,
        band_size_ko=band_trunc,
        search_window_seq=search_window,
    )
# =========================
# Output Generators
# =========================

def write_design_report(design: KnockoutDesign, out_path: str) -> None:
    with open(out_path, "w") as f:
        sg = design.sgrna
        f.write(f"=== TruncYALI Truncation: {design.gene.gene_id} ===\n\n")
        f.write(f"[Target Details]\n")
        f.write(f"Gene ID: {design.gene.gene_id}\n")
        f.write(f"Method: {sg.method} (Score: {sg.score})\n")
        f.write(f"Cut Pos: {sg.cut_pos} ({sg.strand})\n")
        f.write(f"Protospacer: {sg.protospacer}\n\n")

        fwd_oligo = "TTCGATTCCGGGTCGGCGCAGGTTG" + sg.protospacer + "GTTTTA"
        rev_oligo = "GCTCTAAAAC" + revcomp(sg.protospacer) + "CAACCTGCGCCGACCCGGAAT"
        f.write(f"[gRNA Cloning Oligos]\nFwd: {fwd_oligo}\nRev: {rev_oligo}\n\n")

        # Insert the CHOPCHOP 200 bp window sequence section
        f.write(f"[CHOPCHOP 200 bp window (~first third of DELETED promoter)]\n")
        f.write("Paste this sequence into CHOPCHOP (FASTA input) if you want to verify the guides manually.\n")
        f.write(f"{design.search_window_seq}\n\n")

        f.write(f"[Donor Fragment]\n>TRUNC_{design.gene.gene_id}_Donor\n{design.fragment_seq}\n\n")

        f.write("[Verification Primers]\n")
        f.write("Primers are outside homology arms and detect a size shift corresponding to the deleted 5′ segment.\n")
        for p in design.primers:
            f.write(f"{p.name}: {p.seq} (Len: {p.length})\n")

        f.write(f"\n[Expected PCR Bands]\n")
        f.write(f"Unedited/Parental: ~{design.band_size_unedited} bp\n")
        f.write(f"Truncation allele: ~{design.band_size_ko} bp (shorter by ~{design.band_size_unedited - design.band_size_ko} bp)\n")


def write_genbank(design: KnockoutDesign, gb_path: str) -> None:
    record = SeqRecord(
        Seq(design.fragment_seq),
        id=f"TRUNC_{design.gene.gene_id}_donor",
        name=f"TRUNC_{design.gene.gene_id}",
        description=f"TruncYALI donor for 5' truncation of {design.gene.gene_id}. Expected bands: unedited~{design.band_size_unedited}bp, trunc~{design.band_size_ko}bp."
    )
    today = datetime.today().strftime("%d-%b-%Y").upper()
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "linear"
    record.annotations["data_file_division"] = "UNC"
    record.annotations["date"] = today
    record.annotations["accessions"] = [f"TRUNC_{design.gene.gene_id}_DONOR"]
    record.annotations["source"] = "synthetic DNA construct"
    record.annotations["organism"] = "synthetic construct"
    record.annotations["comment"] = "Designed with TruncYALI."

    spacer_len = len(design.sgrna.protospacer)
    modules = [
        (len(GA_LEFT), "misc_feature", "GA_left", "#808080"),
        (spacer_len, "misc_feature", f"sgRNA_{design.gene.gene_id}", "#FF0000"),
        (len(TRACR_SCAFFOLD), "misc_feature", "tracrRNA_scaffold", "#0000FF"),
        (len(TER_RPR1), "terminator", "TerRPR1", "#800080"),
        (len(design.upstream_arm), "misc_recomb", f"UP_{design.gene.gene_id}", "#008000"),
        (len(design.downstream_arm), "misc_recomb", f"DOWN_{design.gene.gene_id}", "#006400"),
        (len(GA_RIGHT), "misc_feature", "GA_right", "#808080"),
    ]

    pos = 0
    for length, f_type, label, color in modules:
        feature = SeqFeature(
            FeatureLocation(pos, pos + length),
            type=f_type,
            qualifiers={
                "label": [label],
                "ApEinfo_fwdcolor": [color],
                "ApEinfo_revcolor": [color],
                "note": [f"TruncYALI Module: {label}"]
            },
        )
        record.features.append(feature)
        pos += length

    SeqIO.write(record, gb_path, "genbank")


def write_png_map(design: KnockoutDesign, png_path: str) -> None:
    # TruncYALI donor map diagram
    if not HAS_DNA_VIEWER: return

    L_ga_left = len(GA_LEFT)
    L_spacer = 20
    L_tracr = len(TRACR_SCAFFOLD)
    L_term = len(TER_RPR1)
    L_up = len(design.upstream_arm)
    L_down = len(design.downstream_arm)
    L_ga_right = len(GA_RIGHT)

    features = []
    i = 0
    map_defs = [
        (L_ga_left, "#808080", "GA_left"),
        (L_spacer, "#FF0000", f"sgRNA_{design.gene.gene_id}"),
        (L_tracr, "#0000FF", "Scaffold"),
        (L_term, "#800080", "Term"),
        (L_up, "#008000", "UP Arm"),
        (L_down, "#006400", "DOWN Arm"),
        (L_ga_right, "#808080", "GA_right")
    ]

    for length, color, label in map_defs:
        features.append(GraphicFeature(start=i, end=i+length, strand=0, color=color, label=label))
        i += length

    record = GraphicRecord(sequence=design.fragment_seq, features=features)
    ax, _ = record.plot(figure_width=10)
    ax.figure.tight_layout()
    ax.figure.savefig(png_path, dpi=300)
    plt.close(ax.figure)

# =========================
# [NEW] Gel Simulation
# =========================

def write_gel_simulation(design: KnockoutDesign, gel_path: str) -> None:
    """
    Creates a simulated agarose gel image using pydna.
    Lanes: Ladder | Unedited (WT) | Knockout (KO)
    """
    if not HAS_PYDNA:
        return

    # Create dummy sequences with physical properties (Dseqrecord)
    wt_band = Dseqrecord("A" * design.band_size_unedited)
    wt_band.name = "Unedited"
    
    ko_band = Dseqrecord("A" * design.band_size_ko)
    ko_band.name = "KO"

    # Setup Lanes
    lanes = [
        GeneRuler_1kb,
        [wt_band],
        [ko_band]
    ]

    try:
        gel_img = gel(lanes)
        gel_img.save(gel_path)
    except Exception as e:
        print(f"Could not save gel image: {e}")


# =========================
# Main
# =========================

def main():
    parser = argparse.ArgumentParser(description="TruncYALI Designer (Calculated Bands for 5' Truncations)")
    parser.add_argument("--gene_id", required=True, help="Gene ID (e.g. YALI1_B12345g)")
    parser.add_argument("--genome", default="W29_genome.fasta", help="Genome FASTA")
    parser.add_argument("--gff", default="W29_annotation.gff", help="Annotation GFF")
    parser.add_argument("--out", default=None, help="Output Base Name")
    parser.add_argument("--trunc_len", type=int, default=50, help="Length (bp) to truncate at the 5' end of the ORF (default: 50 bp)")
    args = parser.parse_args()

    ensure_reference_files(args.genome, args.gff)
    genome = load_genome_fasta(args.genome)
    genes = load_annotation_gff(args.gff)

    if args.gene_id not in genes:
        raise KeyError(f"Gene {args.gene_id} not found in GFF.")

    print(f"Designing 5' truncation for {args.gene_id} (TruncYALI)...")
    design = design_truncation(genome, genes[args.gene_id], trunc_len=args.trunc_len)

    base = args.out or f"{args.gene_id}_trunc_cassette"
    if base.endswith(".txt"): base = base[:-4]

    # Write all outputs
    write_design_report(design, f"{base}.txt")
    write_genbank(design, f"{base}.gb")
    write_png_map(design, f"{base}.png")

    # New Gel Simulation
    write_gel_simulation(design, f"{base}_gel.png")

    print(f"Design Complete (TruncYALI) for {args.gene_id}")
    print(f" - Method:             {design.sgrna.method}")
    print(f" - Score:              {design.sgrna.score}/100")
    print(f" - Unedited Band:      ~{design.band_size_unedited} bp")
    print(f" - Truncation Band:    ~{design.band_size_ko} bp")
    print(f" - Report:             {base}.txt")

    if HAS_PYDNA:
        print(f" - Gel Check:          {base}_gel.png")

if __name__ == "__main__":
    main()