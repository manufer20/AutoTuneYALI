#!/usr/bin/env python3

"""
TuneYALI-style promoter insertion cassette designer for Yarrowia lipolytica (W29/CLIB89).

What TuneYALI builds
--------------------
This module designs a short-arm HDR donor in the TuneYALI cassette layout:

  GA_LEFT — gRNA (20nt) — tracrRNA scaffold — TerRPR1 — UP arm — [GG INSERT] — DOWN arm — GA_RIGHT

For TuneYALI, the homology arms are anchored around a fixed promoter insertion point
upstream of the annotated start codon:
  - UP arm:   (-tune_len-162) .. (-tune_len)   (default tune_len = 500)
  - DOWN arm: 0 .. +162 (first 162 bp of the ORF)

Between arms, a Golden-Gate placeholder is inserted:
  GG_INSERT = GACTGAAGAGCTCTTCA

Neighbor-aware safety rules
--------------------------
We re-use the same neighbour-aware intergenic safety logic as TruncYALI:
  * On '+' strand, shrink the usable upstream promoter so the design does not encroach
    into upstream annotated features (buffers configurable).
  * On '-' strand, the promoter is downstream in genomic coordinates; shrink the usable
    downstream promoter to avoid downstream annotated features.

If there is not enough *safe* promoter sequence to support tune_len=500 (i.e. you
cannot place the UP arm at -662..-500), TuneYALI aborts by default.
Two opt-in behaviors are provided:
  - --use_available : reduce tune_len to the maximum allowed by the safe intergenic space
  - --force         : ignore neighbour-based shrinking (still bounded by contig edges)

Guide selection
---------------
The guide is selected inside the promoter region that will be replaced (the proximal
segment from -tune_len to 0), prioritizing a cut near the UP-arm side of the deleted
segment. CHOPCHOP is used when available; local Doench-2014 scoring is used as fallback.
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

# Promoter truncation rules
PROMOTER_SPAN = 1000         # up to 1 kb upstream treated as "promoter"
BUFFER_SAME_STRAND = 100     # keep ≥100 bp after 3' end of same-strand upstream gene
BUFFER_OPP_STRAND = 500      # keep ≥500 bp after upstream "start" of opposite strand gene
GG_INSERT = "GACTGAAGAGCTCTTCA"  # Golden Gate placeholder between UP and DOWN arms


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


# -------------------------
# Global gene index
# -------------------------

ALL_GENES: Dict[str, GeneInfo] = {}
GENES_BY_CHROM: Dict[str, List[GeneInfo]] = {}


def index_genes_by_chrom(genes: Dict[str, GeneInfo]) -> Dict[str, List[GeneInfo]]:
    """
    Build a per-chromosome index of GeneInfo objects, sorted by start coordinate.
    Used to find upstream neighbours when deciding how far we can truncate.
    """
    by_chrom: Dict[str, List[GeneInfo]] = {}
    for g in genes.values():
        by_chrom.setdefault(g.chrom, []).append(g)
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda x: x.start)
    return by_chrom


def adjust_promoter_start_for_neighbors_plus(
    gene: GeneInfo,
    prom_start: int,
    prom_end: int
) -> int:
    """
    For '+' strand targets:
      - promoter window is [prom_start, prom_end) in 0-based genomic coordinates.
      - We look for the last gene upstream (end <= gene.start) on the same chrom.
      - We then enforce:
          del_start >= last_same_end + 100
          del_start >= last_opp_end  + 500
      by shifting prom_start up to satisfy these constraints.
    """
    if prom_end <= prom_start:
        return prom_start

    neighbors = GENES_BY_CHROM.get(gene.chrom, [])
    if not neighbors:
        return prom_start

    nearest_end_same: Optional[int] = None
    nearest_end_opp: Optional[int] = None

    for g2 in neighbors:
        if g2.gene_id == gene.gene_id:
            continue
        # Only fully upstream features (end <= target start)
        if g2.end <= gene.start:
            if g2.strand == gene.strand:
                if nearest_end_same is None or g2.end > nearest_end_same:
                    nearest_end_same = g2.end
            else:
                if nearest_end_opp is None or g2.end > nearest_end_opp:
                    nearest_end_opp = g2.end

    safe_start = prom_start
    # g2.end is 1-based inclusive; the first base after the gene is at index g2.end (0-based).
    if nearest_end_same is not None:
        safe_start = max(safe_start, nearest_end_same + BUFFER_SAME_STRAND)
    if nearest_end_opp is not None:
        safe_start = max(safe_start, nearest_end_opp + BUFFER_OPP_STRAND)

    if safe_start > prom_start and safe_start < prom_end:
        usable = prom_end - safe_start
        print(
            f"[TruncYALI] Adjusted promoter window for {gene.gene_id} (+ strand) "
            f"to avoid upstream features. Old start: {prom_start}, new start: {safe_start} "
            f"(usable promoter: {usable} bp)."
        )

    # Clamp so that there's still *some* promoter left; the truncation logic itself
    # will enforce >= 50 bp for trunc_len.
    if safe_start >= prom_end:
        # Nothing left; caller will detect prom_len <= 0 and error out.
        return prom_end - 1

    return safe_start


# =========================
# Neighbor-aware shrinking for '-' strand promoter (downstream)
# =========================

def adjust_promoter_end_for_neighbors_minus(
    gene: GeneInfo,
    prom_start: int,
    prom_end: int,
) -> int:
    """\
    For '-' strand targets:
      - promoter window is [prom_start, prom_end) in 0-based genomic coordinates.
      - promoter lies *downstream* of the ORF (increasing genomic coordinate).
      - We look for the nearest downstream gene (start >= gene.end) on the same chrom.
      - We then enforce that the promoter window does not extend too close to that
        downstream feature:
          prom_end <= next_start0 - BUFFER_SAME_STRAND   (same-strand downstream gene)
          prom_end <= next_start0 - BUFFER_OPP_STRAND    (opposite-strand downstream gene)
        where next_start0 is the downstream gene's 0-based start index.

    This is a conservative, annotation-based safety rule intended to avoid deleting
    into a neighbour's promoter/ORF when intergenic regions are short.
    """
    if prom_end <= prom_start:
        return prom_end

    neighbors = GENES_BY_CHROM.get(gene.chrom, [])
    if not neighbors:
        return prom_end

    nearest_start_same: Optional[int] = None
    nearest_start_opp: Optional[int] = None

    for g2 in neighbors:
        if g2.gene_id == gene.gene_id:
            continue
        # Only downstream features on the promoter side for '-' strand targets.
        if g2.start >= gene.end:
            if g2.strand == gene.strand:
                if nearest_start_same is None or g2.start < nearest_start_same:
                    nearest_start_same = g2.start
            else:
                if nearest_start_opp is None or g2.start < nearest_start_opp:
                    nearest_start_opp = g2.start

    safe_end = prom_end

    # g2.start is 1-based inclusive; its 0-based start index is (g2.start - 1).
    if nearest_start_same is not None:
        limit_same = (nearest_start_same - 1) - BUFFER_SAME_STRAND
        safe_end = min(safe_end, max(prom_start, limit_same))

    if nearest_start_opp is not None:
        limit_opp = (nearest_start_opp - 1) - BUFFER_OPP_STRAND
        safe_end = min(safe_end, max(prom_start, limit_opp))

    if safe_end < prom_end and safe_end > prom_start:
        usable = safe_end - prom_start
        print(
            f"[TruncYALI] Adjusted promoter window for {gene.gene_id} (- strand) "
            f"to avoid downstream features. Old end: {prom_end}, new end: {safe_end} "
            f"(usable promoter: {usable} bp)."
        )

    # If bounds collapse, return prom_start so caller can error out cleanly.
    if safe_end <= prom_start:
        return prom_start

    return safe_end


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
        # ============================
        try:
            driver.set_window_size(1400, 1000)
            driver.execute_script("window.scrollTo(0, 0);")
            time.sleep(1.0)
            png_name = f"{gene_label}_chopchop.png" if gene_label else "chopchop_results.png"
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

    five_prime_center = max(1, min(len(chrom_seq), five_prime_center))

    center_idx = five_prime_center - 1
    search_start = max(0, center_idx - window_half)
    search_end = min(len(chrom_seq), center_idx + window_half)
    search_seq = chrom_seq[search_start:search_end]
    
    # 2. TRY CHOPCHOP ROBOT
    guide_seq, robot_score = get_chopchop_guide(search_seq, gene.gene_id)
    
    if guide_seq:
        idx_fwd = chrom_seq.find(guide_seq)
        orf_start_0 = gene.start - 1
        orf_end_0 = gene.end

        if idx_fwd != -1:
            if (
                orf_start_0 <= idx_fwd and
                idx_fwd + 20 <= orf_end_0 and
                idx_fwd >= orf_start_0 + 50 and
                chrom_seq[idx_fwd + 20: idx_fwd + 23].endswith("GG")
            ):
                cut_pos = idx_fwd + 18
                return SgRNA(
                    guide_seq,
                    chrom_seq[idx_fwd + 20: idx_fwd + 23],
                    gene.chrom,
                    cut_pos,
                    "+",
                    robot_score,
                    "CHOPCHOP-Web",
                ), search_seq

        rc = revcomp(guide_seq)
        idx_rev = chrom_seq.find(rc)
        if idx_rev != -1:
            first50_threshold = gene.end - 50
            if (
                orf_start_0 <= idx_rev and
                idx_rev + 20 <= orf_end_0 and
                idx_rev + 20 <= first50_threshold and
                idx_rev >= 3 and
                chrom_seq[idx_rev - 3: idx_rev].startswith("CC")
            ):
                cut_pos = idx_rev + 4
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
    
    half = window_size // 2
    center = int(len(orf_seq) * 0.3)
    start_idx = max(0, center - half)
    end_idx = min(len(orf_seq), start_idx + window_size)
    search_region = orf_seq[start_idx:end_idx]
    
    candidates: List[Tuple[float, SgRNA]] = []
    for i in range(len(search_region) - 23):
        j = start_idx + i
        if j < 50:
            continue

        protospacer = orf_seq[j:j+20]
        pam = orf_seq[j+20:j+23]
        
        if pam[1:] != "GG":
            continue
        if "TTTT" in protospacer:
            continue
        
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


def design_tuneyali(
    genome: Dict[str, str],
    gene: GeneInfo,
    tune_len: int = 500,
    arm_len: int = 162,
    force: bool = False,
    use_available: bool = False,
) -> KnockoutDesign:
    """Design a TuneYALI promoter insertion cassette.

    Default behaviour:
      - Replace the proximal promoter segment of length `tune_len` (default 500 bp)
        immediately upstream of the ATG with a Golden-Gate placeholder.
      - Require enough *safe* promoter sequence to place the UP arm at
        (-tune_len-162 .. -tune_len). If not available, raise ValueError.

    Optional behaviours:
      - use_available=True: reduce tune_len to the maximum allowed by the safe space.
      - force=True: ignore neighbour-aware shrinking (still respects contig boundaries).

    Notes on coordinates:
      - '+' strand: promoter is upstream (decreasing genomic coordinate).
      - '-' strand: promoter is downstream (increasing genomic coordinate).
    """
    chrom_seq = genome[gene.chrom]
    chrom_len = len(chrom_seq)

    if tune_len <= 0:
        raise ValueError("tune_len must be a positive integer.")

    # -----------------------------
    # Determine the usable promoter window with neighbour-aware safety margins
    # -----------------------------
    if gene.strand == '+':
        # 0-based index of first coding base
        tss0 = gene.start - 1

        # Nominal promoter span (we only *replace* the proximal tune_len bp, but we still
        # keep this cap to avoid extending too far if annotation is messy)
        prom_start0 = max(0, tss0 - PROMOTER_SPAN)
        prom_end0 = tss0

        # -----------------------------
        # Neighbour-aware safety should constrain the *REPLACED* segment (-tune_len..0),
        # not the UP homology arm. The UP arm can lie within the "reserved" distal
        # intergenic/promoter space because it is not edited.
        # -----------------------------
        min_del_start = prom_start0

        if GENES_BY_CHROM and not force:
            neighbors = GENES_BY_CHROM.get(gene.chrom, [])
            nearest_end_same: Optional[int] = None
            nearest_end_opp: Optional[int] = None

            for g2 in neighbors:
                if g2.gene_id == gene.gene_id:
                    continue
                # only fully upstream features (end <= target start)
                if g2.end <= gene.start:
                    if g2.strand == gene.strand:
                        if nearest_end_same is None or g2.end > nearest_end_same:
                            nearest_end_same = g2.end
                    else:
                        if nearest_end_opp is None or g2.end > nearest_end_opp:
                            nearest_end_opp = g2.end

            if nearest_end_same is not None:
                min_del_start = max(min_del_start, nearest_end_same + BUFFER_SAME_STRAND)
            if nearest_end_opp is not None:
                min_del_start = max(min_del_start, nearest_end_opp + BUFFER_OPP_STRAND)

        # Ensure the replacement window fits inside the safe intergenic space
        # del_start is the first base we REPLACE (0-based)
        del_start = tss0 - tune_len
        if del_start < min_del_start:
            if use_available:
                new_tune = tss0 - min_del_start
                if new_tune <= 0:
                    raise ValueError(
                        f"TuneYALI cannot be designed for {gene.gene_id}: safe replaceable promoter "
                        f"length is {tss0 - min_del_start} bp (<=0) after neighbour buffers."
                    )
                print(
                    f"[TuneYALI] Safe replaceable promoter upstream of {gene.gene_id} is {tss0 - min_del_start} bp; "
                    f"shrinking tune_len to {new_tune} bp (--use_available)."
                )
                tune_len = new_tune
                del_start = tss0 - tune_len
            else:
                raise ValueError(
                    f"TuneYALI cannot be designed safely for {gene.gene_id}: replaceable promoter length "
                    f"upstream of ATG is {tss0 - min_del_start} bp but needs at least tune_len={tune_len} bp. "
                    f"Use --use_available to shrink automatically, or --force to ignore neighbour constraints."
                )

        # Now place the UP arm immediately upstream of del_start
        up_start = del_start - arm_len
        if up_start < 0:
            if use_available:
                new_tune = max(1, tss0 - arm_len)
                if new_tune < tune_len:
                    print(
                        f"[TuneYALI] UP arm would extend beyond contig start for {gene.gene_id}; "
                        f"shrinking tune_len from {tune_len} to {new_tune} bp (--use_available)."
                    )
                    tune_len = new_tune
                    del_start = tss0 - tune_len
                    up_start = del_start - arm_len
            else:
                raise ValueError(
                    f"TuneYALI cannot be designed for {gene.gene_id}: UP arm would extend beyond contig start. "
                    f"Use --use_available to shrink tune_len, or --force with a smaller --tune_len."
                )

        up_end = del_start
        up_seq = chrom_seq[up_start:up_end]

        down_start = tss0
        down_end = min(chrom_len, tss0 + arm_len)
        down_seq = chrom_seq[down_start:down_end]

        # Deleted/replaced promoter segment: [del_start .. tss0)
        del_end = tss0

        hdr_left = up_start
        hdr_right = down_end

    else:
        # '-' strand: promoter lies downstream of the ORF in genomic coordinates.
        # prom_start0 is the first base immediately upstream of the ATG in gene coordinates.
        # In plus-strand coordinates, that is the first base after the ORF: index gene.end.
        prom_start0 = min(chrom_len, gene.end)
        prom_end0 = min(chrom_len, prom_start0 + PROMOTER_SPAN)
        if GENES_BY_CHROM and not force:
            prom_end0 = adjust_promoter_end_for_neighbors_minus(gene, prom_start0, prom_end0)

        prom_len = max(0, prom_end0 - prom_start0)
        required = tune_len + arm_len
        if prom_len < required:
            if use_available:
                tune_len = prom_len - arm_len
                if tune_len <= 0:
                    raise ValueError(
                        f"TuneYALI cannot be designed for {gene.gene_id}: safe promoter length "
                        f"is {prom_len} bp, which is insufficient even for the UP arm."
                    )
                print(
                    f"[TuneYALI] Safe promoter downstream of {gene.gene_id} is {prom_len} bp; "
                    f"shrinking tune_len to {tune_len} bp (--use_available)."
                )
            else:
                raise ValueError(
                    f"TuneYALI cannot be designed safely for {gene.gene_id}: safe promoter length "
                    f"is {prom_len} bp but needs at least {required} bp for tune_len={tune_len}. "
                    f"Use --use_available to shrink automatically, or --force to ignore neighbour constraints."
                )

        # UP arm is distal (further from ATG) -> starts at prom_start0 + tune_len
        up_start_gen = prom_start0 + tune_len
        up_end_gen = min(chrom_len, up_start_gen + arm_len)
        if up_end_gen - up_start_gen != arm_len:
            raise ValueError(
                f"TuneYALI cannot be designed for {gene.gene_id}: UP arm would extend beyond contig end."
            )
        up_seq_gen = chrom_seq[up_start_gen:up_end_gen]
        up_seq = revcomp(up_seq_gen)

        # DOWN arm is 0..+162 of ORF in gene orientation (i.e. last 162 bp of ORF in genomic plus coords)
        down_end_gen = gene.end  # 0-based slice end immediately after the ORF
        down_start_gen = max(0, down_end_gen - arm_len)
        down_seq_gen = chrom_seq[down_start_gen:down_end_gen]
        if len(down_seq_gen) != arm_len:
            raise ValueError(
                f"TuneYALI cannot be designed for {gene.gene_id}: ORF is too close to contig start for DOWN arm."
            )
        down_seq = revcomp(down_seq_gen)

        # Deleted/replaced promoter segment is proximal: [prom_start0 .. prom_start0+tune_len)
        del_start = prom_start0
        del_end = prom_start0 + tune_len

        hdr_left = down_start_gen
        hdr_right = up_end_gen

    deleted_len = max(0, del_end - del_start)

    # -----------------------------
    # sgRNA selection in the promoter segment that will be replaced
    # -----------------------------
    prom_window = 200
    window_half = prom_window // 2

    if deleted_len <= 0:
        raise ValueError(f"Deleted promoter segment has non-positive length for {gene.gene_id}.")

    if gene.strand == "+":
        # UP edge of deleted segment is del_start (more upstream)
        target_center = del_start + deleted_len // 3
    else:
        # For '-', promoter is to the right; UP edge corresponds to del_end-1
        target_center = (del_end - 1) - deleted_len // 3

    target_center = max(del_start, min(del_end - 1, target_center))
    search_start = max(del_start, target_center - window_half)
    search_end = min(del_end, target_center + window_half)
    search_window = chrom_seq[search_start:search_end]

    sgrna = None
    guide_seq, robot_score = get_chopchop_guide(search_window, gene.gene_id + "_tune")

    if guide_seq:
        idx_fwd = chrom_seq.find(guide_seq, search_start, search_end)
        if idx_fwd != -1 and chrom_seq[idx_fwd + 20: idx_fwd + 23].endswith("GG"):
            cut_pos = idx_fwd + 18
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
            rc = revcomp(guide_seq)
            idx_rev = chrom_seq.find(rc, search_start, search_end)
            if (
                idx_rev != -1
                and idx_rev >= 3
                and chrom_seq[idx_rev - 3: idx_rev].startswith("CC")
            ):
                cut_pos = idx_rev + 4
                sgrna = SgRNA(
                    guide_seq,
                    "NGG",
                    gene.chrom,
                    cut_pos,
                    "-",
                    robot_score,
                    "CHOPCHOP-Web",
                )

    if sgrna is None:
        print("    [TuneYALI] CHOPCHOP promoter search failed; falling back to ORF-based guide search.")
        sgrna, fallback_window = find_sgrna_in_orf(genome, gene)
        if not sgrna:
            raise ValueError(f"No suitable sgRNA found for {gene.gene_id}")
        if not search_window:
            search_window = fallback_window

    # -----------------------------
    # Primer design outside the HDR span
    # -----------------------------
    search_window_len = 200
    buffer = 10

    p_left_end = max(0, hdr_left - buffer)
    p_left_start = max(0, p_left_end - search_window_len)
    left_search_seq = chrom_seq[p_left_start:p_left_end]

    p_right_start = min(chrom_len, hdr_right + buffer)
    p_right_end = min(chrom_len, p_right_start + search_window_len)
    right_search_seq = chrom_seq[p_right_start:p_right_end]

    fw_seq, fw_rel_pos = find_best_primer_with_pos(left_search_seq, 'fw')
    rv_seq, rv_rel_pos = find_best_primer_with_pos(right_search_seq, 'rv')

    fw_genomic_start = p_left_start + fw_rel_pos
    rv_genomic_end = p_right_start + rv_rel_pos + len(rv_seq)

    band_unedited = rv_genomic_end - fw_genomic_start
    # Edited size depends on what is inserted via Golden Gate. Here we report the
    # "empty" placeholder insertion case (GG_INSERT only).
    band_edited = band_unedited - deleted_len + len(GG_INSERT)

    fragment_seq = (
        GA_LEFT
        + sgrna.protospacer.upper()
        + TRACR_SCAFFOLD
        + TER_RPR1
        + up_seq
        + GG_INSERT
        + down_seq
        + GA_RIGHT
    )

    primers = [
        PrimerInfo("TUNE_VERIFY_FWD", fw_seq, fw_genomic_start, len(fw_seq)),
        PrimerInfo("TUNE_VERIFY_REV", rv_seq, rv_genomic_end, len(rv_seq)),
    ]

    # For reporting
    design = KnockoutDesign(
        gene=gene,
        sgrna=sgrna,
        upstream_arm=up_seq,
        downstream_arm=down_seq,
        fragment_seq=fragment_seq,
        primers=primers,
        band_size_unedited=band_unedited,
        band_size_ko=band_edited,
        search_window_seq=search_window,
    )
    return design

# =========================
# Output Generators
# =========================

def write_design_report(design: KnockoutDesign, out_path: str) -> None:
    with open(out_path, "w") as f:
        sg = design.sgrna
        f.write(f"=== TuneYALI Promoter-Insertion Cassette: {design.gene.gene_id} ===\n\n")
        f.write(f"[Target Details]\n")
        f.write(f"Gene ID: {design.gene.gene_id}\n")
        f.write(f"Method: {sg.method} (Score: {sg.score})\n")
        f.write(f"Cut Pos: {sg.cut_pos} ({sg.strand})\n")
        f.write(f"Protospacer: {sg.protospacer}\n\n")

        fwd_oligo = "TTCGATTCCGGGTCGGCGCAGGTTG" + sg.protospacer + "GTTTTA"
        rev_oligo = "GCTCTAAAAC" + revcomp(sg.protospacer) + "CAACCTGCGCCGACCCGGAAT"
        f.write(f"[gRNA Cloning Oligos]\nFwd: {fwd_oligo}\nRev: {rev_oligo}\n\n")

        f.write(f"[CHOPCHOP 200 bp window (~first third of DELETED promoter)]\n")
        f.write("Paste this sequence into CHOPCHOP (FASTA input) if you want to verify the guides manually.\n")
        f.write(f"{design.search_window_seq}\n\n")

        f.write(f"[Donor Fragment]\n>TUNE_{design.gene.gene_id}_Donor\n{design.fragment_seq}\n\n")

        f.write("[Verification Primers]\n")
        f.write("Primers are outside homology arms and detect a size shift corresponding to the deleted 5′ segment.\n")
        for p in design.primers:
            f.write(f"{p.name}: {p.seq} (Len: {p.length})\n")

        f.write(f"\n[Expected PCR Bands]\n")
        f.write(f"Unedited/Parental: ~{design.band_size_unedited} bp\n")
        f.write(
            f"Edited (placeholder only): ~{design.band_size_ko} bp "
            f"(delta: {design.band_size_ko - design.band_size_unedited} bp)\n"
        )
        f.write(
            "NOTE: In TuneYALI, the final edited band size depends on the length of the promoter "
            "inserted by Golden Gate. The value above assumes only the GG placeholder is present.\n"
        )


def write_genbank(design: KnockoutDesign, gb_path: str) -> None:
    record = SeqRecord(
        Seq(design.fragment_seq),
        id=f"TUNE_{design.gene.gene_id}_donor",
        name=f"TUNE_{design.gene.gene_id}",
        description=f"TuneYALI donor for promoter insertion of {design.gene.gene_id}. Expected bands: unedited~{design.band_size_unedited}bp, edited~{design.band_size_ko}bp."
    )
    today = datetime.today().strftime("%d-%b-%Y").upper()
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "linear"
    record.annotations["data_file_division"] = "UNC"
    record.annotations["date"] = today
    record.annotations["accessions"] = [f"TUNE_{design.gene.gene_id}_DONOR"]
    record.annotations["source"] = "synthetic DNA construct"
    record.annotations["organism"] = "synthetic construct"
    record.annotations["comment"] = "Designed with TuneYALI."

    spacer_len = len(design.sgrna.protospacer)
    modules = [
        (len(GA_LEFT), "misc_feature", "GA_left", "#808080"),
        (spacer_len, "misc_feature", f"sgRNA_{design.gene.gene_id}", "#FF0000"),
        (len(TRACR_SCAFFOLD), "misc_feature", "tracrRNA_scaffold", "#0000FF"),
        (len(TER_RPR1), "terminator", "TerRPR1", "#800080"),
        (len(design.upstream_arm), "misc_recomb", f"UP_{design.gene.gene_id}", "#008000"),
        (len(GG_INSERT), "misc_feature", "GG_INSERT", "#FFA500"),
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
                "note": [f"TuneYALI Module: {label}"]
            },
        )
        record.features.append(feature)
        pos += length

    SeqIO.write(record, gb_path, "genbank")


def write_png_map(design: KnockoutDesign, png_path: str) -> None:
    if not HAS_DNA_VIEWER: return

    L_ga_left = len(GA_LEFT)
    L_spacer = 20
    L_tracr = len(TRACR_SCAFFOLD)
    L_term = len(TER_RPR1)
    L_up = len(design.upstream_arm)
    L_insert = len(GG_INSERT)
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
        (L_insert, "#FFA500", "GG_INSERT"),
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

    wt_band = Dseqrecord("A" * design.band_size_unedited)
    wt_band.name = "Unedited"
    
    ko_band = Dseqrecord("A" * design.band_size_ko)
    ko_band.name = "KO"

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
    global BUFFER_SAME_STRAND, BUFFER_OPP_STRAND
    parser = argparse.ArgumentParser(description="TuneYALI Designer (Golden-Gate promoter insertion cassette)")
    parser.add_argument("--gene_id", required=True, help="Gene ID (e.g. YALI1_B12345g)")
    parser.add_argument("--genome", default="W29_genome.fasta", help="Genome FASTA")
    parser.add_argument("--gff", default="W29_annotation.gff", help="Annotation GFF")
    parser.add_argument("--out", default=None, help="Output Base Name")
    parser.add_argument(
        "--tune_len",
        type=int,
        default=500,
        help="Length (bp) of proximal promoter to REPLACE upstream of the ATG (default: 500).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Ignore neighbour-aware safety shrinking (still bounded by contig edges).",
    )
    parser.add_argument(
        "--use_available",
        action="store_true",
        help="If not enough safe promoter exists for tune_len, shrink tune_len to the maximum allowed.",
    )
    parser.add_argument(
        "--buffer_same_strand",
        type=int,
        default=BUFFER_SAME_STRAND,
        help="Safety buffer (bp) from a same-strand neighbor (default: 100).",
    )
    parser.add_argument(
        "--buffer_opp_strand",
        type=int,
        default=BUFFER_OPP_STRAND,
        help="Safety buffer (bp) from an opposite-strand neighbor (default: 500).",
    )
    args = parser.parse_args()
    # Allow overriding neighbor-safety buffers from CLI
    if args.buffer_same_strand < 0 or args.buffer_opp_strand < 0:
        raise ValueError("buffer values must be non-negative.")
    BUFFER_SAME_STRAND = args.buffer_same_strand
    BUFFER_OPP_STRAND = args.buffer_opp_strand

    ensure_reference_files(args.genome, args.gff)
    genome = load_genome_fasta(args.genome)
    genes = load_annotation_gff(args.gff)

    if args.gene_id not in genes:
        raise KeyError(f"Gene {args.gene_id} not found in GFF.")

    # Build global indices for neighbour-aware truncation
    global ALL_GENES, GENES_BY_CHROM
    ALL_GENES = genes
    GENES_BY_CHROM = index_genes_by_chrom(genes)

    print(f"Designing TuneYALI promoter insertion cassette for {args.gene_id}...")
    design = design_tuneyali(
        genome,
        genes[args.gene_id],
        tune_len=args.tune_len,
        force=args.force,
        use_available=args.use_available,
    )

    base = args.out or f"{args.gene_id}_tune_cassette"
    if base.endswith(".txt"): base = base[:-4]

    write_design_report(design, f"{base}.txt")
    write_genbank(design, f"{base}.gb")
    write_png_map(design, f"{base}.png")
    write_gel_simulation(design, f"{base}_gel.png")

    print(f"Design Complete (TuneYALI) for {args.gene_id}")
    print(f" - Method:             {design.sgrna.method}")
    print(f" - Score:              {design.sgrna.score}/100")
    print(f" - Unedited Band:      ~{design.band_size_unedited} bp")
    print(f" - Edited Band (placeholder):    ~{design.band_size_ko} bp")
    print(f" - Report:             {base}.txt")

    if HAS_PYDNA:
        print(f" - Gel Check:          {base}_gel.png")

if __name__ == "__main__":
    main()
