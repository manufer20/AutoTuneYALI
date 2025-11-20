#!/usr/bin/env python3

"""
TuneYALI-style cassette designer for Yarrowia lipolytica (W29).
UPDATED: 
- Calculates PCR Band sizes (Unedited vs. KO).
- Specific naming in GenBank features.
- "Out-Out" Verification Primers.
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
import matplotlib.pyplot as plt

# Optional: DNA Features Viewer
try:
    from dna_features_viewer import GraphicFeature, GraphicRecord
    HAS_DNA_VIEWER = True
except ImportError:
    HAS_DNA_VIEWER = False


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


# =========================
# Reference Files
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
            
            attrs = {k: v for k, v in [x.split("=", 1) for x in cols[8].split(";") if "=" in x]}
            gene_id = None
            for key in ("ID", "locus_tag", "gene", "Name"):
                if attrs.get(key, "").startswith("YALI1_"):
                    gene_id = attrs[key]
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

def calculate_heuristic_score(seq: str) -> int:
    seq = seq.upper()
    score = 50
    gc = gc_content(seq)
    
    if 0.4 <= gc <= 0.6: score += 20
    elif gc < 0.3 or gc > 0.8: score -= 30
    
    if seq[-1] == 'G': score += 15        
    elif seq[-1] == 'T': score -= 10      
    
    if "TTT" in seq: score -= 20          
    if gc_content(seq[-8:]) < 0.25: score -= 15
        
    return max(0, min(100, score))

def find_sgrna_in_orf(genome: Dict[str, str], gene: GeneInfo, window_size: int = 200) -> Optional[SgRNA]:
    chrom_seq = genome[gene.chrom]
    if gene.strand == '+':
        orf_seq = chrom_seq[gene.start - 1 : gene.end]
        offset = gene.start - 1
        strand = '+'
    else:
        orf_seq = revcomp(chrom_seq[gene.start - 1 : gene.end])
        offset = gene.end - 1
        strand = '-'

    if len(orf_seq) < 23: return None
    
    half = window_size // 2
    center = len(orf_seq) // 2
    start_idx = max(0, center - half)
    end_idx = min(len(orf_seq), start_idx + window_size)
    search_region = orf_seq[start_idx:end_idx]
    
    candidates = []
    for i in range(len(search_region) - 23):
        j = start_idx + i
        protospacer = orf_seq[j:j+20]
        pam = orf_seq[j+20:j+23]
        
        if pam[1:] != "GG": continue
        if "TTTT" in protospacer: continue
        
        score = calculate_heuristic_score(protospacer)
        
        if gene.strand == '+':
            cut_pos = offset + j + 20 + 1 - 3
        else:
            cut_pos = gene.start + (len(orf_seq) - (j + 20 + 3)) + 3
            
        candidates.append((score, SgRNA(protospacer, pam, gene.chrom, cut_pos, strand)))

    if not candidates: return None
    candidates.sort(key=lambda x: x[0], reverse=True)
    return candidates[0][1]

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

def design_knockout(genome: Dict[str, str], gene: GeneInfo, arm_len: int = 162) -> KnockoutDesign:
    """
    Designs cassette + "Out-Out" Primers + Calculates PCR Band Sizes.
    """
    sgrna = find_sgrna_in_orf(genome, gene)
    if not sgrna: raise ValueError(f"No suitable sgRNA found for {gene.gene_id}")
    
    up, down = extract_homology_arms(genome, gene, arm_len)
    chrom_seq = genome[gene.chrom]
    
    # Genomic Boundaries
    genomic_start_idx = gene.start - 1
    genomic_end_idx = gene.end
    gene_len = genomic_end_idx - genomic_start_idx

    # Search Windows (Outside Arms)
    search_window = 200
    buffer = 10 
    
    # Left Region (Upstream)
    p_left_end = max(0, genomic_start_idx - arm_len - buffer)
    p_left_start = max(0, p_left_end - search_window)
    left_search_seq = chrom_seq[p_left_start : p_left_end]
    
    # Right Region (Downstream)
    p_right_start = min(len(chrom_seq), genomic_end_idx + arm_len + buffer)
    p_right_end = min(len(chrom_seq), p_right_start + search_window)
    right_search_seq = chrom_seq[p_right_start : p_right_end]
    
    # Find Primers & Relative Positions
    fw_seq, fw_rel_pos = find_best_primer_with_pos(left_search_seq, 'fw')
    rv_seq, rv_rel_pos = find_best_primer_with_pos(right_search_seq, 'rv')

    # Calculate Genomic Coordinates of Primers
    # Fwd primer starts at p_left_start + fw_rel_pos
    fw_genomic_start = p_left_start + fw_rel_pos
    
    # Rev primer (reverse complement) starts at p_right_start + rv_rel_pos
    # But for PCR length, we need the 3' end of the reverse primer on the sense strand,
    # which corresponds to the "end" of the binding site.
    rv_genomic_end = p_right_start + rv_rel_pos + len(rv_seq)

    # PCR Math
    # Unedited = Distance between primers in genome
    band_unedited = rv_genomic_end - fw_genomic_start
    
    # KO = Unedited - Gene_Length (Since arms are fused perfectly)
    # Note: This assumes "clean" deletion.
    band_ko = band_unedited - gene_len

    fragment_seq = GA_LEFT + sgrna.protospacer.upper() + TRACR_SCAFFOLD + TER_RPR1 + up + down + GA_RIGHT
    
    primers = [
        PrimerInfo("KO_VERIFY_FWD", fw_seq, fw_genomic_start, len(fw_seq)),
        PrimerInfo("KO_VERIFY_REV", rv_seq, rv_genomic_end, len(rv_seq))
    ]
    
    return KnockoutDesign(gene, sgrna, up, down, fragment_seq, primers, band_unedited, band_ko)


# =========================
# Output Generators
# =========================

def write_design_report(design: KnockoutDesign, out_path: str) -> None:
    with open(out_path, "w") as f:
        sg = design.sgrna
        f.write(f"=== TuneYALI Knockout: {design.gene.gene_id} ===\n\n")
        f.write(f"[Target Details]\n")
        f.write(f"Gene ID: {design.gene.gene_id}\n")
        f.write(f"Cut Pos: {sg.cut_pos} ({sg.strand})\n")
        f.write(f"Protospacer: {sg.protospacer}\n\n")
        
        fwd_oligo = "TTCGATTCCGGGTCGGCGCAGGTTG" + sg.protospacer + "GTTTTA"
        rev_oligo = "GCTCTAAAAC" + revcomp(sg.protospacer) + "CAACCTGCGCCGACCCGGAAT"
        f.write(f"[gRNA Cloning Oligos]\nFwd: {fwd_oligo}\nRev: {rev_oligo}\n\n")
        
        f.write(f"[Donor Fragment]\n>KO_{design.gene.gene_id}_Donor\n{design.fragment_seq}\n\n")
        
        f.write("[Verification Primers (Out-Out)]\n")
        f.write("Primers are outside homology arms. Plasmid will show NO BAND.\n")
        for p in design.primers:
            f.write(f"{p.name}: {p.seq} (Len: {p.length})\n")
            
        f.write(f"\n[Expected PCR Bands]\n")
        f.write(f"Unedited/Parental: ~{design.band_size_unedited} bp\n")
        f.write(f"Successful KO:     ~{design.band_size_ko} bp\n")
        f.write(f"Donor Plasmid:     NO BAND (Primers do not bind)\n")


def write_genbank(design: KnockoutDesign, gb_path: str) -> None:
    record = SeqRecord(
        Seq(design.fragment_seq),
        id=f"KO_{design.gene.gene_id}_donor",
        name=f"KO_{design.gene.gene_id}",
        description=f"TuneYALI KO {design.gene.gene_id}. Check Primers: Unedited~{design.band_size_unedited}bp, KO~{design.band_size_ko}bp."
    )
    today = datetime.today().strftime("%d-%b-%Y").upper()   # e.g. 19-NOV-2025
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "linear"
    record.annotations["data_file_division"] = "UNC"        # “unclassified”
    record.annotations["date"] = today
    record.annotations["accessions"] = [f"KO_{design.gene.gene_id}_DONOR"]
    record.annotations["source"] = "synthetic DNA construct"
    record.annotations["organism"] = "synthetic construct"
    record.annotations["comment"] = (
        "Designed with KOYALI (TuneYALI-style knockout cassette designer "
        "for Yarrowia lipolytica W29)."
    )

    # Updated Features with Specific Names
    modules = [
        (len(GA_LEFT), "misc_feature", "GA_left", "#808080"),
        (20, "misc_feature", f"sgRNA_{design.gene.gene_id}", "#FF0000"),
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

def main():
    parser = argparse.ArgumentParser(description="TuneYALI Designer (Calculated Bands)")
    parser.add_argument("--gene_id", required=True, help="Gene ID (e.g. YALI1_B12345g)")
    parser.add_argument("--genome", default="W29_genome.fasta", help="Genome FASTA")
    parser.add_argument("--gff", default="W29_annotation.gff", help="Annotation GFF")
    parser.add_argument("--out", default=None, help="Output Base Name")
    args = parser.parse_args()

    ensure_reference_files(args.genome, args.gff)
    genome = load_genome_fasta(args.genome)
    genes = load_annotation_gff(args.gff)

    if args.gene_id not in genes:
        raise KeyError(f"Gene {args.gene_id} not found in GFF.")

    design = design_knockout(genome, genes[args.gene_id])
    
    base = args.out or f"{args.gene_id}_cassette"
    if base.endswith(".txt"): base = base[:-4]

    write_design_report(design, f"{base}.txt")
    write_genbank(design, f"{base}.gb")
    write_png_map(design, f"{base}.png")
    
    print(f"Design Complete for {args.gene_id}")
    print(f" - Unedited Band: ~{design.band_size_unedited} bp")
    print(f" - KO Band:       ~{design.band_size_ko} bp")
    print(f" - Report:        {base}.txt")

if __name__ == "__main__":
    main()