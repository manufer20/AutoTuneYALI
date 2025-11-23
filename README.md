# AutoTuneYALI

Automated design of CRISPR/HDR edits for *Yarrowia lipolytica*.
Given a gene (or coordinates), AutoTuneYALI proposes guides, builds HDR
donors, and outputs primers for assembly and colony-PCR validation.

- **KOYALI** — clean ORF deletions (start→stop) via CRISPR + HDR  
- **TruncYALI** — promoter swaps/truncations for expression attenuation  
- **TuneYALI** — *(planned)* small promoter-tuning libraries

> Defaults target **W29/CLIB89 (RefSeq GCF_001761485.1)** and work well in
> Δku70/NHEJ-deficient strains. Any FASTA/GenBank reference can be used.

---

## How it works

1. Parse the target ORF (± flanks) from the reference.
2. Propose CRISPR cut sites near the intended junctions.
3. Build HDR donors: `up-arm → scar/swap → down-arm` (typically 0.9–1.2 kb).
4. Design primers for donor assembly and colony-PCR (WT vs edited).

---

## What you get

- `guides.tsv` — candidate gRNAs (coord, PAM, strand)
- `donor.fasta` — HDR donor(s) for the edit
- `primers.tsv` — assembly + QC primer pairs with expected sizes
- `summary.json` — chosen sites, spans, parameters
- `design/` — per-target artifacts (maps/plots if emitted)

---

## Quickstart

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt  # or: pip install biopython pandas
# Obtain the W29/CLIB89 reference FASTA/GBK first.

# Clean ORF knockout
python koyali.py --gene DGA1 --genome W29.fasta \
  --left-arm 1000 --right-arm 1000 --scar TTAA \
  --out results/DGA1

# Promoter truncation/swap
python truncyali.py --gene ERG9 --genome W29.fasta \
  --new-promoter Ptef-mini.fasta --left-arm 900 --right-arm 900 \
  --out results/ERG9_promoter_swap

# Help
python koyali.py --help
python truncyali.py --help
