# AutoTuneYALI

Automated design of CRISPR/HDR edits for *Yarrowia lipolytica* (W29/CLIB89).  
Given a **YALI1_ gene ID**, AutoTuneYALI:

- picks a sgRNA (CHOPCHOP or local scoring),
- builds a TuneYALI-style donor cassette with short homology arms, and
- designs outer–outer primers for colony PCR (WT vs edited band).

This repository: <https://github.com/manufer20/AutoTuneYALI>

Current modules:

- **`koyali.py` (KOYALI)** – clean ORF deletions (start→stop codon)  
- **`truncyali.py` (TruncYALI)** – 5′ promoter truncations (e.g. keep last 50 bp)  
- **`tuneyali.py` (TuneYALI)** – *planned* promoter-tuning libraries (not yet stable)

> Defaults target **W29/CLIB89 (RefSeq GCF_001761485.1)** and are intended for
> Δku70 / NHEJ-deficient strains with genomic Cas9. Any FASTA+GFF pair can be used
> if it contains standard `YALI1_…` gene IDs.

---

## Installation

Create an environment (conda or venv) and install the dependencies:

```bash
python -m venv .venv
source .venv/bin/activate  # or .venv\Scripts\activate on Windows

pip install biopython requests matplotlib selenium webdriver-manager \
            dna-features-viewer pydna
