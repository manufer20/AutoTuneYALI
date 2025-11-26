# AutoTuneYALI

Automated design of CRISPR/Cas9 + HDR editing cassettes for *Yarrowia lipolytica* W29/CLIB89.
**Reference genome/annotation:** RefSeq GCF_001761485.1

**AutoTuneYALI generates (per target locus):**
* A TuneYALI-style cassette sequence (GA_LEFT → gRNA → scaffold → TerRPR1 → arms → GA_RIGHT)
* Verification primers (out–out PCR) + expected band sizes
* Optional cassette-map PNG (`dna-features-viewer`)
* Optional simulated gel PNG (`pydna`)

### Modules
* **KOYALI:** Clean, markerless ORF deletion (start→stop) using short HDR arms.
* **TruncYALI:** Promoter truncation (neighbor-aware), keeping a proximal promoter remnant (default: 500 bp) while deleting distal promoter sequence.
* **TuneYALI:** Promoter-swap placeholder (Golden Gate): UP arm (-662..-500), insert `GACTGAAGAGCTCTTCA`, DOWN arm (0..+162).

**Assumption:** Editing in a *Yarrowia* chassis with genomic Cas9 (often integrated at KU70) and a compatible CRISPR backbone using the GA_LEFT/GA_RIGHT overlaps.

---

## Outputs

Each module typically writes:
* `<base>.txt` — guide, cut site, arms, primers, expected PCR bands, CHOPCHOP window (for manual verification)
* `<base>.gb` — GenBank of the cassette (annotated features)
* `<base>.png` — cassette map (optional; requires `dna-features-viewer`)
* `<base>_gel.png` — simulated gel (optional; requires `pydna`)

### Reference files (W29)

By default, scripts will look for these in the current directory:
* `W29_genome.fasta`
* `W29_annotation.gff`

If they are missing, the scripts attempt an automatic download from NCBI Datasets.

---

## Installation

Everything below installs dependencies using conda (no `pip install ...`) and **you need to have Chrome installed in your computer**.

### Windows (fresh install)

**1) Install Miniforge**
Install Miniforge3 for Windows (x64). During installation, allow it to initialize conda.

**2) Install Git for Windows**
Option A (recommended; winget):
```powershell
winget install --id Git.Git -e
```
Option B: install “Git for Windows” using the official installer.

**3) Clone the repo**
```powershell
cd $HOME\Desktop
git clone [https://github.com/manufer20/AutoTuneYALI.git](https://github.com/manufer20/AutoTuneYALI.git)
cd AutoTuneYALI
```

**4) Create environment**
```powershell
conda create -n autotuneyali -y python=3.11
conda activate autotuneyali
```

**5) Install dependencies (choose minimal vs complete)**

*Windows minimal (recommended): no pydna and no dna-features-viewer*
```powershell
conda install -y -c conda-forge biopython requests pandas matplotlib selenium webdriver-manager
```

*Optional extras (only if you want cassette map / gel simulation outputs):*
```powershell
conda install -y -c conda-forge dna-features-viewer
conda install -y -c conda-forge pydna
```

**6) Quick sanity check**
```powershell
python -c "import Bio,requests,pandas,matplotlib; print('OK')"
```

### macOS (fresh install)

**1) Install Miniforge**
Install Miniforge3 for macOS (Apple Silicon or Intel).

**2) Git (usually already present)**
If git is missing:
```bash
xcode-select --install
```

**3) Clone the repo**
```bash
cd ~/Desktop
git clone [https://github.com/manufer20/AutoTuneYALI.git](https://github.com/manufer20/AutoTuneYALI.git)
cd AutoTuneYALI
```

**4) Create environment**
```bash
conda create -n autotuneyali -y python=3.11
conda activate autotuneyali
```

**5) Install dependencies (complete)**
```bash
conda install -y -c conda-forge biopython requests pandas matplotlib selenium webdriver-manager dna-features-viewer pydna
```

**6) Quick sanity check**
```bash
python -c "import Bio,requests,pandas,matplotlib; print('OK')"
```

---

## Running the tools

**General notes:**
* Run from the repository folder (`AutoTuneYALI/`).
* If you don’t pass `--genome` / `--gff`, defaults are: `W29_genome.fasta` and `W29_annotation.gff`.
* On first run, the script may download the W29 references automatically.

### KOYALI (clean ORF deletion)

**Design intent:**
* Chooses an sgRNA in the 5′ third of the ORF (also avoids cuts too close to start codon when possible).
* Builds HDR arms:
    * UP: ~160 bp immediately upstream of the start codon.
    * DOWN: ~160 bp immediately downstream of the stop codon.
* Yields a markerless start→stop deletion after HDR.
* Designs out–out primers and predicts WT vs KO band sizes.

**Windows example (YALI1_A09949g):**
```powershell
cd $HOME\Desktop\AutoTuneYALI
conda activate autotuneyali
python .\koyali.py --gene_id YALI1_A09949g
```

**macOS example (YALI1_A09949g):**
```bash
cd ~/Desktop/AutoTuneYALI
conda activate autotuneyali
python ./koyali.py --gene_id YALI1_A09949g
```

**Optional: set output base name:**
```bash
python ./koyali.py --gene_id YALI1_A09949g --out YALI1_A09949g_KO
```

### TruncYALI (promoter truncation; neighbor-aware)

**Design intent:**
* Deletes distal promoter sequence upstream of the ATG while keeping a proximal remnant.
* Default keep (proximal promoter retained): 500 bp.
* Enforces conservative neighbor-aware safety rules so truncations do not encroach into adjacent annotated features.
* DOWN arm crosses into the ORF: (kept promoter X bp) + (162−X bp ORF).

**Basic usage:**
```bash
python ./truncyali.py --gene_id YALI1_A09949g
```

**Overrides (if supported):**
* `--prom_len` (default 500): minimum promoter to keep.
* `--term_len` (default 100): minimum space reserved near the neighbor/terminator side.

**Example:**
```bash
python ./truncyali.py --gene_id YALI1_A09949g --prom_len 300 --term_len 150
```

### TuneYALI (Golden Gate promoter insertion placeholder)

**Design intent:**
* Creates a cassette with:
    * UP arm: -662..-500 (162 bp segment)
    * INSERT: `GACTGAAGAGCTCTTCA`
    * DOWN arm: 0..+162 (first 162 bp of ORF)
* Uses the same neighbor-aware safety logic as TruncYALI.
* By default, if there is not enough safe intergenic space, the design aborts.

**Optional behaviors (if implemented):**
* `--use_available`: reduce the replaceable promoter length to the maximum allowed.
* `--force`: ignore neighbor buffers (still respects contig boundaries).

**Basic usage:**
```bash
python ./tuneyali.py --gene_id YALI1_A09949g
```

---

## Troubleshooting

**CHOPCHOP automation fails**
If Selenium automation fails (or is blocked), the scripts fall back to local guide scoring.

**Optional outputs not generated**
If you don’t install the optional packages:
* Cassette maps (`.png`) will be skipped without `dna-features-viewer`.
* Gel simulations (`_gel.png`) will be skipped without `pydna`.

**“Gene not found in GFF”**
Confirm the gene ID matches the W29 GFF locus tag format (e.g., `YALI1_XXXXXXg`) and you are using the correct annotation file.

---

## Citation

If you use this tool in academic work, cite the repository:

https://github.com/manufer20/AutoTuneYALI
