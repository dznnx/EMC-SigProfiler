# EMC-SigProfiler (`emcsp`)

EMC-SigProfiler is a lightweight Python tool and library built on top of the **SigProfiler** suite (`SigProfilerMatrixGenerator`, `SigProfilerAssignment`, `SigProfilerPlotting`) to perform mutational signature analysis directly from simple genomic variant TSV files or VCFs. It formats output results into clean TSVs and highly customizable, fully offline-ready HTML reports.

<img width="683" height="913" alt="image" src="https://github.com/user-attachments/assets/7d0ef86a-447f-42ef-93cf-9295727b0549" />


## Installation

### Prerequisites
- **Python**: `>= 3.12`
- **System dependencies**: To embed profile plots in the HTML reports, you need `pdftoppm` (part of the `poppler-utils` package on Ubuntu/Debian) to convert PDF plots to images:
  ```bash
  # Debian/Ubuntu
  sudo apt-get install poppler-utils

  # macOS
  brew install poppler
  ```

### Install from source
Clone the repository and install using `uv` (recommended) or `pip`:
```bash
# Using uv
uv pip install .

# Or using pip
pip install .
```

---

## Input File Format

The input must be a tab-separated values (TSV) file with at least the following columns (`CHROM`, `POS`, `REF`, `ALT`):

```tsv
CHROM	POS	REF	ALT
chr1	123456	A	G
chr2	789012	C	T
chrX	456789	G	A
```

---

## Command Line Interface (CLI) Usage

The tool installs a command line utility named `emcsp`.

### 1. Install Reference Genomes
You must install a reference genome (e.g., `GRCh37` or `GRCh38`) before running signature analysis.
```bash
emcsp install --reference GRCh37
```

**Options:**

- `-r, --reference TEXT`: Reference genome (e.g. `GRCh37` or `GRCh38`) [default: GRCh37]
- `--rsync/--no-rsync`: Use rsync for download [default: rsync]
- `--bash/--no-bash`: Use bash for archive extraction [default: bash]

### 2. Analyze Mutational Signatures
Run signature extraction on an input TSV file.
```bash
emcsp analyze input.tsv --sample-name SampleA --output ./results --format both
```
**Options:**

- `-s, --sample-name TEXT`: Sample name for the analysis [required]
- `-o, --output PATH`: Output directory [default: current working directory]
- `-r, --reference TEXT`: Reference genome (`GRCh37` or `GRCh38`) [default: GRCh37]
- `-f, --format [tsv|html|both]`: Output format [default: tsv]
- `-c, --context-type TEXT`: Mutation context type (e.g., `96`, `192`, `1536`) [default: 96]
- `--plot/--no-plot`: Generate plots for mutational matrices
- `--make-plots/--no-make-plots`: Generate signature assignment plots
- `--cosmic-version FLOAT`: COSMIC signatures version [default: 3.4]
- `--exome/--no-exome`: Enable exome-only analysis [default: no-exome]

---

## Python Library Usage

EMC-SigProfiler can also be imported directly into your Python scripts.

### Example
```python
from pathlib import Path
from emcSP import SampleConfig, analyze_from_tsv

# 1. Define analysis configuration
cfg = SampleConfig(
    sample_name="MySample",
    output_dir=Path("./results"),
    reference="GRCh37",
    output_format="both",       # options: 'tsv', 'html', 'both'
    context_type="96",
    cosmic_version=3.4,
    make_plots=False,
)

# 2. Run analysis
analyze_from_tsv(
    input_tsv=Path("input.tsv"),
    cfg=cfg,
    logger=print  # Accepts any callable logger
)
```

---

## Outputs

Analysis generates the following files in the specified output directory:

1. **TSV Report (`<sample-name>_signatures.tsv`)**:
   Contains mutational signatures sorted by their relative contribution:
   ```tsv
   Signature	Relative Contribution	Mutation Count	Etiology
   SBS5		0.667			20		Unknown (Clock-like)
   SBS1		0.333			10		Deamination of 5-methylcytosine
   ```

2. **HTML Report (`<sample-name>_signatures.html`)**:
   An offline-ready interactive report. Features:
   - Base64-encoded original mutational profile and signature decomposition (reconstruction) plots (requires `pdftoppm` in your environment) embedded directly in the file so the HTML is self-contained.
   - Clinical reconstruction metrics table including **Cosine Similarity**, **Pearson Correlation**, **KL Divergence**, and **Total Mutational Burden** comparison.

The metrics in the HTML output use the following thresholds for 'success', 'warning' or 'error' indications:
- cosine similarity: `success >=0.9`, `warning >=0.8`, `error <0.8`
- correlation: `success >=0.9`, `warning >=0.8`, `error <0.8`
- kl divergence: `success <=0.1`, `warning <=0.2`, `error >0.2`
- total mutations: `low <50`, `moderate 50-200`, `high 200-500`, `very high >500`
