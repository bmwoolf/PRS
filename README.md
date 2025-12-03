# PRS Calculation Pipeline

This pipeline calculates polygenic risk scores (PRS) for 9 diseases using your whole genome sequencing data and the Nucleus Origin reference.

## Setup

### Prerequisites

Install `uv` (fast Python package installer):
```bash
# macOS/Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Or via pip
pip install uv

# Or via homebrew
brew install uv
```

### Installation

**Option 1: Using uv sync (recommended)**
```bash
# Create virtual environment and install all dependencies from pyproject.toml
uv sync

# Activate virtual environment
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

**Option 2: Manual setup**
```bash
# Create virtual environment
uv venv

# Activate virtual environment
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
uv pip install -r requirements.txt
```

## Quick Start

```bash
# Update DEFAULT_VCF_PATH in calculate_prs.py or use --vcf flag

# Test with small sample first (recommended)
python calculate_prs.py --test

# Run full calculation (will take time if VCF lacks rsIDs)
python calculate_prs.py

# Or specify custom file paths
python calculate_prs.py --nucleus path/to/nucleus.tsv --vcf path/to/genome.vcf.gz
```

## Input Files

1. **Nucleus Reference**: `nucleus/NUCLEUS_ORIGIN_V1/NUCLEUS_ORIGIN_V1.tsv`
   - Contains rsIDs, effect alleles (A1), and effect sizes for 9 diseases

2. **Your Genome VCF**: `your_genome.vcf.gz` (or similar)
   - Your personal whole genome sequencing data in VCF format

## Matching SNPs

The script matches SNPs between the reference and your genome in two ways:

1. **By rsID** (preferred): If your VCF has rsIDs in the ID column
2. **By position** (fallback): If rsIDs are missing, the script fetches chromosome/position mappings from Ensembl API

**Important Notes**:
- If your VCF doesn't have rsIDs, position-based matching for 2M+ SNPs via API can take **several hours** on the first run
- The script caches position mappings in `rsid_positions_cache.json` - subsequent runs will be much faster
- For faster results, consider annotating your VCF with rsIDs first (see below)
- You can test with `--test` flag to process only 1000 SNPs first

### Option: Annotate VCF with rsIDs (Recommended)

If your VCF lacks rsIDs, you can annotate it using bcftools:

```bash
# Download dbSNP VCF (one-time setup)
# You'll need a dbSNP VCF file, e.g., from:
# https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151/VCF/

# Annotate your VCF (example - adjust paths as needed)
bcftools annotate \
  -a dbsnp.vcf.gz \
  -c ID \
  -O z \
  -o annotated.vcf.gz \
  your_genome.vcf.gz

# Then index it
bcftools index annotated.vcf.gz

# Update the script to use annotated.vcf.gz
```

## Output Files

- `prs_results.csv`: PRS scores for all 9 diseases (see `prs_results.example.csv` for format)
- `prs_summary.txt`: Summary statistics and SNP counts (see `prs_summary.example.txt` for format)

**Note**: Actual result files are excluded from git via `.gitignore` to protect PII. Example files are provided as templates.

## Diseases Analyzed

1. ALZHEIMERS
2. BREASTCANCER
3. CORONARYARTERYDISEASE
4. ENDOMETRIOSIS
5. HYPERTENSION
6. PROSTATECANCER
7. RHEUMATOIDARTHRITIS
8. TYPE1DIABETES
9. TYPE2DIABETES

## Formula

For each disease D:
```
PRS_D = sum(beta_{i,D} * G_i)
```

where:
- `beta_{i,D}` = effect size from Nucleus file for SNP i and disease D
- `G_i` = genotype dosage (0, 1, or 2) aligned to A1 allele
