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

### Step 1: Create Mapping File (Required if VCF lacks rsIDs)

First, download a dbSNP VCF file (one-time setup):

```bash
# Download dbSNP VCF from NCBI (choose appropriate build, e.g., GRCh38)
# Example for common variants:
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151/VCF/00-common_all.vcf.gz

# Or download chromosome-specific files for smaller downloads:
# wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151/VCF/chr1.vcf.gz
# ... (repeat for other chromosomes)
```

Create the mapping file:

```bash
# Activate environment
source .venv/bin/activate

# Create mapping file (this may take 10-30 minutes depending on dbSNP size)
python scripts/create_mapping_bcftools.py \
  --nucleus <PATH_TO_NUCLEUS_TSV> \
  --dbsnp <PATH_TO_DBSNP_VCF> \
  --output data/mappings/rsid_mapping.tsv

# For testing with limited rsIDs:
python scripts/create_mapping_bcftools.py \
  --nucleus <PATH_TO_NUCLEUS_TSV> \
  --dbsnp <PATH_TO_DBSNP_VCF> \
  --output data/mappings/test_mapping.tsv \
  --limit 1000
```

### Step 2: Calculate PRS

```bash
# Test with small sample first (recommended)
python scripts/calculate_prs.py \
  --nucleus <PATH_TO_NUCLEUS_TSV> \
  --vcf <PATH_TO_YOUR_VCF> \
  --mapping data/mappings/rsid_mapping.tsv \
  --output-dir data/results \
  --test

# Run full calculation
python scripts/calculate_prs.py \
  --nucleus <PATH_TO_NUCLEUS_TSV> \
  --vcf <PATH_TO_YOUR_VCF> \
  --mapping data/mappings/rsid_mapping.tsv \
  --output-dir data/results
```

## Input Files

1. **Nucleus Reference**: `nucleus/NUCLEUS_ORIGIN_V1/NUCLEUS_ORIGIN_V1.tsv`
   - Contains rsIDs, effect alleles (A1), and effect sizes for 9 diseases

2. **Your Genome VCF**: `your_genome.vcf.gz` (or similar)
   - Your personal whole genome sequencing data in VCF format

## Matching SNPs

The script matches SNPs between the reference and your genome in two ways:

1. **By rsID** (preferred): If your VCF has rsIDs in the ID column
2. **By position** (required if VCF lacks rsIDs): Uses a mapping file to match by coordinates

**Important Notes**:
- **Mapping file required**: You must provide a mapping file with `--mapping` option for position-based matching
  - Format: TSV with columns: rsID, chr, pos, ref, alt (header optional)
  - Create using: `python scripts/create_mapping_bcftools.py --nucleus <PATH_TO_NUCLEUS_TSV> --dbsnp <PATH_TO_DBSNP_VCF> --output data/mappings/rsid_mapping.tsv`
- Not all SNPs in the Nucleus reference will be present in your VCF - this is normal
- You can test with `--test` flag to process only 1000 SNPs first

### Creating a Mapping File

The mapping file maps rsIDs to chromosome positions. Create it once and reuse:

**Requirements:**
- bcftools installed (`conda install -c bioconda bcftools`)
- dbSNP VCF file (download from NCBI)

**Steps:**

1. **Download dbSNP VCF** (if not already done):
   ```bash
   # Common variants (recommended, ~2-3 GB):
   wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151/VCF/00-common_all.vcf.gz
   
   # Or chromosome-specific (smaller, combine later):
   wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151/VCF/chr1.vcf.gz
   # ... download other chromosomes as needed
   ```

2. **Create mapping file**:
   ```bash
   python scripts/create_mapping_bcftools.py \
     --nucleus <PATH_TO_NUCLEUS_TSV> \
     --dbsnp <PATH_TO_DBSNP_VCF> \
     --output data/mappings/rsid_mapping.tsv
   ```
   
   This will:
   - Extract all rsIDs from the Nucleus TSV
   - Query the dbSNP VCF for their positions
   - Create a TSV mapping file: `rsID, chr, pos, ref, alt`
   - Take 10-30 minutes depending on dbSNP file size

3. **Use the mapping file** with `calculate_prs.py` (see Quick Start above)

## Output Files

Results are saved to `data/results/`:
- `prs_results.csv`: PRS scores for all 9 diseases (see `prs_results.example.csv` for format)
- `prs_snp_level.csv`: SNP-level details used in calculation
- `prs_summary.txt`: Summary statistics and SNP counts (see `prs_summary.example.txt` for format)

**Note**: Actual result files are excluded from git via `.gitignore` to protect PII. Example files are provided as templates.

## Project Structure

```
PRS/
├── scripts/              # Python scripts
│   ├── calculate_prs.py
│   └── create_mapping_bcftools.py
├── data/
│   ├── mappings/         # rsID mapping files (PII)
│   └── results/          # PRS results (PII)
├── logs/                 # Log files
├── README.md
└── pyproject.toml
```

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
