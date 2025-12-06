# AGENTS.md

Quick reference for AI agents working on this PRS calculation pipeline.

## Project Overview

Calculates polygenic risk scores (PRS) for 9 diseases using whole-genome VCF data and the Nucleus Origin reference dataset.

**Formula**: `PRS_D = Σ(β_i,D × G_i)` where β is effect size and G is genotype dosage (0, 1, or 2).

## Setup

1. **Install uv** (Python package manager):
   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```

2. **Install dependencies**:
   ```bash
   uv sync  # Creates .venv and installs from pyproject.toml
   source .venv/bin/activate
   ```

3. **External dependency**: `bcftools` (required for mapping creation)
   ```bash
   conda install -c bioconda bcftools
   ```

## Project Structure

```
PRS/
├── scripts/
│   ├── calculate_prs.py          # Main PRS calculation script
│   └── create_mapping_bcftools.py # Creates rsID→coordinate mapping
├── data/
│   ├── mappings/                  # rsID mapping files (generated)
│   └── results/                   # PRS output files (generated)
└── pyproject.toml                 # Dependencies: pysam, numpy
```

## Key Scripts

### `scripts/create_mapping_bcftools.py`
- **Purpose**: Creates rsID-to-coordinate mapping from dbSNP VCF
- **Inputs**: Nucleus TSV, dbSNP VCF file
- **Output**: TSV mapping file (rsID, chr, pos, ref, alt)
- **Why needed**: User's VCF may only have coordinates, not rsIDs

### `scripts/calculate_prs.py`
- **Purpose**: Calculates PRS scores for 9 diseases
- **Inputs**: Nucleus TSV, personal VCF, mapping file
- **Outputs**: `prs_results.csv`, `prs_snp_level.csv`, `prs_summary.txt`
- **Matching strategy**: 
  1. Try rsID matching first (if VCF has rsIDs)
  2. Fall back to position-based matching (requires mapping file)

## Key Technical Details

- **Genome build**: GRCh38 (both VCF and dbSNP should match)
- **Allele alignment**: Handles strand flips via reverse complement
- **Genotype format**: `pysam` returns GT as tuple `(0, 1)`, not string `"0/1"`
- **Mapping file format**: TSV with columns: `rsID, chr, pos, ref, alt` (header optional)
- **No Ensembl API**: All rsID mapping done via `bcftools` + dbSNP VCF

## Common Workflows

**First-time setup**:
```bash
# 1. Create mapping file (one-time, takes 10-30 min)
python scripts/create_mapping_bcftools.py \
  --nucleus <PATH_TO_NUCLEUS_TSV> \
  --dbsnp <PATH_TO_DBSNP_VCF> \
  --output data/mappings/rsid_mapping.tsv

# 2. Calculate PRS
python scripts/calculate_prs.py \
  --nucleus <PATH_TO_NUCLEUS_TSV> \
  --vcf <PATH_TO_YOUR_VCF> \
  --mapping data/mappings/rsid_mapping.tsv \
  --output-dir data/results
```

**Testing**:
```bash
# Test with first 1000 SNPs
python scripts/calculate_prs.py --test --nucleus ... --vcf ... --mapping ...
```

## Important Notes

- **PII files**: Results in `data/results/` and `data/mappings/` are gitignored
- **Mapping file required**: `--mapping` is required (no API fallback)
- **dbSNP download**: User must download dbSNP VCF separately (~15GB for 00-All.vcf.gz)
- **Memory efficient**: Mapping creation uses `bcftools view -i ID=@file` to filter before querying

## Diseases Analyzed

ALZHEIMERS, BREASTCANCER, CORONARYARTERYDISEASE, ENDOMETRIOSIS, HYPERTENSION, PROSTATECANCER, RHEUMATOIDARTHRITIS, TYPE1DIABETES, TYPE2DIABETES

