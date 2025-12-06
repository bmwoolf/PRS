#!/usr/bin/env python3
"""
PRS (Polygenic Risk Score) Calculation Pipeline

This script calculates polygenic risk scores for 9 diseases using whole genome
sequencing data and the Nucleus Origin reference.

Usage:
    python calculate_prs.py [--test] [--nucleus PATH] [--vcf PATH] [--mapping PATH] [--output-dir DIR]

Options:
    --test          Process only first 1000 SNPs (for testing)
    --nucleus PATH  Path to Nucleus reference TSV file
    --vcf PATH      Path to your genome VCF file
    --mapping PATH  Path to rsID-to-coordinate mapping file (optional, TSV format: rsID,chr,pos,ref,alt)
    --output-dir DIR Output directory for results (default: current directory)
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, NamedTuple

try:
    import pysam
    import numpy as np
except ImportError:
    print("Error: Required packages not installed.")
    print("Please install using uv:")
    print("  uv pip install pysam numpy")
    print("Or using pip:")
    print("  pip install pysam numpy")
    sys.exit(1)


# Default paths
DEFAULT_NUCLEUS_PATH = "nucleus/NUCLEUS_ORIGIN_V1/NUCLEUS_ORIGIN_V1.tsv"
DEFAULT_VCF_PATH = "your_genome.vcf.gz"

# Diseases analyzed
DISEASES = [
    "ALZHEIMERS",
    "BREASTCANCER",
    "CORONARYARTERYDISEASE",
    "ENDOMETRIOSIS",
    "HYPERTENSION",
    "PROSTATECANCER",
    "RHEUMATOIDARTHRITIS",
    "TYPE1DIABETES",
    "TYPE2DIABETES",
]


class NucleusSNP(NamedTuple):
    """Represents a SNP from the Nucleus reference."""
    rsid: str
    a1: str  # Effect allele
    effect_sizes: Dict[str, float]  # Disease -> effect size


class MatchedSNP(NamedTuple):
    """Represents a matched SNP with genotype information."""
    rsid: str
    chrom: str
    pos: int
    ref: str
    alt: str
    a1: str  # Effect allele from Nucleus
    genotype: int  # 0, 1, or 2 (number of A1 alleles)
    flipped: bool  # Whether alleles were flipped




def load_mapping_file(mapping_path: str) -> Dict[str, Dict]:
    """
    Load rsID-to-coordinate mapping from TSV file.
    
    Expected format: rsID,chr,pos,ref,alt (header optional)
    Returns dict: rsid -> {chrom, pos, ref, alt}
    """
    mapping = {}
    
    open_func = gzip.open if mapping_path.endswith('.gz') else open
    mode = 'rt' if mapping_path.endswith('.gz') else 'r'
    
    with open_func(mapping_path, mode) as f:
        reader = csv.reader(f, delimiter='\t')
        
        # Try to detect header
        first_row = next(reader, None)
        if first_row is None:
            return mapping
        
        # Check if first row is header
        if first_row[0].upper() in ['RSID', 'SNP', 'ID']:
            pass  # Skip header
        else:
            # First row is data, process it
            if len(first_row) >= 4:
                rsid = first_row[0].strip()
                chrom = first_row[1].strip()
                pos = int(first_row[2].strip())
                ref = first_row[3].strip() if len(first_row) > 3 else ""
                alt = first_row[4].strip() if len(first_row) > 4 else ""
                mapping[rsid] = {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt}
        
        # Process remaining rows
        for row in reader:
            if len(row) >= 4:
                rsid = row[0].strip()
                chrom = row[1].strip()
                pos = int(row[2].strip())
                ref = row[3].strip() if len(row) > 3 else ""
                alt = row[4].strip() if len(row) > 4 else ""
                mapping[rsid] = {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt}
    
    return mapping


def load_nucleus_reference(path: str) -> Dict[str, NucleusSNP]:
    """
    Load Nucleus reference file.
    
    Expected format: TSV with columns: SNP, A1, and effect sizes for each disease.
    Returns dict: rsid -> NucleusSNP
    """
    nucleus_data = {}
    
    open_func = gzip.open if path.endswith('.gz') else open
    mode = 'rt' if path.endswith('.gz') else 'r'
    
    with open_func(path, mode) as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            rsid = row['SNP'].strip()
            a1 = row['A1'].strip()
            
            effect_sizes = {}
            for disease in DISEASES:
                try:
                    effect_sizes[disease] = float(row[disease])
                except (ValueError, KeyError):
                    effect_sizes[disease] = 0.0
            
            nucleus_data[rsid] = NucleusSNP(rsid=rsid, a1=a1, effect_sizes=effect_sizes)
    
    return nucleus_data


def normalize_chrom(chrom: str) -> str:
    """Normalize chromosome name (e.g., '1' -> 'chr1', 'chr1' -> 'chr1')."""
    if chrom.startswith('chr'):
        return chrom
    return f"chr{chrom}"


def get_genotype_dosage(gt: str, ref: str, alt: str, a1: str) -> Tuple[Optional[int], bool]:
    """
    Calculate genotype dosage (0, 1, or 2) for A1 allele.
    
    Returns: (dosage, flipped) where flipped indicates if alleles were swapped.
    Handles phased (|) and unphased (/) genotypes.
    """
    # Parse genotype
    if '|' in gt:
        alleles = gt.split('|')
    elif '/' in gt:
        alleles = gt.split('/')
    else:
        return None, False
    
    # Handle missing data
    if '.' in alleles:
        return None, False
    
    try:
        allele_indices = [int(a) for a in alleles]
    except ValueError:
        return None, False
    
    # Get actual alleles: [REF, ALT1, ALT2, ...]
    alt_alleles = alt.split(',') if alt else []
    vcf_alleles = [ref] + alt_alleles
    
    # Check if A1 matches REF
    if a1 == ref:
        dosage = sum(1 for idx in allele_indices if idx == 0)
        return dosage, False
    
    # Check if A1 matches any ALT allele
    for idx in allele_indices:
        if idx > 0 and idx - 1 < len(alt_alleles):
            if alt_alleles[idx - 1] == a1:
                # Count how many times this ALT appears in genotype
                dosage = sum(1 for i in allele_indices if i > 0 and i - 1 < len(alt_alleles) and alt_alleles[i - 1] == a1)
                return dosage, False
    
    # Try reverse complement (flip) - handle strand differences
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    a1_flipped = ''.join(complement.get(b, b) for b in a1)
    ref_flipped = ''.join(complement.get(b, b) for b in ref)
    alt_flipped_list = [''.join(complement.get(b, b) for b in a) for a in alt_alleles]
    
    # Check if flipped A1 matches REF
    if a1_flipped == ref_flipped:
        dosage = sum(1 for idx in allele_indices if idx == 0)
        return dosage, True
    
    # Check if flipped A1 matches any ALT
    for idx in allele_indices:
        if idx > 0 and idx - 1 < len(alt_flipped_list):
            if alt_flipped_list[idx - 1] == a1_flipped:
                dosage = sum(1 for i in allele_indices if i > 0 and i - 1 < len(alt_flipped_list) and alt_flipped_list[i - 1] == a1_flipped)
                return dosage, True
    
    # No match found
    return None, False


def match_snps_by_rsid(vcf_path: str, nucleus_data: Dict[str, NucleusSNP]) -> List[MatchedSNP]:
    """
    Match SNPs between VCF and Nucleus reference by rsID.
    
    Returns list of MatchedSNP objects.
    """
    matched = []
    vcf = pysam.VariantFile(vcf_path)
    
    for record in vcf:
        # Check if ID column has rsID
        if not record.id or record.id == '.':
            continue
        
        rsids = record.id.split(';')
        for rsid in rsids:
            if rsid.startswith('rs') and rsid in nucleus_data:
                nucleus_snp = nucleus_data[rsid]
                
                # Get genotype
                sample = list(record.samples.values())[0]
                gt = sample.get('GT', None)
                if gt is None:
                    continue
                
                ref = record.ref
                alt = ','.join(record.alts) if record.alts else ''
                
                dosage, flipped = get_genotype_dosage(str(gt), ref, alt, nucleus_snp.a1)
                
                if dosage is not None:
                    chrom = normalize_chrom(record.chrom)
                    matched.append(MatchedSNP(
                        rsid=rsid,
                        chrom=chrom,
                        pos=record.pos,
                        ref=ref,
                        alt=alt,
                        a1=nucleus_snp.a1,
                        genotype=dosage,
                        flipped=flipped
                    ))
    
    vcf.close()
    return matched


def match_snps_by_position(
    vcf_path: str,
    nucleus_data: Dict[str, NucleusSNP],
    mapping: Dict[str, Dict],
    test_limit: Optional[int] = None
) -> Tuple[List[MatchedSNP], int, int]:
    """
    Match SNPs between VCF and Nucleus reference by chromosome/position.
    
    Requires mapping file with rsID -> coordinate mappings.
    Returns: (matched_snps, matched_count, dropped_count)
    """
    matched = []
    dropped = 0
    
    # Build position index from VCF
    print("Indexing VCF by position...")
    vcf = pysam.VariantFile(vcf_path)
    vcf_index = {}  # (chrom, pos) -> record
    
    for record in vcf:
        chrom = normalize_chrom(record.chrom)
        key = (chrom, record.pos)
        vcf_index[key] = record
    
    vcf.close()
    print(f"Indexed {len(vcf_index):,} variants from VCF")
    
    # Process Nucleus SNPs
    nucleus_items = list(nucleus_data.items())
    if test_limit:
        nucleus_items = nucleus_items[:test_limit]
    
    print(f"Processing {len(nucleus_items):,} SNPs from Nucleus reference...")
    
    if not mapping:
        print("Error: Mapping file is required for position-based matching.")
        return [], 0, 0
    
    for idx, (rsid, nucleus_snp) in enumerate(nucleus_items):
        if (idx + 1) % 10000 == 0:
            print(f"  Processed {idx + 1:,} / {len(nucleus_items):,} SNPs... (matched: {len(matched)}, dropped: {dropped})")
        
        # Get position from mapping file
        pos_info = mapping.get(rsid)
        
        if not pos_info:
            continue
        
        chrom = normalize_chrom(pos_info['chrom'])
        pos = int(pos_info['pos'])  # Ensure it's an int
        
        # Look up in VCF - try both with and without 'chr' prefix
        key = (chrom, pos)
        record = None
        
        if key in vcf_index:
            record = vcf_index[key]
        else:
            # Try alternative chromosome format
            if chrom.startswith('chr'):
                # Try without 'chr'
                key_alt = (chrom[3:], pos)
                if key_alt in vcf_index:
                    record = vcf_index[key_alt]
            else:
                # Try with 'chr'
                key_alt = (f"chr{chrom}", pos)
                if key_alt in vcf_index:
                    record = vcf_index[key_alt]
        
        if not record:
            continue
        
        # Get genotype
        sample = list(record.samples.values())[0]
        gt = sample.get('GT', None)
        if gt is None:
            dropped += 1
            continue
        
        ref = record.ref
        alt = ','.join(record.alts) if record.alts else ''
        
        dosage, flipped = get_genotype_dosage(str(gt), ref, alt, nucleus_snp.a1)
        
        if dosage is not None:
            matched.append(MatchedSNP(
                rsid=rsid,
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                a1=nucleus_snp.a1,
                genotype=dosage,
                flipped=flipped
            ))
        else:
            dropped += 1
    
    return matched, len(matched), dropped


def calculate_prs(matched_snps: List[MatchedSNP], nucleus_data: Dict[str, NucleusSNP]) -> Dict[str, float]:
    """
    Calculate PRS for each disease.
    
    Formula: PRS_D = sum(beta_{i,D} * G_i)
    where beta_{i,D} is effect size and G_i is genotype dosage (0, 1, or 2).
    """
    prs_scores = {disease: 0.0 for disease in DISEASES}
    
    for matched in matched_snps:
        if matched.rsid not in nucleus_data:
            continue
        
        nucleus_snp = nucleus_data[matched.rsid]
        
        for disease in DISEASES:
            beta = nucleus_snp.effect_sizes[disease]
            prs_scores[disease] += beta * matched.genotype
    
    return prs_scores


def save_results(
    prs_scores: Dict[str, float],
    output_path: str = "prs_results.csv"
):
    """Save PRS results to CSV file."""
    with open(output_path, 'w') as f:
        f.write("Disease,PRS\n")
        for disease in DISEASES:
            f.write(f"{disease},{prs_scores[disease]:.15f}\n")


def save_snp_level_results(
    matched_snps: List[MatchedSNP],
    nucleus_data: Dict[str, NucleusSNP],
    output_path: str = "prs_snp_level.csv"
):
    """Save SNP-level results used in PRS calculation."""
    with open(output_path, 'w') as f:
        writer = csv.writer(f)
        # Header
        header = ["rsID", "chrom", "pos", "ref", "alt", "A1", "genotype", "flipped"] + DISEASES
        writer.writerow(header)
        
        # Data rows
        for matched in matched_snps:
            if matched.rsid not in nucleus_data:
                continue
            
            nucleus_snp = nucleus_data[matched.rsid]
            row = [
                matched.rsid,
                matched.chrom,
                matched.pos,
                matched.ref,
                matched.alt,
                matched.a1,
                matched.genotype,
                "Yes" if matched.flipped else "No"
            ]
            
            # Add contribution to each disease
            for disease in DISEASES:
                beta = nucleus_snp.effect_sizes[disease]
                contribution = beta * matched.genotype
                row.append(f"{contribution:.15f}")
            
            writer.writerow(row)


def save_summary(
    nucleus_count: int,
    matched_count: int,
    dropped_count: int,
    prs_scores: Dict[str, float],
    output_path: str = "prs_summary.txt"
):
    """Save summary statistics to text file."""
    with open(output_path, 'w') as f:
        f.write("PRS Calculation Summary\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Total SNPs in Nucleus reference: {nucleus_count:,}\n")
        f.write(f"SNPs matched and used for PRS calculation: {matched_count:,}\n")
        f.write(f"SNPs dropped (allele mismatch or missing): {dropped_count:,}\n")
        f.write(f"Match rate: {matched_count / nucleus_count * 100:.2f}%\n\n")
        f.write("PRS Results:\n")
        f.write("-" * 70 + "\n")
        f.write(f"{'Disease':<35} {'PRS':>20}\n")
        f.write("-" * 70 + "\n")
        for disease in DISEASES:
            f.write(f"{disease:<35} {prs_scores[disease]:>20.15f}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate PRS scores from genome VCF",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python calculate_prs.py --nucleus /path/to/NUCLEUS_ORIGIN_V1.tsv --vcf /path/to/genome.vcf.gz
  
  # With mapping file (faster, no API calls)
  python calculate_prs.py --nucleus /path/to/NUCLEUS_ORIGIN_V1.tsv --vcf /path/to/genome.vcf.gz --mapping /path/to/mapping.tsv
  
  # Test with first 1000 SNPs
  python calculate_prs.py --nucleus /path/to/NUCLEUS_ORIGIN_V1.tsv --vcf /path/to/genome.vcf.gz --test
        """
    )
    parser.add_argument("--test", action="store_true", help="Process only first 1000 SNPs (for testing)")
    parser.add_argument("--nucleus", type=str, required=True, help="Path to Nucleus reference TSV file")
    parser.add_argument("--vcf", type=str, required=True, help="Path to genome VCF file")
    parser.add_argument("--mapping", type=str, required=True, help="Path to rsID-to-coordinate mapping file (TSV: rsID,chr,pos,ref,alt)")
    parser.add_argument("--output-dir", type=str, default=".", help="Output directory for results")
    
    args = parser.parse_args()
    
    # Validate input files
    if not Path(args.nucleus).exists():
        print(f"Error: Nucleus reference file not found: {args.nucleus}")
        sys.exit(1)
    
    if not Path(args.vcf).exists():
        print(f"Error: VCF file not found: {args.vcf}")
        sys.exit(1)
    
    if args.mapping and not Path(args.mapping).exists():
        print(f"Error: Mapping file not found: {args.mapping}")
        sys.exit(1)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("PRS Calculation Pipeline")
    print("=" * 70)
    print()
    
    print("Loading Nucleus reference...")
    nucleus_data = load_nucleus_reference(args.nucleus)
    print(f"Loaded {len(nucleus_data):,} SNPs from Nucleus reference")
    print()
    
    # Load mapping file (required)
    print(f"Loading mapping file: {args.mapping}")
    mapping = load_mapping_file(args.mapping)
    print(f"Loaded {len(mapping):,} rsID mappings")
    print()
    
    # Match SNPs
    print("Matching SNPs...")
    matched_snps = match_snps_by_rsid(args.vcf, nucleus_data)
    
    if not matched_snps:
        print("No rsIDs found in VCF, using position-based matching...")
        test_limit = 1000 if args.test else None
        matched_snps, matched_count, dropped_count = match_snps_by_position(
            args.vcf,
            nucleus_data,
            mapping=mapping,
            test_limit=test_limit
        )
        print(f"Matched {matched_count:,} SNPs, dropped {dropped_count:,} SNPs")
    else:
        print(f"Matched {len(matched_snps):,} SNPs by rsID")
        matched_count = len(matched_snps)
        dropped_count = 0
    print()
    
    if not matched_snps:
        print("Error: No SNPs matched. Please check your VCF and mapping files.")
        sys.exit(1)
    
    print("Calculating PRS scores...")
    prs_scores = calculate_prs(matched_snps, nucleus_data)
    print()
    
    print("Saving results...")
    results_path = output_dir / "prs_results.csv"
    snp_level_path = output_dir / "prs_snp_level.csv"
    summary_path = output_dir / "prs_summary.txt"
    
    save_results(prs_scores, str(results_path))
    save_snp_level_results(matched_snps, nucleus_data, str(snp_level_path))
    save_summary(
        nucleus_count=len(nucleus_data),
        matched_count=matched_count,
        dropped_count=dropped_count,
        prs_scores=prs_scores,
        output_path=str(summary_path)
    )
    
    print(f"Results saved to:")
    print(f"  - Disease-level PRS: {results_path}")
    print(f"  - SNP-level details: {snp_level_path}")
    print(f"  - Summary: {summary_path}")
    print()
    print("Done!")


if __name__ == "__main__":
    main()
