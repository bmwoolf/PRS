#!/usr/bin/env python3
"""
Create rsID-to-coordinate mapping file using bcftools + dbSNP VCF.

Much faster than API calls - extracts mappings directly from dbSNP VCF.
"""

import argparse
import csv
import subprocess
import sys
from pathlib import Path


def extract_rsids_from_nucleus(nucleus_path: str, limit: int = None) -> list:
    """Extract rsIDs from Nucleus TSV."""
    rsids = []
    open_func = open
    if nucleus_path.endswith('.gz'):
        import gzip
        open_func = gzip.open
    
    mode = 'rt' if nucleus_path.endswith('.gz') else 'r'
    
    with open_func(nucleus_path, mode) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for i, row in enumerate(reader):
            if limit and i >= limit:
                break
            rsid = row.get('SNP', '').strip()
            if rsid.startswith('rs'):
                rsids.append(rsid)
    
    return rsids


def create_rsid_list_file(rsids: list, output_file: str):
    """Create a file with one rsID per line for bcftools query."""
    with open(output_file, 'w') as f:
        for rsid in rsids:
            f.write(f"{rsid}\n")


def extract_mappings_from_dbsnp(dbsnp_vcf: str, rsid_list_file: str, output_tsv: str):
    """
    Use bcftools to extract rsID mappings from dbSNP VCF.
    
    Creates TSV with: rsID, chr, pos, ref, alt
    """
    print(f"Querying dbSNP VCF: {dbsnp_vcf}")
    rsid_count = sum(1 for _ in open(rsid_list_file))
    print(f"Extracting mappings for {rsid_count:,} rsIDs...")
    
    # Check if VCF is indexed, create index if needed
    index_file = dbsnp_vcf + '.tbi' if dbsnp_vcf.endswith('.gz') else dbsnp_vcf + '.csi'
    if not Path(index_file).exists():
        print("Indexing dbSNP VCF (this may take a few minutes)...")
        index_cmd = ['bcftools', 'index', '-t', dbsnp_vcf]
        try:
            subprocess.run(index_cmd, check=True, capture_output=True)
            print("Index created successfully")
        except subprocess.CalledProcessError as e:
            print(f"Warning: Could not create index: {e.stderr.decode()}", file=sys.stderr)
            print("Continuing without index (may be slower)...")
    
    # Load rsIDs into a set for fast lookup
    rsid_set = set()
    with open(rsid_list_file, 'r') as f:
        for line in f:
            rsid = line.strip()
            if rsid:
                rsid_set.add(rsid)
    
    # Use bcftools query to extract all variants, then filter in Python
    # This is more reliable than trying to filter in bcftools
    print("Scanning dbSNP VCF for matching rsIDs...")
    cmd = [
        'bcftools', 'query',
        '-f', '%ID\t%CHROM\t%POS\t%REF\t%ALT\n',
        dbsnp_vcf
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Write to TSV file, filtering for our rsIDs
        matches = 0
        with open(output_tsv, 'w') as f:
            f.write("rsID\tchr\tpos\tref\talt\n")
            for line in result.stdout.strip().split('\n'):
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 5:
                    rsid_field = parts[0]
                    # Handle multiple rsIDs (semicolon-separated)
                    rsids_in_field = [r.strip() for r in rsid_field.split(';') if r.strip().startswith('rs')]
                    
                    # Check if any rsID in this variant matches our list
                    for rsid in rsids_in_field:
                        if rsid in rsid_set:
                            f.write(f"{rsid}\t{parts[1]}\t{parts[2]}\t{parts[3]}\t{parts[4]}\n")
                            matches += 1
                            rsid_set.discard(rsid)  # Remove to avoid duplicates
                            break  # Only write once per variant
        
        print(f"Extracted {matches:,} mappings saved to {output_tsv}")
        if len(rsid_set) > 0:
            print(f"Warning: {len(rsid_set):,} rsIDs not found in dbSNP VCF")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running bcftools: {e.stderr}", file=sys.stderr)
        if "not indexed" in e.stderr.decode().lower():
            print("\nTry indexing the VCF first: bcftools index -t", dbsnp_vcf, file=sys.stderr)
        return False
    except FileNotFoundError:
        print("Error: bcftools not found. Install with: conda install -c bioconda bcftools", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Create rsID mapping file using bcftools + dbSNP VCF",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  # Download dbSNP VCF first (one-time):
  # wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151/VCF/00-common_all.vcf.gz
  
  # Create mapping:
  python create_mapping_bcftools.py \\
    --nucleus NUCLEUS_ORIGIN_V1.tsv \\
    --dbsnp 00-common_all.vcf.gz \\
    --output rsid_mapping.tsv
        """
    )
    parser.add_argument("--nucleus", required=True, help="Path to Nucleus TSV")
    parser.add_argument("--dbsnp", required=True, help="Path to dbSNP VCF file (can be gzipped)")
    parser.add_argument("--output", default="rsid_mapping.tsv", help="Output mapping TSV file")
    parser.add_argument("--limit", type=int, help="Limit number of rsIDs (for testing)")
    parser.add_argument("--rsid-list", default="rsid_list.txt", help="Temporary file for rsID list")
    
    args = parser.parse_args()
    
    if not Path(args.dbsnp).exists():
        print(f"Error: dbSNP VCF not found: {args.dbsnp}")
        print("\nDownload dbSNP VCF from:")
        print("  https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151/VCF/")
        print("  (e.g., 00-common_all.vcf.gz)")
        sys.exit(1)
    
    # Extract rsIDs from Nucleus
    print(f"Loading rsIDs from {args.nucleus}...")
    rsids = extract_rsids_from_nucleus(args.nucleus, args.limit)
    print(f"Found {len(rsids):,} rsIDs")
    
    # Create rsID list file
    print(f"Creating rsID list file...")
    create_rsid_list_file(rsids, args.rsid_list)
    
    # Extract mappings using bcftools
    success = extract_mappings_from_dbsnp(args.dbsnp, args.rsid_list, args.output)
    
    # Clean up
    if Path(args.rsid_list).exists():
        Path(args.rsid_list).unlink()
    
    if success:
        # Count results
        count = sum(1 for _ in open(args.output)) - 1  # Subtract header
        print(f"\nDone! Created mapping file with {count:,} entries")
        print(f"Use with: python calculate_prs.py --mapping {args.output} ...")
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()

