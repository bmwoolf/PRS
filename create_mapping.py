#!/usr/bin/env python3
"""
Create rsID-to-coordinate mapping file from Nucleus TSV.

This script queries Ensembl API to get chromosome/position for each rsID
in the Nucleus reference, then saves to a TSV mapping file.

Usage:
    python create_mapping.py --nucleus PATH --output mapping.tsv [--limit N]
"""

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from typing import Dict, Optional

try:
    import requests
except ImportError:
    print("Error: requests not installed. Run: uv pip install requests")
    sys.exit(1)

ENSEMBL_API_BASE = "https://rest.ensembl.org"
RATE_LIMIT_DELAY = 1.0 / 15  # 15 requests per second


def load_cache(cache_file: str = "rsid_positions_cache.json") -> Dict:
    """Load existing cache."""
    if Path(cache_file).exists():
        with open(cache_file, 'r') as f:
            return json.load(f)
    return {}


def save_cache(cache: Dict, cache_file: str = "rsid_positions_cache.json"):
    """Save cache to file."""
    with open(cache_file, 'w') as f:
        json.dump(cache, f, indent=2)


def fetch_rsid_position(rsid: str, cache: Dict) -> Optional[Dict]:
    """Fetch position from Ensembl API."""
    if rsid in cache:
        return cache[rsid]
    
    try:
        url = f"{ENSEMBL_API_BASE}/variation/human/{rsid}"
        params = {"content-type": "application/json"}
        
        time.sleep(RATE_LIMIT_DELAY)
        response = requests.get(url, params=params, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            mappings = data.get("mappings", [])
            
            for mapping in mappings:
                if mapping.get("assembly_name") == "GRCh38":
                    chrom = mapping.get("seq_region_name", "")
                    pos = mapping.get("start")
                    allele_string = mapping.get("allele_string", "")
                    
                    if "/" in allele_string:
                        alleles = allele_string.split("/")
                        ref = alleles[0]
                        alt = "/".join(alleles[1:])
                    else:
                        ref = ""
                        alt = ""
                    
                    result = {
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt
                    }
                    cache[rsid] = result
                    return result
        
        return None
    except Exception:
        return None


def load_nucleus_rsids(nucleus_path: str, limit: Optional[int] = None) -> list:
    """Load rsIDs from Nucleus TSV."""
    rsids = []
    
    open_func = open
    if nucleus_path.endswith('.gz'):
        import gzip
        open_func = gzip.open
    
    mode = 'rt' if nucleus_path.endswith('.gz') else 'r'
    
    with open_func(nucleus_path, mode) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rsid = row.get('SNP', '').strip()
            if rsid.startswith('rs'):
                rsids.append(rsid)
                if limit and len(rsids) >= limit:
                    break
    
    return rsids


def main():
    parser = argparse.ArgumentParser(
        description="Create rsID-to-coordinate mapping file"
    )
    parser.add_argument("--nucleus", type=str, required=True,
                       help="Path to Nucleus TSV file")
    parser.add_argument("--output", type=str, default="rsid_mapping.tsv",
                       help="Output mapping file path")
    parser.add_argument("--limit", type=int, default=None,
                       help="Limit number of rsIDs to process (for testing)")
    parser.add_argument("--cache", type=str, default="rsid_positions_cache.json",
                       help="Cache file path")
    
    args = parser.parse_args()
    
    if not Path(args.nucleus).exists():
        print(f"Error: Nucleus file not found: {args.nucleus}")
        sys.exit(1)
    
    print("Loading rsIDs from Nucleus reference...")
    rsids = load_nucleus_rsids(args.nucleus, limit=args.limit)
    print(f"Found {len(rsids):,} rsIDs")
    
    print("Loading cache...")
    cache = load_cache(args.cache)
    print(f"Cache contains {len(cache):,} entries")
    
    print("Fetching positions from Ensembl API...")
    print("(This may take a while - positions are cached for future runs)")
    
    mapping_data = []
    fetched = 0
    cached = 0
    
    for idx, rsid in enumerate(rsids):
        if (idx + 1) % 100 == 0:
            print(f"  Processed {idx + 1:,} / {len(rsids):,} ({fetched} fetched, {cached} cached)")
        
        pos_info = fetch_rsid_position(rsid, cache)
        
        if pos_info:
            if rsid in cache and idx < len(rsids) - 100:  # Was already cached
                cached += 1
            else:
                fetched += 1
            
            mapping_data.append({
                'rsid': rsid,
                'chr': pos_info['chrom'],
                'pos': pos_info['pos'],
                'ref': pos_info['ref'],
                'alt': pos_info['alt']
            })
    
    # Save cache periodically
    save_cache(cache, args.cache)
    
    print(f"\nFetched {fetched} new positions, used {cached} from cache")
    print(f"Total mappings: {len(mapping_data):,}")
    
    print(f"\nSaving mapping to {args.output}...")
    with open(args.output, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['rsID', 'chr', 'pos', 'ref', 'alt'])
        for item in mapping_data:
            writer.writerow([
                item['rsid'],
                item['chr'],
                item['pos'],
                item['ref'],
                item['alt']
            ])
    
    print(f"Done! Mapping file saved to {args.output}")
    print(f"\nYou can now use this mapping file with:")
    print(f"  python calculate_prs.py --nucleus ... --vcf ... --mapping {args.output}")


if __name__ == "__main__":
    main()

