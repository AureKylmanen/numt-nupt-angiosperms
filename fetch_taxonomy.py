#!/usr/bin/env python3
"""
fetch_taxonomy.py
-----------------
Given a tab-separated file with a 'species' column (underscore-separated names,
e.g. Acacia_confusa), fetches the taxonomic lineage from NCBI and writes a
mapping file:  species → order (or the nearest named rank if order is absent).

Usage:
    python fetch_taxonomy.py \
        --input  numt_counts.tsv \
        --output species_taxonomy.tsv \
        --email  your@email.com          # required by NCBI Entrez

Output columns:  species  order  family  class  phylum
"""

import argparse
import time
import sys
import pandas as pd
from Bio import Entrez

# ── ranks we want to capture (in priority order for the "display clade") ──────
RANKS_WANTED = ["order", "family", "class", "phylum", "kingdom"]


def fetch_lineage(species_name: str, retries: int = 3) -> dict:
    """Return a dict of rank→name for a single species."""
    query = species_name.replace("_", " ")
    for attempt in range(retries):
        try:
            handle = Entrez.esearch(db="taxonomy", term=query, retmax=1)
            record = Entrez.read(handle)
            handle.close()
            ids = record.get("IdList", [])
            if not ids:
                # Try genus only
                genus = query.split()[0]
                handle = Entrez.esearch(db="taxonomy", term=genus, retmax=1)
                record = Entrez.read(handle)
                handle.close()
                ids = record.get("IdList", [])
            if not ids:
                return {}

            handle = Entrez.efetch(db="taxonomy", id=ids[0], retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            lineage = {}
            for item in records[0].get("LineageEx", []):
                rank = item.get("Rank", "").lower()
                name = item.get("ScientificName", "")
                if rank in RANKS_WANTED:
                    lineage[rank] = name
            return lineage

        except Exception as e:
            print(f"  Attempt {attempt+1} failed for '{species_name}': {e}",
                  file=sys.stderr)
            time.sleep(2 ** attempt)   # exponential back-off
    return {}


def main():
    parser = argparse.ArgumentParser(description="Fetch NCBI taxonomy for species list.")
    parser.add_argument("--input",  required=True,  help="TSV with a 'species' column")
    parser.add_argument("--output", required=True,  help="Output TSV path")
    parser.add_argument("--email",  required=True,  help="Your e-mail (required by NCBI)")
    parser.add_argument("--delay",  type=float, default=0.4,
                        help="Seconds between requests (default 0.4; use 0.1 with API key)")
    parser.add_argument("--api-key", default=None,
                        help="NCBI API key (optional, raises rate limit to 10 req/s)")
    args = parser.parse_args()

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    df = pd.read_csv(args.input, sep="\t")
    if "species" not in df.columns:
        sys.exit("ERROR: input file must have a column named 'species'")

    species_list = df["species"].unique().tolist()
    print(f"Fetching taxonomy for {len(species_list)} species …")

    rows = []
    for i, sp in enumerate(species_list, 1):
        print(f"  [{i}/{len(species_list)}] {sp}", end="  ", flush=True)
        lineage = fetch_lineage(sp)
        row = {"species": sp}
        for rank in RANKS_WANTED:
            row[rank] = lineage.get(rank, "Unknown")
        rows.append(row)
        print(row.get("order", "?"))
        time.sleep(args.delay)

    out = pd.DataFrame(rows, columns=["species"] + RANKS_WANTED)
    out.to_csv(args.output, sep="\t", index=False)
    print(f"\nDone. Taxonomy saved to: {args.output}")


if __name__ == "__main__":
    main()
