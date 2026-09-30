#!/usr/bin/env python3

from Bio import SeqIO
from pathlib import Path

base = "/vol/storage/5T_disk/organelles/downloads"


plastid_gb = Path(base) / "plastid_gb"

plastid_out = Path(base) / "plastid_genes"

plastid_out.mkdir(exist_ok=True)

def extract_cds(gb_dir, out_dir):
    for gb_file in gb_dir.glob("*.gb*"):
        species_name = gb_file.stem  # keep full assembly name
        species_dir = out_dir / species_name
        species_dir.mkdir(exist_ok=True)

        for record in SeqIO.parse(gb_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    gene_name = feature.qualifiers.get("gene", ["unknown"])[0]
                    seq = feature.extract(record.seq)  # handles splicing automatically
                    out_file = species_dir / f"{gene_name}.fasta"
                    out_file.write_text(f">{gene_name}\n{seq}\n")

# Run for mito and plastid
extract_cds(plastid_gb, plastid_out)
