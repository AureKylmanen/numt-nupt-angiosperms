#!/bin/bash

# ============================================================
# Script Name: numt_nupt_detection.sh
# Description: Detect NUMTs and NUPTs
# Author: Aure Kylmänen
# Adapted from Yu-Chi Chen vertebrate numt script
# Created: 26.02.2026
# Requires get_bed2.R (same directory)
# Usage: THREADS=24 WORKING_DIR=/path ./numt_nupt_detection.sh combined_dataset.csv
# ============================================================

set -euo pipefail

# ============================
# CONFIG
# ============================
THREADS="${THREADS:-12}"                             # BLAST threads
WORKING_DIR="${WORKING_DIR:-/vol/storage/5T_disk/organelles}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # get_bed2.R lives here
CSV_FILE="${1:-dataset.csv}"  # CSV: SPECIES,ASSEMBLY,PLASTID,MITO
MERGE_DISTANCE=10000   # distance in bp for merging hits (used in get_bed2.R)
MIN_LENGTH=200         # minimum insertion length retained

# ============================
# TOOLS
# ============================

REQUIRED_TOOLS=("blastn" "makeblastdb" "bedtools" "Rscript")
for tool in "${REQUIRED_TOOLS[@]}"; do
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: Tool $tool not found in PATH."
        exit 1
    fi
done

# ============================
# LOOP THROUGH SPECIES CSV
# ============================

while IFS=',' read -r SPECIES ASSEMBLY PLASTID MITO; do
    # Skip header or empty lines
    [[ -z "$ASSEMBLY" || "$ASSEMBLY" == "ASSEMBLY_ACCESSION" ]] && continue

    echo "===================================="
    echo "Processing species: $SPECIES"

   
    SPECIES_DIR="$WORKING_DIR/species/$SPECIES"
    mkdir -p "$SPECIES_DIR" "$SPECIES_DIR/results" "$SPECIES_DIR/numts" "$SPECIES_DIR/nupts"

    # Paths to downloaded files
    GENOME_FILE=$(find "$WORKING_DIR/downloads/genome_fasta" -name "${ASSEMBLY}_*.fna.gz")
    MITO_FILE=$(find "$WORKING_DIR/downloads/mito_fasta" -name "${MITO}.fasta")
    PLASTID_FILE=$(find "$WORKING_DIR/downloads/plastid_fasta" -name "${PLASTID}.fasta")

    # Check that genome and organelle files exist
    if [[ ! -f "$GENOME_FILE" ]]; then
        echo "WARNING: Genome fasta not found for $SPECIES, skipping..."
        continue
    fi
    if [[ ! -f "$MITO_FILE" ]]; then
        echo "WARNING: Mito fasta not found for $SPECIES, skipping NUMT analysis..."
    fi
    if [[ ! -f "$PLASTID_FILE" ]]; then
        echo "WARNING: Plastid fasta not found for $SPECIES, skipping NUPT analysis..."
    fi

    # Uncompress genome if needed
    GENOME_UNCOMPRESSED="${GENOME_FILE%.gz}"
    if [[ ! -f "$GENOME_UNCOMPRESSED" ]]; then
        gunzip -c "$GENOME_FILE" > "$GENOME_UNCOMPRESSED"
    fi

    # Compute total genome size
    GENOME_SIZE=$(grep -v ">" "$GENOME_UNCOMPRESSED" | tr -d '\n' | wc -c)
    echo "Genome size: $GENOME_SIZE bp"
    echo "$GENOME_SIZE" > "$SPECIES_DIR/results/$SPECIES.genome_size.txt"

    # ============================
    # NUMT DETECTION
    # ============================
    if [[ -f "$MITO_FILE" ]]; then
        echo "Running NUMT detection..."
        NUMT_DIR="$SPECIES_DIR/numts"
        makeblastdb -in "$MITO_FILE" -dbtype nucl -parse_seqids  # optional for mito DB
        makeblastdb -in "$GENOME_UNCOMPRESSED" -dbtype nucl -parse_seqids

        blastn -db "$GENOME_UNCOMPRESSED" \
               -query "$MITO_FILE" \
               -outfmt 7 \
               -word_size 20 \
               -num_threads "$THREADS" \
               -out "$NUMT_DIR/numts.blast.out" || echo "WARNING: Blast failed for NUMTs"

        # Clean blast output
        grep -v "#" "$NUMT_DIR/numts.blast.out" | cut -f2,3,4,9,10 > "$NUMT_DIR/numts.clean.out"

        # Generate BED and extract sequences via Rscript + bedtools
        cd "$NUMT_DIR"

        echo "numts.clean.out" > files_list.txt

        Rscript "$SCRIPT_DIR/get_bed2.R" "$MIN_LENGTH" "$MERGE_DISTANCE" || {
            echo "WARNING: Rscript for NUMTs failed"
            cd - >/dev/null
            continue
        }

        cd - >/dev/null

	if [[ ! -f "$NUMT_DIR/numts.clean.out_AllFrags.bed" ]]; then
		echo "No valid NUMTs for $SPECIES, skipping..."
		continue
	fi

        bedtools getfasta -fi "$GENOME_UNCOMPRESSED" \
                          -bed "$NUMT_DIR/numts.clean.out_AllFrags.bed" \
                          -fo "$NUMT_DIR/numts.seq.fa"

        NUMT_COUNT=$(grep -c ">" "$NUMT_DIR/numts.seq.fa")
        echo "$NUMT_COUNT" > "$SPECIES_DIR/results/$SPECIES.numt_count.txt"
        echo "NUMTs detected: $NUMT_COUNT"

        NUMT_MERGED_COUNT=$(wc -l < "$NUMT_DIR/numts.clean.out_AllFrags_merged${MERGE_DISTANCE}.bed" 2>/dev/null || echo "0")
        echo "$NUMT_MERGED_COUNT" > "$SPECIES_DIR/results/${SPECIES}.numt_merged_count.txt"

        # -------------------------------
        # Calculate genome proportions
        # -------------------------------

        UNMERGED_BED="$NUMT_DIR/numts.clean.out_AllFrags.bed"
        MERGED_BED="$NUMT_DIR/numts.clean.out_AllFrags_merged${MERGE_DISTANCE}.bed"

        if [[ -f "$UNMERGED_BED" ]]; then
            UNMERGED_BP=$(awk '{sum += $3-$2} END {print sum}' "$UNMERGED_BED")
            UNMERGED_PROP=$(echo "$UNMERGED_BP $GENOME_SIZE" | awk '{print $1/$2}')
            echo "$UNMERGED_PROP" > "$SPECIES_DIR/results/$SPECIES.numt_unmerged_prop.txt"
        fi

        if [[ -f "$MERGED_BED" ]]; then
            MERGED_BP=$(awk '{sum += $3-$2} END {print sum}' "$MERGED_BED")
            MERGED_PROP=$(echo "$MERGED_BP $GENOME_SIZE" | awk '{print $1/$2}')
            echo "$MERGED_PROP" > "$SPECIES_DIR/results/$SPECIES.numt_merged_prop.txt"
        fi

        awk '{print $1,$2,$3,$3-$2}' OFS="\t" "$MERGED_BED" > \
        "$SPECIES_DIR/results/$SPECIES.numt_merged_locations_lengths.tsv"
    fi

    # ============================
    # NUPT DETECTION
    # ============================
    if [[ -f "$PLASTID_FILE" ]]; then
        echo "Running NUPT detection..."
        NUPT_DIR="$SPECIES_DIR/nupts"
        makeblastdb -in "$PLASTID_FILE" -dbtype nucl -parse_seqids

        blastn -db "$GENOME_UNCOMPRESSED" \
               -query "$PLASTID_FILE" \
               -outfmt 7 \
               -word_size 20 \
               -num_threads "$THREADS" \
               -out "$NUPT_DIR/nupts.blast.out" || echo "WARNING: Blast failed for NUPTs"

        grep -v "#" "$NUPT_DIR/nupts.blast.out" | cut -f2,3,4,9,10 > "$NUPT_DIR/nupts.clean.out"

        cd "$NUPT_DIR"

        echo "nupts.clean.out" > files_list.txt

        Rscript "$SCRIPT_DIR/get_bed2.R" "$MIN_LENGTH" "$MERGE_DISTANCE" || {
            echo "WARNING: Rscript for NUPTs failed"
            cd - >/dev/null
            continue
        }

        cd - >/dev/null

	if [[ ! -f "$NUPT_DIR/nupts.clean.out_AllFrags.bed" ]]; then
		echo "No valid NUPTs for $SPECIES, skipping..."
		continue
	fi

        bedtools getfasta -fi "$GENOME_UNCOMPRESSED" \
                          -bed "$NUPT_DIR/nupts.clean.out_AllFrags.bed" \
                          -fo "$NUPT_DIR/nupts.seq.fa"

        NUPT_COUNT=$(grep -c ">" "$NUPT_DIR/nupts.seq.fa")
        echo "$NUPT_COUNT" > "$SPECIES_DIR/results/$SPECIES.nupt_count.txt"
        echo "NUPTs detected: $NUPT_COUNT"

        # NUPT merged count
        NUPT_MERGED_COUNT=$(wc -l < "$NUPT_DIR/nupts.clean.out_AllFrags_merged${MERGE_DISTANCE}.bed" 2>/dev/null || echo "0")
        echo "$NUPT_MERGED_COUNT" > "$SPECIES_DIR/results/${SPECIES}.nupt_merged_count.txt"

        UNMERGED_BED="$NUPT_DIR/nupts.clean.out_AllFrags.bed"
        MERGED_BED="$NUPT_DIR/nupts.clean.out_AllFrags_merged${MERGE_DISTANCE}.bed"

        if [[ -f "$UNMERGED_BED" ]]; then
            UNMERGED_BP=$(awk '{sum += $3-$2} END {print sum}' "$UNMERGED_BED")
            UNMERGED_PROP=$(echo "$UNMERGED_BP $GENOME_SIZE" | awk '{print $1/$2}')
            echo "$UNMERGED_PROP" > "$SPECIES_DIR/results/$SPECIES.nupt_unmerged_prop.txt"
        fi

        if [[ -f "$MERGED_BED" ]]; then
            MERGED_BP=$(awk '{sum += $3-$2} END {print sum}' "$MERGED_BED")
            MERGED_PROP=$(echo "$MERGED_BP $GENOME_SIZE" | awk '{print $1/$2}')
            echo "$MERGED_PROP" > "$SPECIES_DIR/results/$SPECIES.nupt_merged_prop.txt"
        fi

        awk '{print $1,$2,$3,$3-$2}' OFS="\t" "$UNMERGED_BED" > \
        "$SPECIES_DIR/results/$SPECIES.nupt_unmerged_locations_lengths.tsv"
    fi

    echo "Finished species: $SPECIES"
    echo "===================================="

    # ============================
    # NUMT vs NUPT OVERLAP
    # ============================

    NUMT_MERGED="$SPECIES_DIR/numts/numts.clean.out_AllFrags_merged${MERGE_DISTANCE}.bed"
    NUPT_MERGED="$SPECIES_DIR/nupts/nupts.clean.out_AllFrags_merged${MERGE_DISTANCE}.bed"

    if [[ -f "$NUMT_MERGED" && -f "$NUPT_MERGED" ]]; then

        bedtools intersect -sorted -a "$NUMT_MERGED" -b "$NUPT_MERGED" -wa -wb > \
            "$SPECIES_DIR/results/$SPECIES.numt_nupt_overlap_detailed.tsv"

        # Count overlapping loci
        OVERLAP_COUNT=$(wc -l < \
            "$SPECIES_DIR/results/$SPECIES.numt_nupt_overlap_detailed.tsv")

        echo "$OVERLAP_COUNT" > \
            "$SPECIES_DIR/results/$SPECIES.numt_nupt_overlap_count.txt"

        # Calculate overlap bp
        bedtools intersect -sorted -a "$NUMT_MERGED" -b "$NUPT_MERGED" > \
            "$SPECIES_DIR/results/tmp_overlap.bed"

        OVERLAP_BP=$(awk '{sum += $3-$2} END {print sum}' \
            "$SPECIES_DIR/results/tmp_overlap.bed")

        echo "$OVERLAP_BP" > \
            "$SPECIES_DIR/results/$SPECIES.numt_nupt_overlap_bp.txt"

        rm "$SPECIES_DIR/results/tmp_overlap.bed"

        echo "NUMT-NUPT overlapping loci: $OVERLAP_COUNT"
        echo "Overlap bp: $OVERLAP_BP"
fi

done < "$CSV_FILE"

echo "All NUMT/NUPT analyses completed."
