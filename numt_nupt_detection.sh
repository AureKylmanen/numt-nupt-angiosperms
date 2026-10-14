#!/bin/bash

# ============================================================
# Script Name: numt_nupt_detection.sh
# Description: Detect NUMTs and NUPTs in a nuclear assembly
# Author: Aure Kylmänen
# Adapted from Yu-Chi Chen vertebrate numt script
# Created: 26.02.2026   Last updated: 14.10.2026
#
# Steps per species:
#   1. BLASTN of the plastome and the mitogenome (all chromosomes of a
#      multi-chromosome mitogenome, concatenated) against the nuclear assembly,
#      without a cap on the number of target sequences.
#   2. Hits on assembly sequences annotated as organelle in their FASTA header
#      (mitochondrion / chloroplast / plastid / organelle) are removed: these are
#      organelle genomes deposited with the assembly, not insertions.
#   3. get_bed2.R: hits >= MIN_LENGTH bp, oriented, sorted -> fragment BED.
#   4. Fragments within MERGE_DISTANCE bp on the same sequence are merged
#      (bedtools merge -d) into insertions.
#   5. Fragment sequences, counts, genome proportions, NUMT-NUPT overlap.
#
# Requires get_bed2.R (same directory), blastn, makeblastdb, bedtools, Rscript,
# samtools (fallback for sequence extraction from very large genomes).
#
# Usage: THREADS=24 WORKING_DIR=/path ./numt_nupt_detection.sh combined_dataset.csv
#   CSV columns: SPECIES,ASSEMBLY,PLASTID,MITO   (MITO: accessions separated by ';'
#   for multi-chromosome mitogenomes). A header line is skipped.
# ============================================================

set -uo pipefail

# ============================
# CONFIG
# ============================
THREADS="${THREADS:-12}"
WORKING_DIR="${WORKING_DIR:-${BASE_DIR:-$PWD}}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CSV_FILE="${1:?usage: numt_nupt_detection.sh species.csv}"
MERGE_DISTANCE=10000     # fragments <= this distance apart are one insertion
MIN_LENGTH=200           # minimum fragment length (bp)
MAX_TARGET_SEQS=1000000  # BLAST default (500) truncates detection in fragmented assemblies
ORGANELLE_HEADER_REGEX="mitochondri|chloroplast|plastid|organelle"

DL="$WORKING_DIR/downloads"
DL_MULTI="$WORKING_DIR/multi_organelles/downloads"   # organelle FASTAs of multi-chromosome mitogenomes
export LC_ALL=C

# ============================
# TOOLS
# ============================
for tool in blastn makeblastdb bedtools Rscript samtools; do
    command -v "$tool" &>/dev/null || { echo "ERROR: $tool not found in PATH"; exit 1; }
done

# ============================
# HELPERS
# ============================
bp() { awk '{s += $3 - $2} END {print s + 0}' "$@"; }

# organelle FASTA for an accession: exact name first, then any version suffix
find_fa() {   # $1 = plastid|mito, $2 = accession
    local k=$1 a=$2 f
    for f in "$DL/${k}_fasta/$a.fasta" "$DL_MULTI/${k}_fasta/$a.fasta"; do
        [[ -s $f ]] && { echo "$f"; return; }
    done
    f=$(ls "$DL/${k}_fasta/${a%.*}".[0-9]*.fasta "$DL_MULTI/${k}_fasta/${a%.*}".[0-9]*.fasta 2>/dev/null | head -1)
    [[ -s $f ]] && echo "$f"
}

# sequences for a BED; bedtools getfasta, or samtools faidx when getfasta fails
# (bedtools getfasta segfaults on some multi-Gb genomes). Headers: chrom:start-end (0-based start).
extract_seqs() {   # $1 genome, $2 bed, $3 out.fa
    if bedtools getfasta -fi "$1" -bed "$2" -fo "$3" 2>/dev/null; then return 0; fi
    echo "  bedtools getfasta failed; using samtools faidx"
    [[ -s $1.fai ]] || samtools faidx "$1" || return 1
    awk '{print $1":"$2+1"-"$3}' "$2" > "$3.regions"
    samtools faidx "$1" -r "$3.regions" \
      | awk '/^>/ { if (match($0, /:[0-9]+-[0-9]+$/)) {
                       r = substr($0, RSTART + 1); split(r, p, "-")
                       print substr($0, 1, RSTART) (p[1] - 1) "-" p[2]; next } }
             { print }' > "$3"
    rm -f "$3.regions"
}

# one organelle class against the nuclear assembly
detect_class() {   # $1 numt|nupt, $2 query FASTA
    local c=$1 Q=$2
    local D="$SPECIES_DIR/${c}s" R="$SPECIES_DIR/results"
    local FRAG="$D/${c}s.clean.out_AllFrags.bed"
    local MERGED="$D/${c}s.clean.out_AllFrags_merged${MERGE_DISTANCE}.bed"
    echo "Running ${c^^} detection..."

    blastn -db "$GENOME" -query "$Q" -outfmt 7 -word_size 20 \
           -max_target_seqs "$MAX_TARGET_SEQS" -num_threads "$THREADS" \
           -out "$D/${c}s.blast.out" || { echo "WARNING: BLAST failed for ${c^^}"; return; }

    # subject, identity, length, sstart, send; hits on annotated organelle sequences removed
    grep -v "^#" "$D/${c}s.blast.out" | cut -f2,3,4,9,10 \
      | awk 'FILENAME == ARGV[1] {o[$1]; next} !($1 in o)' "$ORG_LIST" - > "$D/${c}s.clean.out"

    rm -f "$FRAG" "$MERGED"
    ( cd "$D" && echo "${c}s.clean.out" > files_list.txt \
        && Rscript "$SCRIPT_DIR/get_bed2.R" "$MIN_LENGTH" ) || { echo "WARNING: get_bed2.R failed for ${c^^}"; return; }
    [[ -s $FRAG ]] || { echo "No ${c^^}s >= $MIN_LENGTH bp for $SPECIES"; return; }

    cut -f1-3 "$FRAG" | sort -k1,1 -k2,2n | bedtools merge -d "$MERGE_DISTANCE" > "$MERGED"
    extract_seqs "$GENOME" "$FRAG" "$D/${c}s.seq.fa" || { echo "WARNING: sequence extraction failed for ${c^^}"; return; }

    grep -c ">" "$D/${c}s.seq.fa"                       > "$R/$SPECIES.${c}_count.txt"
    wc -l < "$MERGED"                                    > "$R/$SPECIES.${c}_merged_count.txt"
    echo "$(bp "$FRAG") $GENOME_SIZE"   | awk '{print $1/$2}' > "$R/$SPECIES.${c}_unmerged_prop.txt"
    echo "$(bp "$MERGED") $GENOME_SIZE" | awk '{print $1/$2}' > "$R/$SPECIES.${c}_merged_prop.txt"
    # union of fragments (each bp once; e.g. plastid inverted-repeat hits counted once)
    cut -f1-3 "$FRAG" | sort -k1,1 -k2,2n | bedtools merge | bp > "$R/$SPECIES.${c}_union_bp.txt"
    awk -v OFS="\t" '{print $1,$2,$3,$3-$2}' "$FRAG"   > "$R/$SPECIES.${c}_unmerged_locations_lengths.tsv"
    awk -v OFS="\t" '{print $1,$2,$3,$3-$2}' "$MERGED" > "$R/$SPECIES.${c}_merged_locations_lengths.tsv"
    echo "${c^^}s: $(cat "$R/$SPECIES.${c}_count.txt") fragments, $(cat "$R/$SPECIES.${c}_merged_count.txt") merged insertions"
}

# ============================
# LOOP THROUGH SPECIES CSV
# ============================
while IFS=',' read -r SPECIES ASSEMBLY PLASTID MITO <&3; do
    [[ -z "$SPECIES" || "$SPECIES" == "SPECIES" || "$ASSEMBLY" == "ASSEMBLY_ACCESSION" ]] && continue
    echo "===================================="
    echo "Processing species: $SPECIES"

    SPECIES_DIR="$WORKING_DIR/species/$SPECIES"
    mkdir -p "$SPECIES_DIR/results" "$SPECIES_DIR/numts" "$SPECIES_DIR/nupts"
    R="$SPECIES_DIR/results"

    # ---- nuclear assembly ----
    GENOME_GZ=$(find "$DL/genome_fasta" -maxdepth 1 -name "${ASSEMBLY}_*.fna.gz" | head -1)
    GENOME="${GENOME_GZ%.gz}"
    [[ -s $GENOME || -s $GENOME_GZ ]] || { echo "WARNING: genome FASTA not found for $SPECIES, skipping"; continue; }
    [[ -s $GENOME ]] || gunzip -c "$GENOME_GZ" > "$GENOME"
    [[ -s "$GENOME.nsq" || -s "$GENOME.00.nsq" || -s "$GENOME.nal" ]] \
        || makeblastdb -in "$GENOME" -dbtype nucl -parse_seqids > /dev/null

    GENOME_SIZE=$(grep -v ">" "$GENOME" | tr -d '\n' | wc -c)
    echo "$GENOME_SIZE" > "$R/$SPECIES.genome_size.txt"
    echo "Genome size: $GENOME_SIZE bp"

    # ---- organelle sequences deposited with the assembly ----
    ORG_LIST="$R/$SPECIES.organelle_annotated_headers.txt"
    grep "^>" "$GENOME" | grep -iE "$ORGANELLE_HEADER_REGEX" | awk '{print substr($1, 2)}' > "$ORG_LIST"
    [[ -s $ORG_LIST ]] && echo "Excluding $(wc -l < "$ORG_LIST") organelle-annotated sequence(s) of the assembly"

    # ---- queries ----
    PLASTID_FILE=$(find_fa plastid "$PLASTID")
    MITO_FILE="$SPECIES_DIR/numts/mito_query.fasta"; : > "$MITO_FILE"; MISSING=""
    for a in ${MITO//;/ }; do
        f=$(find_fa mito "$a"); [[ -n $f ]] && cat "$f" >> "$MITO_FILE" || MISSING="$MISSING $a"
    done

    if [[ -n $MISSING || ! -s $MITO_FILE ]]; then
        echo "WARNING: mito FASTA missing for $SPECIES ($MISSING), skipping NUMT analysis"
    else
        detect_class numt "$MITO_FILE"
    fi
    if [[ -z $PLASTID_FILE ]]; then
        echo "WARNING: plastid FASTA not found for $SPECIES ($PLASTID), skipping NUPT analysis"
    else
        detect_class nupt "$PLASTID_FILE"
    fi

    # ---- excluded organelle sequences (for the Supplement) ----
    #   sequence, length, bp covered by excluded NUMT hits, bp covered by excluded NUPT hits
    OL="$R/$SPECIES.organelle_like_sequences.tsv"
    printf "sequence\tlength\tnumt_hit_bp\tnupt_hit_bp\n" > "$OL"
    if [[ -s $ORG_LIST ]]; then
        [[ -s "$GENOME.fai" ]] || samtools faidx "$GENOME"
        T=$(mktemp -d)
        for c in numt nupt; do
            : > "$T/$c"
            [[ -s $SPECIES_DIR/${c}s/${c}s.blast.out ]] || continue
            grep -v "^#" "$SPECIES_DIR/${c}s/${c}s.blast.out" \
              | awk -v OFS="\t" -v m=$MIN_LENGTH '$4 >= m {s = $9 < $10 ? $9 : $10; e = $9 < $10 ? $10 : $9; print $2, s, e}' \
              | awk 'FILENAME == ARGV[1] {o[$1]; next} ($1 in o)' "$ORG_LIST" - \
              | sort -k1,1 -k2,2n | bedtools merge \
              | awk -v OFS="\t" '{b[$1] += $3 - $2} END {for (k in b) print k, b[k]}' > "$T/$c"
        done
        awk -v OFS="\t" 'FILENAME == ARGV[1] {o[$1]; next}
                         FILENAME == ARGV[2] {n[$1] = $2; next}
                         FILENAME == ARGV[3] {p[$1] = $2; next}
                         ($1 in o) {print $1, $2, n[$1] + 0, p[$1] + 0}' \
            "$ORG_LIST" "$T/numt" "$T/nupt" "$GENOME.fai" | sort -k1,1 >> "$OL"
        rm -rf "$T"
    fi

    # ============================
    # NUMT vs NUPT OVERLAP
    # ============================
    #   overlap_bp     bp shared by merged NUMT and NUPT insertions (co-location)
    #   overlap_count  NUMT insertions overlapping at least one NUPT insertion
    #   shared.bed     sequence called by both queries on the fragment level
    #                  (plastid-derived regions of the mitogenome, MTPTs); shared_bp its length
    MN="$SPECIES_DIR/numts/numts.clean.out_AllFrags_merged${MERGE_DISTANCE}.bed"
    MP="$SPECIES_DIR/nupts/nupts.clean.out_AllFrags_merged${MERGE_DISTANCE}.bed"
    UN="$SPECIES_DIR/numts/numts.clean.out_AllFrags.bed"
    UP="$SPECIES_DIR/nupts/nupts.clean.out_AllFrags.bed"
    if [[ -s $MN && -s $MP && -s $UN && -s $UP ]]; then
        T=$(mktemp -d)
        sort -k1,1 -k2,2n "$MN" > "$T/mn"; sort -k1,1 -k2,2n "$MP" > "$T/mp"
        cut -f1-3 "$UN" | sort -k1,1 -k2,2n | bedtools merge > "$T/un"
        cut -f1-3 "$UP" | sort -k1,1 -k2,2n | bedtools merge > "$T/up"
        bedtools intersect -sorted -a "$T/mn" -b "$T/mp" -wa -wb > "$R/$SPECIES.numt_nupt_overlap_detailed.tsv"
        bedtools intersect -sorted -a "$T/un" -b "$T/up"         > "$R/$SPECIES.numt_nupt_shared.bed"
        bedtools intersect -sorted -a "$T/mn" -b "$T/mp" | bp    > "$R/$SPECIES.numt_nupt_overlap_bp.txt"
        bedtools intersect -sorted -u -a "$T/mn" -b "$T/mp" | wc -l > "$R/$SPECIES.numt_nupt_overlap_count.txt"
        bp "$R/$SPECIES.numt_nupt_shared.bed" > "$R/$SPECIES.numt_nupt_shared_bp.txt"
        echo "NUMT-NUPT overlap: $(cat "$R/$SPECIES.numt_nupt_overlap_count.txt") NUMT loci, $(cat "$R/$SPECIES.numt_nupt_overlap_bp.txt") bp; shared fragment sequence $(bp "$R/$SPECIES.numt_nupt_shared.bed") bp"
        rm -rf "$T"
    fi

    echo "Finished species: $SPECIES"
done 3< <(tr -d '\r' < "$CSV_FILE")

echo "All NUMT/NUPT analyses completed."
