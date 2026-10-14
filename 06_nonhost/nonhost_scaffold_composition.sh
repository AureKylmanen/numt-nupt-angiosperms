#!/usr/bin/env bash
# ============================================================
# Script Name: nonhost_scaffold_composition.sh
# Description: Organellar composition of the contigs carrying high-confidence
#              non-host dcNUxT candidates (input to nonhost_summary.py, which
#              applies the exclusion rules below).
# Author: Aure Kylmänen
#
# Each candidate contig is searched against the pooled organelle database of its
# class (dc-megablast, e <= 1e-20) and the fraction of its length covered by
# merged hit intervals is recorded. Contigs longer than 300 kb are treated as
# nuclear without testing.
#
# Exclusion rules (Methods):
#   dcNUPT  contig organellar fraction >= 0.95. Plastid candidate contigs were
#           bimodal (0.00-0.05 or 0.95-1.00), so any cut-off between the groups
#           gives the same result.
#   dcNUMT  any detectable organellar coverage (fraction > 0). Mitochondrial
#           candidate contigs formed a continuum (0.28-0.97), reflecting the
#           recombinant structure of plant mitogenomes.
#
# Input : all_nupt_query_categories.tsv, all_numt_query_categories.tsv
#         (written by alien_nupt_detection.sh / alien_numt_detection.sh)
#         dcnuxt_unique.tsv (02_cnuxt_kaks/04_build_funnel_table.py)
# Output: all_{nupt,numt}_query_categories_final.tsv
#             category tables restricted to the final unique dcNUxT set
#             (species + Sequence without .mafft = species + qseqid)
#         organellar_fractions.tsv   class, species, contig, length, fraction
#
# Usage (from the working directory): bash nonhost_scaffold_composition.sh
# ============================================================
set -uo pipefail

BASE_DIR="${BASE_DIR:-$PWD}"
CSV="$BASE_DIR/combined_dataset.csv"
PLASTID_DB="/vol/storage/pt_database/plastid_db"
MITO_DB="/vol/storage/mt_database/mito_db"
THREADS="${THREADS:-6}"
NUCLEAR_LEN=300000
cd "$BASE_DIR" || exit 1
TMP=$(mktemp -d)
UNIQUE=dcnuxt_unique.tsv
[[ -s $UNIQUE ]] || { echo "missing $UNIQUE" >&2; exit 1; }

# restrict the category tables to the final unique dcNUxT set
for set in nupt numt; do
    [[ -s all_${set}_query_categories.tsv ]] || continue
    awk -F'\t' -v c="$set" '
        FNR == 1 { if (NR == 1) for (i = 1; i <= NF; i++) h[$i] = i; else print; next }
        NR == FNR { if (toupper($h["class"]) == toupper(c)) { q = $h["Sequence"]; sub(/\.mafft$/, "", q); k[$h["Species"] "\t" q] = 1 } next }
        ($1 "\t" $2) in k' "$UNIQUE" all_${set}_query_categories.tsv > all_${set}_query_categories_final.tsv
done

for set in nupt numt; do
    cat=all_nupt_query_categories_final.tsv; db=$PLASTID_DB
    [[ $set == numt ]] && { cat=all_numt_query_categories_final.tsv; db=$MITO_DB; }
    [[ -s $cat ]] || { echo "missing $cat" >&2; continue; }
    # contig of each candidate: query ID = <contig>_<start>-<end>_<gene>[_REV]
    awk -F'\t' 'NR > 1 && $3 == "HIGH_CONF_NONHOST" {s = $2; sub(/_[0-9]+-[0-9]+_[^_]+(_REV)?$/, "", s); print $1 "\t" s}' "$cat" \
      | sort -u |
    while IFS=$'\t' read -r sp scaf; do
        a=$(tr -d '\r' < "$CSV" | awk -F',' -v s="$sp" '$1 == s {print $2; exit}')
        [[ -z "$a" ]] && { printf '%s\t%s\t%s\tNA\tNA\n' "$set" "$sp" "$scaf"; continue; }
        g=""
        for c in downloads/genome_fasta/${a}*.fna multi_organelles/downloads/genome_fasta/${a}*.fna; do
            [[ -s "$c" ]] || continue
            [[ -s "$c.fai" ]] || samtools faidx "$c" 2>/dev/null
            awk -v s="$scaf" '$1 == s {f = 1; exit} END {exit !f}' "$c.fai" 2>/dev/null && { g="$c"; break; }
        done
        [[ -z "$g" ]] && { printf '%s\t%s\t%s\tNA\tNA\n' "$set" "$sp" "$scaf"; continue; }
        L=$(awk -v s="$scaf" '$1 == s {print $2; exit}' "$g.fai")
        if (( L > NUCLEAR_LEN )); then printf '%s\t%s\t%s\t%s\t0.00\n' "$set" "$sp" "$scaf" "$L"; continue; fi
        samtools faidx "$g" "$scaf" > "$TMP/contig.fa" 2>/dev/null || continue
        C=$(blastn -task dc-megablast -query "$TMP/contig.fa" -db "$db" -outfmt '6 qstart qend' \
                -max_target_seqs 1 -evalue 1e-20 -num_threads "$THREADS" 2>/dev/null \
            | sort -k1,1n | awk '{if ($1 > e) {s += e - b; b = $1; e = $2} else if ($2 > e) e = $2} END {s += e - b; print s + 0}')
        awk -v st="$set" -v s="$sp" -v x="$scaf" -v l="$L" -v c="$C" \
            'BEGIN {printf "%s\t%s\t%s\t%s\t%.2f\n", st, s, x, l, (l ? c / l : 0)}'
    done
done > organellar_fractions.tsv
rm -rf "$TMP"

echo "wrote organellar_fractions.tsv; next: python3 06_nonhost/nonhost_summary.py"
