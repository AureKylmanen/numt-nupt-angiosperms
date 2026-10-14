#!/bin/bash
# ============================================================
# Script Name: mtpt_coverage.sh
# Description: Fraction of each plastome that is also present in the species'
#              own mitogenome (mitochondrial plastid DNA, MTPT). Sequence of this
#              kind is detected by both the NUMT and the NUPT query, so it
#              explains NUMT-NUPT overlap (numt_nupt_shared.bed).
# Author: Aure Kylmänen
#
# Plastome and mitogenome are the accessions in the species CSV, read from the
# same locations as numt_nupt_detection.sh. Multi-chromosome mitogenomes
# (accessions separated by ';') are concatenated.
#
# blastn plastome vs mitogenome (word_size 20, e <= 1e-10); hit intervals on the
# plastome are merged, so an MTPT in the inverted repeat covers both IR copies.
#
# Output: species/mtpt_coverage.tsv
#   species, plastid_acc, plastid_len, n_mito_seqs, mito_len, mtpt_cov_bp, mtpt_frac,
#   longest_hit_bp, longest_hit_pid, bp_in_hits_ge10kb_ge99.5
#   (the last column: plastome bp in long, near-identical blocks)
#         species/mtpt_errors.log
#
# Usage: WORKING_DIR=/path bash mtpt_coverage.sh combined_dataset.csv
# ============================================================
set -uo pipefail
W="${WORKING_DIR:-${BASE_DIR:-$PWD}}"
CSV="${1:-$W/combined_dataset.csv}"
DL=$W/downloads
MO=$W/multi_organelles/downloads
OUT=$W/species/mtpt_coverage.tsv
LOG=$W/species/mtpt_errors.log
export LC_ALL=C
TMP=$(mktemp -d)

find_fa() {   # $1 = plastid|mito, $2 = accession
    local k=$1 a=$2 f
    for f in "$DL/${k}_fasta/$a.fasta" "$MO/${k}_fasta/$a.fasta"; do [[ -s $f ]] && { echo "$f"; return; }; done
    f=$(ls "$DL/${k}_fasta/${a%.*}".[0-9]*.fasta "$MO/${k}_fasta/${a%.*}".[0-9]*.fasta 2>/dev/null | head -1)
    [[ -s $f ]] && echo "$f"
}

printf "species\tplastid_acc\tplastid_len\tn_mito_seqs\tmito_len\tmtpt_cov_bp\tmtpt_frac\tlongest_hit_bp\tlongest_hit_pid\tbp_in_hits_ge10kb_ge99.5\n" > $OUT
: > $LOG
while IFS=, read -r s asm pl mt <&3; do
    [[ -z $s || $s == SPECIES ]] && continue
    PF=$(find_fa plastid "$pl")
    : > $TMP/mito.fa; miss=""
    for a in ${mt//;/ }; do
        f=$(find_fa mito "$a"); [[ -n $f ]] && cat "$f" >> $TMP/mito.fa || miss="$miss $a"
    done
    if [[ -z $PF || -n $miss || ! -s $TMP/mito.fa ]]; then
        echo "$s: plastid=${PF:-MISSING($pl)} mito_missing:${miss:- none}" >> $LOG
        printf "%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n" "$s" "$pl" >> $OUT; continue
    fi
    plen=$(grep -v '>' "$PF" | tr -d '\n\r' | wc -c)
    nm=$(grep -c '>' $TMP/mito.fa); mlen=$(grep -v '>' $TMP/mito.fa | tr -d '\n\r' | wc -c)
    timeout 900 blastn -query "$PF" -subject $TMP/mito.fa -word_size 20 -evalue 1e-10 \
        -outfmt "6 qseqid qstart qend length pident" > $TMP/hits.tsv 2>> $LOG
    merged_bp() { awk -v OFS='\t' '{a = $2 < $3 ? $2 : $3; b = $2 < $3 ? $3 : $2; print $1, a - 1, b}' \
                  | sort -k1,1 -k2,2n | bedtools merge 2>/dev/null | awk '{s += $3 - $2} END {print s + 0}'; }
    cov=$(merged_bp < $TMP/hits.tsv)
    read lh lp < <(sort -k4,4nr $TMP/hits.tsv | head -1 | awk '{print $4 + 0, $5 + 0}'); lh=${lh:-0}; lp=${lp:-0}
    big=$(awk '$4 >= 10000 && $5 >= 99.5' $TMP/hits.tsv | merged_bp)
    awk -v s="$s" -v p="$pl" -v pl=$plen -v nm=$nm -v ml=$mlen -v c=$cov -v lh=$lh -v lp=$lp -v bg=$big \
        'BEGIN {printf "%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%d\t%.1f\t%d\n", s, p, pl, nm, ml, c, (pl > 0 ? c / pl : 0), lh, lp, bg}' >> $OUT
done 3< <(tr -d '\r' < "$CSV")
rm -rf $TMP
echo "done: $(($(wc -l < $OUT) - 1)) species, $(grep -c . $LOG) problems (see $LOG)"
