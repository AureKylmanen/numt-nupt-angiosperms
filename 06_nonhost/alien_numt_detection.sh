#!/usr/bin/env bash

# =============================================================================
# Script Name: alien_numt_detection.sh
# Description: Detect alien NUMTs for species with downloaded reference genomes using mitogenome database
# Author: Aure Kylmänen
# Created: 20.04.2026
#
# Purpose:
#   Identify candidate non-host (alien) DNA mitochondrial nuclear insertions
#   (dcNUMTs) by comparing each candidate against:
#     1. a pooled mitochondrial-genome BLAST database, and
#     2. every available mitogenome for the candidate's host species.
#
#   Host identity is measured directly from the host mitogenomes rather than
#   inferred from the absence of a host accession in the pooled database.
#
# Requirements:
#   - Bash 4+
#   - NCBI BLAST+ (blastn)
#   - awk, sed, grep, find, sort, cut, tr, wc, xargs
#   - A completed NUMT Ka/Ks analysis
#   - MAFFT alignment files produced by the upstream NUMT pipeline
#   - A BLAST database containing the mitochondrial genomes to be searched
#
# Inputs:
#   - Species metadata CSV with columns:
#       SPECIES,ASSEMBLY,PLASTID,MITO
#   - Per-species NUMT Ka/Ks tables
#   - Per-species .mafft alignment files
#   - Host mitogenome FASTA files
#   - Pooled mitochondrial BLAST database
#
# Outputs:
#   - Combined BLAST hit table
#   - Per-query classification table
#   - Per-species summary table
#   - Per-species intermediate files and logs
# =============================================================================

set -euo pipefail

############################################
# CONFIG
############################################

BASE_DIR="${BASE_DIR:-$PWD}"   # working directory; override with BASE_DIR=...
export BASE_DIR
SP_DIR="$BASE_DIR/species"
SPECIES_LIST="$BASE_DIR/combined_dataset.csv"        # SPECIES,ASSEMBLY,PLASTID,MITO

# Pooled mitochondrial genome BLAST database.
MITO_DB="/vol/storage/mt_database/mito_db"

# Host mitogenome FASTAs. Primary tree uses versioned filenames, the multi-mito
# tree unversioned -- both are searched, both forms tried.
HOST_MITO_DIR="$BASE_DIR/downloads/mito_fasta"
HOST_MITO_DIR2="$BASE_DIR/multi_organelles/downloads/mito_fasta"

THREADS=26
KS_THRESHOLD=0.1          # Ks floor defining dcNUMTs
IDENTITY_THRESHOLD=95     # Identity threshold used for classification (%)
HIGH_CONF_DIFF=5          # Non-host identity margin required for high confidence (pp)
MIN_QCOV=80               # Minimum query coverage required for classification
MIN_SEQ_LEN=200           # Minimum candidate sequence length in bp
BLAST_TASK="dc-megablast" # SAME task for both searches
EVALUE="1e-10"
MAX_TARGETS=500

RESUME="${RESUME:-0}"

# Upstream NUMT Ka/Ks analysis settings.
RESULTS_SUBDIR="results"
ALIGN_SUBDIR="alignments/mito"
KAKS_BASENAME_FMT="%s_numt_kaks.csv"         # as written by 02_cnuxt_kaks/numt_gene_alignment.sh

COMBINED_OUT="$BASE_DIR/all_numt_blast_hits.tsv"
SUMMARY_OUT="$BASE_DIR/numt_summary_table.tsv"
CATEGORIES_OUT="$BASE_DIR/all_numt_query_categories.tsv"

############################################
# HELPERS
############################################

# Escape a string for safe use as a basic regular-expression fragment.
escape_regex() {
    printf '%s' "$1" | sed -e 's/[][\.*^$(){}?+|\\\/]/\\&/g'
}

init_outputs() {
    [[ -f "$COMBINED_OUT" ]] || printf 'species\tqseqid\tsseqid\tpident\tqcovhsp\tlength\tmismatch\tgapopen\tevalue\tbitscore\tstitle\n' > "$COMBINED_OUT"
    [[ -f "$SUMMARY_OUT" ]] || printf 'species\tstatus\ttotal_candidates\thost\tnonhost_low\thigh_conf_nonhost\tnonhost_host_absent\tno_hit\tn_host_mitogenomes\n' > "$SUMMARY_OUT"
    [[ -f "$CATEGORIES_OUT" ]] || printf 'species\tqseqid\tcategory\tbest_pident\tbest_qcov\tbest_bitscore\tbest_sseqid\tbest_stitle\thost_pident\thost_qcov\tmargin\thost_acc_matched\n' > "$CATEGORIES_OUT"
}

purge_species_rows() {
    local sp="$1" f
    for f in "$COMBINED_OUT" "$SUMMARY_OUT" "$CATEGORIES_OUT"; do
        [[ -f "$f" ]] || continue
        awk -F'\t' -v sp="$sp" 'NR==1 || $1 != sp' "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"
    done
}

emit_summary() {
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$@" >> "$SUMMARY_OUT"
}

# Resolve one accession to a FASTA across both trees, versioned or not.
resolve_mito_fasta() {
    local acc="$1" base="${1%%.*}" d
    for d in "$HOST_MITO_DIR" "$HOST_MITO_DIR2"; do
        [[ -d "$d" ]] || continue
        [[ -s "$d/${acc}.fasta"  ]] && { printf '%s' "$d/${acc}.fasta";  return 0; }
        [[ -s "$d/${base}.fasta" ]] && { printf '%s' "$d/${base}.fasta"; return 0; }
        local c
        c=$(LC_ALL=C find "$d" -maxdepth 1 -name "${base}.[0-9]*.fasta" | sort | head -n1)
        [[ -n "$c" ]] && { printf '%s' "$c"; return 0; }
    done
    return 1
}

############################################
# WORKFLOW
############################################
#
# For each species the pipeline:
#   1. Reads the host mitogenome accessions from the metadata CSV.
#   2. Resolves every host accession to an available FASTA file.
#   3. Selects candidate dcNUMTs using the configured Ka/Ks threshold.
#   4. Extracts candidate sequences from the upstream alignments.
#   5. Searches candidates against the pooled mitochondrial BLAST database.
#   6. Searches candidates independently against every host mitogenome.
#   7. Classifies each query using host identity, non-host identity, coverage,
#      and the configured identity margin.
#   8. Writes per-species and combined summary tables.
#
# Run one species:
#   ./alien_numt_detection.sh "Species_directory"
#
# Run all species:
#   ./alien_numt_detection.sh
#
# Resume a previously started species:
#   RESUME=1 ./alien_numt_detection.sh "Species_directory"
#
# Before publication/use, review the CONFIG block and set paths and thresholds
# appropriate for the dataset and analysis.
#
############################################

############################################
# MAIN PER-SPECIES ROUTINE
############################################

process_species() {
    local SPECIES_DIR="$1"
    local SPECIES
    SPECIES=$(basename "$SPECIES_DIR")

    echo "=============================================="
    echo "Processing species: $SPECIES"

    local CSV_LINE HOST_NAME MITO_FIELD HOST_RE
    CSV_LINE=$(awk -F',' -v sp="$SPECIES" '$1==sp {print; exit}' "$SPECIES_LIST")

    if [[ -z "$CSV_LINE" ]]; then
        echo "[ABORT] $SPECIES not in $SPECIES_LIST -- cannot define host."
        emit_summary "$SPECIES" "NO_CSV_ENTRY" 0 0 0 0 0 0 0
        return
    fi

    HOST_NAME=$(echo "$CSV_LINE" | cut -d',' -f1 | tr '_' ' ')
    MITO_FIELD=$(echo "$CSV_LINE" | cut -d',' -f4 | tr -d '[:space:]')   #

    if [[ -z "$MITO_FIELD" ]]; then
        echo "[ABORT] No mito accession in CSV for $SPECIES."
        emit_summary "$SPECIES" "NO_MITO_ACC" 0 0 0 0 0 0 0
        return
    fi

    # Split semicolon-separated accessions; resolve each to a FASTA.
    local -a HOST_ACCS=() HOST_FASTAS=() MISSING_ACCS=()
    local _a
    while IFS= read -r _a; do
        [[ -n "$_a" ]] || continue
        HOST_ACCS+=("$_a")
        local fa
        if fa=$(resolve_mito_fasta "$_a"); then
            HOST_FASTAS+=("$fa")
        else
            MISSING_ACCS+=("$_a")
        fi
    done < <(printf '%s\n' "$MITO_FIELD" | tr ';' '\n')
    # The process substitution supplies a terminating newline so the final
    # accession is also consumed by read.
    if (( ${#HOST_FASTAS[@]} == 0 )); then
        echo "[ABORT] No host mitogenome FASTA found for: ${HOST_ACCS[*]}"
        emit_summary "$SPECIES" "NO_HOST_FASTA" 0 0 0 0 0 0 "${#HOST_ACCS[@]}"
        return
    fi
    if (( ${#MISSING_ACCS[@]} > 0 )); then
        echo "  [WARN] ${#MISSING_ACCS[@]} of ${#HOST_ACCS[@]} host mitogenomes missing: ${MISSING_ACCS[*]}"
        echo "         Host identity may be underestimated for fragments derived from them."
    fi

    HOST_RE=$(escape_regex "$HOST_NAME")

    local RESULTS_DIR="$SPECIES_DIR/$RESULTS_SUBDIR"
    local ALIGN_DIR="$SPECIES_DIR/$ALIGN_SUBDIR"
    local OUT_DIR="$SPECIES_DIR/dcnumt_pipeline"

    [[ "$RESUME" != "1" ]] && rm -rf "$OUT_DIR"
    mkdir -p "$OUT_DIR/logs"

    if [[ "$RESUME" == "1" && -f "$OUT_DIR/query_categories.tsv" ]]; then
        echo "[SKIP] RESUME=1 and $SPECIES already processed."
        return
    fi

    purge_species_rows "$SPECIES"

    echo "  host organism    : $HOST_NAME"
    echo "  host mitogenomes : ${#HOST_FASTAS[@]} (${HOST_ACCS[*]})"
    echo "  identity floor   : ${IDENTITY_THRESHOLD}%   margin: ${HIGH_CONF_DIFF}pp   qcov: ${MIN_QCOV}%"

    # ------------------------------------------------------------------
    # [1] Ka/Ks table -> dcNUMT ids
    # ------------------------------------------------------------------
    local KAKS_FILE
    KAKS_FILE="$RESULTS_DIR/$(printf "$KAKS_BASENAME_FMT" "$SPECIES")"

    if [[ ! -s "$KAKS_FILE" ]]; then
        echo "[SKIP] No NUMT Ka/Ks file: $KAKS_FILE"
        emit_summary "$SPECIES" "NO_KAKS_FILE" 0 0 0 0 0 0 "${#HOST_FASTAS[@]}"
        return
    fi

    local NEWEST_ALN
    NEWEST_ALN=$(find "$ALIGN_DIR" -maxdepth 1 -name '*.mafft' -newer "$KAKS_FILE" -print -quit 2>/dev/null || true)
    if [[ -d "$ALIGN_DIR" && -z "$NEWEST_ALN" ]]; then
        if find "$ALIGN_DIR" -maxdepth 1 -name '*.mafft' -print -quit 2>/dev/null | grep -q .; then
            echo "[WARN] All .mafft in $ALIGN_DIR predate $KAKS_FILE; alignment files may be older than the Ka/Ks analysis." \
                | tee -a "$OUT_DIR/logs/staleness.log"
        fi
    fi

    echo "[1] Filtering Ks > $KS_THRESHOLD ..."
    awk -F'\t' -v ks="$KS_THRESHOLD" '
        NR==1 { next }
        $4 != "NA" && $4 != "" && $4 != "-" && ($4+0) > (ks+0) { print $1 }
    ' "$KAKS_FILE" | sort -u > "$OUT_DIR/dcnumt_ids.txt"

    local NUM
    NUM=$(wc -l < "$OUT_DIR/dcnumt_ids.txt")
    if [[ "$NUM" -eq 0 ]]; then
        echo "[SKIP] No sequences with Ks > $KS_THRESHOLD"
        emit_summary "$SPECIES" "NO_DCNUMTS" 0 0 0 0 0 0 "${#HOST_FASTAS[@]}"
        return
    fi
    echo "    dcNUMT ids: $NUM"

    # ------------------------------------------------------------------
    # [2] Extract fragment sequences: gene first, fragment second
    # ------------------------------------------------------------------
    echo "[2] Extracting fragment sequences ..."
    : > "$OUT_DIR/dcnumts_all.fa"
    : > "$OUT_DIR/logs/extraction.log"

    local n_extracted=0
    while read -r aln_id; do
        [[ -n "$aln_id" ]] || continue
        local f="$ALIGN_DIR/$aln_id"
        [[ "$f" == *.mafft ]] || f="${f}.mafft"

        if [[ ! -f "$f" ]]; then
            printf 'MISSING_ALIGNMENT\t%s\n' "$aln_id" >> "$OUT_DIR/logs/extraction.log"
            continue
        fi
        local nrec
        nrec=$(grep -c '^>' "$f" || true)
        # Record 1 is the gene, record 2 the fragment. Alignments with a third
        # record (the reverse-orientation fragment) still yield a valid pair.
        if [[ "$nrec" -lt 2 ]]; then
            printf 'TOO_FEW_RECORDS\t%s\t%s\n' "$aln_id" "$nrec" >> "$OUT_DIR/logs/extraction.log"
            continue
        fi
        if [[ "$nrec" -gt 2 ]]; then
            printf 'MULTI_FRAGMENT_USING_RECORD2\t%s\t%s\n' "$aln_id" "$nrec" >> "$OUT_DIR/logs/extraction.log"
        fi

        awk -v qid="$(basename "$f" .mafft)" -v minlen="$MIN_SEQ_LEN" '
            /^>/ { n++; next }
            n==2 { frag = frag $0 }
            END {
                gsub(/[-\r]/, "", frag)
                if (length(frag) >= minlen) { print ">" qid; print toupper(frag) }
            }
        ' "$f" >> "$OUT_DIR/dcnumts_all.fa"
        n_extracted=$((n_extracted+1))
    done < "$OUT_DIR/dcnumt_ids.txt"

    local SEQ_COUNT UNIQ_IDS
    SEQ_COUNT=$(grep -c '^>' "$OUT_DIR/dcnumts_all.fa" || true)
    UNIQ_IDS=$(grep '^>' "$OUT_DIR/dcnumts_all.fa" | sort -u | wc -l || true)

    if [[ "$SEQ_COUNT" -ne "$UNIQ_IDS" ]]; then
        echo "[ABORT] Duplicate query IDs: $SEQ_COUNT sequences, $UNIQ_IDS unique."
        emit_summary "$SPECIES" "DUPLICATE_QUERY_IDS" "$SEQ_COUNT" 0 0 0 0 0 "${#HOST_FASTAS[@]}"
        return
    fi
    if [[ "$SEQ_COUNT" -eq 0 ]]; then
        echo "[SKIP] No sequences >= ${MIN_SEQ_LEN} bp"
        emit_summary "$SPECIES" "NO_SEQS_PASSING_LENGTH" 0 0 0 0 0 0 "${#HOST_FASTAS[@]}"
        return
    fi
    echo "    queries: $SEQ_COUNT (from $n_extracted alignments)"

    # ------------------------------------------------------------------
    # [3] BLAST A: pooled mitochondrial database
    # ------------------------------------------------------------------
    echo "[3] BLAST A: pooled mitochondrial database ..."
    blastn \
        -task "$BLAST_TASK" \
        -query "$OUT_DIR/dcnumts_all.fa" \
        -db "$MITO_DB" \
        -out "$OUT_DIR/blast_db_raw.tsv" \
        -outfmt "6 qseqid sseqid pident qcovhsp length mismatch gapopen evalue bitscore stitle" \
        -max_target_seqs "$MAX_TARGETS" \
        -max_hsps 1 \
        -evalue "$EVALUE" \
        -num_threads "$THREADS"

    # ------------------------------------------------------------------
    # [4] BLAST B: EVERY host mitogenome, best identity per query
    # ------------------------------------------------------------------
    echo "[4] BLAST B: ${#HOST_FASTAS[@]} host mitogenome(s) ..."
    : > "$OUT_DIR/blast_host_raw.tsv"
    local hf hacc
    for hf in "${HOST_FASTAS[@]}"; do
        hacc=$(basename "$hf" .fasta)
        blastn \
            -task "$BLAST_TASK" \
            -query "$OUT_DIR/dcnumts_all.fa" \
            -subject "$hf" \
            -outfmt "6 qseqid sseqid pident qcovhsp length mismatch gapopen evalue bitscore" \
            -max_hsps 1 \
            -evalue "$EVALUE" 2>/dev/null \
        | awk -v a="$hacc" 'BEGIN{OFS="\t"} { print $0, a }' \
        >> "$OUT_DIR/blast_host_raw.tsv" || true
    done

    # Best host HSP per query across ALL host mitogenomes, ranked by bitscore.
    awk -F'\t' -v cov="$MIN_QCOV" 'BEGIN{OFS="\t"}
        ($4+0) >= (cov+0) {
            q=$1
            if (!(q in bs) || ($9+0) > bs[q]) { bs[q]=$9+0; pid[q]=$3+0; qc[q]=$4+0; acc[q]=$10 }
        }
        END { for (q in bs) print q, pid[q], qc[q], bs[q], acc[q] }
    ' "$OUT_DIR/blast_host_raw.tsv" | sort -k1,1 > "$OUT_DIR/host_identity.tsv"

    echo "    queries with a host alignment: $(wc -l < "$OUT_DIR/host_identity.tsv") / $SEQ_COUNT"

    # ------------------------------------------------------------------
    # [5] Tag database hits as host / non-host, matching ANY host accession
    # ------------------------------------------------------------------
    echo "[5] Partitioning database hits ..."
    printf '%s\n' "${HOST_ACCS[@]}" | sed 's/\..*//' | sort -u > "$OUT_DIR/host_accs.txt"

    awk -F'\t' -v hostre="$HOST_RE" -v cov="$MIN_QCOV" '
        BEGIN { OFS="\t" }
        FNR==NR { hacc[$1]=1; next }
        ($4+0) < (cov+0) { next }
        {
            is_host = 0
            for (a in hacc) if (index($2, a) > 0) { is_host = 1; break }
            if (!is_host && hostre != "" && \
                $10 ~ ("(^|[[:space:]|])" hostre "([[:space:]]|$|,|\\.)")) is_host = 1
            print is_host, $0
        }
    ' "$OUT_DIR/host_accs.txt" "$OUT_DIR/blast_db_raw.tsv" > "$OUT_DIR/blast_db_tagged.tsv"

    echo "    database hits attributed to host: $(awk -F'\t' '$1==1' "$OUT_DIR/blast_db_tagged.tsv" | wc -l)"

    awk -F'\t' 'BEGIN{OFS="\t"} { $1=""; sub(/^\t/,""); print }' "$OUT_DIR/blast_db_tagged.tsv" \
        | awk -v sp="$SPECIES" 'BEGIN{OFS="\t"} { print sp, $0 }' >> "$COMBINED_OUT"

    # ------------------------------------------------------------------
    # [6] Classify
    # ------------------------------------------------------------------
    echo "[6] Classifying ..."
    touch "$OUT_DIR/host_identity.tsv" "$OUT_DIR/blast_db_tagged.tsv"

    awk -F'\t' -v sp="$SPECIES" -v id="$IDENTITY_THRESHOLD" -v diff="$HIGH_CONF_DIFF" \
        -v hostfile="$OUT_DIR/host_identity.tsv" \
        -v idsfile="$OUT_DIR/dcnumts_all.fa" \
        -v dbfile="$OUT_DIR/blast_db_tagged.tsv" '
        BEGIN { OFS="\t" }
        FILENAME == hostfile { hpid[$1]=$2+0; hqc[$1]=$3+0; hacc[$1]=$5; hseen[$1]=1; next }
        FILENAME == idsfile  { if ($0 ~ /^>/) allq[substr($0,2)]=1; next }
        FILENAME == dbfile {
            is_host=$1; q=$2; sid=$3; pid=$4+0; qc=$5+0; bs=$10+0; sti=$11
            if (pid < (id+0)) next
            if (is_host) {
                if (!(q in hostbs) || bs > hostbs[q]) { hostbs[q]=bs; hostpid_db[q]=pid }
            } else if (!(q in nhbs) || bs > nhbs[q]) {
                nhbs[q]=bs; nhpid[q]=pid; nhqc[q]=qc; nhsid[q]=sid; nhsti[q]=sti
            }
        }
        END {
            for (q in allq) {
                hp="NA"; hq="NA"; ha="NA"; have_host=0
                if (q in hseen)           { hp=hpid[q]; hq=hqc[q]; ha=hacc[q]; have_host=1 }
                else if (q in hostpid_db) { hp=hostpid_db[q]; ha="db_record";  have_host=1 }

                have_nh = (q in nhbs)
                host_ok = (have_host && hp != "NA" && (hp+0) >= (id+0))

                if (!have_nh && !host_ok) {
                    print sp,q,"NO_HIT","NA","NA","NA","NA","NA",hp,hq,"NA",ha; continue
                }
                if (host_ok && (!have_nh || (hp+0) >= nhpid[q])) {
                    print sp,q,"HOST",hp,hq,"NA","host_mitogenome","host_mitogenome",hp,hq,"NA",ha
                    continue
                }
                if (!have_host) {
                    print sp,q,"NONHOST_HOST_ABSENT",nhpid[q],nhqc[q],nhbs[q],nhsid[q],nhsti[q],"NA","NA","NA","NA"
                    continue
                }
                margin = nhpid[q] - hp
                cat = (margin >= (diff+0)) ? "HIGH_CONF_NONHOST" : "NONHOST_LOW_CONF"
                printf "%s\t%s\t%s\t%.3f\t%.1f\t%.1f\t%s\t%s\t%.3f\t%s\t%.3f\t%s\n", \
                       sp,q,cat,nhpid[q],nhqc[q],nhbs[q],nhsid[q],nhsti[q],hp,hq,margin,ha
            }
        }
    ' "$OUT_DIR/host_identity.tsv" "$OUT_DIR/dcnumts_all.fa" "$OUT_DIR/blast_db_tagged.tsv" \
      > "$OUT_DIR/query_categories.tsv"

    cat "$OUT_DIR/query_categories.tsv" >> "$CATEGORIES_OUT"

    # ------------------------------------------------------------------
    # [7] Summary
    # ------------------------------------------------------------------
    local N_HOST N_LOW N_HC N_ABSENT N_NOHIT TOTAL
    N_HOST=$(awk -F'\t'   '$3=="HOST"'                "$OUT_DIR/query_categories.tsv" | wc -l)
    N_LOW=$(awk -F'\t'    '$3=="NONHOST_LOW_CONF"'    "$OUT_DIR/query_categories.tsv" | wc -l)
    N_HC=$(awk -F'\t'     '$3=="HIGH_CONF_NONHOST"'   "$OUT_DIR/query_categories.tsv" | wc -l)
    N_ABSENT=$(awk -F'\t' '$3=="NONHOST_HOST_ABSENT"' "$OUT_DIR/query_categories.tsv" | wc -l)
    N_NOHIT=$(awk -F'\t'  '$3=="NO_HIT"'              "$OUT_DIR/query_categories.tsv" | wc -l)

    TOTAL=$((N_HOST + N_LOW + N_HC + N_ABSENT + N_NOHIT))
    (( TOTAL != SEQ_COUNT )) && echo "[WARN] categories ($TOTAL) != queries ($SEQ_COUNT)" \
        | tee -a "$OUT_DIR/logs/consistency.log"

    emit_summary "$SPECIES" "OK" "$SEQ_COUNT" "$N_HOST" "$N_LOW" "$N_HC" "$N_ABSENT" "$N_NOHIT" "${#HOST_FASTAS[@]}"

    echo "  total candidates      : $SEQ_COUNT"
    echo "  host                  : $N_HOST"
    echo "  non-host, low conf    : $N_LOW"
    echo "  non-host, HIGH CONF   : $N_HC"
    echo "  non-host, host absent : $N_ABSENT"
    echo "  no hit (<${IDENTITY_THRESHOLD}%)        : $N_NOHIT"
    echo "Done: $SPECIES"
    echo
}

############################################
# MAIN
############################################

if [[ ! -f "${MITO_DB}.nsq" && ! -f "${MITO_DB}.00.nsq" && ! -f "${MITO_DB}.ndb" ]]; then
    echo "ERROR: mitochondrial BLAST database not found at $MITO_DB" >&2
    echo "       Set MITO_DB in the CONFIG block (see)." >&2
    exit 1
fi

init_outputs

echo "=============================================="
echo "Reading dcNUMTs from : */$RESULTS_SUBDIR/$(printf "$KAKS_BASENAME_FMT" '<SPECIES>')"
echo "Reading alignments   : */$ALIGN_SUBDIR/"
echo "Mito database        : $MITO_DB"
_pat=$(printf "$KAKS_BASENAME_FMT" '*')
_n_res=$(find "$SP_DIR" -maxdepth 3 -path "*/$RESULTS_SUBDIR/$_pat" 2>/dev/null | wc -l)
_n_dc=$(find "$SP_DIR" -maxdepth 3 -path "*/$RESULTS_SUBDIR/$_pat" 2>/dev/null \
        | xargs -r awk -F'\t' 'FNR>1 && $4!="NA" && $4!="" && ($4+0)>0.1' 2>/dev/null | wc -l)
echo "Species with a NUMT Ka/Ks table : $_n_res"
echo "Total rows with Ks > 0.1        : $_n_dc"
echo "=============================================="
echo

if [[ $# -eq 1 ]]; then
    process_species "$1"
else
    for dir in "$SP_DIR"/*; do
        [[ -d "$dir/$RESULTS_SUBDIR" ]] || continue
        process_species "$dir"
    done
fi

echo "Pipeline complete."
echo "  Combined BLAST hits : $COMBINED_OUT"
echo "  Per-query categories: $CATEGORIES_OUT"
echo "  Summary table       : $SUMMARY_OUT"
