#!/bin/bash

# ================================================================
# Script Name: alien_nupt_detection.sh
# Description: Detect non-host (alien) dcNUPTs using TWO BLAST searches:
#                (a) against the pooled plastid genome database
#                (b) against the host's OWN plastid genome
#              so that host identity is always measured directly and is
#              never inferred from absence in the database hit list.
# Author:      Aure Kylmanen
# Created:     20.04.2026
# Revised:     24.08.2026
#
# This script identifies candidate divergent plastid-to-nuclear insertions
# (dcNUPTs) from an upstream Ka/Ks analysis and compares each candidate with
# a pooled plastid database and the host species' plastid genome.
#
# The host comparison is performed independently of the pooled database.
# This is important because a host sequence may not occur among the top
# database hits even when the candidate is genuinely host-derived.
#
# Inputs:
#   - Species metadata CSV: SPECIES,ASSEMBLY,PLASTID,MITO
#   - Per-species NUPT Ka/Ks tables
#   - Per-species MAFFT alignment files
#   - Host plastid FASTA files
#   - Pooled plastid BLAST database
#
# Outputs:
#   - Combined BLAST hit table
#   - Per-query classification table
#   - Per-species summary table
#   - Per-species intermediate files and logs
#
# Dependencies:
#   - Bash 4+
#   - NCBI BLAST+ (blastn)
#   - awk, sed, grep, find, sort, cut, tr, wc, xargs
#
# Run one species:
#   ./alien_nupt_detection.sh "Species_directory"
#
# Run all species:
#   ./alien_nupt_detection.sh
#
# Resume an existing species result:
#   RESUME=1 ./alien_nupt_detection.sh "Species_directory"
# ================================================================

set -euo pipefail

############################################
# CONFIG
############################################

BASE_DIR="/vol/storage/5T_disk/organelles"
SP_DIR="$BASE_DIR/species"
SPECIES_LIST="$BASE_DIR/combined_dataset.csv"           # SPECIES,ASSEMBLY,PLASTID,MITO
PLASTID_DB="/vol/storage/pt_database/plastid_db"        # Pooled plastid BLAST database
HOST_PLASTID_DIR="$BASE_DIR/downloads/plastid_fasta"    # Host plastid FASTA files
HOST_PLASTID_DIR2="$BASE_DIR/multi_organelles/downloads/plastid_fasta"

THREADS=26
KS_THRESHOLD=0.1          # Minimum Ks used to select dcNUPTs
IDENTITY_THRESHOLD=95     # Minimum identity used for classification (%)
HIGH_CONF_DIFF=5          # Minimum non-host identity advantage for high confidence (pp)
MIN_QCOV=80               # Minimum query coverage per HSP (%)
MIN_SEQ_LEN=200           # Minimum ungapped candidate length (bp)
BLAST_TASK="dc-megablast" # BLAST task used for both searches
EVALUE="1e-10"
MAX_TARGETS=500

RESUME="${RESUME:-0}"

# Source of the NUPT candidates. This filename pattern must match the
# Ka/Ks table produced by the upstream NUPT analysis.
RESULTS_SUBDIR="results"
ALIGN_SUBDIR="alignments/plastid"
KAKS_BASENAME_FMT="%s_nupt_kaks.csv"           # as written by 02_cnuxt_kaks/nupt_gene_alignment.sh

COMBINED_OUT="$BASE_DIR/all_nupt_blast_hits.tsv"
SUMMARY_OUT="$BASE_DIR/nupt_summary_table.tsv"
CATEGORIES_OUT="$BASE_DIR/all_nupt_query_categories.tsv"

############################################
# HELPERS
############################################

# Escape regex metacharacters so an organism name can safely be used as a
# secondary match against a BLAST subject title.
escape_regex() {
    printf '%s' "$1" | sed -e 's/[][\.*^$(){}?+|\\\/]/\\&/g'
}

# Create project-level output files and headers if they do not already exist.
init_outputs() {
    [[ -f "$COMBINED_OUT" ]] || printf 'species\tqseqid\tsseqid\tpident\tqcovhsp\tlength\tmismatch\tgapopen\tevalue\tbitscore\tstitle\n' > "$COMBINED_OUT"
    [[ -f "$SUMMARY_OUT" ]] || printf 'species\tstatus\ttotal_candidates\thost\tnonhost_low\thigh_conf_nonhost\tnonhost_host_absent\tno_hit\tdb_host_records\n' > "$SUMMARY_OUT"
    [[ -f "$CATEGORIES_OUT" ]] || printf 'species\tqseqid\tcategory\tbest_pident\tbest_qcov\tbest_bitscore\tbest_sseqid\tbest_stitle\thost_pident\thost_qcov\tmargin\n' > "$CATEGORIES_OUT"
}

# Remove existing rows for one species before reprocessing it. This makes
# reruns replace the species' previous results instead of duplicating them.
purge_species_rows() {
    local sp="$1" f
    for f in "$COMBINED_OUT" "$SUMMARY_OUT" "$CATEGORIES_OUT"; do
        [[ -f "$f" ]] || continue
        awk -F'\t' -v sp="$sp" 'NR==1 || $1 != sp' "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"
    done
}

# Append one per-species summary row.
emit_summary() {
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$@" >> "$SUMMARY_OUT"
}

############################################
# MAIN PER-SPECIES ROUTINE
############################################

process_species() {
    local SPECIES_DIR="$1"
    local SPECIES
    SPECIES=$(basename "$SPECIES_DIR")

    echo "=============================================="
    echo "Processing species: $SPECIES"

    # Resolve the host plastid accession and organism name from the metadata.
    local CSV_LINE HOST_ACC HOST_NAME HOST_ACC_BASE HOST_RE
    CSV_LINE=$(awk -F',' -v sp="$SPECIES" '$1==sp {print; exit}' "$SPECIES_LIST")

    if [[ -z "$CSV_LINE" ]]; then
        echo "[ABORT] $SPECIES not found in $SPECIES_LIST -- cannot define host."
        emit_summary "$SPECIES" "NO_CSV_ENTRY" 0 0 0 0 0 0 0
        return
    fi

    HOST_NAME=$(echo "$CSV_LINE" | cut -d',' -f1 | tr '_' ' ')
    HOST_ACC=$(echo "$CSV_LINE" | cut -d',' -f3)
    HOST_ACC_BASE="${HOST_ACC%%.*}"

    if [[ -z "$HOST_ACC" ]]; then
        echo "[ABORT] No plastid accession in CSV for $SPECIES."
        emit_summary "$SPECIES" "NO_PLASTID_ACC" 0 0 0 0 0 0 0
        return
    fi

    # Resolve both versioned and unversioned accession filenames in both
    # configured host FASTA directories. If multiple versioned files exist,
    # stop rather than selecting one arbitrarily.
    local HOST_FASTA="" _d
    for _d in "$HOST_PLASTID_DIR" "$HOST_PLASTID_DIR2"; do
        [[ -d "$_d" ]] || continue

        if [[ -s "$_d/${HOST_ACC}.fasta" ]]; then
            HOST_FASTA="$_d/${HOST_ACC}.fasta"
            break
        fi

        if [[ -s "$_d/${HOST_ACC_BASE}.fasta" ]]; then
            HOST_FASTA="$_d/${HOST_ACC_BASE}.fasta"
            break
        fi

        local -a cands=()
        while IFS= read -r c; do
            cands+=("$c")
        done < <(
            LC_ALL=C find "$_d" -maxdepth 1 -name "${HOST_ACC_BASE}.[0-9]*.fasta" | sort
        )

        if (( ${#cands[@]} == 1 )); then
            HOST_FASTA="${cands[0]}"
            break
        elif (( ${#cands[@]} > 1 )); then
            echo "[ABORT] Ambiguous host plastid FASTA for $HOST_ACC: ${cands[*]}"
            emit_summary "$SPECIES" "AMBIGUOUS_HOST_FASTA" 0 0 0 0 0 0 0
            return
        fi
    done

    if [[ -z "$HOST_FASTA" || ! -s "$HOST_FASTA" ]]; then
        echo "[ABORT] Host plastid FASTA missing for $HOST_ACC"
        emit_summary "$SPECIES" "NO_HOST_FASTA" 0 0 0 0 0 0 0
        return
    fi

    HOST_RE=$(escape_regex "$HOST_NAME")

    local RESULTS_DIR="$SPECIES_DIR/$RESULTS_SUBDIR"
    local ALIGN_DIR="$SPECIES_DIR/$ALIGN_SUBDIR"
    local OUT_DIR="$SPECIES_DIR/dcnupt_pipeline"

    # Clear stale per-species intermediates unless resume mode is explicitly
    # requested.
    if [[ "$RESUME" != "1" ]]; then
        rm -rf "$OUT_DIR"
    fi
    mkdir -p "$OUT_DIR/logs"

    if [[ "$RESUME" == "1" && -f "$OUT_DIR/query_categories.tsv" ]]; then
        echo "[SKIP] RESUME=1 and $SPECIES already processed."
        return
    fi

    purge_species_rows "$SPECIES"

    echo "  host organism : $HOST_NAME"
    echo "  host plastid  : $HOST_ACC"
    echo "  identity floor: ${IDENTITY_THRESHOLD}%   margin: ${HIGH_CONF_DIFF}pp   qcov: ${MIN_QCOV}%"

    # ------------------------------------------------------------------
    # [1] Locate the NUPT Ka/Ks table and select candidates.
    # ------------------------------------------------------------------
    local KAKS_FILE
    KAKS_FILE="$RESULTS_DIR/$(printf "$KAKS_BASENAME_FMT" "$SPECIES")"

    if [[ ! -s "$KAKS_FILE" ]]; then
        echo "[SKIP] No NUPT Ka/Ks file: $KAKS_FILE"
        emit_summary "$SPECIES" "NO_KAKS_FILE" 0 0 0 0 0 0 0
        return
    fi

    # Guard against accidentally supplying a mitochondrial NUMT table.
    case "$(basename "$KAKS_FILE")" in
        *numt*)
            echo "[ABORT] Resolved a NUMT table for a NUPT analysis: $KAKS_FILE"
            emit_summary "$SPECIES" "WRONG_ORGANELLE_TABLE" 0 0 0 0 0 0 0
            return
            ;;
    esac

    # Warn if all available alignments are older than the Ka/Ks table.
    local NEWEST_ALN
    NEWEST_ALN=$(find "$ALIGN_DIR" -maxdepth 1 -name '*.mafft' -newer "$KAKS_FILE" -print -quit 2>/dev/null || true)
    if [[ -d "$ALIGN_DIR" && -z "$NEWEST_ALN" ]]; then
        if find "$ALIGN_DIR" -maxdepth 1 -name '*.mafft' -print -quit 2>/dev/null | grep -q .; then
            echo "[WARN] All .mafft files in $ALIGN_DIR predate $KAKS_FILE." \
                 "Alignments may not correspond to the current Ka/Ks table." \
                 | tee -a "$OUT_DIR/logs/staleness.log"
        fi
    fi

    echo "[1] Filtering Ks > $KS_THRESHOLD ..."
    awk -F'\t' -v ks="$KS_THRESHOLD" '
        NR==1 { next }
        $4 != "NA" && $4 != "" && $4 != "-" && ($4+0) > (ks+0) { print $1 }
    ' "$KAKS_FILE" | sort -u > "$OUT_DIR/dcnupt_ids.txt"

    local NUM
    NUM=$(wc -l < "$OUT_DIR/dcnupt_ids.txt")
    if [[ "$NUM" -eq 0 ]]; then
        echo "[SKIP] No sequences with Ks > $KS_THRESHOLD"
        emit_summary "$SPECIES" "NO_DCNUPTS" 0 0 0 0 0 0 0
        return
    fi
    echo "    dcNUPT ids: $NUM"

    # ------------------------------------------------------------------
    # [2] Extract candidate fragment sequences.
    #
    # The upstream alignment is expected to contain exactly two records:
    # reference sequence first and candidate fragment second. The candidate
    # is therefore the second record. The .mafft basename is used as the
    # query ID so the BLAST result maps directly to the upstream identifier.
    # ------------------------------------------------------------------
    echo "[2] Extracting fragment sequences ..."
    : > "$OUT_DIR/dcnupts_all.fa"
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
        if [[ "$nrec" -ne 2 ]]; then
            printf 'UNEXPECTED_RECORD_COUNT\t%s\t%s\n' "$aln_id" "$nrec" >> "$OUT_DIR/logs/extraction.log"
            continue
        fi

        local qid
        qid=$(basename "$f" .mafft)

        awk -v qid="$qid" -v minlen="$MIN_SEQ_LEN" '
            /^>/ { n++; seq=""; next }
            n==2 { seq = seq $0 }
            END {
                gsub(/[-\r]/, "", seq)
                if (length(seq) >= minlen) {
                    print ">" qid
                    print toupper(seq)
                }
            }
        ' "$f" >> "$OUT_DIR/dcnupts_all.fa"

        n_extracted=$((n_extracted+1))
    done < "$OUT_DIR/dcnupt_ids.txt"

    local SEQ_COUNT UNIQ_IDS
    SEQ_COUNT=$(grep -c '^>' "$OUT_DIR/dcnupts_all.fa" || true)
    UNIQ_IDS=$(grep '^>' "$OUT_DIR/dcnupts_all.fa" | sort -u | wc -l || true)

    if [[ "$SEQ_COUNT" -ne "$UNIQ_IDS" ]]; then
        echo "[ABORT] Duplicate query IDs: $SEQ_COUNT sequences, $UNIQ_IDS unique."
        emit_summary "$SPECIES" "DUPLICATE_QUERY_IDS" "$SEQ_COUNT" 0 0 0 0 0 0
        return
    fi

    if [[ "$SEQ_COUNT" -eq 0 ]]; then
        echo "[SKIP] No sequences >= ${MIN_SEQ_LEN} bp to BLAST"
        emit_summary "$SPECIES" "NO_SEQS_PASSING_LENGTH" 0 0 0 0 0 0 0
        return
    fi

    echo "    queries: $SEQ_COUNT (from $n_extracted alignments)"

    # ------------------------------------------------------------------
    # [3] BLAST A: pooled plastid database.
    # ------------------------------------------------------------------
    echo "[3] BLAST A: pooled plastid database ..."
    blastn \
        -task "$BLAST_TASK" \
        -query "$OUT_DIR/dcnupts_all.fa" \
        -db "$PLASTID_DB" \
        -out "$OUT_DIR/blast_db_raw.tsv" \
        -outfmt "6 qseqid sseqid pident qcovhsp length mismatch gapopen evalue bitscore stitle" \
        -max_target_seqs "$MAX_TARGETS" \
        -max_hsps 1 \
        -evalue "$EVALUE" \
        -num_threads "$THREADS"

    # ------------------------------------------------------------------
    # [4] BLAST B: host's own plastid genome.
    #
    # The same BLAST task and E-value are used as in BLAST A so the host and
    # database comparisons are made under comparable search conditions.
    # ------------------------------------------------------------------
    echo "[4] BLAST B: host plastid genome ($HOST_ACC) ..."
    : > "$OUT_DIR/blast_host_raw.tsv"

    blastn \
        -task "$BLAST_TASK" \
        -query "$OUT_DIR/dcnupts_all.fa" \
        -subject "$HOST_FASTA" \
        -out "$OUT_DIR/blast_host_raw.tsv" \
        -outfmt "6 qseqid sseqid pident qcovhsp length mismatch gapopen evalue bitscore" \
        -max_hsps 1 \
        -evalue "$EVALUE" || true

    # Select the best qualifying host HSP per query by bitscore.
    awk -F'\t' -v cov="$MIN_QCOV" 'BEGIN{OFS="\t"}
        ($4+0) >= (cov+0) {
            q=$1
            if (!(q in bs) || ($9+0) > bs[q]) {
                bs[q]=$9+0
                pid[q]=$3+0
                qc[q]=$4+0
            }
        }
        END {
            for (q in bs) print q, pid[q], qc[q], bs[q]
        }
    ' "$OUT_DIR/blast_host_raw.tsv" | sort -k1,1 > "$OUT_DIR/host_identity.tsv"

    local N_WITH_HOST
    N_WITH_HOST=$(wc -l < "$OUT_DIR/host_identity.tsv")
    echo "    queries with a host alignment: $N_WITH_HOST / $SEQ_COUNT"

    # ------------------------------------------------------------------
    # [5] Tag pooled-database hits as host or non-host.
    #
    # Host membership is determined primarily from the plastid accession,
    # using a versionless accession match. The organism name is retained as
    # a secondary anchored check for database records whose accession naming
    # differs from the metadata.
    # ------------------------------------------------------------------
    echo "[5] Partitioning database hits ..."

    awk -F'\t' -v acc="$HOST_ACC_BASE" -v hostre="$HOST_RE" -v cov="$MIN_QCOV" '
        BEGIN { OFS="\t" }
        ($4+0) < (cov+0) { next }
        {
            is_host = 0
            if (index($2, acc) > 0) is_host = 1
            else if (hostre != "" && $10 ~ ("(^|[[:space:]|])" hostre "([[:space:]]|$|,|\\.)")) is_host = 1
            print is_host, $0
        }
    ' "$OUT_DIR/blast_db_raw.tsv" > "$OUT_DIR/blast_db_tagged.tsv"

    local DB_HOST_RECORDS
    DB_HOST_RECORDS=$(awk -F'\t' '$1==1' "$OUT_DIR/blast_db_tagged.tsv" | wc -l)

    echo "    database hits attributed to host: $DB_HOST_RECORDS"
    if [[ "$DB_HOST_RECORDS" -eq 0 ]]; then
        echo "    NOTE: host plastid not represented among database hits."
        echo "          Host identity comes from the direct host search."
    fi

    # Append all coverage-passing database hits to the combined project table.
    awk -F'\t' 'BEGIN{OFS="\t"} { $1=""; sub(/^\t/,""); print }' "$OUT_DIR/blast_db_tagged.tsv" \
        | awk -v sp="$SPECIES" 'BEGIN{OFS="\t"} { print sp, $0 }' >> "$COMBINED_OUT"

    # ------------------------------------------------------------------
    # [6] Classify each candidate.
    #
    # HOST:
    #   A qualifying host alignment is at least as similar as the best
    #   qualifying non-host hit.
    #
    # HIGH_CONF_NONHOST:
    #   The best non-host identity is at least HIGH_CONF_DIFF percentage
    #   points above the measured host identity.
    #
    # NONHOST_LOW_CONF:
    #   A qualifying non-host hit is present, but the identity margin is
    #   below HIGH_CONF_DIFF.
    #
    # NONHOST_HOST_ABSENT:
    #   A qualifying non-host hit is present, but no qualifying host
    #   alignment was detected in the direct host search.
    #
    # NO_HIT:
    #   Neither a qualifying host nor qualifying non-host hit is available.
    # ------------------------------------------------------------------
    echo "[6] Classifying ..."

    touch "$OUT_DIR/host_identity.tsv" "$OUT_DIR/blast_db_tagged.tsv"

    awk -F'\t' -v sp="$SPECIES" -v id="$IDENTITY_THRESHOLD" -v diff="$HIGH_CONF_DIFF" \
        -v hostfile="$OUT_DIR/host_identity.tsv" \
        -v idsfile="$OUT_DIR/dcnupts_all.fa" \
        -v dbfile="$OUT_DIR/blast_db_tagged.tsv" '
        BEGIN { OFS="\t" }

        FILENAME == hostfile {
            hpid[$1]=$2+0
            hqc[$1]=$3+0
            hseen[$1]=1
            next
        }

        FILENAME == idsfile {
            if ($0 ~ /^>/) allq[substr($0,2)]=1
            next
        }

        FILENAME == dbfile {
            is_host=$1
            q=$2
            sid=$3
            pid=$4+0
            qc=$5+0
            bs=$10+0
            sti=$11

            if (pid < (id+0)) next

            if (is_host) {
                if (!(q in hostbs) || bs > hostbs[q]) {
                    hostbs[q]=bs
                    hostpid_db[q]=pid
                }
            } else if (!(q in nhbs) || bs > nhbs[q]) {
                nhbs[q]=bs
                nhpid[q]=pid
                nhqc[q]=qc
                nhsid[q]=sid
                nhsti[q]=sti
            }
        }

        END {
            for (q in allq) {
                hp="NA"
                hq="NA"
                have_host=0

                if (q in hseen) {
                    hp=hpid[q]
                    hq=hqc[q]
                    have_host=1
                } else if (q in hostpid_db) {
                    hp=hostpid_db[q]
                    have_host=1
                }

                have_nh=(q in nhbs)
                host_ok=(have_host && hp != "NA" && (hp+0) >= (id+0))

                if (!have_nh && !host_ok) {
                    print sp,q,"NO_HIT","NA","NA","NA","NA","NA",hp,hq,"NA"
                    continue
                }

                if (host_ok && (!have_nh || (hp+0) >= nhpid[q])) {
                    printf "%s\t%s\tHOST\t%s\tNA\tNA\thost_plastome\thost_plastome\t%.3f\t%s\tNA\n", \
                           sp,q,hp,hp,hq
                    continue
                }

                if (!have_host) {
                    print sp,q,"NONHOST_HOST_ABSENT",nhpid[q],nhqc[q],nhbs[q],nhsid[q],nhsti[q],"NA","NA","NA"
                    continue
                }

                margin=nhpid[q]-hp
                cat=(margin >= (diff+0)) ? "HIGH_CONF_NONHOST" : "NONHOST_LOW_CONF"

                printf "%s\t%s\t%s\t%.3f\t%.1f\t%.1f\t%s\t%s\t%.3f\t%s\t%.3f\n", \
                       sp,q,cat,nhpid[q],nhqc[q],nhbs[q],nhsid[q],nhsti[q],hp,hq,margin
            }
        }
    ' "$OUT_DIR/host_identity.tsv" "$OUT_DIR/dcnupts_all.fa" "$OUT_DIR/blast_db_tagged.tsv" \
      > "$OUT_DIR/query_categories.tsv"

    cat "$OUT_DIR/query_categories.tsv" >> "$CATEGORIES_OUT"

    # ------------------------------------------------------------------
    # [7] Summarize classifications and check that every query was assigned.
    # ------------------------------------------------------------------
    local N_HOST N_LOW N_HC N_ABSENT N_NOHIT TOTAL
    N_HOST=$(awk -F'\t'   '$3=="HOST"'                "$OUT_DIR/query_categories.tsv" | wc -l)
    N_LOW=$(awk -F'\t'    '$3=="NONHOST_LOW_CONF"'    "$OUT_DIR/query_categories.tsv" | wc -l)
    N_HC=$(awk -F'\t'     '$3=="HIGH_CONF_NONHOST"'   "$OUT_DIR/query_categories.tsv" | wc -l)
    N_ABSENT=$(awk -F'\t' '$3=="NONHOST_HOST_ABSENT"' "$OUT_DIR/query_categories.tsv" | wc -l)
    N_NOHIT=$(awk -F'\t'  '$3=="NO_HIT"'              "$OUT_DIR/query_categories.tsv" | wc -l)

    TOTAL=$((N_HOST + N_LOW + N_HC + N_ABSENT + N_NOHIT))
    if [[ "$TOTAL" -ne "$SEQ_COUNT" ]]; then
        echo "[WARN] Category counts ($TOTAL) != query count ($SEQ_COUNT)." \
             | tee -a "$OUT_DIR/logs/consistency.log"
    fi

    emit_summary "$SPECIES" "OK" "$SEQ_COUNT" "$N_HOST" "$N_LOW" "$N_HC" "$N_ABSENT" "$N_NOHIT" "$DB_HOST_RECORDS"

    echo "  total candidates        : $SEQ_COUNT"
    echo "  host                    : $N_HOST"
    echo "  non-host, low conf      : $N_LOW"
    echo "  non-host, HIGH CONF     : $N_HC"
    echo "  non-host, host absent   : $N_ABSENT"
    echo "  no hit (<${IDENTITY_THRESHOLD}%) : $N_NOHIT"
    echo "Done: $SPECIES"
    echo
}

############################################
# MAIN
############################################

# Confirm that the configured pooled BLAST database exists before processing.
if [[ ! -f "${PLASTID_DB}.nsq" && ! -f "${PLASTID_DB}.00.nsq" && ! -f "${PLASTID_DB}.ndb" ]]; then
    echo "ERROR: plastid BLAST database not found at $PLASTID_DB" >&2
    echo "       Set PLASTID_DB in the CONFIG block." >&2
    exit 1
fi

init_outputs

echo "=============================================="
echo "Reading dcNUPTs from : */$RESULTS_SUBDIR/$(printf "$KAKS_BASENAME_FMT" '<SPECIES>')"
echo "Reading alignments   : */$ALIGN_SUBDIR/"
echo "Plastid database     : $PLASTID_DB"

_n_res=$(find "$SP_DIR" -maxdepth 3 -path "*/$RESULTS_SUBDIR/$(printf "$KAKS_BASENAME_FMT" '*')" 2>/dev/null | wc -l)
_n_dc=$(find "$SP_DIR" -maxdepth 3 -path "*/$RESULTS_SUBDIR/$(printf "$KAKS_BASENAME_FMT" '*')" 2>/dev/null \
        | xargs -r awk -F'\t' 'FNR>1 && $4!="NA" && $4!="" && ($4+0)>0.1' 2>/dev/null | wc -l)

echo "Species with a NUPT Ka/Ks table : $_n_res"
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
