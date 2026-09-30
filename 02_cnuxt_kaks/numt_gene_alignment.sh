#!/bin/bash

# ============================================================
# Script Name: numt_gene_alignment.sh
# Description: NUMT Ka/Ks pipeline for plants with BLAST prefilter
# Adapted from nupt_gene_alignment.sh
# Author: Aure Kylmänen
# Created: 08.06.2026
# ============================================================

set -euo pipefail
IFS=$'\n\t'

############################################################
# CONFIGURATION
############################################################

BASE_DIR="/vol/storage/5T_disk/organelles"
SPECIES_LIST="$1"

THREADS=12
MAFFT_THREADS=4
BLAST_THREADS=4
KAKS_BIN="/vol/storage/software/KaKs_Calculator2.0/bin/Linux/KaKs_Calculator"

GENETIC_CODE_MITO=1   # standard genetic code for plant mitochondria

############################################################
# UTILITIES
############################################################
orf_filter() {
    local seqfile=$1
    local min_orf=$2
    local tmppep=$(mktemp)

    transeq -sequence "$seqfile" -outseq "$tmppep" -frame 6 -auto >/dev/null 2>&1 || {
        rm -f "$tmppep"
        return 1
    }

    local max_aa_len=0
    while read -r line; do
        [[ "$line" =~ ^\> ]] && continue
        for chunk in $(echo "$line" | tr '*' '\n'); do
            len=${#chunk}
            (( len > max_aa_len )) && max_aa_len=$len
        done
    done < "$tmppep"

    rm -f "$tmppep"

    local max_nt=$(( max_aa_len * 3 ))

    if (( max_nt >= min_orf )); then
        return 0
    else
        return 1
    fi
}

############################################################
# ALIGNMENT (forward + reverse)
############################################################
align_pair() {
    local gene=$1
    local frag=$2
    local outdir=$3

    local gene_name frag_name tmp_pair tmp_pair_rev revfrag

    gene_name=$(basename "$gene" .fasta)
    frag_name=$(basename "$frag" .fa)

    mkdir -p "$outdir"

    if [[ ! -s "$gene" ]]; then
        echo "Skipping empty gene file: $gene" >&2
        return 0
    fi

    if [[ ! -s "$frag" ]]; then
        echo "Skipping empty fragment file: $frag" >&2
        return 0
    fi

    tmp_pair=$(mktemp)
    cat "$gene" "$frag" > "$tmp_pair"

    if [[ ! -s "$tmp_pair" ]] || [[ $(grep -c '^>' "$tmp_pair") -lt 2 ]]; then
        echo "Skipping invalid alignment input: $gene + $frag" >&2
        rm -f "$tmp_pair"
        return 0
    fi

    mafft --quiet --thread "$MAFFT_THREADS" "$tmp_pair" \
        > "$outdir/${frag_name}_${gene_name}.mafft" 2>/dev/null || {
            echo "MAFFT failed: $gene + $frag" >&2
            rm -f "$tmp_pair"
            return 0
    }

    revfrag=$(mktemp)
    revseq -sequence "$frag" -outseq "$revfrag" -snucleotide -auto >/dev/null 2>&1

    if [[ ! -s "$revfrag" ]]; then
        echo "Warning: revseq failed for $frag" >&2
        rm -f "$tmp_pair" "$revfrag"
        return 0
    fi

    tmp_pair_rev=$(mktemp)
    cat "$gene" "$revfrag" > "$tmp_pair_rev"

    mafft --quiet --thread "$MAFFT_THREADS" "$tmp_pair_rev" \
        > "$outdir/${frag_name}_${gene_name}_REV.mafft" 2>/dev/null

    rm -f "$tmp_pair" "$tmp_pair_rev" "$revfrag"
}

export -f align_pair
export -f orf_filter

############################################################
# PROCESS into AXT
############################################################
process_alignment() {
    local aln=$1
    local code=$2
    local axt_out=$3

    local tmp
    tmp=$(mktemp)
    seqkit seq -w 0 "$aln" > "$tmp"

    local seq1 seq2
    seq1=$(grep -v '^>' "$tmp" | sed -n '1p')
    seq2=$(grep -v '^>' "$tmp" | sed -n '2p')

    rm -f "$tmp"

    [[ -z "$seq1" || -z "$seq2" ]] && return 0

    local seq1_trim seq2_trim num num1
    seq1_trim=$(echo "$seq1" | sed 's/^-*//;s/-*$//')
    seq2_trim=$(echo "$seq2" | sed 's/^-*//;s/-*$//')

    num=$(echo "$seq1_trim" | tr -cd '-' | wc -c | tr -d ' ')
    num1=$(echo "$seq2_trim" | tr -cd '-' | wc -c | tr -d ' ')

    if [[ "$num" -ne 0 || "$num1" -ne 0 ]]; then
        return 0
    fi

    local l1 l2 len s1 s2
    l1=$(echo "$seq2" | sed 's/-*$//' | tr -cd '-' | wc -c | tr -d ' ')
    l2=$(echo "$seq2" | sed 's/^-*//' | tr -cd '-' | wc -c | tr -d ' ')
    len=${#seq1}

    # FIXED: guard against negative substring length. If the fragment
    # (seq2) is almost entirely gaps, l1+l2 can exceed len, which makes
    # ${seq1:l1:len-l1-l2} a negative-length substring expression -
    # a hard bash error ("substring expression < 0") that can abort
    # the rest of the calling loop under set -e, not just this entry.
    if (( l1 + l2 >= len )); then
        return 0
    fi

    if (( l2 != 0 && l1 != 0 )); then
        s1=${seq1:l1:len-l1-l2}
        s2=${seq2:l1:len-l1-l2}
    elif (( l2 == 0 )); then
        s1=${seq1:l1}
        s2=${seq2:l1}
    else
        s1=${seq1:0:len-l2}
        s2=${seq2:0:len-l2}
    fi

    [[ "$s1" == *"-"* ]] && return 0

    local len_s1 mod
    len_s1=${#s1}
    if (( len_s1 < 3 )); then
        return 0
    fi

    # ---- Codon-frame correction ----
    # FIXED: trim amounts were previously swapped (mod==1 trimmed 2,
    # mod==2 trimmed 1), which never actually produced a length
    # divisible by 3. Correct mapping: mod==1 -> trim 1; mod==2 -> trim 2.
    mod=$(( len_s1 % 3 ))
    if (( mod == 1 && len_s1 >= 1 )); then
        s1=${s1:0:len_s1-1}
        s2=${s2:0:len_s1-1}
    elif (( mod == 2 && len_s1 >= 2 )); then
        s1=${s1:0:len_s1-2}
        s2=${s2:0:len_s1-2}
    fi

    local tmpdna tmppep aa1 aa2
    tmpdna=$(mktemp)
    tmppep=$(mktemp)

    {
        echo ">gene"
        echo "$s1"
        echo ">frag"
        echo "$s2"
    } > "$tmpdna"

    transeq -sequence "$tmpdna" -outseq "$tmppep" \
        -table "$code" -auto >/dev/null 2>&1 || {
        rm -f "$tmpdna" "$tmppep"
        return 0
    }

    aa1=$(grep -A1 '^>gene' "$tmppep" | tail -n1 | tr -d '\n')
    aa2=$(grep -A1 '^>frag' "$tmppep" | tail -n1 | tr -d '\n')

    rm -f "$tmpdna" "$tmppep"

    [[ -z "$aa1" || -z "$aa2" ]] && return 0
    [[ ${#aa1} -ne ${#aa2} ]] && return 0

    # Check for internal stop codons (Chen et al. requirement)
    local gene_prot_trim="${aa1%\*}"
    local frag_prot_trim="${aa2%\*}"
    [[ "$gene_prot_trim" == *"*"* ]] && return 0
    [[ "$frag_prot_trim" == *"*"* ]] && return 0

    local id
    id=$(basename "$aln")
    {
        echo "$id"
        echo "${s1^^}"
        echo "${s2^^}"
        echo
    } >> "$axt_out"
}

export -f process_alignment

############################################################
# BLAST utilities
############################################################
make_mito_db() {
    local mito_dir="$1"
    local all_fasta="$mito_dir/all_mito.fasta"

    local files=()
    for f in "$mito_dir"/*.fasta; do
        [[ -f "$f" && -s "$f" ]] && files+=("$f")
    done

    if (( ${#files[@]} == 0 )); then
        echo "WARNING: no valid mito FASTA files in $mito_dir, skipping DB creation" >&2
        return 1
    fi

    cat "${files[@]}" > "$all_fasta"

    if [[ ! -s "$all_fasta" ]]; then
        echo "WARNING: concatenation produced empty FASTA in $mito_dir" >&2
        return 1
    fi

    if [[ ! -f "$mito_dir/mito_db.nsq" ]]; then
        makeblastdb -in "$all_fasta" -dbtype nucl \
            -out "$mito_dir/mito_db" >/dev/null 2>&1 || {
            echo "ERROR: makeblastdb failed for $mito_dir" >&2
            return 1
        }
    fi
}

blast_filter_frag() {
    local frag_file="$1"
    local mito_dir="$2"
    local out_file="$3"

    blastn -query "$frag_file" \
           -db "$mito_dir/mito_db" \
           -outfmt "6 qseqid sseqid pident qcovs length evalue bitscore" \
           -perc_identity 70 \
           -qcov_hsp_perc 50 \
           -evalue 1e-5 \
           -max_target_seqs 50 \
           -num_threads "$BLAST_THREADS" \
           > "$out_file"
}

############################################################
# MAIN
############################################################
export -f orf_filter align_pair process_alignment
export -f blast_filter_frag

while IFS=',' read -r SPECIES ASSEMBLY PLASTID MITO; do
    [[ "$SPECIES" == "SPECIES" || -z "$SPECIES" ]] && continue
    echo "=============================="
    echo "Processing $SPECIES"

    SPECIES_DIR="$BASE_DIR/species/$SPECIES"
    NUMT_FILE="$SPECIES_DIR/numts/numts.seq.fa"
    MITO_GENES=$(ls -d "$BASE_DIR/downloads/mito_genes/${MITO}"* \
        2>/dev/null | head -n1 || echo "")
    ALIGN_MITO="$SPECIES_DIR/alignments/mito"
    FRAGMENTS_NUMT="$SPECIES_DIR/fragments/numt"
    RESULTS="$SPECIES_DIR/results"

    mkdir -p "$ALIGN_MITO" "$FRAGMENTS_NUMT" "$RESULTS"

    if [[ ! -f "$NUMT_FILE" ]]; then
        echo "WARNING: NUMT file not found for $SPECIES, skipping..."
        continue
    fi

    if [[ ! -d "$MITO_GENES" ]]; then
        echo "ERROR: mito gene directory missing for $SPECIES: $MITO" >&2
        continue
    fi

    if ! make_mito_db "$MITO_GENES"; then
        echo "Skipping $SPECIES due to invalid mito DB" >&2
        continue
    fi

    MIN_ORF=90
    echo "Using ORF cutoff: $MIN_ORF bp for $SPECIES"

    # Split NUMTs to fragments and prefilter by ORF
    awk '/^>/{if(seq){print seq >> out; close(out)}; name=$0; gsub(/^>/,"",name); safe=name; gsub(/[^A-Za-z0-9_.-]/,"_",safe); out="'"$FRAGMENTS_NUMT"'/"safe".fa"; print ">"safe > out; seq=""; next}
         {print >> out}
         END{if(seq){close(out)}}' "$NUMT_FILE"

    for frag in "$FRAGMENTS_NUMT"/*.fa; do
        [[ -f "$frag" ]] || continue
        if ! orf_filter "$frag" "$MIN_ORF"; then
            rm -f "$frag"
        fi
    done

    # Align only candidate genes from BLAST
    for frag in "$FRAGMENTS_NUMT"/*.fa; do
        [[ -f "$frag" ]] || continue
        blast_out="${frag}.blast.out"
        blast_filter_frag "$frag" "$MITO_GENES" "$blast_out"

        # Strip accession prefix from gene IDs (e.g. PX251578.1_cox3 -> cox3)
        mapfile -t genes < <(cut -f2 "$blast_out" | sort -u | \
            sed 's/^[^_]*\.[0-9]*_//')

        for gene_id in "${genes[@]}"; do
            gene_file="$MITO_GENES/${gene_id}.fasta"
            [[ -f "$gene_file" ]] || \
                gene_file=$(find "$MITO_GENES" -maxdepth 1 \
                    -name "${gene_id}*.fasta" | head -n1)
            [[ -f "$gene_file" ]] || continue

            align_pair "$gene_file" "$frag" "$ALIGN_MITO" &
            while (( $(jobs -rp | wc -l) >= THREADS )); do
                sleep 0.1
            done
        done
    done
    wait

    AXT="$RESULTS/numt_input.axt"
    : > "$AXT"

    for aln in "$ALIGN_MITO"/*.mafft; do
        [[ -f "$aln" ]] || continue
        process_alignment "$aln" "$GENETIC_CODE_MITO" "$AXT"
    done

    if [[ -s "$AXT" ]]; then
        echo "Running KaKs_Calculator for $SPECIES"
        "$KAKS_BIN" -i "$AXT" -c 1 \
            -o "$RESULTS/${SPECIES}_numt_kaks.csv" -m LWL
    else
        echo "No valid NUMT alignments for $SPECIES (AXT empty), skipping KaKs."
    fi

    echo "Finished $SPECIES"

done < "$SPECIES_LIST"

echo "Pipeline complete."
