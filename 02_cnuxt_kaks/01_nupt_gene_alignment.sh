#!/bin/bash

# ============================================================
# Script Name: 01_nupt_gene_alignment.sh
# Description: NUPT Ka/Ks pipeline for plants with BLAST prefilter
# Created: 27.3.2026
# Author: Aure Kylmänen
# ============================================================

set -euo pipefail
IFS=$'\n\t'

############################################################
# CONFIGURATION
############################################################

BASE_DIR="${BASE_DIR:-$PWD}"   # working directory; override with BASE_DIR=...
export BASE_DIR
SPECIES_LIST="$1"         # CSV: SPECIES,ASSEMBLY,PLASTID,MITO

THREADS=12                # max concurrent fragment/gene align jobs
MAFFT_THREADS=4            # threads per MAFFT process (tune to machine)
BLAST_THREADS=4
KAKS_BIN="/vol/storage/software/KaKs_Calculator2.0/bin/Linux/KaKs_Calculator"

GENETIC_CODE_PLASTID=11   # bacterial/plant plastid


############################################################
# UTILITIES
############################################################
orf_filter() {
    local seqfile=$1
    local min_orf=$2
    local tmppep=$(mktemp)

    # Translate all 6 frames
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

    # Convert to nucleotide length
    max_nt=$(( max_aa_len * 3 ))

    if (( max_nt >= min_orf )); then
        return 0  # keep fragment
    else
        return 1  # discard
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

    # Skip empty or missing inputs
    if [[ ! -s "$gene" ]]; then
        echo "Skipping empty gene file: $gene" >&2
        return 0
    fi

    if [[ ! -s "$frag" ]]; then
        echo "Skipping empty fragment file: $frag" >&2
        return 0
    fi

    # forward alignment: gene + fragment
    tmp_pair=$(mktemp)
    cat "$gene" "$frag" > "$tmp_pair"

    # Skip malformed combined FASTA
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

    # reverse-complement fragment to a unique temp file
    revfrag=$(mktemp)
    revseq -sequence "$frag" -outseq "$revfrag" -snucleotide -auto >/dev/null 2>&1

    if [[ ! -s "$revfrag" ]]; then
        echo "Warning: revseq failed for $frag" >&2
        rm -f "$tmp_pair" "$revfrag"
        return 0
    fi

    # reverse alignment: gene + reverse fragment
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

    # Re-wrap to 1-line sequences
    local tmp
    tmp=$(mktemp)
    seqkit seq -w 0 "$aln" > "$tmp"

    # Expect exactly two sequences: gene first, fragment second
    local seq1 seq2
    seq1=$(grep -v '^>' "$tmp" | sed -n '1p')
    seq2=$(grep -v '^>' "$tmp" | sed -n '2p')

    rm -f "$tmp"

    # Skip empty / malformed alignments
    [[ -z "$seq1" || -z "$seq2" ]] && return 0

    # Remove leading/trailing gaps to look for internal gaps
    local seq1_trim seq2_trim num num1
    seq1_trim=$(echo "$seq1" | sed 's/^-*//;s/-*$//')
    seq2_trim=$(echo "$seq2" | sed 's/^-*//;s/-*$//')

    num=$(echo "$seq1_trim" | tr -cd '-' | wc -c | tr -d ' ')
    num1=$(echo "$seq2_trim" | tr -cd '-' | wc -c | tr -d ' ')

    # Require no internal gaps in either (vertebrate script required this)
    if [[ "$num" -ne 0 || "$num1" -ne 0 ]]; then
        return 0
    fi

    # Now trim based on fragment coverage, following Python logic:
    # l1: internal gaps before trailing block, l2: internal gaps after leading block
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

    # Reject if gene still has internal gaps after trimming
    [[ "$s1" == *"-"* ]] && return 0

    # Codon-correct to length multiple of 3 by trimming from end of both sequences
    local len_s1 mod
    len_s1=${#s1}
    if (( len_s1 < 3 )); then
        return 0
    fi

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

    # Translation validation: gene & fragment AA lengths must match
    local tmpdna tmppep aa1 aa2
    tmpdna=$(mktemp)
    tmppep=$(mktemp)

    {
        echo ">gene"
        echo "$s1"
        echo ">frag"
        echo "$s2"
    } > "$tmpdna"

    transeq -sequence "$tmpdna" -outseq "$tmppep" -table "$code" -auto >/dev/null 2>&1 || {
        rm -f "$tmpdna" "$tmppep"
        return 0
    }

    aa1=$(grep -A1 '^>gene' "$tmppep" | tail -n1 | tr -d '\n')
    aa2=$(grep -A1 '^>frag' "$tmppep" | tail -n1 | tr -d '\n')

    rm -f "$tmpdna" "$tmppep"

    [[ -z "$aa1" || -z "$aa2" ]] && return 0
    [[ ${#aa1} -ne ${#aa2} ]] && return 0

    # Check for internal stop codons (Chen et al. requirement)
    # Remove legitimate trailing stop codon then check for internal ones
    local gene_prot_trim="${aa1%\*}"
    local frag_prot_trim="${aa2%\*}"
    [[ "$gene_prot_trim" == *"*"* ]] && return 0
    [[ "$frag_prot_trim" == *"*"* ]] && return 0

    # Write to AXT: 3 lines + blank line
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
make_plastid_db() {
    local plastid_dir="$1"
    local all_fasta="$plastid_dir/all_plastid.fasta"

    # collect only NON-empty fasta files
    local files=()
    for f in "$plastid_dir"/*.fasta; do
        [[ "$(basename "$f")" == all_* ]] && continue      # our own concatenation
        [[ -f "$f" && -s "$f" ]] && files+=("$f")
    done

    # if nothing valid → skip completely
    if (( ${#files[@]} == 0 )); then
        echo "WARNING: no valid plastid FASTA files in $plastid_dir, skipping DB creation" >&2
        return 1
    fi

    # rebuild only if needed
    cat "${files[@]}" > "$all_fasta"

    # if still empty → abort
    if [[ ! -s "$all_fasta" ]]; then
        echo "WARNING: concatenation produced empty FASTA in $plastid_dir" >&2
        return 1
    fi

    # rebuilt every run so a changed gene set is never searched with a stale DB
    makeblastdb -in "$all_fasta" -dbtype nucl -out "$plastid_dir/plastid_db" >/dev/null 2>&1 || {
        echo "ERROR: makeblastdb failed for $plastid_dir" >&2
        return 1
    }
}

blast_filter_frag() {
    local frag_file="$1"
    local plastid_dir="$2"
    local out_file="$3"

    blastn -query "$frag_file" \
           -db "$plastid_dir/plastid_db" \
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
    NUPT_FILE="$SPECIES_DIR/nupts/nupts.seq.fa"
    PLASTID_GENES=$(ls -d "$BASE_DIR/downloads/plastid_genes/${PLASTID}"* 2>/dev/null | head -n1 || echo "")
    ALIGN_PLASTID="$SPECIES_DIR/alignments/plastid"
    FRAGMENTS_NUPT="$SPECIES_DIR/fragments/nupt"
    RESULTS="$SPECIES_DIR/results"

    # start clean: fragments, alignments and Ka/Ks output of an earlier run would
    # otherwise be reused (stale fragments from an older detection run, e.g. on
    # organelle-annotated sequences, would end up in the AXT)
    rm -rf "$ALIGN_PLASTID" "$FRAGMENTS_NUPT"
    rm -f "$RESULTS/${SPECIES}_nupt_kaks.csv" "$RESULTS/${SPECIES}_nupt_kaks.csv".pre_* "$RESULTS/nupt_input.axt"
    mkdir -p "$ALIGN_PLASTID" "$FRAGMENTS_NUPT" "$RESULTS"

    if [[ ! -d "$PLASTID_GENES" ]]; then
        echo "ERROR: plastid gene directory missing: $PLASTID_GENES" >&2
        continue
    fi

    if ! make_plastid_db "$PLASTID_GENES"; then
        echo "Skipping $SPECIES due to invalid plastid DB" >&2
        continue
    fi

    MIN_ORF=90
    echo "Using ORF cutoff: $MIN_ORF bp for $SPECIES"

    # split NUPTs to fragments and prefilter by ORF
    awk '/^>/{if(seq){print seq >> out; close(out)}; name=$0; gsub(/^>/,"",name); safe=name; gsub(/[^A-Za-z0-9_.-]/,"_",safe); out="'$FRAGMENTS_NUPT'/"safe".fa"; print ">"safe > out; seq=""; next}
         {print >> out}
         END{if(seq){close(out)}}' "$NUPT_FILE"

    for frag in "$FRAGMENTS_NUPT"/*.fa; do
        [[ -f "$frag" ]] || continue
        if ! orf_filter "$frag" "$MIN_ORF"; then
            rm -f "$frag"
        fi
    done

    # Align only candidate genes from BLAST
    for frag in "$FRAGMENTS_NUPT"/*.fa; do
        [[ -f "$frag" ]] || continue
        blast_out="${frag}.blast.out"
        blast_filter_frag "$frag" "$PLASTID_GENES" "$blast_out"

        mapfile -t genes < <(cut -f2 "$blast_out" | sort -u)

        for gene_id in "${genes[@]}"; do
            # gene_id may include pathless ID - resolve .fasta
            gene_file="$PLASTID_GENES/${gene_id}.fasta"
            [[ -f "$gene_file" ]] || gene_file=$(find "$PLASTID_GENES" -maxdepth 1 -name "${gene_id}*.fasta" | head -n1)
            [[ -f "$gene_file" ]] || continue

            align_pair "$gene_file" "$frag" "$ALIGN_PLASTID" &
            # job control
            while (( $(jobs -rp | wc -l) >= THREADS )); do
                sleep 0.1
            done
        done
    done
    wait

    AXT="$RESULTS/nupt_input.axt"
    : > "$AXT"

    for aln in "$ALIGN_PLASTID"/*.mafft; do
        [[ -f "$aln" ]] || continue
        process_alignment "$aln" "$GENETIC_CODE_PLASTID" "$AXT"
    done

    if [[ -s "$AXT" ]]; then
        echo "Running KaKs_Calculator for $SPECIES"
        "$KAKS_BIN" -i "$AXT" -c 11 -o "$RESULTS/${SPECIES}_nupt_kaks.csv" -m LWL
    else
        echo "No valid alignments for $SPECIES (AXT empty), skipping KaKs."
    fi

    echo "Finished $SPECIES"

done < "$SPECIES_LIST"

echo "Pipeline complete."
