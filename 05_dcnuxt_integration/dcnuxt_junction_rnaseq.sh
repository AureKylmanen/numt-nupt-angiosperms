#!/usr/bin/env bash
# dcnuxt_junction_rnaseq.sh
#
# Test whether a dcNUxT locus is transcribed FROM ITS NUCLEAR POSITION.
#
# The discriminator is junction-spanning reads. A read lying wholly inside
# the insertion is ambiguous: the organellar transcript (atpA, psbA, atp1 ...)
# is abundant and matches the insertion sequence almost perfectly, which is
# why diagnostic-position BLAST keeps returning "plastid, not nuclear".
# A read that starts in unique nuclear flank and continues into the insertion
# can only come from the nuclear locus.
#
# Reports, per locus:
#   j5, j3        reads spanning the 5' and 3' insertion boundaries
#   flank5, flank3  mean depth over 300 bp of flank
#   inside        mean depth inside the insertion (context only - NOT evidence)
#
# A zero is only informative if flank depth is decent. Zero junction reads at
# 40x flank coverage means something; at 2x it means nothing. Both are printed
# so the distinction is never lost.
#
# Requires: sra-tools (fasterq-dump), hisat2, samtools
#
# Usage:
#   ./dcnuxt_junction_rnaseq.sh -s Phtheirospermum_japonicum -r SRR123,SRR124
#   ./dcnuxt_junction_rnaseq.sh -s Ziziphus_jujuba -r SRR999 -t 16 -k
#
#   -s species (must match combined_dataset.csv)
#   -r comma-separated SRA run accessions
#   -t threads (default 8)
#   -m margin in bp a read must overhang each side of a junction (default 10)
#   -q minimum MAPQ (default 20); with organellar decoys in the reference,
#      plastid transcript reads over a NUPT become multi-mappers and drop out
#   -k keep FASTQ and BAM (default: delete FASTQ, keep BAM)

set -euo pipefail

BASE=/vol/storage/5T_disk/organelles
LOCI="$BASE/dcnuxt_overlap_v2_exon_cluster.tsv"
THREADS=8
MARGIN=10
MINQ=20
KEEP=0
WORK="$BASE/rna_validation"

while getopts "s:r:t:m:q:b:l:w:k" opt; do
  case $opt in
    s) SPECIES="$OPTARG" ;;
    r) RUNS="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    m) MARGIN="$OPTARG" ;;
    q) MINQ="$OPTARG" ;;
    b) BASE="$OPTARG" ;;
    l) LOCI="$OPTARG" ;;
    w) WORK="$OPTARG" ;;
    k) KEEP=1 ;;
    *) echo "see header for usage" >&2; exit 1 ;;
  esac
done

: "${SPECIES:?-s species required}"
: "${RUNS:?-r run accessions required}"

for tool in fasterq-dump hisat2 hisat2-build samtools; do
  command -v "$tool" >/dev/null || { echo "ERROR: $tool not in PATH" >&2; exit 1; }
done

mkdir -p "$WORK/$SPECIES"
cd "$WORK/$SPECIES"

# ---------------------------------------------------------------- genome
ACC=$(awk -F, -v s="$SPECIES" '$1==s{print $2}' "$BASE/combined_dataset.csv")
[[ -n "$ACC" ]] || { echo "ERROR: $SPECIES not in combined_dataset.csv" >&2; exit 1; }
GEN=""
for cand in "$BASE"/downloads/genome_fasta/"${ACC}"_*_genomic.fna; do
  [[ -s "$cand" ]] && { GEN="$cand"; break; }
done
[[ -f "$GEN" ]] || { echo "ERROR: no unzipped genome for $ACC" >&2; exit 1; }

# Loci for this species, excluding intronic and organellar-contig calls.
awk -F'\t' -v s="$SPECIES" '
  NR>1 && $1==s && $24!="INTRONIC_ONLY" && $24!="ORGANELLAR_CONTIG" {
    print $4"\t"$5"\t"$6"\t"$3"\t"$12"\t"$24"\t"$2
  }' "$LOCI" | sort -u > loci.tsv

N=$(wc -l < loci.tsv)
(( N > 0 )) || { echo "no testable loci for $SPECIES" >&2; exit 1; }
echo "== $SPECIES: $N loci, genome $(basename "$GEN")"

# Map to the loci-bearing scaffolds PLUS the host plastome and mitogenome.
# The decoys are essential: without them a read from the abundant plastid
# atpA/psbA transcript has no true home in the index, maps to the NUPT, and
# is assigned MAPQ 60 because it is unique *within the restricted reference*.
# That produces sharp-edged coverage blocks over the insertion that look like
# transcription but are not. With the organelles present those reads map to
# their real source, or become multi-mappers with low MAPQ.
cut -f1 loci.tsv | sort -u > scaffolds.txt
[[ -f "${GEN}.fai" ]] || samtools faidx "$GEN"
samtools faidx "$GEN" $(tr '\n' ' ' < scaffolds.txt) > target.fa

PL=$(awk -F, -v s="$SPECIES" '$1==s{print $3}' "$BASE/combined_dataset.csv" || true)
MT=$(awk -F, -v s="$SPECIES" '$1==s{print $4}' "$BASE/combined_dataset.csv" || true)
DECOY=0
for acc in $(echo "${PL:-};${MT:-}" | tr ';' '\n' | grep -v '^$' || true); do
  # ls returns non-zero when ANY operand is missing, which under
  # `set -e` + `pipefail` kills the script even on a successful match.
  # Glob into an array instead.
  f=""
  for pat in "$BASE"/downloads/plastid_fasta/"${acc}"*.fasta \
             "$BASE"/downloads/mito_fasta/"${acc}"*.fasta \
             "$BASE"/multi_organelles/downloads/plastid_fasta/"${acc}"*.fasta \
             "$BASE"/multi_organelles/downloads/mito_fasta/"${acc}"*.fasta \
             "$BASE"/multi_organelles/downloads/plastid_fasta/"${acc%.*}"*.fasta \
             "$BASE"/multi_organelles/downloads/mito_fasta/"${acc%.*}"*.fasta; do
    [[ -s "$pat" ]] && { f="$pat"; break; }
  done
  if [[ -n "$f" ]]; then
    cat "$f" >> target.fa
    DECOY=$((DECOY+1))
    echo "   decoy: $(basename "$f")"
  else
    echo "   WARNING: no FASTA found for organelle accession $acc" >&2
  fi
done
echo "-- reference: $(wc -l < scaffolds.txt) scaffold(s) + $DECOY organellar decoy(s)"
(( DECOY > 0 )) || echo "   WARNING: no organellar decoys found - inside_depth will be inflated" >&2
samtools faidx target.fa

if [[ ! -f target_idx.1.ht2 ]]; then
  echo "-- building hisat2 index over $(wc -l < scaffolds.txt) scaffold(s)"
  hisat2-build -p "$THREADS" target.fa target_idx >/dev/null 2>&1
fi

# ---------------------------------------------------------------- reads
BAMS=()
IFS=',' read -ra RUNARR <<< "$RUNS"
for RUN in "${RUNARR[@]}"; do
  if [[ -f "${RUN}.bam" ]]; then
    echo "-- $RUN: BAM exists, reusing"
    BAMS+=("${RUN}.bam"); continue
  fi
  echo "-- $RUN: downloading"
  fasterq-dump --split-3 -e "$THREADS" -O . "$RUN" >/dev/null

  if [[ -f "${RUN}_1.fastq" && -f "${RUN}_2.fastq" ]]; then
    READS=(-1 "${RUN}_1.fastq" -2 "${RUN}_2.fastq")
  elif [[ -f "${RUN}.fastq" ]]; then
    READS=(-U "${RUN}.fastq")
  else
    echo "   WARNING: no FASTQ produced for $RUN, skipping" >&2; continue
  fi

  echo "-- $RUN: mapping"
  # --no-unal keeps the BAM small; the alignment rate line goes to the log
  hisat2 -p "$THREADS" -x target_idx "${READS[@]}" --no-unal 2> "${RUN}.hisat2.log" \
    | samtools sort -@ "$THREADS" -o "${RUN}.bam" -
  samtools index "${RUN}.bam"
  grep 'overall alignment rate' "${RUN}.hisat2.log" | sed "s/^/   $RUN /"
  (( KEEP )) || rm -f "${RUN}"*.fastq
  BAMS+=("${RUN}.bam")
done

(( ${#BAMS[@]} > 0 )) || { echo "ERROR: no BAMs produced" >&2; exit 1; }

if (( ${#BAMS[@]} > 1 )); then
  samtools merge -f -@ "$THREADS" merged.bam "${BAMS[@]}"
  samtools index merged.bam
  BAM=merged.bam
else
  BAM="${BAMS[0]}"
fi

# ---------------------------------------------------------------- count
OUT="${SPECIES}_junction_counts.tsv"
printf 'species\tkind\tdonor_gene\thost_gene\tverdict\tchrom\tstart\tend\tj5\tj3\tflank5_depth\tflank3_depth\tinside_depth\n' > "$OUT"

mean_depth() {  # chrom start end  (0-based half-open)
  local c=$1 s=$2 e=$3
  (( s < 0 )) && s=0
  samtools depth -a -Q "$MINQ" -r "${c}:$((s+1))-${e}" "$BAM" 2>/dev/null \
    | awk '{t+=$3; n++} END{printf "%.1f", (n? t/n : 0)}'
}

span_count() {  # chrom pos margin -- reads with pos inside an aligned block,
                # overhanging by >= margin on both sides
  local c=$1 p=$2 m=$3
  samtools view -q "$MINQ" "$BAM" "${c}:$((p-m))-$((p+m))" 2>/dev/null \
  | awk -v p="$p" -v m="$m" '
      {
        # Walk the CIGAR. The junction must fall within an aligned (M/=/X)
        # block: a read that splices OVER the insertion (…M<N>M…) covers the
        # coordinate but skips the sequence, and is evidence against, not for.
        pos=$4; cig=$6; s=pos; e=pos-1; ref=pos; ok=0
        while (match(cig, /^[0-9]+[MIDNSHP=X]/)) {
          len = substr(cig, 1, RLENGTH-1) + 0
          op  = substr(cig, RLENGTH, 1)
          if (op ~ /[M=X]/) {
            if (ref <= p-m && ref+len-1 >= p+m) ok=1   # aligned block covers it
            ref += len; e = ref-1
          } else if (op ~ /[DN]/) {
            ref += len; e = ref-1
          }
          cig = substr(cig, RLENGTH+1)
        }
        if (ok && s <= p-m) n++
      }
      END{print n+0}'
}

while IFS=$'\t' read -r CHR S E GENE HOST VERD KIND; do
  j5=$(span_count "$CHR" "$((S+1))" "$MARGIN")
  j3=$(span_count "$CHR" "$E" "$MARGIN")
  f5=$(mean_depth "$CHR" $((S-300)) "$S")
  f3=$(mean_depth "$CHR" "$E" $((E+300)))
  ind=$(mean_depth "$CHR" "$S" "$E")
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$SPECIES" "$KIND" "$GENE" "$HOST" "$VERD" "$CHR" "$S" "$E" \
    "$j5" "$j3" "$f5" "$f3" "$ind" >> "$OUT"
done < loci.tsv

echo
column -t "$OUT" 2>/dev/null || cat "$OUT"
echo
echo "Wrote $WORK/$SPECIES/$OUT"
echo "Reading it: j5/j3 are the evidence. Both non-zero = transcribed from the"
echo "nuclear locus. Zero with high flank depth = not transcribed. Zero with low"
echo "flank depth = untested. inside_depth alone proves nothing - the organellar"
echo "transcript maps there too."
