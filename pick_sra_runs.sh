#!/usr/bin/env bash
# pick_sra_runs.sh
#
# For each species with testable dcNUxT loci, pick up to N paired-end RNA-seq
# runs in a sensible depth window and write them to sra_runs.tsv.
#
# Why the window: runs below ~5M spots give flank depth too low for a zero to
# mean anything, and runs above ~50M spots are usually mislabelled genomic or
# Hi-C data (Phtheirospermum DRR523892-97 were 345-770M spots). Both waste time.
#
# Usage: ./pick_sra_runs.sh [species_list.txt] [runs_per_species]

set -uo pipefail
LIST=${1:-rna_species_testable.txt}
NRUNS=${2:-3}
MIN_SPOTS=5000000
MAX_SPOTS=50000000
OUT=sra_runs.tsv

: > "$OUT"
while read -r sp; do
  [[ -z "$sp" ]] && continue
  name="${sp//_/ }"
  info=""
  for try in 1 2 3; do
    info=$(esearch -db sra -query "\"$name\"[Organism] AND \"RNA-Seq\"[Strategy]" </dev/null 2>/dev/null \
           | efetch -format runinfo 2>/dev/null)
    [[ -n "$info" ]] && break
    sleep 5
  done

  runs=$(printf '%s\n' "$info" | awk -F, -v lo=$MIN_SPOTS -v hi=$MAX_SPOTS '
      NR==1 { for(i=1;i<=NF;i++){h[$i]=i}; next }
      $h["LibraryLayout"]=="PAIRED" && $h["spots"]+0>=lo && $h["spots"]+0<=hi {
        print $h["spots"]"\t"$h["Run"]
      }' | sort -k1,1 -rn | head -"$NRUNS" | cut -f2 | paste -sd,)

  if [[ -z "$runs" ]]; then
    # fall back to single-end if no paired run fits
    runs=$(printf '%s\n' "$info" | awk -F, -v lo=$MIN_SPOTS -v hi=$MAX_SPOTS '
        NR==1 { for(i=1;i<=NF;i++){h[$i]=i}; next }
        $h["spots"]+0>=lo && $h["spots"]+0<=hi { print $h["spots"]"\t"$h["Run"] }' \
      | sort -k1,1 -rn | head -"$NRUNS" | cut -f2 | paste -sd,)
    [[ -n "$runs" ]] && runs="${runs}"  # marked in output below
  fi

  printf "%s\t%s\n" "$sp" "${runs:-NONE}" | tee -a "$OUT"
  sleep 1
done < "$LIST"

echo "---"
awk -F'\t' '$2=="NONE"{n++} $2!="NONE"{y++} END{print y+0" species with runs, "n+0" without"}' "$OUT"
