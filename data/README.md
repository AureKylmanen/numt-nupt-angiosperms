# Merged NUPT and NUMT loci

Nuclear insertions of plastid (NUPT) and mitochondrial (NUMT) DNA detected in the
nuclear assemblies of 537 angiosperm species.

## Files

| File | Contents |
|---|---|
| `nupt_merged_beds.tar.gz` | 537 files, `<Species>_nupt_merged.bed` |
| `numt_merged_beds.tar.gz` | 537 files, `<Species>_numt_merged.bed` |


Species names in file names match the species table used in the manuscript.

## Format

Standard BED3, tab-separated, no header. One row per merged locus:

| Column | Description |
|---|---|
| 1 | Nuclear scaffold/chromosome accession of the species' assembly (e.g. `CM112924.1`) |
| 2 | Locus start (0-based) |
| 3 | Locus end (exclusive) |

Example:

```
CM112924.1	3246182	3246418
CM112924.1	9093750	9117052
```

## How loci were defined

1. Organelle genomes (plastome or mitogenome) were aligned to the nuclear assembly
   with BLASTN (dc-megablast) without a cap on target sequences.
2. Hits on scaffolds annotated as organellar in the assembly were removed.
3. The remaining hit fragments were merged into loci when separated by ≤10 kb.

Because of step 3, a merged locus can include nuclear sequence between the
organelle-derived fragments; locus length is therefore an upper bound on the
organelle-derived sequence it contains. These merged loci are the units used for
the NUPT/NUMT counts and genomic proportions in the manuscript.

Coordinates refer to the nuclear assembly version listed for each species in the
species table. Produced by `numt_nupt_detection.sh` in this repository.
