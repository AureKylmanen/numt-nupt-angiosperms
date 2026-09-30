# NUMTs and NUPTs in angiosperm genomes

Code for detecting nuclear insertions of mitochondrial (NUMT) and plastid (NUPT) DNA across more than 500 angiosperm reference genomes, and for testing how these insertions evolve after transfer: how much of each genome they occupy, whether that content is phylogenetically structured, whether coding insertions are under selection, whether old insertions have been integrated into nuclear genes, and whether any insertions derive from a non-host organelle.

---

## Terminology

| Term | Meaning |
|---|---|
| NUMT / NUPT | Nuclear sequence of mitochondrial / plastid origin (BLASTN hit ≥ 200 bp) |
| NUxT | Either class |
| cNUxT | *Coding* NUxT: contains an ORF ≥ 90 bp and yields a gap-free codon alignment to an organellar protein-coding gene with a calculable Ka/Ks |
| dcNUxT | *Diverged* coding NUxT: cNUxT with Ks > 0.1 against the host organelle gene, i.e. not a recent transfer |
| Non-host ("alien") dcNUxT | dcNUxT whose best organellar match is a species other than the host |

---

## Repository layout

```
01_detection/              NUMT/NUPT detection, counts and genome proportions
02_cnuxt_kaks/             Coding insertions, Ka/Ks, dcNUxT definition, gene-level summaries
03_comparative/            Phylogenetic signal and PGLS on insertion content
04_clustering_selection/   Shared-origin clustering of dcNUxTs and codeml branch tests
05_dcnuxt_integration/     Genomic context of dcNUxTs, host-gene structure, RNA-seq junctions
06_nonhost/                Non-host (alien) dcNUxT detection
figures/                   Plotting scripts not tied to a single stage
```

---

## Requirements

Versions are those stated in the manuscript Methods where given; pin the rest from the
analysis environment (`conda env export --no-builds > environment.yml`).

| Software | Version | Used in |
|---|---|---|
| BLAST+ (`blastn`, `makeblastdb`) | | 01, 02, 04, 06 |
| bedtools | | 01 |
| EMBOSS (`transeq`, `revseq`) | 6.5.7 | 02 |
| seqkit | | 02, 04 |
| MAFFT | 7.490 | 02, 04 |
| KaKs_Calculator | 2.0 | 02, 04 |
| GNU parallel | | 04 |
| IQ-TREE 2 | | 04 |
| pal2nal | v14 | 04 |
| PAML (`codeml`) | 4.10 | 04 |
| HISAT2, samtools, sra-tools | | 05 |
| Python 3 | | Biopython, pandas, numpy, scipy, dendropy, matplotlib |
| R | | ape, nlme, dplyr, phytools, geiger, ggplot2 |

---
## Data layout and configuration

All scripts assume one working directory (default `/vol/storage/5T_disk/organelles`,
set as `BASE_DIR` / `WORKING_DIR` / `BASE` at the top of each script, or via
`NUXT_BASE` / `NUXT_OUT` for the R scripts in `03_comparative`).

The species manifest is a CSV with a header row:

```
SPECIES,ASSEMBLY,PLASTID,MITO
Arabidopsis_thaliana,GCA_000001735.2,NC_000932.1,NC_037304.1
```

`SPECIES` is the underscore-separated binomial used as directory name and tree tip
label; `ASSEMBLY` is the GenBank/RefSeq assembly accession; `PLASTID` and `MITO` are
organelle genome accessions.

```
<BASE>/
├── combined_dataset.csv            species manifest
├── downloads/
│   ├── genome_fasta/               <ASSEMBLY>_*.fna.gz (+ GFF, assembly reports)
│   ├── plastid_fasta/  mito_fasta/ <ACCESSION>.fasta
│   ├── plastid_gb/                 GenBank records for plastid gene extraction
│   └── plastid_genes/  mito_genes/ <ACCESSION>/<gene>.fasta (one CDS per file)
└── species/<SPECIES>/
    ├── numts/  nupts/              BLAST output, fragment BEDs, insertion FASTA
    ├── fragments/  alignments/     per-fragment FASTA and MAFFT alignments (stage 02)
    └── results/                    per-species counts, proportions, Ka/Ks tables
```
---
 
## 01_detection — NUMT/NUPT detection
 
| Script | Purpose |
|---|---|
| `numt_nupt_detection.sh` | BLASTN of the mitogenome and plastome against the nuclear assembly, fragment extraction, counts, genome proportions, NUMT–NUPT overlap |
| `get_bed2.R` | Called by the above: filters hits, orients coordinates, writes fragment and merged BEDs |
 
```bash
THREADS=24 WORKING_DIR=/vol/storage/5T_disk/organelles \
  ./01_detection/numt_nupt_detection.sh combined_dataset.csv
```
 
Parameters: `blastn -word_size 20`; hits shorter than 200 bp are discarded
(`MIN_LENGTH`); hits on the same scaffold less than 10 kb apart are merged into one
insertion (`MERGE_DISTANCE`). Proportions are insertion bp / assembly bp, both for
individual fragments ("unmerged") and merged insertions.
 
Outputs in `species/<SPECIES>/results/`: `*.numt_count.txt`, `*.numt_merged_count.txt`,
`*.numt_unmerged_prop.txt`, `*.numt_merged_prop.txt` (and NUPT equivalents),
`*.genome_size.txt`, `*.numt_nupt_overlap_{count,bp}.txt`. Insertion sequences go to
`numts/numts.seq.fa` and `nupts/nupts.seq.fa`.
 
---
 
## 02_cnuxt_kaks — coding insertions and Ka/Ks
 
Run in this order:
 
| # | Script | Purpose |
|---|---|---|
| 1 | `gb_to_fasta.py` | Extract plastid CDS from GenBank records into `downloads/plastid_genes/<ACC>/<gene>.fasta` |
| 2 | `nupt_gene_alignment.sh` / `numt_gene_alignment.sh` | ORF filter (≥ 90 bp in any of six frames) → BLASTN pre-filter against the host's organelle genes (≥ 70 % identity, ≥ 50 % query cover, e ≤ 1e-5) → MAFFT alignment of each fragment and its reverse complement to each candidate gene → gap-free, in-frame, stop-free codon alignments → KaKs_Calculator (LWL) |
| 3 | `merge_kaks.py` | Combine per-species tables into `all_species_{nupt,numt}_kaks.tsv`: adds `Species`, makes column names identifier-safe (`Ka/Ks` → `Ka_Ks`, `P-Value(Fisher)` → `P_Value_Fisher`), keeps rows with a numeric Ka/Ks |
| 4 | `bh_correction_dcnuxt.py` | Keep Ks > 0.1 (dcNUxTs); Benjamini–Hochberg correction of the Fisher exact test p-values, separately for NUPTs and NUMTs |
| 5 | `mannwhitney_kaks_comparison.py` | Compare dcNUPT and dcNUMT Ka/Ks distributions (Mann–Whitney U, rank-biserial r) |
| 6 | `build_gene_lengths.py` | Reference length of each organellar gene, averaged across accessions |
| 7 | `gene_level_summary.py` | Per gene: cNUxT and dcNUxT counts, species breadth, fraction old (ret%), fraction under purifying selection (sel%) |
| 8 | `gene_predictors.py` | Gene length vs transfer and retention; length-normalised transfer vs retention; NUPT vs NUMT comparison (Fisher r-to-z) |
 
```bash
./02_cnuxt_kaks/nupt_gene_alignment.sh combined_dataset.csv
./02_cnuxt_kaks/numt_gene_alignment.sh combined_dataset.csv
python3 02_cnuxt_kaks/merge_kaks.py nupt      # species/ and multi_organelles/species/
python3 02_cnuxt_kaks/merge_kaks.py numt
python3 02_cnuxt_kaks/bh_correction_dcnuxt.py all_species_nupt_kaks.tsv all_species_numt_kaks.tsv
python3 02_cnuxt_kaks/mannwhitney_kaks_comparison.py dcnupt_final_bh.tsv dcnumt_final_bh.tsv
python3 02_cnuxt_kaks/build_gene_lengths.py            # run from <BASE>
python3 02_cnuxt_kaks/gene_level_summary.py --latex
python3 02_cnuxt_kaks/gene_predictors.py --min-cnuxt 5 --min-dcnuxt 3
```
 
**Genetic codes.** Plant plastids use NCBI translation table 11 and plant mitochondria
the standard code, table 1 (C-to-U RNA editing rather than codon reassignment).
KaKs_Calculator and pal2nal are run with `-c 11` / `-codontable 11` for NUPTs and `1`
for NUMTs. codeml is run with `icode = 0` for both, because table 11 assigns the same
amino acid to every codon as table 1 and differs only in permitted start codons, so
synonymous/non-synonymous site classification is unaffected.
 
---
 
## 03_comparative — phylogenetic signal and PGLS
 
Insertion proportions are strongly right-skewed and are log10-transformed before any
Brownian-motion-based analysis. The tree is a TimeTree backbone with missing species
grafted to their closest relative (`tree/species_names_grafted.nwk`), made ultrametric
by extending terminal branches.
 
| Script | Purpose | Main output (`$NUXT_OUT`) |
|---|---|---|
| `common.R` | Shared data loading, organelle-flag exclusions, λ-profiled PGLS | — |
| `01_phylo_signal.R` | Blomberg's K (tip randomisation) and Pagel's λ (LRT vs 0 and vs 1) for NUMT, NUPT, NUPT:NUMT ratio and scaffold N50 | `phylo_signal_table.tsv` |
| `02_pgls_models.R` | PGLS (`nlme::gls` + `corPagel`, λ profiled on a 0.02 grid because free estimation fails to converge for several models) with 95 % profile CI and OLS alongside | `pgls_table.tsv`, `lambda_profiles.tsv`, `pgls_fits.rds`, `excluded_species.tsv` |
| `03_plot_pgls_panels.R` | Scatter panels with PGLS fit and 95 % band | `Figures/*.pdf` |
| `tech_effect_test.R` | Effect of sequencing technology (long-read / hybrid / short-read) on insertion content: rank tests, OLS, PGLS, organelle-scaffold sensitivity, within-genus comparison | `tech_effect/` |
| `assembly_quality.R` | Whether scaffold N50 and technology explain insertion content and the NUMT–NUPT difference in λ | stdout |
 
```bash
Rscript 03_comparative/01_phylo_signal.R      # PHYLOSIG_NSIM=10000 for finer K P-values
Rscript 03_comparative/02_pgls_models.R
Rscript 03_comparative/03_plot_pgls_panels.R
Rscript 03_comparative/tech_effect_test.R final_species_summary_v2.tsv tree/assembly_stats.tsv
```
 
Species whose organelle genome was flagged are excluded only from models of the
affected class (`exclude_numt` / `exclude_nupt` in `organelle_lengths_final.tsv`).
 
---
 
## 04_clustering_selection — shared-origin clusters and codeml
 
Separate NUPT and NUMT runs, writing to `clustering_{nupt,numt}/` and
`selection_{nupt,numt}/`. The Python steps take `--organelle nupt|numt` (default
`nupt`).
 
| # | Script | Purpose |
|---|---|---|
| 1 | `01_cluster_dcnupts.sh` / `01_cluster_dcnumts.sh` | Collect dcNUxTs per gene (identical copies within a species collapsed), pairwise Ks between all dcNUxTs of a gene, link pairs with Ks < 0.1 where neither copy is ≥ 98 % identical to any organelle genome in the pooled database, connected components → clusters and singletons |
| 2 | `02_prepare_background_alignments_{nupt,numt}.sh` | Background protein alignment per gene across all accessions (sequences < 80 % of median length removed); add each cluster/singleton with `mafft --add` |
| 3 | `03_build_trees_iqtree.py` | IQ-TREE 2 (LG, `-fast`) tree per combined alignment |
| 4 | `04_rebuild_subtrees_and_alignments.py` | Resolve polytomies, climb up to 3 nodes above the insertion MRCA (≤ 50 taxa), prune, codon alignment with pal2nal |
| 5 | `05_run_codeml_stage.py` | Label foreground branches, run codeml two-ratio branch model vs ω_fg = 1, LRT (χ², df = 1), BH correction, write `results/selection_results.tsv` |
 
```bash
bash   04_clustering_selection/01_cluster_dcnupts.sh
bash   04_clustering_selection/02_prepare_background_alignments_nupt.sh
python3 04_clustering_selection/03_build_trees_iqtree.py --organelle nupt --jobs 7 --threads 4
python3 04_clustering_selection/04_rebuild_subtrees_and_alignments.py --organelle nupt --levels 3 --max-taxa 50
python3 04_clustering_selection/05_run_codeml_stage.py --organelle nupt --threads 5
 
# NUMTs
bash   04_clustering_selection/01_cluster_dcnumts.sh
bash   04_clustering_selection/02_prepare_background_alignments_numt.sh
python3 04_clustering_selection/03_build_trees_iqtree.py --organelle numt
python3 04_clustering_selection/04_rebuild_subtrees_and_alignments.py --organelle numt
python3 04_clustering_selection/05_run_codeml_stage.py --organelle numt
```
 
Branch labelling follows Chen et al.'s `numt_dnds.7.py`: for clusters, all branches
below the MRCA of the insertions are foreground and the stem (integration) branch is
background; for singletons the terminal branch is foreground. Foreground ω ≥ 10
indicates dS ≈ 0 and is reported as uninformative. `05_run_codeml_stage.py` skips
finished runs unless `--force` is given, so re-running it only re-parses results.
 
---
 
## 05_dcnuxt_integration — genomic context and transcription
 
| # | Script | Purpose |
|---|---|---|
| 1 | `dcnuxt_annotation_overlap.py` | Map dcNUxTs back to genome coordinates, classify against the GFF (CDS > UTR/exon > intron > pseudogene > intergenic), compare with a within-scaffold permutation null. `--audit-only` checks GFF availability and seqid concordance first |
| 2 | `dcnuxt_exon_cluster.py` | For loci inside annotated genes: host-gene exon structure, fraction of the gene formed by the insertion, whether the annotation merely restates the donor gene; verdicts CANDIDATE / CIRCULAR / SINGLE_EXON / INTRONIC_ONLY / ORGANELLAR_CONTIG; writes RNA-seq target BED/FASTA |
| 3 | `dcnuxt_junction_rnaseq.sh` | HISAT2 mapping of public RNA-seq (with organellar decoys) and counts of reads spanning the nuclear-flank/insertion junctions — the only reads that can come from the nuclear copy |
 
```bash
python3 05_dcnuxt_integration/dcnuxt_annotation_overlap.py \
    --base /vol/storage/5T_disk/organelles \
    --gff-dir /vol/storage/5T_disk/organelles/downloads/genome_fasta \
    --perms 1000 --out-prefix dcnuxt_overlap
python3 05_dcnuxt_integration/dcnuxt_exon_cluster.py \
    --loci dcnuxt_overlap_loci.tsv --csv combined_dataset.csv \
    --gff-dir downloads/genome_fasta --genome-dir downloads/genome_fasta \
    --flank 300 --out-prefix dcnuxt_overlap
./05_dcnuxt_integration/dcnuxt_junction_rnaseq.sh -s Juglans_regia -r SRRxxxxxxx -t 16
```
 
---
 
## 06_nonhost — non-host (alien) dcNUxTs
 
`alien_nupt_detection.sh` and `alien_numt_detection.sh` BLAST each dcNUxT against
(a) the pooled organelle genome database and (b) every available organelle genome of
the host species, so host identity is measured directly rather than inferred from its
absence among top database hits. Each query is classified by whether its best match is
host or non-host and by the identity margin between them.
 
```bash
./06_nonhost/alien_nupt_detection.sh                    # all species
./06_nonhost/alien_nupt_detection.sh Species_name       # one species
RESUME=1 ./06_nonhost/alien_numt_detection.sh           # resume
```
 
Plant mitogenomes frequently carry genuinely foreign DNA via horizontal transfer,
especially in parasitic lineages, so a non-host best hit for a NUMT is not by itself
evidence of mislabelling or contamination.
 
---
 
## figures
 
| Script | Output |
|---|---|
| `plot_ks_distribution.R` | Ks histogram of cNUxTs with the dcNUxT threshold: `Rscript figures/plot_ks_distribution.R nupt` (or `numt`) |
| `fetch_taxonomy.py` | Order/family lookup from NCBI Taxonomy for clade labels |
 
Figures of the PGLS panels are written by `03_comparative/03_plot_pgls_panels.R`.
 
---
 
## References
 
- Chen Y.-C. *et al.* (2025) *Genome Research*. Predecessor vertebrate NUMT pipeline.
- Wang D. *et al.* (2010) KaKs_Calculator 2.0. *Genomics Proteomics Bioinformatics* 8:77–80.
- Yang Z. (2007) PAML 4. *Mol Biol Evol* 24:1586–1591.
- Katoh K., Standley D.M. (2013) MAFFT v7. *Mol Biol Evol* 30:772–780.
- Suyama M. *et al.* (2006) PAL2NAL. *Nucleic Acids Res* 34:W609–W612.
- Minh B.Q. *et al.* (2020) IQ-TREE 2. *Mol Biol Evol* 37:1530–1534.
- Revell L.J. (2012) phytools. *Methods Ecol Evol* 3:217–223.
- NCBI genetic codes: <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
