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
00_data/                   Organelle accession audit, gene sets, assembly statistics
01_detection/              NUMT/NUPT detection, counts and genome proportions
02_cnuxt_kaks/             Coding insertions, Ka/Ks, dcNUxT definition, gene-level summaries
03_comparative/            Phylogenetic signal and PGLS on insertion content
04_clustering_selection/   Shared-origin clustering of dcNUxTs and codeml branch tests
05_dcnuxt_integration/     Genomic context of dcNUxTs, host-gene structure, RNA-seq junctions
06_nonhost/                Non-host (alien) dcNUxT detection
tables/                    Supplementary tables
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
| samtools | | 01, 05, 06 |
| HISAT2, sra-tools | | 05 |
| Python 3 | | Biopython, pandas, numpy, scipy, dendropy, matplotlib |
| R | | ape, nlme, dplyr, phytools, geiger, ggplot2 |
| NCBI `datasets` CLI | | 03 (`assembly_technology.py --fetch-tech`) |

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
│   ├── plastid_gb/  mito_gb/       GenBank records for gene extraction
│   ├── plastid_genes/  mito_genes/ <ACCESSION>/<gene>.fasta (one CDS per file)
│   └── mito_genes_merged/          <SPECIES>/<gene>.fasta (multi-chromosome mitogenomes)
├── multi_organelles/downloads/     organelle FASTAs of multi-chromosome mitogenomes
└── species/<SPECIES>/
    ├── numts/  nupts/              BLAST output, fragment BEDs, insertion FASTA
    ├── fragments/  alignments/     per-fragment FASTA and MAFFT alignments (stage 02)
    └── results/                    per-species counts, proportions, Ka/Ks tables
```
---

## 01_detection — NUMT/NUPT detection
 
| Script | Purpose |
|---|---|
| `numt_nupt_detection.sh` | BLASTN of the mitogenome and plastome against the nuclear assembly, removal of hits on organelle sequences deposited with the assembly, fragment extraction, merging, counts, genome proportions, NUMT–NUPT overlap |
| `get_bed2.R` | Called by the above: keeps hits ≥ 200 bp, orients and sorts them, writes the fragment BED |
| `mtpt_coverage.sh` | Fraction of each plastome present in the species' own mitogenome (MTPTs), which explains sequence called by both queries |
| `build_species_summary.py` | One row per species from the per-species result files → `final_species_summary_v2.tsv`, the table all later stages read. Excludes the two species with partial nuclear assemblies (*Tilia amurensis*, *Tripidium arundinaceum*). Rerun after stage 02 step 3 to add cNUxT/dcNUxT counts |
 
```bash
THREADS=24 WORKING_DIR=/vol/storage/5T_disk/organelles \
  ./01_detection/numt_nupt_detection.sh combined_dataset.csv
./01_detection/mtpt_coverage.sh combined_dataset.csv
python3 01_detection/build_species_summary.py combined_dataset.csv
```
 
Parameters:
- `blastn -word_size 20 -max_target_seqs 1000000`. The BLAST default of 500 target
  sequences truncates detection in fragmented assemblies, where insertions are spread
  over thousands of scaffolds.
- Multi-chromosome mitogenomes (`MITO` accessions separated by `;`) are searched as one
  concatenated query.
- Hits on assembly sequences whose FASTA header annotates them as organelle
  (mitochondrion, chloroplast, plastid, organelle) are removed. These are the species'
  organelle genomes deposited with the nuclear assembly, not insertions. The removed
  sequences are listed per species in `results/<SPECIES>.organelle_like_sequences.tsv`.
- Hits shorter than 200 bp are discarded (`MIN_LENGTH`). Fragments on the same sequence
  up to 10 kb apart are merged into one insertion with `bedtools merge -d 10000`
  (`MERGE_DISTANCE`).
- Proportions are insertion bp / assembly bp, both for individual fragments
  ("unmerged") and merged insertions.
 
Outputs in `species/<SPECIES>/results/`: `*.numt_count.txt`, `*.numt_merged_count.txt`,
`*.numt_unmerged_prop.txt`, `*.numt_merged_prop.txt`, `*.numt_{merged,unmerged}_locations_lengths.tsv`
(and NUPT equivalents), `*.{numt,nupt}_union_bp.txt`, `*.genome_size.txt`, `*.organelle_annotated_headers.txt`, `*.organelle_like_sequences.tsv`, and the
NUMT–NUPT overlap: `*.numt_nupt_overlap_bp.txt` (bp shared by merged insertions),
`*.numt_nupt_overlap_count.txt` (NUMT insertions overlapping a NUPT insertion) and
`*.numt_nupt_shared.bed` (sequence called by both queries at fragment level). Insertion
sequences go to `numts/numts.seq.fa` and `nupts/nupts.seq.fa`. `mtpt_coverage.sh` writes
`species/mtpt_coverage.tsv`.
 
---
 
## 02_cnuxt_kaks — coding insertions and Ka/Ks
 
Run in this order:
 
| # | Script | Purpose |
|---|---|---|
| 1 | `01_nupt_gene_alignment.sh` / `01_numt_gene_alignment.sh` | ORF filter (≥ 90 bp in any of six frames) → BLASTN pre-filter against the host's organelle genes (≥ 70 % identity, ≥ 50 % query cover, e ≤ 1e-5) → MAFFT alignment of each fragment and its reverse complement to each candidate gene → gap-free, in-frame, stop-free codon alignments → KaKs_Calculator (LWL). Multi-chromosome mitogenomes use the merged gene set |
| 2 | `02_dedup_kaks_rows.py` | Collapse Ka/Ks rows for the same locus (≥ 80 % reciprocal overlap) to one row; dry run by default, `--apply` rewrites (backup `.pre_dedup`) |
| 3 | `03_combine_kaks.py` | Combine per-species tables into `all_species_{nupt,numt}_kaks.tsv`; species with < 10 annotated genes are left out of coding analyses (`coding_species_status.tsv`) |
| 4 | `04_build_funnel_table.py` | Table 1 (NUxT → cNUxT → dcNUxT → unique → purifying); unique set = identical alignment fragment within species × gene collapsed → `dcnuxt_unique.tsv` |
| 5 | `05_bh_correction_dcnuxt.py` | Benjamini–Hochberg correction of the Fisher exact test p-values on the unique dcNUxTs, separately for NUPTs and NUMTs |
| 6 | `06_mannwhitney_kaks_comparison.py` | Compare dcNUPT and dcNUMT Ka/Ks distributions (Mann–Whitney U, rank-biserial r) |
| 7 | `07_build_gene_lengths.py` | Reference length of each organellar gene, averaged across accessions |
| 8 | `08_gene_level_summary.py` | Per gene: cNUxT and dcNUxT counts, species breadth, fraction old (ret%), fraction under purifying selection (sel%) |
| 9 | `09_gene_predictors.py` | Gene length vs transfer and retention; length-normalised transfer vs retention; NUPT vs NUMT comparison (Fisher r-to-z) |

Unnamed ORFs (`orf…`, `hypothetical`, `ycf`-style locus tags without a gene name) are
counted in Table 1 but left out of gene-level analyses (steps 8–9), clustering and
branch models, because they cannot be matched across species.
 
```bash
./02_cnuxt_kaks/01_nupt_gene_alignment.sh combined_dataset.csv
./02_cnuxt_kaks/01_numt_gene_alignment.sh combined_dataset.csv
python3 02_cnuxt_kaks/02_dedup_kaks_rows.py combined_dataset.csv --apply
python3 02_cnuxt_kaks/03_combine_kaks.py combined_dataset.csv
python3 01_detection/build_species_summary.py combined_dataset.csv    # adds cNUxT/dcNUxT counts
python3 02_cnuxt_kaks/04_build_funnel_table.py final_species_summary_v2.tsv
python3 02_cnuxt_kaks/05_bh_correction_dcnuxt.py all_species_nupt_kaks.tsv all_species_numt_kaks.tsv
python3 02_cnuxt_kaks/06_mannwhitney_kaks_comparison.py dcnupt_final_bh.tsv dcnumt_final_bh.tsv
python3 02_cnuxt_kaks/07_build_gene_lengths.py
python3 02_cnuxt_kaks/08_gene_level_summary.py --latex
python3 02_cnuxt_kaks/09_gene_predictors.py --min-cnuxt 5 --min-dcnuxt 3
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
| `common.R` | Shared data loading (species summary, tree, parasite flag), λ-profiled PGLS | — |
| `01_phylo_signal.R` | Blomberg's K (tip randomisation) and Pagel's λ (LRT vs 0 and vs 1) for NUMT, NUPT, NUPT:NUMT ratio and scaffold N50 | `phylo_signal_table.tsv` |
| `02_pgls_models.R` | PGLS (`nlme::gls` + `corPagel`, λ profiled on a 0.02 grid because free estimation fails to converge for several models) with 95 % profile CI and OLS alongside | `pgls_table.tsv`, `lambda_profiles.tsv`, `pgls_fits.rds` |
| `03_plot_pgls_panels.R` | Scatter panels with PGLS fit and 95 % band | `Figures/*.pdf` |
| `tech_effect_test.R` | Effect of sequencing technology (long-read / hybrid / short-read) on insertion content: rank tests, OLS, PGLS, organelle-scaffold sensitivity, within-genus comparison | `tech_effect/` |
| `assembly_quality.R` | Whether scaffold N50 and technology explain insertion content and the NUMT–NUPT difference in λ | stdout |
| `raw_scaffold_composition.py` | Organellar fraction of every scaffold carrying merged insertions; flags species whose content may be inflated by organellar contigs retained in the assembly | `raw_scaffold_composition.tsv` |
| `assembly_technology.py` | Contiguity from `.fai` indices and declared sequencing technology from NCBI (`--fetch-tech`) | `assembly_technology.tsv` |
| `tech_insertion_metrics.py` | Per species: high-identity hits and long merged insertions, for the technology test | `insertion_metrics.tsv` |
| `mtpt_sensitivity.R` | PGLS rerun with sequence called by both queries (MTPT-explained) removed from NUMT and from NUPT content | `mtpt_sensitivity/` |
 
```bash
Rscript 03_comparative/01_phylo_signal.R      # PHYLOSIG_NSIM=10000 for finer K P-values
Rscript 03_comparative/02_pgls_models.R
Rscript 03_comparative/03_plot_pgls_panels.R
python3 03_comparative/raw_scaffold_composition.py
python3 03_comparative/assembly_technology.py --fetch-tech
python3 03_comparative/tech_insertion_metrics.py
Rscript 03_comparative/tech_effect_test.R final_species_summary_v2.tsv assembly_technology.tsv
Rscript 03_comparative/assembly_quality.R
Rscript 03_comparative/mtpt_sensitivity.R
```

Set `OMP_NUM_THREADS=1` / `OPENBLAS_NUM_THREADS=1` for the R scripts; multithreaded
BLAS slows the many small PGLS fits. Parasitic genera are flagged in `common.R`
(`PARASITE_GENERA`).
 
---
 
## 04_clustering_selection — shared-origin clusters and codeml
 
Separate NUPT and NUMT runs, writing to `clustering_{nupt,numt}/` and
`selection_{nupt,numt}/`. The Python steps take `--organelle nupt|numt` (default
`nupt`).
 
| # | Script | Purpose |
|---|---|---|
| 1 | `01_cluster_dcnupts.sh` / `01_cluster_dcnumts.sh` | Collect named-gene dcNUxTs per gene (identical copies within a species collapsed; dcNUMTs annotated as *atpA* pooled with *atp1*), pairwise Ks between all dcNUxTs of a gene, link pairs with Ks < 0.1 where neither copy is ≥ 98 % identical to any organelle genome in the pooled database, connected components → clusters and singletons |
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
finished runs unless `--force` is given, so re-running it only re-parses results. Each codeml run has a 1 h time limit.
 
---
 
## 05_dcnuxt_integration — genomic context and transcription
 
| # | Script | Purpose |
|---|---|---|
| 1 | `dcnuxt_annotation_overlap.py` | Map the unique dcNUxTs (`--final-set`) back to genome coordinates, classify against the GFF (CDS > UTR/exon > intron > pseudogene > intergenic) and compare with a within-scaffold permutation null (per-replicate counts, two-sided empirical P). Loci that cannot be placed on their scaffold are left out of both observed and expected. `--audit-only` checks GFF availability and seqid concordance first |
| 2 | `dcnuxt_exon_cluster.py` | For loci inside annotated genes: host-gene exon structure, fraction of the gene formed by the insertion, whether the annotation merely restates the donor gene; verdicts CANDIDATE / CIRCULAR / SINGLE_EXON / INTRONIC_ONLY / ORGANELLAR_CONTIG; writes RNA-seq target BED/FASTA |
| 3 | `dcnuxt_annotation_overlap.py --exclude-organellar` | Sensitivity rerun without loci on ORGANELLAR_CONTIG scaffolds |
| 4 | `dcnuxt_junction_rnaseq.sh` | HISAT2 mapping of public RNA-seq (with organellar decoys) and counts of reads spanning the nuclear-flank/insertion junctions — the only reads that can come from the nuclear copy |

```bash
python3 05_dcnuxt_integration/dcnuxt_annotation_overlap.py \
    --base /vol/storage/5T_disk/organelles --gff-dir downloads/genome_fasta \
    --final-set dcnuxt_unique.tsv --perms 1000 --out-prefix dcnuxt_overlap
python3 05_dcnuxt_integration/dcnuxt_exon_cluster.py \
    --loci dcnuxt_overlap_loci.tsv --csv combined_dataset.csv \
    --gff-dir downloads/genome_fasta --genome-dir downloads/genome_fasta \
    --flank 300 --out-prefix dcnuxt_overlap
python3 05_dcnuxt_integration/dcnuxt_annotation_overlap.py \
    --base /vol/storage/5T_disk/organelles --gff-dir downloads/genome_fasta \
    --final-set dcnuxt_unique.tsv --exclude-organellar dcnuxt_overlap_exon_cluster.tsv \
    --perms 1000 --out-prefix dcnuxt_overlap_noorg
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
./06_nonhost/nonhost_scaffold_composition.sh
```

`nonhost_scaffold_composition.sh` restricts the category tables to the unique dcNUxT
set (`dcnuxt_unique.tsv`) and measures the organellar fraction of each contig carrying
a high-confidence non-host candidate (dc-megablast against the pooled database of the
same class; contigs > 300 kb treated as nuclear). 
 
Plant mitogenomes frequently carry genuinely foreign DNA via horizontal transfer,
especially in parasitic lineages, so a non-host best hit for a NUMT is not by itself
evidence of mislabelling or contamination.
 

## figures
 
| Script | Output |
|---|---|
| `plot_ks_distribution.R` | Ks histogram of cNUxTs with the dcNUxT threshold: `Rscript figures/plot_ks_distribution.R nupt` (or `numt`) |
| `fetch_taxonomy.py` | Order/family lookup from NCBI Taxonomy for clade labels |
| `plot_proportion_fan.py` | Fan phylogenetic tree coloured by log10(proportion) |
| `make_colorbar.py` |  Generate a standalone colourbar PDF to go with the fan plots. |
Figures of the PGLS panels are written by `03_comparative/03_plot_pgls_panels.R`.
 
---
 
## Run order

00_data → 01_detection (detection, MTPT coverage, species summary) → 02 steps 1–3 →
species summary again → 02 steps 4–9 → 03_comparative → 04_clustering_selection →
05_dcnuxt_integration → 06_nonhost → tables, figures.

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
