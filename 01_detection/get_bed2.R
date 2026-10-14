# get_bed2.R
# Called by numt_nupt_detection.sh from inside species/<SPECIES>/{numts,nupts}/.
# Adapted from Yu-Chi Chen's vertebrate NUMT pipeline.
#
# Input : files listed in files_list.txt; each is BLAST tabular output reduced to
#         subject, identity, alignment length, sstart, send (cut -f2,3,4,9,10).
# Output: <file>_AllFrags.bed  hits >= min_length bp, start < end, sorted by
#         sequence and start. Coordinates are the BLAST subject coordinates as
#         reported (start column = smaller BLAST coordinate), as in the original
#         pipeline; fragment IDs downstream (<seq>_<start>-<end>) depend on this.
#
# Merging of nearby fragments into insertions is done afterwards with
# `bedtools merge -d`. 
#
# Usage: Rscript get_bed2.R <min_length>

options(scipen = 999)   # never write coordinates in scientific notation
args <- commandArgs(trailingOnly = TRUE)
ins_length <- as.numeric(args[1])

files <- read.table("files_list.txt", header = FALSE, stringsAsFactors = FALSE)

for (inp in seq_len(nrow(files))) {
  f <- files[inp, 1]
  if (file.info(f)$size == 0) next
  data <- read.table(f, header = FALSE, stringsAsFactors = FALSE,
                     colClasses = c("character", "numeric", "numeric", "numeric", "numeric"),
                     comment.char = "", quote = "")
  dat <- data[data[, 3] >= ins_length, ]
  if (nrow(dat) == 0) next

  # orient reverse-strand hits so that start < end
  s <- pmin(dat[, 4], dat[, 5])
  e <- pmax(dat[, 4], dat[, 5])
  bed <- data.frame(chrom = dat[, 1], start = s, end = e, stringsAsFactors = FALSE)
  bed <- bed[order(bed$chrom, bed$start), ]

  write.table(bed, paste0(f, "_AllFrags.bed"),
              quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
}
