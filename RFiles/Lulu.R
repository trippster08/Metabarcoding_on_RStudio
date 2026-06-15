## File Housekeeping ===========================================================
# Install lulu if you haven't already. You only need to do this once.
install_github("tobiasgf/lulu")

# Load all R packages you may need, if necessary
library(rBLAST)
library(Biostrings)
library(lulu)
library(digest)
library(dplyr)
library(tibble)
## Make Lulu matchlist =========================================================
dir.create("ref/repseq_db")

rBLAST::makeblastdb(
  "data/results/PROJECTNAME_rep-seq.fas",
  db_name = "ref/repseq_db/repseq_db",
  dbtype = "nucl"
)

rBLAST::blast(db = "ref/rep_seqs/refseqdb")

# Make a DNAStringSet object from our representative sequences
sequences_dna <- BioStrings::DNAStringSet(setNames(
  repseq_nochim_md5_asv$ASV,
  repseq_nochim_md5_asv$md5
))

# You can also get this from the fasta file we downloaded earlier.
sequences_fasta <- BioStrings::readDNAStringSet("data/results/PROJECTNAME_rep-seq.fas")

# They make the same thing.
head(sequences_rep_seq)
head(sequences_fasta)


lulu_blast <- rBLAST::predict(
  refseqdb,
  sequences.dna,
  outfmt = "6 qseqid sseqid pident",
  BLAST_args = "-perc_identity 85 -qcov_hsp_perc 80"
)

lulu_matchlist <- lulu_blast %>%
  dplyr::select(qseqid, sseqid, pident)


## Run Lulu Analysis ===========================================================
# We need to have our representative sequences (the sequences we are going to
# blast, Now we have to reformat our representative-sequence table to be a named vector
View(repseq_nochim_md5_asv)

seqtab_nochim_transpose_md5_lulu <- seqtab.nochim.transpose.md5 %>%
  tibble::column_to_rownames(var = "ASV")


curated_asv <- lulu::lulu(
  seqtab_nochim_transpose_md5_lulu,
  lulu_matchlist,
  minimum_match = 84,
  minimum_relative_cooccurence = 0.95,
  minimum_ratio_type = "min",
  minimum_ratio = 1
)

# Get the feature table out of this object.
feattab_lulu <- curated_asv$curated_table

# Get the representative sequences from the feature table and add md5 hash
# (which) we have to make anew.
repseq_lulu <- feattab_lulu$

repseq_lulu_md5 <- c()
for (i in seq_along(repseq_lulu)) {
  repseq_lulu_md5[i] <- digest::digest(
    repseq_lulu[i],
    serialize = FALSE,
    algo = "md5"
  )
}


# How long did this take to run? We can check that:

print(curated_asv$runtime)

# Check how many ASVs/OTUs lulu has counted as "valid"

valids <- curated_asv$curated_count
print(valids)

# Check how many ASVs/OTUs lulu regarded as errors and discarded:

errors <- curated_asv$discarded_count
print(errors)

# Of course, the number of valid ASVs plus error ASVs should equal the total number
# of ASVs you started with:

total <- sum(valids, errors)
print(total)


# So you can do some simple math to figure out proportion of ASVs lulu has marked as
# erroneous and/or valid:

prop_error <- errors / total
print(prop_error)

prop_valid <- valids / total
print(prop_valid)

save(
  repseq_db,
  lulu_blast,
  lulu_matchlist,
  seqtab_nochim_transpose_md5_lulu,
  curated_asv,
  feattab_lulu,
  file = "data/working/lulu.RData"
)
