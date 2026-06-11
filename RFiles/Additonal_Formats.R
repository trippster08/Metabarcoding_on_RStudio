# 4 FORMAT AND EXPORT DADA2 RESULTS ############################################

# This section of the gives lots of options for creating and exporting various
# ways of looking at the DADA2 results. This section is not necessary if you
# just want to use your results directly in taxonomic assignment or phylip or
# another program. Each format below has a description that mentions why you may
# want to use a particular format.

## Load Libraries ==============================================================
# Load all R packages you may need if necessary.
library(dada2)
library(digest)
library(tibble)
library(tidyr)
library(digest)
library(seqinr)

## Create and Export Sequence-List Table =======================================
# This creates a table containing three columns: sample name, ASV, and read
# count. Each row is a separate sample/ASV combination. This is a tidy table,
# which means each column contains a single variable and each rowh a single
# observation. This is a good table format for storage of DADA2 results because
# it can be easily concatenated with other sequence-list tables in Excel or any
# text-editing software (unlike the sequence-table), yet it still contains all
# the information needed from our trimming and denoising steps.

# You will also need to make a Sequence-List table if you want to export a
# feature-to-fasta file later. A feature-to-fasta file contains every
# combination of ASV and sample, and each sequence is named with the sample
# name, md5 hash and number of reads of that ASV in that sample. It's a good
# way to look at ASV distributions in a phylogenetic tree.

# Convert the sequence-table from your DADA2 results into a tibble,
# converting row names to column 1, labeled "sample". A tibble is a more
# versatile data.frame, but it does not have row headings
# (among other differences, see https://tibble.tidyverse.org/). We'll need
# this to be a tibble for the next step.
seqtab_nochim_tb <- tibble::as_tibble(seqtab_nochim, rownames = "sample")

# The sequence-table has a column with sample names, and N columns of ASV's
# containing count values. We want all the count data to be in a single column,
# so we use a tidyr command called "pivot_longer" to make the table "tall",
# which means the table goes from 17 x 2811 to 47770 x 3 for example
# (47770 = 2810 x 17. 2810 instead of 2811 because the first column of the
# original table contains sample names, not counts). This makes the table tidier
# (meaning that each column is now a true variable).
seqtab_nochim_tall <- seqtab_nochim_tb %>%
  tidyr::pivot_longer(
    !sample,
    names_to = "ASV",
    values_to = "count"
  )
# Look at your new table.
head(seqtab_nochim_tall)
# Look at the dimensions of this table.
dim(seqtab_nochim_tall)

# Remove rows with sequence counts = 0. This removes any sample in which a
# particular ASV was not found.
seqtab_nochim_tall_nozero <- subset(seqtab_nochim_tall, count != 0)
# Look at the dimensions of this table to compare to the previous.
dim(seqtab_nochim_tall_nozero)

# Export Sequence-List Table. This includes 3 columns: sample name, ASV
# sequence, and read count for that sample name ASV combo. This is a tidy table
# (each column a variable), and is a good way to include all the possible data
# involved.
write.table(
  seqtab_nochim_tall_nozero,
  file = "data/results/PROJECTNAME_SeqList_Tall.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

## Create and Export feature-to-fasta ==========================================
# This creates a fasta file containing all the ASV's for each sample. Each ASV
# will be labeled with the sample name, ASV hash, and number of reads of that
# ASV in that sample. This was derived from a python script from Matt Kweskin
# called featuretofasta.py (hence the name).

# Save the ASV sequences from the sequence-list table
# (seqtab_nochim_tall_nozera) as a new list.
repseq_tall <- seqtab_nochim_tall_nozero$ASV
# Look at the top of your new list
head(repseq_tall)

# Convert the sequences into md5 hashs, as we did earlier. md5 hashs are
# consistent across jobs, meaning identical sequences from different projects or
# being converted by different programs will result in the same hash (i.e.
# hashs here will match hashs above)
repseq_tall_md5 <- c()
for (i in seq_along(repseq_tall)) {
  repseq_tall_md5[i] <- digest::digest(
    repseq_tall[i],
    serialize = FALSE,
    algo = "md5"
  )
}
# Examine hashed sequences to make sure conversion worked
head(repseq_tall_md5)

# Attach the ASV hashes as a column (called "md5") to the tall table. The
# table should now have 4 columns, and each row of the "md5" column should
# be a md5 hash of its respective ASV.
seqtab_nochim_tall_nozero_md5 <- cbind(
  seqtab_nochim_tall_nozero,
  md5 = repseq_tall_md5
)
# Check to make sure new table contains "feature" column
colnames(seqtab_nochim_tall_nozero_md5)

# Create a new column in this table that contains "sample", "feature", and
# "count", concatenated. This is the heading for each sequence in the fasta file
# created by Matt Kweskin's script "featuretofasta.py"
seqtab_nochim_tall_nozero_md5_name <- tidyr::unite(
  seqtab_nochim_tall_nozero_md5,
  sample_feature_count,
  c(sample, md5, count),
  sep = "_",
  remove = FALSE
)
# Check to make sure the new names are correct, for at least the first 4 rows.
# And check column headings if you like.
head(seqtab_nochim_tall_nozero_md5_name)[1:4, 1]
colnames(seqtab_nochim_tall_nozero_md5_name)

# Reorder columns if you want, it doesn't really matter for the next step, but
# it may matter in how you look at the table. You could even remove the "sample"
# "feature" and "count" columns since they are not needed anymore in the table
# (by using [,c(1,3)]), but it would only be for ease of visualizing the table,
# and wouldn't change the next "write.fasta" step. You can also export this
# table with reordered columns as a sequence-list table with an additional
# column of ASV hashes. Follow the directions for exporting the sequence-list
# table.
seqtab_nochim_tall_nozero_md5_name_reorder <- seqtab_nochim_tall_nozero_md5_name[, c(
  1,
  3,
  2,
  5,
  4
)]
# Check reordering
colnames(seqtab_nochim_tall_nozero_md5_name_reorder)

# Create a fasta-formatted file of each row sequence (i.e. ASV), with a heading
# of "sample_feature_count".
seqinr::write.fasta(
  sequences = as.list(seqtab_nochim_tall_nozero_md5_name$ASV),
  names = seqtab_nochim_tall_nozero_md5_name$sample_feature_count,
  open = "w",
  as.string = FALSE,
  file.out = "data/results/PROJECTNAME_sample-feature-count.fas"
)
