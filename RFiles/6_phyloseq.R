# 7 PHYLOSEQ ###################################################################
# Phyloseq relies on a phyloseq object for its operations. Phyloseq objects
#  consist of 5 components, but not all are required, only ones that are needed
# for the desired operations. The 5 components are: otu_table, sample_data,
# tax_table, phy_tree, and refseq (not sure why the last doesn't get an
# underscore). We will go into more details about the components below, when we
# configure them.

## File Housekeeping ===========================================================
# Load all R packages you may need, if necessary

library(dada2)
library(digest)
library(phyloseq)
library(dplyr)
library(tibble)
library(readr)
library(seqinr)
library(ape)
library(DECIPHER)
library(ade4)
library(Biostrings)
library(bioseq)
library(phangorn)
# Set up your working directory. If you created your new project in the
# directory you want as your working directory (or came directory from the
# previous step in the pipeline), you don't need to do this, and
# skip to the next RStudio command. If you need to set your working directory,
# substitute your own path for the one below.
setwd("/Users/smithb/Dropbox (Smithsonian)/Projects_Metabarcoding/PROJECTNAME")

# We run this on a gene-specific basis, so you will do this for each gene you
# wish to analyze.
## GENE 1
## Prepare Components to be Imported Into Phyloseq =============================

# Load your R objects from the DADA2 analysis.
load("data/working/9_output_GENE.RData")
# Load your taxonomy results
load("data/working/4_Assign_Taxonomy_GENE.RData")
### sequence-table -------------------------------------------------------------


# Make phyloseq otu_table from the sequence-table (columns of ASV md5 hashes,
# rows of samples). If you want to use a feature-table (columns of samples,
# rows of ASV/OTUs) instead, use "taxa_are_rows = TRUE"
OTU_md5 <- phyloseq::otu_table(seqtab_nochim_md5, taxa_are_rows = FALSE)

### tax_table ------------------------------------------------------------------
# Our current dada2 taxonomy table has a taxonomy table and a bootstrap table.
# Phyloseq only needs the taxonomy portion. Also, this table has ASVs, while
# we would prefer md5 hashes. We'll do a left join to replace ASVs with md5
# hashes. We will also convert to a tibble for easier manipulation.

# Make a new taxonomy-only table, and replace the current rownames (ASVs) with
# md5 hashes.
taxonomy_tax_md5 <- tibble::as_tibble(taxonomy$tax, rownames = "ASV") %>%
  dplyr::left_join(
    tibble::as_tibble(repseq_nochim_md5_asv),
    by = "ASV"
  ) %>%
  dplyr::relocate(md5, .before = 1) %>%
  tibble::column_to_rownames("md5")

# Make phyloseq tax-table from our taxonomy-only table. Phyloseq requires this
# in matrix format, so we'll convert.
TAX_md5 <- phyloseq::tax_table(as.matrix(taxonomy_tax_md5))

### sample_data ----------------------------------------------------------------
# Import metadata (here as a tab-delimited text file, see examples for
# formatting) and convert to a sample_data object.

# Import your metadata file. I usually use a tab-delimited file (sep = "\t"),
# but you can also use a comma delimited metadata file (sep = ","). Your sample
# names need to match the sample names in the sequencing run. You may need
# to define some columns whose type may be incorrectly determined (some
# columns may contain only numbers and be recognized as a numeric column,
# but the values are actually discreet, and should be characters instead)
meta <- readr::read_tsv(
  "metadata.tsv",
  col_types = cols(
    replicate = col_character(),
    depth = col_character()
  ),
  show_col_types = TRUE
)

# Also, you may have some variables that are numbers but you want to keep as
# characters. They contain discreet variables, but because they are numbers,
# R reads them as continuous. For example, we may have two filter sizes,
# 22 and 45. While these are numbers, "filter_size" is not a continuous
# variable, we only have two discreet sizes. To ensure that R recognizes these
# appropriately (as characters and not numbers), I have added an optional
# "colClasses = c()" argument, which defines any column (or multiple columns)
# as a particular data type. This is important when looking at plots downstream.
# To check data type for a all columns in your table, use "str(metadata)".
metadata <- phyloseq::sample_data(read.delim(
  "data/working/PROJECTNAME_metadata.tsv", 
  sep = "\t", 
  header = TRUE,
  colClasses = c(water_replicate = "character", filter_size = "character"),
  row.names = "sample_name"
))

# Make sure that the metadata has the same samples (and only the
# samples, no extras) as the sequence table. We will do a right
# join with sample_names_filtered (which contains the names of the samples
# that have made it through the pipeline) to make sure these are the same.
meta_filtered <- meta %>%
  filter(., SeqID %in% sample_names_filtered) %>%
  tibble::column_to_rownames("SeqID")


# Look at data type of all the columns of the table.
str(metadata)
# Make a phyloseq sample_data file from "metadata"
SAMPLE <- phyloseq::sample_data(meta_filtered)


### refseq ---------------------------------------------------------------------
# The refseq phyloseq-class item must contain sequences of equal length, which
# in most cases means it needs to be aligned first. We will align using
# DECIPHER.
# We already have a list of ASVs in "repseq", which we obtained using a DADA2
# command called getSequences. getSequences extracts sequences from a DADA2
# object, which in the case of "repseq" was "seqtab.nochim". However, the
# ASVs in "repseq" are not named. Since we have been using md5 hashes as ASV
# names up to this point, we should do the same here, using the list we already
# made called "repseq.md5".

# If you have not made "repseq" or "repseq.md5" go to
# "4 Metabarcoding_R_Pipeline_RStudio_FormatandExportFiles.txt", section "Create
# And USE md5 Hash".

# Make a new list of ASVs from the representative sequences, and add md5 hashes
# as names.
sequences <- repseq
names(sequences) <- repseq.md5

# Convert these sequences into a DNAString, which is the format of sequences
# used by DECIPHER, and many other phylogenetic programs in R.
sequences_dna <- Biostrings::DNAStringSet(sequences)

# Align using DECIPHER. DECIPHER "Performs profile-to-profile alignment of
# multiple unaligned sequences following a guide tree" (from the manual). We do
# not give a preliminary guide tree, so one is automatically created. For each
# iteration, a new guide tree is created based on the previous alignment, and
# the sequences are realigned. For each refinement, portions of the sequences
# are realigned to the original, and the best alignment is kept. useStructures
# probably should be FALSE if you are using COI or another protein-coding
# region. If you are using RNA, then "useStructures=TRUE" may give a better
# alignment. However, since our sequences are in DNA format, you first have
# to convert sequences.dna into an RNAStringSet. Alignments can take a long time
# if you have lots of ASVs. Use more refinements and/or interations to get a
# "better" alignment, but increasing these will take more computing time.

# Make an alignment from your DNA or RNA data, and change "useStructures" 
# accordingly.
alignment_dna <- DECIPHER::AlignSeqs(
  sequences_dna,
  guideTree = NULL,
  anchor = NA,
  gapOpening = c(-15, -10),
  gapExtension = c(-3, -2),
  terminalGap = c(-15, -10),
  iterations = 2,
  refinements = 3,
  processors = NULL,
  useStructures = FALSE
)
# Look at a brief "summary" of the alignment. This shows the alignment length,
# the first 5 and last 5 ASVs and the first 50 and last 50 bps of the alignment.
alignment_dna

# You can also look at the entire alignment in your browser.
DECIPHER::BrowseSeqs(alignment_dna)

# Create a reference sequence (refseq) object from the alignment. This contains
# the ASV sequences, using the md5 hashes as names.
REFSEQ_md5 <- Biostrings::DNAStringSet(alignment_dna, use.names = TRUE)
# Look at the refseq object, just to make sure it worked
head(REFSEQ_md5)

### phy_tree -------------------------------------------------------------------
# We can create a phylogenetic (or at least, a phenetic) tree using the
# alignment we just created using the program ape.

# For ape, the aligned sequences must be in binary format (DNAbin, which reduces
# the size of large datasets), so we first convert the DNAstring alignment.
alignment_dnabin <- bioseq::as_DNAbin(alignment_dna)

#Create pairwise distance matrix in ape. There are many different models to use.
# Here we are using the Tamura Nei 93 distance measure. Turning the resulting
# distances into a matrix (using as.matrix = TRUE) results in a table of
# pairwise distances. Using "as.matrix = FALSE" results in a
# distances are meaningless to look at, but either type can be used for
# tree-building).
pairwise_tn93 <- ape::dist.dna(
  alignment_dnabin,
  model = "tn93",
  as.matrix = TRUE
)

# NOTE: If your sequences are highly divergent, pairwise distance will not be
# calculatable by dist.dna, and you must use another method. I rarely have
# a distance matrix that does not contain any NaN's (the result for pairs
# without a distance). ape has a alternative tree-building command for each
# method that is meant to deal with a some NaNs. However, if there are too many
# NaNs, tree-building will not work well. In this case, we can obtain
# maximum-likelihood distances using the program phangorn, which seems to be
# able to obtain distances even from highly divergent sequences.

# Check the number of NaNs in dist.dna, and the proportion of all distance. I
# don't know how many NaNs are too many, but if there are more than a few, I
# would prefer to be safe and use ml distances.
length(is.nan(pairwise_tn93))
length(pairwise_tn93) / length(is.nan(pairwise_tn93))

# Obtain pairwise distances with the "dist.ml" command. It currently only has
# two models of evolution: JC69 and F81.
# First, "dist.ml" requires the alignment to be in the phyDat format, which can
# be converted from the dnabin format (but not the DNAStringSet format) using
# "as.phyDat".
alignment.phyDat <- phangorn::as.phyDat(alignment.dnabin)

pairwise.ml <- phangorn::dist.ml(
    alignment.phyDat,
    model = "F81"
)

# Check the number of NaNs in the ml distances.
length(is.nan(pairwise.ml))

# Make an improved neighbor-joining tree out of our pairwise distance matrix.
# ape has other tree-building phenetic methods to use as well, such as nj or
# upgm. If there are any NaNs in your data, use "bionjs", otherwise us "bionj".
# If you used a different distance model (ml, K80, F84, etc) replace "tn93"
# with the model used.
tree.tn93.bionj <- ape::bionj(pairwise.tn93)

# Look at the tree. Of course, if there are a lot of ASV's, the tree is pretty
# much indecipherable, even with lots of editing using ggplot2. We will look at
# the tree in greater detail from the phyloseq object, below.
plot(tree.tn93.bionj)

# Create a phyloseq tree object from our neihbor-joining tree.
TREE.md5 <-  phyloseq::phy_tree(tree.tn93.bionj)

### phyloseq object ------------------------------------------------------------
# We use the components created above to create a phyloseq object. For the
# otu_table, we need to tell phyloseq the orientation of the table (Samples as
# rows vs. ASV's, here labelled "taxa", as rows). Remember, the default output
# of dada2 is taxa as rows. You can use our transposed table (feature-table),
# and make "taxa_are_rows = TRUE", but it didn't work as well for me, so I
# suggest you stick with the default format.

physeq <- phyloseq::phyloseq(
  OTU.md5,
  TAX.md5,
  SAMPLE.md5,
  REFSEQ.md5,
  TREE.md5
)


# There are lots of parameters of the phyloseq object you can look at. Part of
# the reason for looking at these is to make sure the values are what you
# expect.
phyloseq::ntaxa(physeq)
phyloseq::nsamples(physeq)
phyloseq::sample_names(physeq)
phyloseq::sample_variables(physeq)
phyloseq::otu_table(physeq)[1:5,1:5]
phyloseq::tax_table(physeq)[1:5,1:5]
phyloseq::phy_tree(physeq)
phyloseq::taxa_names(physeq)[1:20]

# We can also look at the tree in greater detail, including adding labels that
# show factor variables for each branch. "color", "shape", and "size" are all
# terminal branch labels that can show variable values.  "base.spacing" and
# "min.abundance" are for making these labels easier to see when there are a 
# lot.
phyloseq::plot_tree(
  physeq,
  ladderize = "left",
  color = "factor1",
  shape = "factor2",
  size = "abundance",
  base.spacing = 0.03,
  min.abundance = 2,
  label.tips = "Class"
)



phyloseq::plot_richness(physeq, x="factor1", measures=c("Shannon", "Fisher"), color = "factor2") 
ord.euclidean <- phyloseq::ordinate(physeq, "MDS", "euclidean")
phyloseq::plot_ordination(physeq, ord.euclidean, type = "samples", color = "factor1")
phyloseq::plot_ordination(physeq, ord.euclidean, type = "samples", color = "factor2")
