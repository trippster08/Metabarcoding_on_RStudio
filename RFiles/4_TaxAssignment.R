# 8 ASSIGN TAXONOMY ############################################################
# This is a work in progress and in the prelimiary stages of development, please
# contact me for any help or questions.

## Load Libraries = ============================================================
# Load all R packages you may need, if necessary

library(dada2)
library(digest)
library(rBLAST)
library(tidyverse)
library(seqinr)
library(R.utils)

## Get Reference Database ======================================================
# If you have your own reference database, you can either enter the path when
# we assign taxonomy below, or move it into the ref/ directory in the main
# project directory

# To use a supplied reference library, download into the ref folder. However, I
# strongly recommend using your own reference.

# Use this link for a full reference database based on the Midori COI database
# (www.reference-midori.info)
reference_url <- "https://www.dropbox.com/scl/fi/7ba1096klgr11n70f8zkx/MIDORI2_UNIQ_NUC_GB260_CO1_DADA2_noInsect.fasta.gz?rlkey=ehr1f9de5x1llr05u89ko07qi&dl=1"
# Use this link for a much reduced reference database. This database only
# contains a single sequence for each Genus in the original Midori pipeline. It
# will not give you great taxonomic assignment, and is really only useful for
# demonstrative purposes (to see if it works) or for a very rough assignment.
reference_url <- "https://www.dropbox.com/scl/fi/1d2ibge46ytbpqoxxry5w/midori_COI_genus_dada2.fasta.gz?rlkey=9lknfo5e03ndzfg86qrml0rg2&dl=1"
# Specify destination file and download
ref_gzip <- sub("\\.fasta\\.gz.*$", ".fasta.gz", basename(reference_url))
download.file(reference_url, paste0("ref/", ref_gzip), mode = "wb")

# Specify decompressed file and unzip
reference_fasta <- sub("\\.gz$", "", paste0("ref/", ref_gzip))
R.utils::gunzip(paste0("ref/", ref_gzip), destname = reference_fasta, remove = TRUE)

## Assign Taxonomy With DADA2 ==================================================
# We first assign taxonomy using the assignTaxonomy function in DADA2.
# Use whatever table you obtained
# for your earlier analyses (e.g. seqtab.nochim). This command can take a long
# time for large reference databases and large numbers of ASVs. Your reference
# file will have a different path and name then here, please make sure to use
# the correct path and name.
# From the dada2 manual: "assignTaxonomy(...) expects a training fasta file (or
# compressed fasta file) in which the taxonomy corresponding to each sequence is
# encoded in the id line in the following fashion (the second sequence is not
# assigned down to level 6)." Note that all levels above the lowest level must
# have a space, even if there is no name for that level (I've added a third
# example to demonstrate).
#>Level1;Level2;Level3;Level4;Level5;Level6;
#ACCTAGAAAGTCGTAGATCGAAGTTGAAGCATCGCCCGATGATCGTCTGAAGCTGTAGCATGAGTCGATTTTCACATTCAGGGATACCATAGGATAC
#>Level1;Level2;Level3;Level4;Level5;
#CGCTAGAAAGTCGTAGAAGGCTCGGAGGTTTGAAGCATCGCCCGATGGGATCTCGTTGCTGTAGCATGAGTACGGACATTCAGGGATCATAGGATAC"
#>Level1;Level2;;Level4;Level5;Level6
#CGCTAGAAAGTCGTAGAAGGCTCGGAGGTTTGAAGCATCGCCCGATGGGATCTCGTTGCTGTAGCATGAGTACGGACATTCAGGGATCATAGGATAC"

# We want to define the taxonomic levels of our reference database.
# taxLevels defines what taxonomic rank each of the levels shown in the above
# example represents.
tax_levels <- c(
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "species"
)
# Assign taxonomy
# If you used your own database, change the path after "seqtab_nochim" to your
# reference database.
# tryRC determines whether to also include the reverse
# complement of each sequence. outputBootstraps results in a second table with
# bootstrap support values for each taxonomic assignment. minBoot gives the
# minimum bootstrap required to assign taxonomy.
taxonomy <- dada2::assignTaxonomy(
  seqtab_nochim,
  reference_fasta,
  taxLevels = tax_levels,
  tryRC = FALSE,
  minBoot = 50,
  outputBootstraps = TRUE,
  multithread = TRUE,
  verbose = TRUE
)
# Save this object, it took a long time to get.
save(taxonomy, file = "data/working/tax_rdp.Rdata")

## Examine and Manipulate Taxonomy =============================================
# Look at the taxonomic assignments
View(taxonomy$tax)
View(taxonomy$boot)

# You can check to see all the unique values exist in each column
unique(taxonomy$tax[, "Phylum"])
unique(taxonomy$tax[, "Class"])
unique(taxonomy$tax[, "Order"])
unique(taxonomy$tax[, "Family"])
# You can also get see how many times each unique value exists in that column
# using table
table(taxonomy$tax[, "Phylum"])

### Combine taxonomy and bootstrap tables --------------------------------------
# Join the two tables using an inner-join with dplyr (it shouldn't matter here
# what kind of join you use since the two tables should have the exact same
# number of rows and row headings (actually, now column 1)). I amend bootstrap
# column names with "_boot" (e.g. the bootstrap column for genus would be
# "Genus_boot"). I also add the md5 hash, and rearrange the columns
taxonomy_rdp <- dplyr::inner_join(
  tibble::as_tibble(taxonomy$tax, rownames = "ASV"),
  tibble::as_tibble(taxonomy$boot, rownames = "ASV"),
  by = "ASV",
  suffix = c("", "_boot")
) %>%
  dplyr::mutate(md5 = repseq_nochim_md5_asv$md5) %>%
  dplyr::select(
    md5,
    ASV,
    Phylum,
    Phylum_boot,
    Class,
    Class_boot,
    Order,
    Order_boot,
    Family,
    Family_boot,
    Genus,
    Genus_boot,
    species,
    species_boot
  )
# Look at this table
View(taxonomy_rdp)

# Save these objects
save(
  taxonomy,
  taxonomy_rdp,
  file = "data/working/tax_rdp.RData"
)
# Export this table as a .tsv file. I name it with Project Name,
# the reference library used, and taxonomy.

write.table(
  taxonomy_rdp,
  file = paste0("data/results/", project_name, "taxonomy.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

## Assign Taxonomy With BLAST+ =================================================

# We can also assign taxonomy using BLAST. Here we will use the program rBLAST
# to identify our ASVs. rBLAST allows you to connect directly to the NCBI
# server, or use a locally saved reference database (in BLAST format)
# One of the reasons I'm using rBLAST is that it has a command to make a
# BLAST-formatted database from a fasta file.

# We can either use the previous DADA2-formatted Midori reference database, or 
# your own personal database, or use the NCBI genbank database installed in
# Hydra.

# First, we need to make sure BLAST+ is installed on your computer, or is in
# your path if you are using Hydra and it's already installed version.
# Directions for installing BLAST+ locally are coming soon. rBLAST also has
# directions for installing BLAST+ on their GitHub page:
# https://github.com/mhahsler/rBLAST/blob/devel/INSTALL. These instructions
# also include how to ensure the path is set correctly.

# If you are running BLAST+ on Hydra, add the blast+ program binary to the
# path so it can be run with rBLAST. This needs to be run every time you open
# this project, just like loading packages. This path will be automatically
# updated to use the latest version of the database
Sys.setenv(
  PATH = paste(Sys.getenv("PATH"),
  "/share/apps/bioinformatics/blast/2.15.0/bin",
  sep = .Platform$path.sep)
)


### DADA2 formatted Reference Library ------------------------------------------
# If we want to use the DADA2 formatted library from the DAADA2 taxonomic
# assignment step, we first need to make a BLAST database from the fasta file
# we downloaded earlier. We use the makeblastdb function from rBLAST to do this

# The first line is the name of the reference database object - change to the
# database you are using. The second line is the name you want to give your
# BLAST database. Currently, this is set up to create a directory in ref/ with
# the name you chose, and will create all the necessary database files in that
# directory, all using the chosen name also. I use the same directory name and
# file name, you do not have to.
rBLAST::rmakeblastdb(
  reference_fasta,
  db_name = "ref/midori_COI/midori_COI",
  dbtype = "nucl"
)

# Next we create a BLAST database object for this database. This command opens
# the database for running BLAST via predict. Give
# the relative path to the database, and include the name you gave the database
# in the previous step. I.e. db should be the same here as db_name is above
midori_coi_db <- rBLAST::blast(db = "ref/midori_COI/midori_COI", type = "blastn")

# We need to have our representative sequences (the sequences we are going to
# blast). We have to reformat our representative-sequence table to be a named
# vector. Here is our representative sequence table from DADA2.
View(repseq_nochim_md5_asv)

# Make a DNAStringSet object from our representative sequences
sequences_dna <- Biostrings::DNAStringSet(setNames(
  repseq_nochim_md5_asv$ASV,
  repseq_nochim_md5_asv$md5
))
# Look at this object
sequences_dna
# You can also get this from the fasta file we downloaded earlier.
sequences_fasta <- Biostrings::readDNAStringSet("data/results/PROJECTNAME_rep-seq_GENE.fas")

# They make the same thing.
head(sequences_dna)
head(sequences_fasta)

# Finally, we blast our representative sequences against the database we created
tax_blast <- rBLAST::predict(
  midori_coi_db,
  sequences_dna,
  outfmt = "6 qseqid sacc staxids sscinames scomnames qcovs pident",
  BLAST_args = "-perc_identity 85 -max_target_seqs 1 -qcov_hsp_perc 80"
)
View(tax_blast)

# Now we split the taxonomy (sseqid in the table) into a column for each
# taxonomic level
taxonomy_blast <- tax_blast %>%
  tidyr::separate(
    sseqid,
    into = tax_levels,
    sep = ";",
    remove = FALSE
  )
View(taxonomy_blast)

### Hydra NCBI Reference Library -----------------------------------------------
# Hydra has a local copy of the NCBI Genbank database installed. It has the
# nt database, but also a smaller mitochondrial-only database called "mito". We
# will use this mito database for our BLAST searches.

# Create a BLAST database object of the hydra blast mito database, this opens
# it for blasting later via predict
blast_mito_db <- rBLAST::blast(db = "/scratch/dbs/blast/v5/mito")

# We need to have our representative sequences (the sequences we are going to
# blast). We have to reformat our representative-sequence table to be a named
# vector. Here is our representative sequence table from DADA2.
View(repseq_nochim_md5_asv)

# Make a DNAStringSet object from our representative sequences
sequences_dna <- Biostrings::DNAStringSet(setNames(
  repseq_nochim_md5_asv$ANML$ASV,
  repseq_nochim_md5_asv$ANML$md5
))
sequences_dna

# Set the outputs we want from blast. You can add your own.
# qseqid = ASV hast of query sequence
# pident = % identity of query sequence with reference
# saccver = ncbi accession number for reference sequence
# staxids = ncbi taxonomic identifior of reference sequence
fmt <- "qseqid pident saccver staxids"
# Finally, we blast our representative sequences against the database we created
tax_blast <- rBLAST::predict(
  blast_mito_db,
  sequences_dna,
  custom_format = fmt,
  BLAST_args = "-perc_identity 85 -max_target_seqs 1 -qcov_hsp_perc 80"
)
# Make taxid a character instead of a number, for later use (plus, it's not a
# number, it's an identifier)
tax_blast$staxids <- as.character(tax_blast$staxids)
# Let's also change column headings to make more sense
tax_blast <- setNames(tax_blast, c("ASV", "%Identity", "Accession", "taxid"))
View(tax_blast)

# Currently, the Hydra NCBI mito database does not have taxonomic information
# included, so we cannot get classifications directly, only Taxonomic IDs
# (taxIDs). We therefore need to use taxize to get the taxonomic classifications
# from the taxids. This can take some time.

# Use taxize to get classifications from taxid (staxids). Outputs a list of
# taxids, each with a dataframe containing taxid, rank, and name for that rank.
taxonomy <- taxize::classification(tax_blast$taxid, db = "ncbi")

### Make Classification Table --------------------------------------------------
# I want a single table of classifications for each taxid, not a list of
# dataframes. We will then add these classifications to the tax_blast table.
# This is a bit involved, unfortunately

# Keep only unique taxids from taxonomy (remove duplicates).
taxonomy_unique <- taxonomy[!duplicated(names(taxonomy))]
# Make a vector of unique taxids.
taxid_unique <- names(taxonomy_unique)

# Save the classification ranks that we want to keep as a vector. You can save
# whatever ranks you want, but make sure they match the rank names used by NCBI.
tax_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

# Make a tibble with the number of rows equal to the number of unique taxids,
# and the number of columns equal to the number of ranks
taxonomy_table <- tibble::as_tibble(
  matrix("", nrow = length(taxid_unique), ncol = length(ranks)),
  .name_repair = "minimal"
)
# Name the columns from the ranks vector.
colnames(taxonomy_table) <- tax_ranks
# Add a first column to the tibble as the taxids.
taxonomy_table <- tibble::add_column(taxonomy_table, taxid = taxid_unique, .before = 1)

# For each item in the taxonomy_unique list, and for each rank in each list
# item, take the value for that item and rank and save it in the tibble at the
# row for that item and under the column for that rank. If the rank does not
# exist for that item, then use a blank. This fills the blank tibble we just
# made.
for (i in seq_along(taxonomy_unique)) {
  df <- taxonomy_unique[[i]]
  for (r in ranks) {
    taxonomy_table[i,r] <- ifelse(any(df$rank == r),
                                  df$name[df$rank == r],
                                  "")
  }
}

### Add Classification to Taxonomy Table ---------------------------------------
# Left join the taxonomy_table to tax_blast.
tax_blast_class <- dplyr::left_join(
  tax_blast,
  taxonomy_table,
  by = taxid
)
View(tax_blast_class)

# Export this table as a .tsv file. I name it with Project Name,
# the reference library used, and gene.
write.table(
  tax_blast_class,
  file = "data/results/GENE/PROJECT_NCBI_blast_taxonomy_GENE.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


### Combine DADA2 and BLAST Taxonomy Results -----------------------------------
# We can also combine the DADA2 and BLAST taxonomic assignment results into a
# single table to get a direct comparison of the two methods.
taxonomy_rdp_blast <- dplyr::left_join(
  taxonomy_rdp,
  tax_blast_class,
  by = ASV
)


# Export this table as a .tsv file.
write.table(
  taxonomy,
  file = "data/results/GENE/PROJECTNAME_MIDORI_NCBI_BLAST_rdp_taxonomy_GENE.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


# Save all objects to this point.
save.image(file = "data/working/4_Assign_Taxonomy.RData")