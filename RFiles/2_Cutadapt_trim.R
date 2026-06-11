# 2 RUN CUTADAPT ###############################################################

## Load Libraries ==============================================================
# Load all R packages you may need, if not coming directly from
# "1_Metabarcoding_R_Pipeline_ComputerPrep".

library(tidyverse)
library(ShortRead)
library(R.utils)

## File Housekeeping ===========================================================

# check to make sure you still have the correct working directory
getwd()
# Load the objects from the previous step, if you are not coming directly from
# that step.
load("data/working/1_RStudioPrep.RData")

# Make a list of all the files in your "data/raw" folder for trimming.
reads_to_trim <- list.files("data/raw", pattern = "\\.fastq\\.gz$")
head(reads_to_trim)
# Separate files by read direction (R1,R2), and save each
reads_to_trim_F <- reads_to_trim[stringr::str_detect(reads_to_trim, "_R1(?=[_.])")]
reads_to_trim_R <- reads_to_trim[stringr::str_detect(reads_to_trim, "_R2(?=[_.])")]

# Look to ensure that there are the same number of F and R reads
length(reads_to_trim_F)
length(reads_to_trim_R)

# Get the names of all the samples in "data/raw". The name is everything
# before the illumina barcode
sample_names_raw <- sapply(
  strsplit(
    basename(
      reads_to_trim[stringr::str_detect(reads_to_trim, "_R1(?=[_.])")]
    ),
    "_S\\d{1,3}_"
  ),
  `[`,
  1
)
head(sample_names_raw)

# Count the number of reads in each sample.
sequence_counts_raw <- sapply(
  paste(
    path_to_raw_reads,
    reads_to_trim[stringr::str_detect(reads_to_trim, "_R1(?=[_.])")],
    sep = "/"
  ),
  function(file) {
    fastq_data <- ShortRead::readFastq(file)
    length(fastq_data)
  }
)

# Name these counts with your sample names
names(sequence_counts_raw) <- sample_names_raw
# Look at the raw read counts
sequence_counts_raw

# Define the path to your primer definition fasta file, if you have more than
# one potential primer to trim. This path will be different for each user.

# At LAB, we use both Nextera and iTru sequencing primers. Currently, our Truseq
# primers include spacers between the sequencing primer and the amplicon primer
# (our Nextera primers do not contain these spacers). Make sure the primer
# definition file includes primers with these spacers attached, otherwise all
# reads will be discared as untrimmed. We have primer definition files for the
# standard COI, 12S MiFish, and 18S_V4 iTru primer (with spacers) and for COI
# and 12S MiFish nextera primers (without spacers). Our primer definition files
# will be downloaded when you download the pipeline. To use a pair of these
# files, just replace "PRIMERF" or "PRIMERR" with the name of the actual primer
# file in the directory "primers" in the paths defined below. If you have
# different primers to remove, open an existing primer file and follow the
# formatting in creating your new primer files.

# If your reads are short, and there is potential for readthrough, you need to
# tell cutadapt to look for primers on the 3' end of each read, as well. These
# primers will be ther reverse complement of the normal primers. They also will
# not be anchored, so the files don't need to include any spacers, and if they
# are not found, the read will still be kept.

# If you know that there will not be any readthrough, you can remove the two
# paths to the RC primers, and two entire lines from the cutadapt arguments:
# following and including "-a" and following and including "-A". Also remove
# "-n 2", because you don't need to run cutadapt twice, since each read will
# only have one primer.

# Create the path to the primer files that will be used by cutadapt.
PrimerF <- file.path(path_to_primers, "active/PrimerF.fas")
PrimerR <- file.path(path_to_primers, "active/PrimerR.fas")
PrimerF_RC <- file.path(path_to_primers, "active/PrimerF_RC.fas")
PrimerR_RC <- file.path(path_to_primers, "active/PrimerR_RC.fas")
# Create empty primer files to be used by cutadapt. These will be filled
# momentarily
file.create(PrimerF, PrimerR, PrimerF_RC, PrimerR_RC)

# Validate the gene names you have entered previously have available primers
invalid_genes <- setdiff(genes, available_primers)
if (length(invalid_genes) > 0) {
  stop(paste("Invalid gene names:", paste(invalid_genes, collapse = ", ")))
}

# Populate primer files to be used by cutadapt. This loop goes through each gene
# finds the primer files (F and R) for that gene, and adds it to the cutadapt
# primer file. We do this only for F and R primers if there is no chance of
# read-through, and additionally for F_RC and R_RC primers if a gene is chosen
# that may have read-through (MiFish, for example).

RC_found <- FALSE
for (gene in genes) {
  cat(
    readLines(file.path(path_to_primers, paste0(gene, "-F.fas"))),
    file = PrimerF,
    sep = "\n",
    append = TRUE
  )
  cat(
    readLines(file.path(path_to_primers, paste0(gene, "-R.fas"))),
    file = PrimerR,
    sep = "\n",
    append = TRUE
  )

  if (gene %in% RC_primers) {
    cat(
      readLines(file.path(path_to_primers, paste0(gene, "-F_RC.fas"))),
      file = PrimerF_RC,
      sep = "\n",
      append = TRUE
    )
    cat(
      readLines(file.path(path_to_primers, paste0(gene, "-R_RC.fas"))),
      file = PrimerR_RC,
      sep = "\n",
      append = TRUE
    )
    RC_found <- TRUE
  }
}


## Run Cutadapt ================================================================
# Save the path to the cutadapt binary.

# Use this path if you are running on the RStudio server.
cutadapt_binary <- "/share/apps/bioinformatics/cutadapt/5.0/bin/cutadapt"

# If you are running locally and installed cutadapt yourself, change the path
# below to your installed cutadapt binary.
cutadapt_binary <- "/home/macdonaldk/.conda/envs/cutadapt/bin/cutadapt"

# The following for loop runs cutadapt on paired samples, one pair at a time.
# It will first determine whether read-through may occur. If not, it will
# only trim 5' primers, and will discard pairs of reads for which

# fmt: skip
if (!RC_found) {
  # Run cutadapt, only removing 5' primers
  cat("\nYour reads should not have potential read-through, so we are only removing primers from the 5' end of each read")
  for (i in seq_along(sample_names_raw)) {
    system2(
      cutadapt_binary,
      args = c(
        "-e 0.2 --discard-untrimmed --minimum-length 30 --cores=8",
        "-g", paste0("file:", PrimerF),
        "-G", paste0("file:", PrimerR),
        "-o", paste0(
          "data/working/trimmed_sequences/{name}/",
          sample_names_raw[i],
          "_trimmed_R1.fastq"
        ), "-p", paste0(
          "data/working/trimmed_sequences/{name}/",
          sample_names_raw[i],
          "_trimmed_R2.fastq"
        ),
        paste0("data/raw/", reads_to_trim_F[i]),
        paste0("data/raw/", reads_to_trim_R[i])
      )
    )
  }
} else {
  # Run cutadapt, first removing 3' primers, then run again removing 5' primers
  cat("\nYour reads have potential read-through, so we will attempt to trim primers from 3' read ends before 5'")
  for (i in seq_along(sample_names_raw)) {
    system2(
      cutadapt_binary,
      args = c(
        "-e 0.2 --cores=8 -O 6",
        "-a", paste0("file:", PrimerR_RC),
        "-A", paste0("file:", PrimerF_RC),
        "-o", paste0(
          "data/raw/fastq/",
          sample_names_raw[i],
          "_trimmed_R1.fastq"
        ), "-p", paste0(
          "data/raw/fastq/",
          sample_names_raw[i],
          "_trimmed_R2.fastq"
        ),
        paste0("data/raw/", reads_to_trim_F[i]),
        paste0("data/raw/", reads_to_trim_R[i])
      )
    )
    system2(
      cutadapt_binary,
      args = c(
        "-e 0.2 --discard-untrimmed --minimum-length 30 --cores=8",
        "-g", paste0("file:", PrimerF),
        "-G", paste0("file:", PrimerR),
        "-o", paste0(
          "data/working/trimmed_sequences/{name}/",
          sample_names_raw[i],
          "_trimmed_R1.fastq"
        ),
        "-p", paste0(
          "data/working/trimmed_sequences/{name}/",
          sample_names_raw[i],
          "_trimmed_R2.fastq"
        ),
        paste0("data/raw/fastq/", sample_names_raw[i], "_trimmed_R1.fastq"),
        paste0("data/raw/fastq/", sample_names_raw[i], "_trimmed_R2.fastq")
      )
    )
  }
  file.remove(file.path("data/raw/fastq", list.files("data/raw/fastq")))
}

# Cutadapt has an issue where it does not correctly compress large read files
# so we have to do it here ourselves. Unfortunately, this takes some time.

# Loop through each gene-specific folder
for (dir in path_to_trimmed) {
  # List all .fastq files in the folder
  fastq_files <- list.files(dir, pattern = "\\.fastq$", full.names = TRUE)

  # Compress each file
  for (file in fastq_files) {
    gzip(file, overwrite = TRUE, compression = 2, remove = TRUE)
  }
}

# Save all objects in case you need to stop here.
save.image(file = "data/working/2_trim.RData")

## Parameter Descriptions ------------------------------------------------------
# We are including our default parameters for cutadapt. You can change these
# parameters if you have prefer others.

# -e 0.2 allows an error rate of 0.2 (20% of primer basepairs can me wrong)

# --minimum-length 30 removes all reads that are not at least 30 bp. However,
# as currently implemented in cutadapt, this does not always work correctly,
# and sometimes it removes the sequence of the reads, but not the name, leaving
# empty reads. We deal with this later in the pipeline.

# -O 3 This is the minimum number of basepairs that the primer must overlap the
# read to be counted and removed. Three is the default. This only affects the
# 3' primer. The 5' primer is anchored, which means it must be found in its
# entirety, or the read is removed

# -n 2 This is the number of times to run cutadapt on each sample. You may need
# to run cutadapt twice if you have read-through and need to remove primers from
# the 3' end as well as the 5' end. Cutadapt will only remove one primer from
# a sample each time it's run.

# --cores=0 tells cutadapt how many cores to use while trimming. 0 sets cutadapt
# to automatically detect the number of cores.

# -g and -G are the paths to the 5' primers with spacers
# -a and -A are the paths to the 3' RC primers. If you are certain there is no
# possibility of read-through, you can omit -a and -A

# -o is the output for trimmed R1 reads
# -p is the output for the trimmed paired R2 reads

# The final line contains the two paired input files
