# 4 DADA2 ######################################################################
# We use Dada2 to filter and trim reads, estimate error rates and use these
# estimates to denoise reads, merge paired reads, and remove chimeric sequences

## Load Libraries ==============================================================
# Load all R packages you may need if not coming directly from the previous
# step.
library(dada2)
library(digest)
library(tidyverse)
library(seqinr)
library(ShortRead)

## File Housekeeping ===========================================================

# Set up your working directory. If you created your new project in the
# directory you want as your working directory (or came directory from the
# previous step in the pipeline), you don't need to do this, and
# skip to the next RStudio command. If you need to set your working directory,
# substitute your own path for the one below.
# check to make sure you still have the correct working directory
getwd()
# Load the objects from the previous step, if you are not coming directly from
# that step.
load("data/working/2_trim.RData")
### Trimmed Reads --------------------------------------------------------------
# Make a list of all gene-specific trimmed reads (it's a list of 3 gene-specific
# vectors containing trimmed read names), with each gene-specific item
# containing a list of Forward(R1) and a list of Reverse(R2) reads.
trimmed_reads <- setNames(
  lapply(genes, function(gene) {
    list(
      F = sort(list.files(
        path_to_trimmed[[gene]],
        pattern = "_R1.fastq.gz",
        full.names = TRUE
      )),
      R = sort(list.files(
        path_to_trimmed[[gene]],
        pattern = "_R2.fastq.gz",
        full.names = TRUE
      ))
    )
  }),
  genes
)

# Make list of gene-specific sample names for trimmed reads, and populate list
# with names from gene-specific trimmed read files.
sample_names_trimmed <- setNames(
  lapply(genes, function(gene) {
    files <- trimmed_reads[[gene]]$F
    if (length(files) == 0 || is.null(files)) {
      return(character(0))
    }
    sapply(
      strsplit(basename(files), "_trimmed"),
      `[`,
      1
    )
  }),
  genes
)

### Remove empty sample files --------------------------------------------------
# This saves the R1 fastq for the sample file only if both the R1 and R2 sample
# files have reads.

# Create a list to store non-empty trimmed reads
actual_trimmed_reads <- setNames(vector("list", length(genes)), genes)
# Create a lis to store trimmed read counts
sequence_counts_trimmed <- setNames(vector("list", length(genes)), genes)

# For each gene, create a list of trimmed forward and trimmed reverse reads.
# These are ephemeral for each gene, and are not saved after this loop.
for (gene in genes) {
  # Count reads in R1 files
  sequence_counts_trimmed[[gene]] <- sapply(
    trimmed_reads[[gene]]$F,
    function(file) {
      fq <- ShortRead::readFastq(file)
      length(fq)
    }
  )

  # Name the read counts using sample_names_trimmed and store them
  names(sequence_counts_trimmed[[gene]]) <- sample_names_trimmed[[gene]]

  # Keep only files with at least 1 read
  valid_indices <- which(sequence_counts_trimmed[[gene]] > 0)

  actual_trimmed_reads[[gene]] <- list(
    F = trimmed_reads[[gene]]$F[valid_indices],
    R = trimmed_reads[[gene]]$R[valid_indices]
  )
}

# Remove read files with zero reads for each gene
for (gene in genes) {
  original_samples <- sample_names_trimmed[[gene]]
  kept_files <- actual_trimmed_reads[[gene]]$F

  # Extract sample names from kept files
  kept_samples <- if (length(kept_files) > 0) {
    sapply(strsplit(basename(kept_files), "_trimmed"), `[`, 1)
  } else {
    character(0)
  }

  # Reset names for trimmed without removed samples
  sample_names_trimmed <- setNames(
    lapply(genes, function(gene) {
      files <- actual_trimmed_reads[[gene]]$F
      if (length(files) == 0 || is.null(files)) {
        return(character(0))
      }

      sapply(
        strsplit(basename(files), "_trimmed"),
        `[`,
        1
      )
    }),
    genes
  )

  # Find filtered-out samples
  removed_samples <- setdiff(original_samples, kept_samples)

  # Print results
  if (length(removed_samples) > 0) {
    cat(
      "\nAfter trimming,",
      length(removed_samples),
      "samples had zero reads for",
      gene,
      "and were removed:\n"
    )
    print(removed_samples)
  } else {
    cat("\nAll samples contained reads for", gene, "\n")
  }
}

# Count the number of samples remaining for each gene, and print
for (gene in genes) {
  nsamples <- length(sample_names_trimmed[[gene]])
  cat("\nWe will analyze", nsamples, "samples for", gene, "\n")
}


## Make Quality Plots ==========================================================

# This creates DADA2 quality plots for each gene. It aggregates all samples and
# creates a single plot. These are the same quality plots as Qiime2, but
# maybe slightly easier to interpret.

# For these plots, the green line is the mean quality score at that position,
# the orange lines are the quartiles (solid for median, dashed for 25% and 75%)
# and the red line represents the proportion of reads existing at that position.

# Create gene-specific list for the quality plots that will be saved and printed
quality_plots_better <- setNames(vector("list", length(genes)), genes)

# Make quality plots
for (gene in genes) {
  # For each gene create a list of direction (F and R) for storing quality plots
  quality_plots <- list(F = NULL, R = NULL)
  quality_plots_better[[gene]] <- list(F = NULL, R = NULL)
  # For each direction for each gene, create quality plots
  for (direction in c("F", "R")) {
    reads <- actual_trimmed_reads[[gene]][[direction]]
    sample_names <- sample_names_trimmed[[gene]]
    quality_plots <- ggplot2::plotQualityProfile(
      reads[1:length(sample_names)],
      aggregate = TRUE
    )
    plot_build <- ggplot2::ggplot_build(quality_plots)
    max_x <- plot_build$layout$panel_params[[1]]$x.range[2]
    # Make the quality plots easier to interpret by changing the x-axis scale,
    # creating vertical lines every 10 bp to better determine quality scores at
    # length, add name of gene to plot
    quality_plots_better[[gene]][[direction]] <- quality_plots +
      scale_x_continuous(
        limits = c(0, max_x),
        breaks = seq(0, max_x, 10)
      ) +
      geom_vline(
        xintercept = seq(0, max_x, 10),
        color = "blue",
        linewidth = 0.25
      ) +
      annotate(
        "text",
        x = 30,
        y = 2,
        label = paste0(gene, " ", direction),
        size = 6
      )

    # Save the plot as a PDF
    ggplot2::gggsave(
      filename = file.path(
        path_to_results[[gene]],
        paste0(
          "plots/",
          project_name,
          "_",
          gene,
          "_qualplot",
          direction,
          ".pdf"
        )
      ),
      plot = quality_plots_better[[gene]][[direction]],
      width = 9,
      height = 9
    )
  }
}

save.image(file = "data/working/3_qual.RData")

## Trim And Filter Reads =======================================================
# We will next use the DADA2 command filterandTrim to remove poor quality reads
# and truncate the 3' end of reads to reduce error inclusion.

# Enter your genes and truncation value for each gene. Follow the pattern shown
# (gene1, R1_truncation_length_gene1, R2_truncation_length_gene1, gene2 etc).
# Make sure to include quotation marks around gene names, but DO NOT put
# truncation numbers in quotation marks (they will be seen as words, not
# numbers).
gene_trunc <- c("gene1", 250, 250, "gene2", 250, 250, "gene3", 220, 220)

# Now we need to make a gene-specific list of truncation values
# Create the empty list that will be populated below
truncation_list <- list()

# Loop through the vector gene_trunc so that each truncation value will be put
# in it's correct place in the list of gene-specific truncation values
for (i in seq(1, length(gene_trunc), by = 3)) {
  gene <- gene_trunc[i]
  trunc_vals <- as.integer(gene_trunc[(i + 1):(i + 2)])
  truncation_list[[gene]] <- trunc_vals
}

# This creates file paths for the reads that will be quality filtered with dada2
# in the next step and a list for storing filtered reads.
filtered_reads <- setNames(vector("list", length(genes)), genes)
# These paths tell DADA2 where to save the filtered reads. These directories
# have not been created yet, but DADA2 will create them for us.
path_to_filtered <- setNames(
  lapply(genes, function(gene) {
    paste0("data/working/filtered_sequences/", gene)
  }),
  genes
)

# For each gene, create a list for forward and reverse read files in
# filtered_reads for DADA2 to populate, then create a list of sample names for
# these read files from trimmed sample names, and add these names to the
# filtered_reads list.
for (gene in genes) {
  sample_names_filtered <- sample_names_trimmed[[gene]]
  filtered_reads[[gene]] <- list(F = NULL, R = NULL)
  for (direction in c("F", "R")) {
    filtered_reads[[gene]][[direction]] <- setNames(
      file.path(
        path_to_filtered[[gene]],
        paste0(
          sample_names_filtered,
          "_filt_",
          direction,
          ".fastq.gz"
        )
      ),
      sample_names_filtered
    )
  }
}

# This filters all reads depending upon the quality (as assigned by the user)
# and trims the ends off the reads for all samples as determined by the quality
# plots. It also saves into "out" the number of reads in each sample that were
# filtered, and how many remain after filtering

# On Windows set multithread=FALSE. If you get errors while running this, change
# to multithread = FALSE, because "error messages and tracking are not handled
# gracefully when using the multithreading functionality".

# "truncLen=c(X,Y)" is how you tell Dada2 where to truncate all forward (X) and
# reverse (Y) reads. Using "0" means reads will not be truncated.
# maxEE sets how many expected errors are allowed before a read is filtered out.

# The amount to truncate is a common question, and very unsettled. I usually
# truncate at the point just shorter than where the red line (proportion of
# reads) in the quality plot reaches 100%.

# This filters all reads depending upon the quality (as assigned by the user)
# and trims the ends off the reads for all samples as determined by the quality
# plots. It also saves into "out" the number of reads in each sample that were
# filtered, and how many remain after filtering

# The amount to truncate is a common question, and very unsettled. I usually
# truncate at the point just shorter than where the red line (proportion of
# reads) in the quality plot reaches 100%.

# Most pipelines use a low maxEE (maximum number of expected errors), but I tend
# to relax this value (from 0,0 to 4,4) because it increases the number of reads
# that are kept, and Dada2 incorporates quality scores in its error models, so
# keeping poorer-quality reads does not adversely effect the results, except in
# very low quality reads. However, increasing maxEE does increase computational
# time.
for (gene in genes) {
  dada2::dfilterAndTrim(
    actual_trimmed_reads[[gene]]$F,
    filtered_reads[[gene]]$F,
    actual_trimmed_reads[[gene]]$R,
    filtered_reads[[gene]]$R,
    truncLen = c(truncation_list[[gene]]),
    maxN = 0,
    maxEE = c(4, 4),
    rm.phix = TRUE,
    truncQ = 2,
    compress = TRUE,
    multithread = TRUE,
    verbose = TRUE
  )
}

# Create a new sample_names_filtered, since some samples no longer have
# reads after filtering, and therefore no longer exist in the directory
for (gene in genes) {
  sample_names_filtered <- setNames(
    lapply(genes, function(gene) {
      files <- list.files(path_to_filtered[[gene]], pattern = "filt_F")
      sapply(strsplit(basename(files), "_filt_F"), `[`, 1)
    }),
    genes
  )
}

# Update filtered_reads, since some samples no longer have
# reads after filtering, and therefore no longer exist in the directory
for (gene in genes) {
  for (direction in c("F", "R")) {
    F_filtered <- list.files(
      path_to_filtered[[gene]],
      pattern = "_filt_F",
      full.names = TRUE
    )
    R_filtered <- list.files(
      path_to_filtered[[gene]],
      pattern = "_filt_R",
      full.names = TRUE
    )

    filtered_reads[[gene]] <- list(
      F = setNames(
        list.files(
          path_to_filtered[[gene]],
          pattern = "filt_F",
          full.names = TRUE
        ),
        sapply(strsplit(basename(F_filtered), "_filt_F"), `[`, 1)
      ),
      R = setNames(
        list.files(
          path_to_filtered[[gene]],
          pattern = "filt_R",
          full.names = TRUE
        ),
        sapply(strsplit(basename(R_filtered), "_filt_R"), `[`, 1)
      )
    )
  }
}

# Create a list to contain the gene-specific filtered read counts
sequence_counts_filtered <- setNames(vector("list", length(genes)), genes)

# For each gene, count the number of forward reads, and add names to the read
# counts
for (gene in genes) {
  sequence_counts_filtered[[gene]] <- sapply(
    filtered_reads[[gene]]$F,
    function(file) {
      fq <- ShortRead::readFastq(file)
      length(fq)
    }
  )
  sample_names <- sapply(
    strsplit(basename(filtered_reads[[gene]]$F), "_filt_F"),
    `[`,
    1
  )
  names(sequence_counts_filtered[[gene]]) <- sample_names
}

# Report removed samples
for (gene in genes) {
  removed <- setdiff(
    sample_names_trimmed[[gene]],
    sample_names_filtered[[gene]]
  )
  cat("\nHere are the samples removed after filtering for", gene, "\n")
  print(removed)
}

# Save all the objects created to this point in this section
save.image(file = "data/working/4_filter.RData")

## Estimate Error Rates ========================================================

# Here we use a portion of the data to determine error rates. These error rates
# will be used in the next (denoising) step to narrow down the sequences to a
# reduced and corrected set of unique sequences

# Create a list to contain the gene-specific error rates
errors <- setNames(vector("list", length(genes)), genes)

# Make a loop to determine errors for each gene. First add read direction to
# each gene in the list, then model errors
for (gene in genes) {
  cat("\nModelling error rates and creating plots for", gene, "\n")
  errors[[gene]] <- list(F = NULL, R = NULL)
  for (direction in c("F", "R")) {
    errors[[gene]][[direction]] <- dada2::learnErrors(
      filtered_reads[[gene]][[direction]],
      nbases = 1e+08,
      errorEstimationFunction = loessErrfun,
      multithread = TRUE,
      randomize = FALSE,
      MAX_CONSIST = 10,
      OMEGA_C = 0,
      qualityType = "Auto",
      verbose = FALSE
    )
  }
}

# We can visualize the estimated error rates to make sure they don't look too
# crazy. The red lines are error rates expected under the "...nominal definition
# of the Q-score." The black dots are "...observed error rates for each
# consensus quality score." The black line shows the "...estimated error rates
# after convergence of the machine-learning algorithm." I think the main things
# to look at here are to make sure that each black line is a good fit to the
# observed error rates, and that estimated error rates decrease with increased
# quality. Note: They almost never match. What do we do about that?...? I don't
# know.

# First we make a list to store the gene-specific error plots
error_plots <- setNames(vector("list", length(genes)), genes)
# Make a loop to create error plots for each gene. First add read direction to
# each gene in the list, then make plots, then export plots as a pdf
for (gene in genes) {
  error_plots[[gene]] <- list(F = NULL, R = NULL)

  for (direction in c("F", "R")) {
    err <- errors[[gene]][[direction]]
    error_plots[[gene]][[direction]] <- dada2::plotErrors(err, nominalQ = TRUE)
    ggplot2::ggsave(
      filename = file.path(
        path_to_results[[gene]],
        paste0(
          "plots/",
          project_name,
          "_",
          gene,
          "_errorplot_",
          direction,
          ".pdf"
        )
      ),
      plot = error_plots[[gene]][[direction]],
      width = 9,
      height = 9
    )
  }
}

# Save all the objects created to this point in this section
save.image(file = "data/working/5_error.RData")

## Denoise =====================================================================
# This applies the "core sample inference algorithm" (i.e. denoising) in dada2
# to get corrected unique sequences. The two main inputs are the first, which is
# the filtered sequences (filtered_F), and "err =" which is the error file from
# learnErrors (effF).

# First make a list to store the gene-specific denoised sequences
denoised <- setNames(vector("list", length(genes)), genes)
# Loop through all genes, denoising reads, after adding read direction to the
# list.
for (gene in genes) {
  cat("\nDenoising reads for", gene, "\n")
  denoised[[gene]] <- list(F = NULL, R = NULL)
  for (direction in c("F", "R")) {
    denoised[[gene]][[direction]] <- dada2::dada(
      filtered_reads[[gene]][[direction]],
      err = errors[[gene]][[direction]],
      errorEstimationFunction = loessErrfun,
      selfConsist = FALSE,
      pool = FALSE,
      multithread = TRUE,
      verbose = TRUE
    )
  }
  cat("\nDenoising is complete for", gene, "\n")
}

# Save all the objects created to this point in this section
save.image(file = "data/working/6_denoise.RData")

## Merge Paired Sequences ======================================================

# Here we merge the paired reads. mergePairs calls for the forward denoising
# result (denoised), then the forward filtered and truncated reads
# (filtered_reads), then the same for the reverse reads.

# You can change the minimum overlap (minOverlap), and the number of mismatches
# that are allowed in the overlap region (maxMismatch). Default values are
# shown.

# mergePairs results in a data.frame from each sample that contains a row for
# "...each unique pairing of forward/reverse denoised sequences." The data.frame
# also contains multiple columns describing data for each unique merged
# sequence.

# First, make a list to hold the gene-specific merged reads.
merged_reads <- setNames(vector("list", length(genes)), genes)
# Loop through each gene, creating gene-specific merged sequences
for (gene in genes) {
  cat("\nMerging forward and reverse reads for", gene, "\n")
  merged_reads[[gene]] <- dada2::mergePairs(
    denoised[[gene]]$F,
    filtered_reads[[gene]]$F,
    denoised[[gene]]$R,
    filtered_reads[[gene]]$R,
    minOverlap = 12,
    maxMismatch = 0,
    verbose = TRUE
  )
  cat("\nMerging is complete for", gene)
}
# Save all the objects created to this point in this section
save.image(file = "data/working/7_merge.RData")

## Create Sequence-Table =======================================================

# Now we make a sequence-table containing columns for unique sequence (ASV),
# rows for each sample, and numbers of reads for each ASV/sample combination as
# the body of the table. This table is in the reverse orientation of the feature
# table output by Qiime2 (which uses samples as columns), but we easily
# transpose this table if needed later (and we will later). I think for now, I
# will use "sequence-table" for the table with columns of sequences, and
# "feature-table" for tables with columns of samples.

# First make a list to hold the gene-specific sequence-table
seqtab <- setNames(vector("list", length(genes)), genes)

# Loop through each gene, making gene-specific sequence-tables
for (gene in genes) {
  seqtab[[gene]] <- dada2::makeSequenceTable(merged_reads[[gene]])

  # This describes the dimensions of the gene-specific table just made
  # First the number of samples
  cat(
    "\nThis is the number of samples for your",
    gene,
    "Sequence-Table:",
    length(rownames(seqtab[[gene]]))
  )
  # Then the number of ASVs
  cat(
    "\nThis is the number of ASVs for your",
    gene,
    "Sequence-Table:",
    length(colnames(seqtab[[gene]]))
  )
}

## Remove Chimeric Sequences ===================================================

# Here we remove chimera sequences.

# Make lists that will be populated later by gene-specific data
seqtab_nochim <- setNames(vector("list", length(genes)), genes)
repseq_chimera_md5 <- setNames(vector("list", length(genes)), genes)

# Loop through each gene, remove chimeras from the gene-specific sequence-table
# and make a fasta of removed chimeric sequences
for (gene in genes) {
  cat(
    "\nRemoving chimeric sequences and creating new sequence table for",
    gene,
    "\n"
  )
  # remove chimeric sequences from seqtab and save in seqtab_nochim
  seqtab_nochim[[gene]] <- dada2::removeBimeraDenovo(
    seqtab[[gene]],
    method = "consensus",
    multithread = TRUE,
    verbose = TRUE
  )

  # Print to the log the number of ASVs remaining in the table.
  cat(
    "\nThis is the number of ASVs for your chimera-free",
    gene,
    "Sequence-Table:",
    length(colnames(seqtab_nochim[[gene]])),
    "\n"
  )

  # Make a list of the ASVs that are considered chimeras.
  chimeras_list <- dada2::isBimeraDenovoTable(
    seqtab[[gene]],
    multithread = TRUE,
    verbose = TRUE
  )

  # This makes a new vector containing all the ASV's (unique sequences) returned
  # by dada2.
  repseq_all <- dada2::getSequences(seqtab[[gene]])
  # Get a list of just chimera ASVs by filtering all sequences by chimera_list
  repseq_chimera <- repseq_all[chimeras_list]
  # Make and add md5 hash to the repseq_chimera
  repseq_chimera_md5[[gene]] <- c()
  for (i in seq_along(repseq_chimera)) {
    repseq_chimera_md5[[gene]][i] <- digest::digest(
      repseq_chimera[i],
      serialize = FALSE,
      algo = "md5"
    )
  }

  # Export chimeric sequences as fastas
  seqinr::write.fasta(
    sequences = as.list(repseq_chimera),
    names = repseq_chimera_md5[[gene]],
    open = "w",
    as.string = FALSE,
    file.out = file.path(
      path_to_results[[gene]],
      paste0(
        "additional_results/",
        project_name,
        "_",
        gene,
        "_rep-seq_chimeras.fas"
      )
    )
  )
}

# Save all the objects created to this point in this section
save.image(file = "data/working/8_chimera.RData")

## Examine Sequence Lengths ====================================================

# This shows the length of the representative sequences (ASV's). Typically,
# there are a lot of much longer and much shorter sequences.
# Create empty list for length tables
seq_length_table <- setNames(vector("list", length(genes)), genes)

# Loop through genes, creating a sequence length table for each
for (gene in genes) {
  # Count the number of bp for each gene-specific ASV from the sequence-tables
  seq_length_table[[gene]] <- table(nchar(dada2::getSequences(seqtab_nochim[[gene]])))

  # Export this table as a .tsv
  write.table(
    seq_length_table[[gene]],
    file = file.path(
      path_to_results[[gene]],
      paste0(
        "additional_results/",
        project_name,
        "_ASV_lengths_table_",
        gene,
        ".tsv"
      )
    ),
    quote = FALSE,
    sep = "\t",
    row.names = TRUE,
    col.names = NA
  )
}

## Track Reads Through Dada2 Process ===========================================

# Here, we look at how many reads made it through each step. This is similar to
# the stats table that we look at in Qiime2 or "track" from the DADA2 tutorial.
# This is a good quick way to see if something is wrong (i.e. only a small
# proportion make it through).

# Make lists that will be populated later by gene-specific data
sequence_counts_postfiltered <- setNames(vector("list", length(genes)), genes)
track_reads <- setNames(vector("list", length(genes)), genes)

# Make a table for the post-filtered samples, including denoised,
# merged, and non-chimeric read counts
for (gene in genes) {
  getN <- function(x) sum(dada2::getUniques(x))
  sequence_counts_postfiltered[[gene]] <- tibble::tibble(
    Sample_ID = sample_names_filtered[[gene]],
    Denoised_Reads_F = sapply(denoised[[gene]]$F, getN),
    Denoised_Reads_R = sapply(denoised[[gene]]$R, getN),
    Merged_Reads = sapply(merged_reads[[gene]], getN),
    Non_Chimeras = as.integer(rowSums(seqtab_nochim[[gene]]))
  )

  # Then we are going to add the post-filtered read count data to the three
  # count data objects we already have (raw, trimmed, filtered).
  track_reads[[gene]] <- tibble::tibble(
    Sample_ID = names(sequence_counts_raw),
    Raw_Reads = as.numeric(sequence_counts_raw),
  ) %>%
    dplyr::left_join(
      tibble::tibble(
        sequence_counts_trimmed[[gene]],
        Sample_ID = names(sequence_counts_trimmed[[gene]]),
        Trimmed_Reads = as.numeric(sequence_counts_trimmed[[gene]])
      ),
      join_by(Sample_ID)
    ) %>%
    dplyr::left_join(
      tibble::tibble(
        sequence_counts_filtered[[gene]],
        Sample_ID = names(sequence_counts_filtered[[gene]]),
        Filtered_Reads = as.numeric(sequence_counts_filtered[[gene]])
      ),
      join_by(Sample_ID)
    ) %>%
    dplyr::left_join(
      sequence_counts_postfiltered[[gene]],
      join_by(Sample_ID)
    ) %>%
    dplyr::mutate(Proportion_Trimmed_Passed = Non_Chimeras / Trimmed_Reads) %>%
    dplyr::mutate(Proportion_Raw_Passed = Trimmed_Reads / Raw_Reads) %>%
    dplyr::select(
      Sample_ID,
      Raw_Reads,
      Trimmed_Reads,
      Filtered_Reads,
      Denoised_Reads_F,
      Denoised_Reads_R,
      Merged_Reads,
      Non_Chimeras,
      Proportion_Trimmed_Kept,
      Proportion_Gene
    )
}

# Export this table as a .tsv
for (gene in genes) {
  write.table(
    track_reads[[gene]],
    file = file.path(
      path_to_results[[gene]],
      paste0(
        "additional_results/",
        project_name,
        "_track_reads_",
        gene,
        ".tsv"
      )
    ),
    quote = FALSE,
    sep = "\t",
    row.names = TRUE,
    col.names = NA
  )
}

## Export Sequence-Table and Rep-Seq-Table =====================================
# This exports a sequence-table: columns of ASV's, rows of samples, and
# values = number of reads.

# If you have mulitple Miseqruns for the same project that will need to be
# combined for further analyses, you may want to name this file
# "PROJECTNAME_MISEQRUN_sequence-table.tsv" to differentiate different runs.
# In "5 Metabarcoding_R_Pipeline_RStudio_ImportCombine" we'll show how to
# combine data from separate runs for analyses.

# NOTE!!!
# The Sequence-Table we have now (seqtab_nochim) is very unwieldy, since each
# column name is an entire ASV. Instead, we will convert ASVs using the md5
# encryption model, creating a 32bit representative "hash" of each ASV. Every
# hash is essentially unique to the ASV it is representing. We would then
# replace the ASVs in the column headings with their representative md5 hash.
# However, having an ASV hash as a column heading requires the creation of a
# Representative Sequence list, which tells us which hash represents which ASV.

### Create And Use md5 Hash ----------------------------------------------------
# To create a Sequence list with md5 hash instead of ASVs, we first need to
# create a list of md5 hash's of all ASV's.

# This makes a new vector containing all the ASV's (unique sequences) returned
# by dada2. We are going to use this list to create md5 hashes.

repseq_nochim <- setNames(vector("list", length(genes)), genes)
repseq_nochim_md5 <- setNames(vector("list", length(genes)), genes)
seqtab_nochim_md5 <- setNames(vector("list", length(genes)), genes)
repseq_nochim_md5_asv <- setNames(vector("list", length(genes)), genes)
feattab_nochim_md5 <- setNames(vector("list", length(genes)), genes)
# Use the program digest (in a For Loop) to create a new vector containing the
# unique md5 hashes of the representative sequences (ASV's). This results in
# identical feature names to those assigned in Qiime2.
for (gene in genes) {
  repseq_nochim[[gene]] <- dada2::getSequences(seqtab_nochim[[gene]])

  repseq_nochim_md5[[gene]] <- c()
  for (i in seq_along(repseq_nochim[[gene]])) {
    repseq_nochim_md5[[gene]][i] <- digest::digest(
      repseq_nochim[[gene]][i],
      serialize = FALSE,
      algo = "md5"
    )
  }

  # Add md5 hash to the sequence-table from the DADA2 analysis.

  seqtab_nochim_md5[[gene]] <- seqtab_nochim[[gene]]
  colnames(seqtab_nochim_md5[[gene]]) <- repseq_nochim_md5[[gene]]

  # Export this sequence table with column headings as md5 hashs instead of ASV
  # sequences
  write.table(
    seqtab_nochim_md5[[gene]],
    file = file.path(
      path_to_results[[gene]],
      paste0(
        "/",
        project_name,
        "_sequence-table-md5_",
        gene,
        ".tsv"
      )
    ),
    quote = FALSE,
    sep = "\t",
    row.names = TRUE,
    col.names = NA
  )

  # Create an md5/ASV table, with each row as an ASV and it's representative md5
  # hash.
  repseq_nochim_md5_asv[[gene]] <- tibble::tibble(
    md5 = repseq_nochim_md5[[gene]],
    ASV = repseq_nochim[[gene]]
  )
  # This exports all the ASVs and their respective md5 hashes as a two-column
  # table.
  write.table(
    repseq_nochim_md5_asv[[gene]],
    file = file.path(
      path_to_results[[gene]],
      paste0(
        "/",
        project_name,
        "_representative-seq-md5-table_",
        gene,
        ".tsv"
      )
    ),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
}

## Create and Export Feature-Table ===========================================
# This creates and exports a feature-table: row of ASV's (shown as a md5 hash
# instead of sequence), columns of samples, and values = number of reads. With
# this table you will also need a file that relates each ASV to it's
# representative md5 hash. We download this in the next section.

for (gene in genes) {
  # Transpose the sequence-table, and convert the result into a tibble.
  feattab_nochim_md5[[gene]] <- t(seqtab_nochim_md5[[gene]])
  # Export as tsv
  write.table(
    feattab_nochim_md5[[gene]],
    file = file.path(
      path_to_results[[gene]],
      paste0(
        "/",
        project_name,
        "_feature-table_md5_",
        gene,
        ".tsv"
      )
    ),
    quote = FALSE,
    sep = "\t",
    row.names = TRUE,
    col.names = NA
  )
}

## Export Representative Sequences fasta =======================================
# Here we export our our representative sequences, either as a fasta (with the
# md5 hash as the ASV name), or as a table with ASV and md5 hash as columns.

for (gene in genes) {
  # This exports all the ASVs in fasta format, with ASV hash as the sequence
  # name. This is analogous to the representative sequence output in Qiime2.
  write.fasta(
    sequences = as.list(repseq_nochim_md5_asv[[gene]]$ASV),
    names = repseq_nochim_md5_asv[[gene]]$md5,
    open = "w",
    as.string = FALSE,
    file.out = file.path(
      path_to_results[[gene]],
      paste0(
        "/",
        project_name,
        "_rep-seq_",
        gene,
        ".fas"
      )
    )
  )
}

## Create feature-to-fasta =====================================================
# This creates a fasta file containing all the ASV's for each sample. Each ASV
# will be labeled with the sample name, ASV hash, and number of reads of that
# ASV in that sample. This was derived from a python script from Matt Kweskin
# called featuretofasta.py (hence the name).

### Create Sequence-List
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
# The sequence-table has a column with sample names, and N columns of ASV's
# containing count values. We want all the count data to be in a single column,
# so we use a tidyr command called "pivot_longer" to make the table "tall",
# which means the table goes from 17 x 2811 to 47770 x 3 for example
# (47770 = 2810 x 17. 2810 instead of 2811 because the first column of the
# original table contains sample names, not counts). This makes the table tidier
# (meaning that each column is now a true variable).
seqtab_nochim_tall <- setNames(vector("list", length(genes)), genes)
repseq_tall <- setNames(vector("list", length(genes)), genes)
repseq_tall_md5 <- setNames(vector("list", length(genes)), genes)
seqtab_nochim_tall_md5 <- setNames(vector("list", length(genes)), genes)

for (gene in genes) {
  seqtab_nochim_tall[[gene]] <- tibble::as_tibble(
    seqtab_nochim[[gene]],
    rownames = "sample"
  ) %>%
    tidyr::pivot_longer(
      !sample,
      names_to = "ASV",
      values_to = "count"
    ) %>%
    subset(count != 0)

  # Save the ASV sequences from the sequence-list table
  # (seqtab_nochim_tall_nozera) as a new list.
  repseq_tall[[gene]] <- seqtab_nochim_tall[[gene]]$ASV

  # Convert the sequences into md5 hashs, as we did earlier. md5 hashs are
  # consistent across jobs, meaning identical sequences from different projects
  # or being converted by different programs will result in the same hash (i.e.
  # hashs here will match hashs above)
  repseq_tall_md5[[gene]] <- c()
  for (i in seq_along(repseq_tall[[gene]])) {
    repseq_tall_md5[[gene]][i] <- digest::digest(
      repseq_tall[[gene]][i],
      serialize = FALSE,
      algo = "md5"
    )
  }

  # Attach the ASV hashes as a column (called "md5") to the tall table. The
  # table should now have 4 columns, and each row of the "md5" column should
  # be a md5 hash of its respective ASV.
  # Also create a new column in this table that contains "sample", "feature",
  # and "count", concatenated. This is the heading for each sequence in the
  # fasta file created by Matt Kweskin's script "featuretofasta.py"
  seqtab_nochim_tall_md5[[gene]] <- seqtab_nochim_tall[[gene]] %>%
    dplyr::mutate(md5 = repseq_tall_md5[[gene]]) %>%
    dplyr::mutate(sample_md5_count = paste(sample, md5, count, sep = "_")) %>%
    dplyr::select(sample, md5, count, sample_md5_count, ASV)
  ### Export feature-to-fasta fas file for each gene ---------------------------
  # Create a fasta-formatted file of each row sequence (i.e. ASV), with a
  # heading of "sample_feature_count".
  seqinr::write.fasta(
    sequences = as.list(seqtab_nochim_tall_md5[[gene]]$ASV),
    names = seqtab_nochim_tall_md5[[gene]]$sample_md5_count,
    open = "w",
    as.string = FALSE,
    file.out = file.path(
      path_to_results[[gene]],
      paste0(
        "additional_results/",
        project_name,
        "_feature-to-fasta_",
        gene,
        ".fas"
      )
    )
  )
}
save.image(file = "data/working/9_output.RData")
