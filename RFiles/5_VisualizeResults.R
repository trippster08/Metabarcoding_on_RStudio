library(vegan)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tibble)
library(scales)

# First, lets look at read counts for each sample
read_count <- rowSums(seqtab_nochim_md5)
View(read_count)

# Create a dataframe of the number of reads for each sample
read_count <- tibble::enframe(rowSums(seqtab_nochim_md5))
View(read_count)

read_count <- tibble::enframe(rowSums(seqtab_nochim_md5)) %>%
  rename(
    Sample_ID = name,
    reads = value
  )
View(read_count)
# We can plot this out in a bar plot
read_count_plot <- ggplot2::ggplot(
  read_count,
  aes(x = Sample_ID, y = reads)
) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_discrete(
    labels = read_count$Sample_ID,
    guide = guide_axis(angle = 90)
  )
read_count_plot

# We can also look at the number of ASVs for each sample, first by creating a
# dataframe of the number of ASVs for each sample
asv_count <- tibble::enframe(apply(
  seqtab_nochim_md5,
  1,
  function(row) sum(row != 0)
)) %>%
  rename(
    Sample_ID = name,
    ASVs = value
  )
# and plotting this out in a bar plot
asv_count_plot <- ggplot2::ggplot(
  asv_count,
  aes(x = Sample_ID, y = ASVs)
) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_discrete(labels = asv_count$Sample_ID, guide = guide_axis(angle = 90))
asv_count_plot

# We can also plot number of reads for each sample versus number of ASVs for
# each sample.
# First, we need to join our asv count table with our read count table
read_count_asv_count <- dplyr::left_join(
  read_count,
  asv_count,
  by = join_by(Sample_ID)
)
# Then we can plot these two variables
read_count_asv_count_plot <- ggplot2::ggplot(
  read_count_asv_count,
  aes(x = reads, y = ASVs)
) +
  geom_point() +
  scale_x_continuous(labels = scales::comma)
read_count_asv_count_plot


# We can also rarefy the data to get the expected number of ASVs for each
# sample if they had the same number of reads as the smallest sample.
# Firest, we need to determine which sample has the least number of reads.
raremin <- min(rowSums(seqtab_nochim_md5))
# Then we can use this value and the vegan command rarefy to calculate the
# expected number of ASVs for each sample if all had the same number of reads
# as the smallest sample.
asv_count_rarefied <- tibble::enframe(vegan::rarefy(seqtab_nochim_md5, raremin)) %>%
  rename(
    Sample_ID = name,
    expected_ASVs = value
  )
asv_rarefied_plot <- ggplot2::ggplot(
  asv_count_rarefied,
  aes(x = Sample_ID, y = expected_ASVs)
) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_discrete(
    labels = asv_count_rarefied$Sample_ID,
    guide = guide_axis(angle = 90)
  )
asv_rarefied_plot

# Now lets add this expected_ASV column to our table already containg both the
# reads counts and ASV counts.

read_count_asv_count_expected_asv <- dplyr::left_join(
  read_count_asv_count,
  asv_count_rarefied,
  by = join_by(Sample_ID)
)
# Then we can plot these actual ASV vs expected ASV
asv_count_expected_asv_plot <- ggplot2::ggplot(
  reads_count_asv_count_expected_asv,
  aes(x = reads, y = ASVs)
) +
  geom_point() +
  scale_x_continuous(labels = scales::comma) +
  geom_abline(intercept = 0, slope = 1)

asv_count_expected_asv_plot

# We can look at some basic diversity measures by sample. We first create both
# Simpson  and Shannon-Weaver diversity indices.
simpson <- vegan::diversity(seqtab_nochim_md5, index = "simpson")
shannon <- vegan::diversity(seqtab_nochim_md5, index = "shannon")

# Make a quick histogram of each
hist(simpson)
hist(shannon)


# We can make a dataframe of Shannon-Weaver index measures per sample
shannon_sample <- tibble::enframe(shannon) %>%
  rename(
    Sample_ID = name,
    shannon = value
  )
# And plot them
shannon_plot <- ggplot2::ggplot(
  shannon_sample,
  aes(x = Sample_ID, y = shannon)
) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_discrete(
    labels = shannon_sample$Sample_ID,
    guide = guide_axis(angle = 90)
  )
shannon_plot

# We can do the same for the Simpson index
simpson_sample <- tibble::enframe(simpson) %>%
  rename(
    Sample_ID = name,
    simpson = value
  )
simpson_plot <- ggplot2::ggplot(
  simpson_sample,
  aes(x = Sample_ID, y = simpson)
) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_discrete(
    labels = simpson_sample$Sample_ID,
    guide = guide_axis(angle = 90)
  )
simpson_plot


# We can also look at some basic variables for all these graphs. Say we want to
# look at depth in ASV counts.
# We need metadata, so we have to import your metadata.
meta <- read.delim(
  "dataset1.tsv",
  header = TRUE,
  sep = "\t"
)
# However, you may have some data that looks like a continuous variable that is
# actually a discreete variable (such as filter size or depth). Check to see how
# read.delim interpreted your data
str(meta)
# If you want to change the data type of some columns, you can add the arguement
# "colClasses" to read.delim
meta <- read.delim(
  "dataset1.tsv",
  header = TRUE,
  sep = "\t",
  colClasses = c(depth_ft = "character", arms_num = "character")
)
str(meta)
View(meta)

# Your metadata may have samples that are not on this run, so you can perform
# a left_join to add the metadata only to the samples that you are analyzing.
# Lets add the metadata to our dataframe containing read counts, ASV counts, and
# expected ASVs
read_count_asv_count_expected_asv_meta <- dplyr::left_join(
  read_count_asv_count_expected_asv,
  meta,
  by = join_by(Sample_ID)
)
View(read_count_asv_count_expected_asv_meta)


# Lets go back to our original read_count and asv_count plots, but use the new
# dataframe we just created that has both counts and expected ASV, but lets
# also add a component, coloring by depth.
read_count_plot
read_count_plot <- ggplot2::ggplot(
  read_count_asv_count_expected_asv_meta,
  aes(x = Sample_ID, y = reads, fill = depth_ft)
) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_discrete(
    labels = read_count_asv_count_expected_asv_meta$Sample_ID,
    guide = guide_axis(angle = 90)
  )
read_count_plot

asv_count_plot
asv_count_plot <- ggplot2::ggplot(
  read_count_asv_count_expected_asv_meta,
  aes(x = Sample_ID, y = ASVs, fill = depth_ft)
) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_discrete(
    labels = read_count_asv_count_expected_asv_meta$Sample_ID,
    guide = guide_axis(angle = 90)
  )
asv_count_plot

asv_rarefiled_plot
asv_rarefied_plot <- ggplot2::ggplot(
  read_count_asv_count_expected_asv_meta,
  aes(x = Sample_ID, y = expected_ASVs, fill = depth_ft)
) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_discrete(
    labels = read_count_asv_count_expected_asv_meta$Sample_ID,
    guide = guide_axis(angle = 90)
  )
asv_rarefied_plot

# We can also look at the plot of reads vs ASVs colored by a variable and using
# different shapes for different fractions

read_count_asv_count_meta_plot <- ggplot2::ggplot(
  read_count_asv_count_expected_asv_meta,
  aes(x = reads, y = ASVs)
) +
  geom_point(
    size = 4,
    aes(fill = depth_ft, color = depth_ft, shape = fraction)
  ) +
  scale_x_continuous(labels = scales::comma)
read_count_asv_count_meta_plot

## Rarefaction Curves ==========================================================
# Next we want to look at some rarefaction curves. vegan can give you a
# rarefaction curve on it's own, but your ability to modify it is limited
vegan::rarecurve(seqtab_nochim_md5, step = 1000)

# Instead we use the arguement "tidy = TRUE", and rarecurve gives you the
# rarefied curvews in tabular form, so you can make your own figure.
rarecurve_df <- vegan::rarecurve(
  seqtab_nochim_md5,
  step = 500,
  tidy = TRUE
) %>%
  rename(
    Sample_ID = Site,
    ASVs = Species,
    reads = Sample
  )

rarecurve_df <- vegan::rarecurve(
  seqtab_nochim_md5,
  step = 500,
  tidy = TRUE
) %>%
  dplyr::rename(
    Sample_ID = Site,
    ASVs = Species,
    reads = Sample
  )
head(rarecurve_df)

# Lets add the metadata to this, so we can look at treatment affects
rarecurve_df_meta <- dplyr::left_join(
  rarecurve_df,
  meta,
  join_by(Sample_ID)
)
head(rarecurve_df_meta)

# You can look at the different variables and see how many different values
# there are
unique(rarecurve_df_meta$Sample_ID)
unique(rarecurve_df_meta$fraction)
unique(rarecurve_df_meta$depth_ft)
unique(rarecurve_df_meta$retrieval_year)
# You can also see how many insances of each value there are
table(rarecurve_df_meta$fraction)
table(rarecurve_df_meta$depth_ft)

# Make a line plot of this data.
rarecurve_df_meta_plot <- ggplot2::ggplot(rarecurve_df_meta) +
  geom_line(
    aes(
      x = reads,
      y = ASVs,
      color = fraction,
      linetype = depth_ft,
    ),
    linewidth = 0.75
  ) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash"))
rarecurve_df_meta_plot

# We need to group by Sample_ID, so each sample will form a separate line.
rarecurve_df_meta_plot <- ggplot2::ggplot(rarecurve_df_meta) +
  geom_line(
    aes(
      x = reads,
      y = ASVs,
      color = fraction,
      linetype = depth_ft,
      group = Sample_ID
    ),
    linewidth = 0.75
  ) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash"))
rarecurve_df_meta_plot

# Add an upper limit to the x-axis (reads) to see the expected number of ASVs
# found in each sample with read depth equal to the sample with the least
# number of reads (which should equal the values found in asv.count.rarefied).
rarecurve_df_meta_plot <- ggplot2::ggplot(rarecurve_df_meta) +
  geom_line(
    aes(
      x = reads,
      y = ASVs,
      color = fraction,
      linetype = depth_ft,
      group = Sample_ID
    ),
    linewidth = 0.75
  ) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
  scale_x_continuous(limits = c(0, raremin))
rarecurve_df_meta_plot

## Non-Metric Multidimensional Scaling =========================================
# Now lets move to some ways to look at how similar your samples are usingj
# ordination techniques.
# We are going to use NMDS (Non-Metric Multidimensional Scaling).
# First we need to run a NMDS analysis on our sequence-table.
nmds <- vegan::metaMDS(seqtab_nochim_md5, k = 2, distance = "bray")
# Look at your stress value, it is important in telling you the "goodness of
# fit" of your ordination. You want this value below 0.2.

# You can use vegan to plot this, but like the rarefaction plot, it is limited
# in what you can do, so I convert this to nmds scores that are plotable. A lot
# of the plots in this section was greatly helped by these websites:
# https://jkzorz.github.io/ and
# https://eddatascienceees.github.io/tutorial-rayrr13/
nmds_scores <- tibble::as_tibble(scores(nmds)$sites, rownames = "Sample_ID")
# Add your metadata to this new table
nmds_scores_meta <- dplyr::left_join(
  nmds_scores,
  meta,
  join_by(Sample_ID)
)
head(nmds_scores_meta)

# Plot this table, differentiating two of your variables of interest by color
# and symbol shape.
nmds_scores_meta_plot <- ggplot2::ggplot(nmds_scores_meta, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 4, aes(fill = depth_ft, color = depth_ft, shape = fraction))
nmds_scores_meta_plot

# We can do some very basic statistics to see if these factors are significantly
# different.
# Well use a ANOSIM (Analysis of Similarities) test for two different variables.
# First, looking at depth
anosim_depth <- vegan::anosim(
  x = seqtab_nochim_md5,
  grouping = read_count_asv_count_expected_asv_meta$depth_ft,
  permutations = 9999,
  distance = "bray"
)
anosim_depth
# Then, looking at fraction
anosim_fraction <- vegan::anosim(
  x = seqtab_nochim_md5,
  grouping = read_count_asv_count_expected_asv_meta$fraction,
  permutations = 9999,
  distance = "bray"
)
anosim_fraction
