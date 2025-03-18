# SET UP #################################################


# Load/install packages

if (!require("pacman")) # pacman = package manager
  install.packages("pacman")

pacman::p_load(
  cowplot,
  dplyr,
  GGally,
  ggplot2,
  ggthemes,
  ggvis,
  gridExtra,
  httr,
  lubridate,
  patchwork,
  pheatmap,
  plotly,
  rio,
  rmarkdown,
  shiny,
  stringr,
  tidyr,
  tidyverse
)


# Import data - ** MAKE SURE TO CHANGE THE FILE DIRECTORY PATHS :) **

tumIndelData <- import("~/Downloads/pnet-dataviz/MIP data - tumor indels.csv") %>%
  rename_with( ~ gsub("[ -]", "_", .x)) # Remove annoying spaces/hyphen in "Depth of variant-supporting bases"

tumCovData <- import("~/Downloads/pnet-dataviz/MIP data - tumor coverage.csv")


# PROCESS DATA #################################################


# Reshape tumor coverage data to long format for easier merging

tumCovData_byMIPname <- tumCovData %>%
  pivot_longer(
    cols = -c(Chr, Start, End, MIPname),
    names_to = "SampleName",
    values_to = "Coverage"
  ) %>%
  group_by(MIPname, SampleName) %>%
  summarise(TotalCoverage = sum(Coverage, na.rm = TRUE),
            .groups = "drop")


# Merge tumor indel and tumor coverage tables based on MIPname and SampleName
# Compute square root variant frequency

nonAdj_mergedData <- tumIndelData %>%
  left_join(tumCovData_byMIPname,
            by = c("In50bp" = "MIPname", "SampleName" = "SampleName")) %>%
  mutate(SqrtVarFreq = sqrt((Depth_of_variant_supporting_bases / TotalCoverage) * 100))


# ADJUST MERGED DATA TABLE FOR HETEROZYGOUS MUTATIONS:

# Identify variants with a frequency ~1/2 of the highest square root variant frequency for
# SampleName-In50bp pairs with multiple variants, within +/- the leniency factor
# times the highest square root variant frequency

LENIENCY_FACTOR <- 0.15 # ** THIS CAN BE ADJUSTED **

hetAdj_mergedData <- nonAdj_mergedData %>%
  group_by(SampleName, In50bp) %>%
  mutate(
    max_SqrtVarFreq = max(SqrtVarFreq, na.rm = TRUE),
    lower_bound = (max_SqrtVarFreq / 2) - (LENIENCY_FACTOR * max_SqrtVarFreq),
    # CHANGE
    upper_bound = (max_SqrtVarFreq / 2) + (LENIENCY_FACTOR * max_SqrtVarFreq)
  ) %>%
  mutate(
    Depth_of_variant_supporting_bases = ifelse(
      SqrtVarFreq >= lower_bound & SqrtVarFreq <= upper_bound,
      0,
      Depth_of_variant_supporting_bases
    )
  ) %>%
  ungroup()


# ADJUST MERGED DATA TABLE FOR NON-FRAMESHIFT MUTATIONS:

# Identify non-frameshift mutations by checking if the REF and ALT sequences
# differ by a factor of 3 base pairs, set their depth of variant supporting bases to 0

nonFSAdj_mergedData <- nonAdj_mergedData %>%
  mutate(
    Depth_of_variant_supporting_bases = ifelse(
      abs(nchar(ALT) - nchar(REF)) %% 3 == 0,
      0,
      Depth_of_variant_supporting_bases
    )
  )


# Extract gene names from MIP names

nonAdj_mergedData <- nonAdj_mergedData %>%
  mutate(GeneName = sub("_.*", "", In50bp))  

hetAdj_mergedData <- hetAdj_mergedData %>%
  mutate(GeneName = sub("_.*", "", In50bp))  

nonFSAdj_mergedData <- nonFSAdj_mergedData %>%
  mutate(GeneName = sub("_.*", "", In50bp))

tumCovData_byMIPname <- tumCovData_byMIPname %>%
  mutate(GeneName = sub("_.*", "", MIPname))


# Compute total depth of variant supporting bases

nonAdj_totVarDepth <- nonAdj_mergedData %>%
  group_by(SampleName, GeneName) %>%
  summarise(
    TOTAL_Depth_of_variant_supporting_bases = sum(Depth_of_variant_supporting_bases, na.rm = TRUE)
  ) %>%
  ungroup()

hetAdj_totVarDepth <- hetAdj_mergedData %>%
  group_by(SampleName, GeneName) %>%
  summarise(
    TOTAL_Depth_of_variant_supporting_bases = sum(Depth_of_variant_supporting_bases, na.rm = TRUE)
  ) %>%
  ungroup()

nonFSAdj_totVarDepth <- nonFSAdj_mergedData %>%
  group_by(SampleName, GeneName) %>%
  summarise(
    TOTAL_Depth_of_variant_supporting_bases = sum(Depth_of_variant_supporting_bases, na.rm = TRUE)
  ) %>%
  ungroup()


# Compute grand total coverage for each gene

totCov <- tumCovData_byMIPname %>%
  group_by(SampleName, GeneName) %>%
  summarise(GRAND_TotalCoverage = sum(TotalCoverage, na.rm = TRUE)) %>%
  ungroup()


# Merge (adjusted) tumor indel and tumor coverage tables by gene name

nonAdj_finalData <- nonAdj_totVarDepth %>%
  left_join(totCov, by = c("SampleName", "GeneName")) %>%
  mutate(TOTAL_nonAdj_SqrtVarFreq = sqrt((
    TOTAL_Depth_of_variant_supporting_bases / GRAND_TotalCoverage
  ) * 100
  ))


hetAdj_finalData <- hetAdj_totVarDepth %>%
  left_join(totCov, by = c("SampleName", "GeneName")) %>%
  mutate(TOTAL_hetAdj_SqrtVarFreq = sqrt((
    TOTAL_Depth_of_variant_supporting_bases / GRAND_TotalCoverage
  ) * 100
  ))


nonFSAdj_finalData <- nonFSAdj_totVarDepth %>%
  left_join(totCov, by = c("SampleName", "GeneName")) %>%
  mutate(TOTAL_nonFSAdj_SqrtVarFreq = sqrt((
    TOTAL_Depth_of_variant_supporting_bases / GRAND_TotalCoverage
  ) * 100
  ))


# Pivot final data tables into heat map-friendly matrices

nonAdj_heatmapMatrix <- nonAdj_finalData %>%
  select(GeneName, SampleName, TOTAL_nonAdj_SqrtVarFreq) %>%
  pivot_wider(names_from = SampleName,
              values_from = TOTAL_nonAdj_SqrtVarFreq,
              values_fill = 0) %>%
  column_to_rownames(var = "GeneName") %>%
  as.matrix()

hetAdj_heatmapMatrix <- hetAdj_finalData %>%
  select(GeneName, SampleName, TOTAL_hetAdj_SqrtVarFreq) %>%
  pivot_wider(names_from = SampleName,
              values_from = TOTAL_hetAdj_SqrtVarFreq,
              values_fill = 0) %>%
  column_to_rownames(var = "GeneName") %>%
  as.matrix()

nonFSAdj_heatmapMatrix <- nonFSAdj_finalData %>%
  select(GeneName, SampleName, TOTAL_nonFSAdj_SqrtVarFreq) %>%
  pivot_wider(names_from = SampleName,
              values_from = TOTAL_nonFSAdj_SqrtVarFreq,
              values_fill = 0) %>%
  column_to_rownames(var = "GeneName") %>%
  as.matrix()

# Compute column-wise (sample) and row-wise (gene) mean variant frequencies

nonAdj_colMeans <- colMeans(nonAdj_heatmapMatrix)
nonAdj_rowMeans <- rowMeans(nonAdj_heatmapMatrix)

hetAdj_colMeans <- colMeans(hetAdj_heatmapMatrix)
hetAdj_rowMeans <- rowMeans(hetAdj_heatmapMatrix)

nonFSAdj_colMeans <- colMeans(nonFSAdj_heatmapMatrix)
nonFSAdj_rowMeans <- rowMeans(nonFSAdj_heatmapMatrix)


# Order row-wise mean variant frequencies in increasing/decreasing order (output = integer vector)
# Re-order genes accordingly in a new character vector

nonAdj_rowOrder_inc <- order(nonAdj_rowMeans, decreasing = FALSE)
nonAdj_rowOrder_dec <- order(nonAdj_rowMeans, decreasing = TRUE)
nonAdj_ordGeneNames <- rownames(nonAdj_heatmapMatrix)[nonAdj_rowOrder_inc]

hetAdj_rowOrder_inc <- order(hetAdj_rowMeans, decreasing = FALSE)
hetAdj_rowOrder_dec <- order(hetAdj_rowMeans, decreasing = TRUE)
hetAdj_ordGeneNames <- rownames(hetAdj_heatmapMatrix)[hetAdj_rowOrder_inc]

nonFSAdj_rowOrder_inc <- order(nonFSAdj_rowMeans, decreasing = FALSE)
nonFSAdj_rowOrder_dec <- order(nonFSAdj_rowMeans, decreasing = TRUE)
nonFSAdj_ordGeneNames <- rownames(nonFSAdj_heatmapMatrix)[nonFSAdj_rowOrder_inc]


# Convert gene vs. square root variant frequency data from heat map matrix to a data frame format,
# which will be used to create the mean square root variant frequency bar chart.

nonAdj_rowMeans_df <- data.frame(GeneName = rownames(nonAdj_heatmapMatrix),
                                 MeanSqrtFreq = nonAdj_rowMeans)

hetAdj_rowMeans_df <- data.frame(GeneName = rownames(hetAdj_heatmapMatrix),
                                 MeanSqrtFreq = hetAdj_rowMeans)

nonFSAdj_rowMeans_df <- data.frame(GeneName = rownames(nonFSAdj_heatmapMatrix),
                                 MeanSqrtFreq = nonFSAdj_rowMeans)

# Re-order genes in the data frame

nonAdj_rowMeans_df$GeneName <- factor(nonAdj_rowMeans_df$GeneName, levels = nonAdj_ordGeneNames)

hetAdj_rowMeans_df$GeneName <- factor(hetAdj_rowMeans_df$GeneName, levels = hetAdj_ordGeneNames)

nonFSAdj_rowMeans_df$GeneName <- factor(nonFSAdj_rowMeans_df$GeneName, levels = nonFSAdj_ordGeneNames)


# Re-order heat map matrix rows (genes) by mean square root variant frequency, decreasing top to bottom

ord_nonAdj_heatmapMatrix <- nonAdj_heatmapMatrix[nonAdj_rowOrder_dec, ]

ord_hetAdj_heatmapMatrix <- hetAdj_heatmapMatrix[hetAdj_rowOrder_dec, ]

ord_nonFSAdj_heatmapMatrix <- nonFSAdj_heatmapMatrix[nonFSAdj_rowOrder_dec, ]



# CREATE HEATMAPS #################################################

# The heatmaps below contain a title and a column-wise dendogram.
# The rows (genes) are ordered by mean square root variant frequency, decreasing top to bottom.
# The columns (samples) are ordered by hierarchical clustering as generated by pheatmap.

# They are used to obtain the column-wise (sample) hierachical clustering as generated by pheatmap,
# which will be used to generate the correspondingly ordered mean square root variant frequency bar chart.

nonAdj_heatmap <- pheatmap(
  ord_nonAdj_heatmapMatrix,
  cluster_rows = FALSE, # ** HIDE row-wise (gene) hierachical clustering **
  cluster_cols = TRUE, # ** SHOW column-wise (sample) hierarchical clustering **
  color = colorRampPalette(
    c(
      "steelblue",
      "green",
      "lawngreen",
      "khaki1",
      "darkgoldenrod1",
      "orangered",
      "red"
    )
  )(100),
  scale = "none",
  main = "Square Root Variant Frequency Heatmap",   # ** SHOW title **
  show_rownames = TRUE,                             # ** SHOW row labels **
  show_colnames = TRUE,                             # ** SHOW column labels **
  legend = TRUE,                                    # ** SHOW legend **
  border_color = "white",                           # ** SHOW gridlines **
  angle_col = 315,
)

hetAdj_heatmap <- pheatmap(
  ord_hetAdj_heatmapMatrix,
  cluster_rows = FALSE, # ** HIDE row-wise (gene) hierachical clustering **
  cluster_cols = TRUE, # ** SHOW column-wise (sample) hierarchical clustering **
  color = colorRampPalette(
    c(
      "steelblue",
      "green",
      "lawngreen",
      "khaki1",
      "darkgoldenrod1",
      "orangered",
      "red"
    )
  )(100),
  scale = "none",
  main = "Het. Mut. Adj. Sqrt. Var. Freq. Heatmap",     # ** SHOW title **
  show_rownames = TRUE,                                 # ** SHOW row labels **
  show_colnames = TRUE,                                 # ** SHOW column labels **
  legend = TRUE,                                        # ** SHOW legend **
  border_color = "white",                               # ** SHOW gridlines **
  angle_col = 315,
)

nonFSAdj_heatmap <- pheatmap(
  ord_nonFSAdj_heatmapMatrix,
  cluster_rows = FALSE, # ** HIDE row-wise (gene) hierachical clustering **
  cluster_cols = TRUE, # ** SHOW column-wise (sample) hierarchical clustering **
  color = colorRampPalette(
    c(
      "steelblue",
      "green",
      "lawngreen",
      "khaki1",
      "darkgoldenrod1",
      "orangered",
      "red"
    )
  )(100),
  scale = "none",
  main = "Non-FS Mut. Adj. Sqrt. Var. Freq. Heatmap",   # ** SHOW title **
  show_rownames = TRUE,                                 # ** SHOW row labels **
  show_colnames = TRUE,                                 # ** SHOW column labels **
  legend = TRUE,                                        # ** SHOW legend **
  border_color = "white",                               # ** SHOW gridlines **
  angle_col = 315,
)



# Extract column-wise (sample) hierarchical clustering order from heat maps as integer vectors
# Re-order samples accordingly in a new character vector

nonAdj_colOrder <- nonAdj_heatmap$tree_col$order
nonAdj_ordSampleNames <- colnames(nonAdj_heatmapMatrix)[nonAdj_colOrder]

hetAdj_colOrder <- hetAdj_heatmap$tree_col$order
hetAdj_ordSampleNames <- colnames(hetAdj_heatmapMatrix)[hetAdj_colOrder]

nonFSAdj_colOrder <- nonFSAdj_heatmap$tree_col$order
nonFSAdj_ordSampleNames <- colnames(nonFSAdj_heatmapMatrix)[nonFSAdj_colOrder]


# Convert sample vs. square root variant frequency data from heat map matrix to a data frame format,
# which will be used to create the mean square root variant frequency bar chart.

nonAdj_colMeans_df <- data.frame(SampleName = colnames(nonAdj_heatmapMatrix),
                                 MeanSqrtFreq = nonAdj_colMeans)

hetAdj_colMeans_df <- data.frame(SampleName = colnames(hetAdj_heatmapMatrix),
                                 MeanSqrtFreq = hetAdj_colMeans)

nonFSAdj_colMeans_df <- data.frame(SampleName = colnames(nonFSAdj_heatmapMatrix),
                                 MeanSqrtFreq = nonFSAdj_colMeans)


# Re-order samples in the data frame as a factor

nonAdj_colMeans_df$SampleName <- factor(nonAdj_colMeans_df$SampleName, levels = nonAdj_ordSampleNames)

hetAdj_colMeans_df$SampleName <- factor(hetAdj_colMeans_df$SampleName, levels = hetAdj_ordSampleNames)

nonFSAdj_colMeans_df$SampleName <- factor(nonFSAdj_colMeans_df$SampleName, levels = nonFSAdj_ordSampleNames)


# Re-order heat map matrix rows (genes) by mean square root variant frequency, decreasing top to bottom
# AND re-order columns (samples) by heirarchical clustering

ord_nonAdj_heatmapMatrix <- nonAdj_heatmapMatrix[nonAdj_rowOrder_dec, nonAdj_colOrder]

ord_hetAdj_heatmapMatrix <- hetAdj_heatmapMatrix[hetAdj_rowOrder_dec , hetAdj_colOrder]

ord_nonFSAdj_heatmapMatrix <- nonFSAdj_heatmapMatrix[nonFSAdj_rowOrder_dec , nonFSAdj_colOrder]


# Create upright bar charts for column-wise (sample) mean square root variant frequency

nonAdj_colMeans_barChart <- ggplot(nonAdj_colMeans_df, aes(x = SampleName, y = MeanSqrtFreq)) +
  geom_bar(stat = "identity", fill = "steelblue4") +
  scale_y_continuous(position = "right") +
  labs(y = "Mean Sqrt.\nVar. Freq. %") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 7, hjust = 0.5), # ** ADJUST IF NEEDED **
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 10, r = 51, b = 0, l = 4), # ** ADJUST IF NEEDED **
  )

hetAdj_colMeans_barChart <- ggplot(hetAdj_colMeans_df, aes(x = SampleName, y = MeanSqrtFreq)) +
  geom_bar(stat = "identity", fill = "steelblue4") +
  scale_y_continuous(position = "right") +
  labs(y = "Mean Sqrt.\nVar. Freq. %") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 7, hjust = 0.5), # ** ADJUST IF NEEDED **
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 10, r = 51, b = 0, l = 4), # ** ADJUST IF NEEDED **
  )

nonFSAdj_colMeans_barChart <- ggplot(nonFSAdj_colMeans_df, aes(x = SampleName, y = MeanSqrtFreq)) +
  geom_bar(stat = "identity", fill = "steelblue4") +
  scale_y_continuous(position = "right") +
  labs(y = "Mean Sqrt.\nVar. Freq. %") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 7, hjust = 0.5), # ** ADJUST IF NEEDED **
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 10, r = 51, b = 0, l = 4), # ** ADJUST IF NEEDED **
  )

# Create upright bar charts for row-wise (sample) mean square root variant frequency

nonAdj_rowMeans_barChart <- ggplot(nonAdj_rowMeans_df, aes(x = MeanSqrtFreq, y = GeneName)) +
  geom_bar(stat = "identity", fill = "maroon") +
  scale_x_continuous(position = "bottom") +
  labs(x = "Mean Sqrt.\nVar. Freq. %") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 7, hjust = 0.5), # ** ADJUST IF NEEDED **
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 4, r = 0, b = 79, l = 10), # ** ADJUST IF NEEDED **
  )

hetAdj_rowMeans_barChart <- ggplot(hetAdj_rowMeans_df, aes(x = MeanSqrtFreq, y = GeneName)) +
  geom_bar(stat = "identity", fill = "maroon") +
  scale_x_continuous(position = "bottom") +
  labs(x = "Mean Sqrt.\nVar. Freq. %") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 7, hjust = 0.5), # ** ADJUST IF NEEDED **
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 4, r = 0, b = 79, l = 10), # ** ADJUST IF NEEDED **
  )

nonFSAdj_rowMeans_barChart <- ggplot(nonFSAdj_rowMeans_df, aes(x = MeanSqrtFreq, y = GeneName)) +
  geom_bar(stat = "identity", fill = "maroon") +
  scale_x_continuous(position = "bottom") +
  labs(x = "Mean Sqrt.\nVar. Freq. %") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 7, hjust = 0.5), # ** ADJUST IF NEEDED **
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 4, r = 0, b = 79, l = 10), # ** ADJUST IF NEEDED **
  )

nonAdj_rowMeans_barChart_bigMar <- ggplot(nonAdj_rowMeans_df, aes(x = MeanSqrtFreq, y = GeneName)) +
  geom_bar(stat = "identity", fill = "maroon") +
  scale_x_continuous(position = "bottom") +
  labs(x = "Mean Sqrt.\nVar. Freq. %") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 7, hjust = 0.5), # ** ADJUST IF NEEDED **
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 68, r = 0, b = 79, l = 10), # ** ADJUST IF NEEDED **
  )

hetAdj_rowMeans_barChart_bigMar <- ggplot(hetAdj_rowMeans_df, aes(x = MeanSqrtFreq, y = GeneName)) +
  geom_bar(stat = "identity", fill = "maroon") +
  scale_x_continuous(position = "bottom") +
  labs(x = "Mean Sqrt.\nVar. Freq. %") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 7, hjust = 0.5), # ** ADJUST IF NEEDED **
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 68, r = 0, b = 79, l = 10), # ** ADJUST IF NEEDED **
  )

nonFSAdj_rowMeans_barChart_bigMar <- ggplot(nonFSAdj_rowMeans_df, aes(x = MeanSqrtFreq, y = GeneName)) +
  geom_bar(stat = "identity", fill = "maroon") +
  scale_x_continuous(position = "bottom") +
  labs(x = "Mean Sqrt.\nVar. Freq. %") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 7, hjust = 0.5), # ** ADJUST IF NEEDED **
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 68, r = 0, b = 79, l = 10), # ** ADJUST IF NEEDED **
  )


# The heatmaps below do NOT contain a title or a column-wise dendogram.
# The rows (genes) are ordered by mean square root variant frequency, decreasing top to bottom.
# The columns (samples) are ordered by hierarchical clustering as generated by pheatmap.

naked_nonAdj_heatmap <- pheatmap(
  ord_nonAdj_heatmapMatrix,
  cluster_rows = FALSE, # ** HIDE row-wise (gene) hierachical clustering **
  cluster_cols = FALSE, # ** HIDE column-wise (sample) hierarchical clustering **
  color = colorRampPalette(
    c(
      "steelblue",
      "green",
      "lawngreen",
      "khaki1",
      "darkgoldenrod1",
      "orangered",
      "red"
    )
  )(100),
  scale = "none",
  main = NA,                            # ** HIDE title **
  show_rownames = TRUE,                 # ** SHOW row labels **
  show_colnames = TRUE,                 # ** SHOW column labels **
  angle_col = 315,
  legend = TRUE,                        # ** SHOW legend **
  border_color = "white",               # ** SHOW gridlines **
  plot.margin = margin(0, 0, 0, 0),     # ** ADJUST IF NEEDED **
)

naked_hetAdj_heatmap <- pheatmap(
  ord_hetAdj_heatmapMatrix,
  cluster_rows = FALSE, # ** HIDE row-wise (gene) hierachical clustering **
  cluster_cols = FALSE, # ** HIDE column-wise (sample) hierarchical clustering **
  color = colorRampPalette(
    c(
      "steelblue",
      "green",
      "lawngreen",
      "khaki1",
      "darkgoldenrod1",
      "orangered",
      "red"
    )
  )(100),
  scale = "none",
  main = NA,                            # ** HIDE title **
  show_rownames = TRUE,                 # ** SHOW row labels **
  show_colnames = TRUE,                 # ** SHOW column labels **
  angle_col = 315,
  legend = TRUE,                        # ** SHOW legend **
  border_color = "white",               # ** SHOW gridlines **
  plot.margin = margin(0, 0, 0, 0),     # ** ADJUST IF NEEDED **
)

naked_nonFSAdj_heatmap <- pheatmap(
  ord_nonFSAdj_heatmapMatrix,
  cluster_rows = FALSE, # ** HIDE row-wise (gene) hierachical clustering **
  cluster_cols = FALSE, # ** HIDE column-wise (sample) hierarchical clustering **
  color = colorRampPalette(
    c(
      "steelblue",
      "green",
      "lawngreen",
      "khaki1",
      "darkgoldenrod1",
      "orangered",
      "red"
    )
  )(100),
  scale = "none",
  main = NA,                            # ** HIDE title **
  show_rownames = TRUE,                 # ** SHOW row labels **
  show_colnames = TRUE,                 # ** SHOW column labels **
  angle_col = 315,
  legend = TRUE,                        # ** SHOW legend **
  border_color = "white",               # **SHOW gridlines **
  plot.margin = margin(0, 0, 0, 0),     # ** ADJUST IF NEEDED **
)

BLANK_PLOT <- ggplot() +
  theme_void()






# DISPLAY HEAT MAPS #################################################

# ** DISPLAY INDIVIDUAL HEAT MAP W/ DENDOGRAM + MEAN SQUARE ROOT VARIANT FREQUENCY BAR CHARTS ** 

grid.arrange(
  nonAdj_rowMeans_barChart_bigMar, nonAdj_heatmap$gtable,
  ncol = 2,
  nrow = 1,
  heights = c(10),
  widths = c(1, 10)
)

grid.arrange(
  hetAdj_rowMeans_barChart_bigMar, hetAdj_heatmap$gtable,
  ncol = 2,
  nrow = 1,
  heights = c(10),
  widths = c(1, 10)
)

grid.arrange(
  nonFSAdj_rowMeans_barChart_bigMar, nonFSAdj_heatmap$gtable,
  ncol = 2,
  nrow = 1,
  heights = c(10),
  widths = c(1, 10)
)

# ** DISPLAY ALL HEAT MAPS W/ DENDOGRAMS + MEAN SQUARE ROOT VARIANT FREQUENCY BAR CHARTS ** 

grid.arrange(
  nonAdj_rowMeans_barChart_bigMar, nonAdj_heatmap$gtable, hetAdj_rowMeans_barChart_bigMar, hetAdj_heatmap$gtable, nonFSAdj_rowMeans_barChart_bigMar, nonFSAdj_heatmap$gtable,
  ncol = 6,
  nrow = 1,
  heights = c(10),
  widths = c(1, 10, 1, 10, 1, 10)
)


# ** DISPLAY INDIVIDUAL NAKED HEAT MAP + MEAN SQUARE ROOT VARIANT FREQUENCY BAR CHARTS ** 

grid.arrange(
  BLANK_PLOT, nonAdj_colMeans_barChart,
  nonAdj_rowMeans_barChart, naked_nonAdj_heatmap$gtable,
  ncol = 2,
  nrow = 2,
  heights = c(1, 10),
  widths = c(1, 10)
)

grid.arrange(
  BLANK_PLOT, hetAdj_colMeans_barChart,
  hetAdj_rowMeans_barChart, naked_hetAdj_heatmap$gtable,
  ncol = 2,
  nrow = 2,
  heights = c(1, 10),
  widths = c(1, 10)
)

grid.arrange(
  BLANK_PLOT, nonFSAdj_colMeans_barChart,
  nonFSAdj_rowMeans_barChart, naked_nonFSAdj_heatmap$gtable,
  ncol = 2,
  nrow = 2,
  heights = c(1, 10),
  widths = c(1, 10)
)

# ** DISPLAY ALL INDIVIDUAL NAKED HEAT MAPS + MEAN SQUARE ROOT VARIANT FREQUENCY BAR CHARTS ** 

grid.arrange(
  BLANK_PLOT, nonAdj_colMeans_barChart, BLANK_PLOT, hetAdj_colMeans_barChart, BLANK_PLOT, nonFSAdj_colMeans_barChart,
  nonAdj_rowMeans_barChart, naked_nonAdj_heatmap$gtable, hetAdj_rowMeans_barChart, naked_hetAdj_heatmap$gtable, nonFSAdj_rowMeans_barChart, naked_nonFSAdj_heatmap$gtable,
  ncol = 6,
  nrow = 2,
  heights = c(1, 10),
  widths = c(1, 10, 1, 10, 1, 10)
)




# THE END :) #################################################