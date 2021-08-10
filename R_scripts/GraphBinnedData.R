library(data.table)
library(ggplot2)

parseBinData = function(binnedCountsFilePath, binnedColorDomainsFilePath = NA, backgroundBinCountsFilePath = NA) {

  # Read in the data and derive a name for it
  binnedCountsFilePath = choose.files(multi = FALSE, caption = "Main Binned File")
  binnedCountsTable = fread(binnedCountsFilePath)

  # Create a bin start column from the bin range column.
  binnedCountsTable[,Bin_Start := vapply(strsplit(`Bin_Start-End`, '-'),
                                         function(x) as.integer(x[1]), FUN.VALUE = integer(1))]

  # Create a column for counts normalized by total counts
  binnedCountsTable[Feature_Counts == 0, Feature_Counts := 1]# Pseudo-counts!
  binnedCountsTable[,Counts_Per_Million := Feature_Counts / sum(Feature_Counts) * 10^6]

  # Add a column for the majority color domain if the binned color domains file path was given.
  if (!is.na(binnedColorDomainsFilePath)) {
    binnedChromatinDomainColorsFilePath = choose.files(multi = FALSE, caption = "Bin Color Domains File")
    binnedChromatinDomainColorsTable = fread(binnedChromatinDomainColorsFilePath)
    binnedCountsTable[,Majority_Domain_Color := binnedChromatinDomainColorsTable$Domain_Color]
  }

  # Add in a complementary (background) data set (If we have the relevant file).
  if (!is.na(backgroundBinCountsFilePath)) {
    backgroundBinCountsFilePath = choose.files(multi = FALSE, caption = "Background Binned File")
    backgroundBinCountsTable = fread(backgroundBinCountsFilePath)
    backgroundBinCountsDataName = strsplit(basename(backgroundBinCountsFilePath),'.', fixed = T)[[1]][1]

    backgroundBinCountsTable[,Bin_Start := vapply(strsplit(`Bin_Start-End`, '-'),
                                           function(x) as.integer(x[1]), FUN.VALUE = integer(1))]


    backgroundBinCountsTable[Feature_Counts == 0, Feature_Counts := 1]# Pseudo-counts!
    backgroundBinCountsTable[,Counts_Per_Million := Feature_Counts / sum(Feature_Counts) * 10^6]

    # Create a log ratio column in the main counts file using the complementary data set.
    binnedCountsTable[, Log_Ratio := log(Counts_Per_Million/backgroundBinCountsTable$Counts_Per_Million,2)]
  }

  return(binnedCountsTable)

}

# Create a graph to bin in each chromosome for one data set.
binInChromosomes = function(binnedCountsTable, baseTitle = "Binned Counts:") {
  for (chromosome in unique(binnedCountsTable$Chromosome)) {

    # Skip those funky chromosome fragments...
    if (grepl('_',chromosome)) next

    relevantCountsTable = binnedCountsTable[Chromosome == chromosome]

    barplot(relevantCountsTable[, Feature_Counts],
            names.arg = relevantCountsTable[, Bin_Start],
            xlab = "Chromosome Position", ylab = "Counts", main = paste(baseTitle,chromosome))

  }
}

# Default comparison is damage and repair.
firstDataSet = "Repair"
compDataSet = "Damage"

# Create a graph for each chromosome comparing across two data sets
binInChromosomesAcrossDataSets = function(binnedCountsTable, compBinCountsTable, title = "Binned Counts:") {
  for (chromosome in unique(binnedCountsTable$Chromosome)) {

    # Skip those funky chromosome fragments...
    if (grepl('_',chromosome)) next

    # Plot the first data set...
    relevantCountsTable = binnedCountsTable[Chromosome == chromosome]

    barplot(relevantCountsTable[, Counts_Per_Million],
            names.arg = relevantCountsTable[, Bin_Start], col = "#ca0020", border = "#ca0020",
            xlab = "Chromosome Position", ylab = "Counts Per Million Total", main = paste(title, chromosome))

    # Then overlay the second data set!
    par(new = T)

    compRelevantCountsTable = compBinCountsTable[Chromosome == chromosome]

    barplot(compRelevantCountsTable[, Counts_Per_Million], yaxt = 'n',
            names.arg = compRelevantCountsTable[, Bin_Start], col = "#0571b0", border = "#0571b0")

    # What if we keep track of where they intersect?
    par(new = T)
    barplot(pmin(relevantCountsTable$Counts_Per_Million, compRelevantCountsTable$Counts_Per_Million), yaxt = 'n',
            names.arg = compRelevantCountsTable[, Bin_Start], col = "#683968", border = "#683968")


    # Finally, add a legend
    legend("topleft", c(firstDataSet, compDataSet), col=c("#ca0020", "#0571b0"),
           lwd=c(3,3), cex = 0.8, bg = "white")

  }
}

# Create a graph for each chromosome that uses the log ratio between the two data sets and
# colors based on majority chromatin domain.
yAxisLabel = "Repair to Damage log Ratio"
colorPlot = TRUE
ylim = NULL
for (chromosome in list("chrX", "chrY", c("chr2L","chr2R"), c("chr3L","chr3R"), "chr4")) {

  plot = ggplot(binnedCountsTable[Chromosome %in% chromosome],
                aes(Bin_Start, Log_Ratio)) +
    geom_bar(stat = "identity", color = "Black", fill = "Black") +
    labs(title = title, x = "Chromosome Position", y = yAxisLabel) +
    facet_grid(~Chromosome, space = "free", scales = "free") +
    coord_cartesian(ylim = ylim) +
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 15), axis.text = element_text(size = 12),
          strip.text = element_text(size = 15), panel.spacing = unit(0.1, "lines"), legend.position = "none")

  if (colorPlot) {
    plot = plot + geom_bar(aes(color = Majority_Domain_Color, fill = Majority_Domain_Color), stat = "identity") +
      scale_fill_manual(values = c("BLACK" = "black", "BLUE" = "blue", "GREEN" = "green",
                                   "RED" = "red", "YELLOW" = "gold", "GRAY" = "gray")) +
      scale_color_manual(values = c("BLACK" = "black", "BLUE" = "blue", "GREEN" = "green",
                                    "RED" = "red", "YELLOW" = "gold", "GRAY" = "gray"))
  }

  print(plot)

}

title = "PLACEHOLDER TITLE"
ylim = NULL

# Get the median for each of the color domains.
colorLogRatioMedians = binnedCountsTable[Majority_Domain_Color != "GRAY", .(Median = median(Log_Ratio)),
                                         Majority_Domain_Color]

# Plot as scatter.
ggplot(binnedCountsTable[Majority_Domain_Color != "GRAY"],
       aes(Majority_Domain_Color, Log_Ratio, color = Majority_Domain_Color)) +
  scale_color_manual(name = "Medians",
                     labels = setNames(round(colorLogRatioMedians$Median,4),
                                       colorLogRatioMedians$Majority_Domain_Color),
                     values = c("BLACK" = "black", "BLUE" = "blue", "GREEN" = "green",
                                "RED" = "red", "YELLOW" = "gold")) +
  geom_jitter(width = 0.2, shape = 1, size = 2) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 2, colour = "purple") +
  labs(title = title, x = xAxislabel, y = yAxisLabel) +
  coord_cartesian(ylim = ylim) +
  theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 15), axis.title.x = element_blank(),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14))

# Plot as boxes (without outliers)
ggplot(binnedCountsTable[Majority_Domain_Color != "GRAY"],
       aes(Majority_Domain_Color, Log_Ratio, color = Majority_Domain_Color)) +
  geom_boxplot(outlier.shape = NA) +
  scale_color_manual(name = "Medians",
                     labels = setNames(round(colorLogRatioMedians$Median,4),
                                       colorLogRatioMedians$Majority_Domain_Color),
                     values = c("BLACK" = "black", "BLUE" = "blue", "GREEN" = "green",
                                "RED" = "red", "YELLOW" = "gold")) +
  labs(title = title, y = yAxisLabel) +
  coord_cartesian(ylim = ylim) +
  theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 15), axis.title.x = element_blank(),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14))

# Use the following commands (modified as necessary) to create a table of medians by time points and color domains.
masterMedianTable = data.table()

# Read in a new binned counts file and get the median log ratio values for each color (except GRAY).
# Depends on the background (compBinCounts) and color domains files already being read in.
timepoint = NA # Make sure to set this manually for each time point!
binnedCountsFilePath = choose.files(multi = FALSE)
binnedCountsTable = fread(binnedCountsFilePath)
binnedCountsTable[,Counts_Per_Million := Feature_Counts / sum(Feature_Counts)]
binnedCountsTable[,Majority_Domain_Color := binnedChromatinDomainColorsTable$Domain_Color]
compBinCountsTable[,Counts_Per_Million := Feature_Counts / sum(Feature_Counts)]
binnedCountsTable[, Log_Ratio := log(Counts_Per_Million/compBinCountsTable$Counts_Per_Million,2)]

newMedianTable = binnedCountsTable[Majority_Domain_Color != "GRAY", .(Median = median(Log_Ratio), Time = timepoint),
                                   Majority_Domain_Color]

# I think this is really inefficient for a lot of time points, but I think it should do alright for a few.
masterMedianTable = rbindlist(list(masterMedianTable,newMedianTable))

# Graph a line for each color.
ggplot(masterMedianTable, aes(Time, Median, color = Majority_Domain_Color)) +
  scale_color_manual(values = c("BLACK" = "black", "BLUE" = "blue", "GREEN" = "green",
                                "RED" = "red", "YELLOW" = "gold"), guide = FALSE) +
  labs(title = title, x = xAxisLabel, y = yAxisLabel) +
  geom_line() + geom_point() +
  theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title = element_text(size = 15))
