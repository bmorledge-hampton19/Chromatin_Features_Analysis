library(data.table)
library(ggplot2)

# Given the path to a binned counts file (and potentially paths to files containing info on the domains
# or background binned counts) return a binned counts data.table.
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


# Create a graph for each chromosome comparing across two data sets
# NOTE: Kinda old.  Uses weird base-R plotting.
binInChromosomesAcrossDataSets = function(binnedCountsTable, compBinCountsTable, title = "Binned Counts:",
                                          firstDataSet = "Repair", compDataSet = "Damage") {
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

# Create a graph for each (given) chromosome set in chromosomeSets, that uses the log ratio
# between the two data sets and colors based on majority chromatin domain.
# AKA: The whole plotting enchilada with ggplot
# (Could probably make more modular in the future is some features are not desired.)
# chromosomeSets example: list("chrX", "chrY", c("chr2L","chr2R"), c("chr3L","chr3R"), "chr4")
createggplotBinPlots = function(binnedCountsTable, chromosomeSets, yAxisLabel = "Log Ratio",
                                title = "", ylim = NULL, colorPlot = TRUE) {

  for (chromosomeSet in chromosomeSets) {

    plot = ggplot(binnedCountsTable[Chromosome %in% chromosomeSet],
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
        scale_fill_identity() + scale_color_identity()
        #scale_fill_manual(values = c("BLACK" = "black", "BLUE" = "blue", "GREEN" = "green",
        #                             "RED" = "red", "YELLOW" = "gold", "GRAY" = "gray")) +
        #scale_color_manual(values = c("BLACK" = "black", "BLUE" = "blue", "GREEN" = "green",
        #                              "RED" = "red", "YELLOW" = "gold", "GRAY" = "gray"))
    }

    print(plot)

  }
}


# Plot the distribution of log ratio values for each individual color domain throughout the genome.
# Plotted as either a scatter plot or box plot based on given "pseudo constant" values.
SCATTER = 1
BOXPLOT = 2
plotLogRatioDistribution = function(binnedCountsTable, plotType, title = "",
                                    yAxisLabel = "Log Ratio", ylim = NULL) {

  # Deprecated?
  #
  # Get the median for each of the color domains (Not really used right now though...)
  # colorLogRatioMedians = binnedCountsTable[Majority_Domain_Color != "GRAY", .(Median = median(Log_Ratio)),
  #                                         Majority_Domain_Color]
  #
  # Alternative to scale_color_identity with legend for median values.
  # scale_color_manual(name = "Medians",
  #                    labels = setNames(round(colorLogRatioMedians$Median,4),
  #                                      colorLogRatioMedians$Majority_Domain_Color),
  #                    values = c("BLACK" = "black", "BLUE" = "blue", "GREEN" = "green",
  #                               "RED" = "red", "YELLOW" = "gold"))

  # Basic plotting framework.
  thisPlot = ggplot(binnedCountsTable[Majority_Domain_Color != "GRAY"],
                aes(Majority_Domain_Color, Log_Ratio, color = Majority_Domain_Color)) +
    scale_color_identity()

  # Plot as scatter
  if (plotType == SCATTER) {
    thisPlot = thisPlot + geom_jitter(width = 0.2, height = 0, shape = 1, size = 2) +
      stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 2, colour = "purple")
  }
  # Plot as boxes (without outliers)
  else if (plotType == BOXPLOT) {
  thisPlot = thisPlot + geom_boxplot(outlier.shape = NA)
  }
  else stop("Invalid \"plotType\" argument given.")

  # Finishing touches
  thisPlot = thisPlot + labs(title = title, y = yAxisLabel) +
    coord_cartesian(ylim = ylim) +
    theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title = element_text(size = 15),
          axis.text.x = element_text(size = 15), axis.title.x = element_blank(),
          legend.text = element_text(size = 12), legend.title = element_text(size = 14))

  # Go! (Do we need this?)
  thisPlot

}


# Graphs median log_ratio counts for different domains across different time points.
# Requires two parallel lists of binned counts tables and time points.
plotLogRatioMediansOverTime = function(binnedCountsTables, timePoints, title = "",
                                       xAxisLabel = "Timepoint", yAxisLabel = "Log Ratio", ylim = NULL) {

  masterMedianTable = rbindlist(mapply(getMedianTable, binnedCountsTables, timePoints))

  # Graph a line for each color.
  ggplot(masterMedianTable, aes(Time, Median, color = Majority_Domain_Color)) +
    scale_color_identity() +
    labs(title = title, x = xAxisLabel, y = yAxisLabel) +
    geom_line() + geom_point() +
    coord_cartesian(ylim = ylim) +
    theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title = element_text(size = 15))

}


# A function for producing a table of median log_ratio values with the domain color and timepoint included as columns.
getMedianTable = function(binnedCountsTable, timePoint) {

  return(binnedCountsTable[Majority_Domain_Color != "GRAY", .(Median = median(Log_Ratio), Time = timePoint),
                           Majority_Domain_Color])

}
