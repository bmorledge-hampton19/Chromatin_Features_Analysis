library(data.table)
library(ggplot2)

# Given the path to a binned counts file (and potentially paths to files containing info on the domains
# or background binned counts) return a binned counts data.table.
parseBinData = function(binnedCountsFilePath, binnedColorDomainsFilePath = NA, backgroundBinCountsFilePath = NA) {

  # Read in the data
  binnedCountsTable = fread(binnedCountsFilePath)

  # Create a bin start column from the bin range column.
  binnedCountsTable[,Bin_Start := vapply(strsplit(`Bin_Start-End`, '-'),
                                         function(x) as.integer(x[1]), FUN.VALUE = integer(1))]

  # Create a column for counts normalized by total counts
  binnedCountsTable[Feature_Counts == 0, Feature_Counts := 1]# Pseudo-counts!
  binnedCountsTable[,Counts_Per_Million := Feature_Counts / sum(Feature_Counts) * 10^6]

  # Add a column for the majority color domain if the binned color domains file path was given.
  if (!is.na(binnedColorDomainsFilePath)) {
    binnedColorDomainsTable = fread(binnedColorDomainsFilePath)
    binnedCountsTable[,Majority_Domain_Color := binnedColorDomainsTable$Domain_Color]
  }

  # Add in a complementary (background) data set (If we have the relevant file).
  if (!is.na(backgroundBinCountsFilePath)) {
    backgroundBinCountsTable = fread(backgroundBinCountsFilePath)

    backgroundBinCountsTable[,Bin_Start := vapply(strsplit(`Bin_Start-End`, '-'),
                                           function(x) as.integer(x[1]), FUN.VALUE = integer(1))]

    backgroundBinCountsTable[Feature_Counts == 0, Feature_Counts := 1]# Pseudo-counts!
    backgroundBinCountsTable[,Counts_Per_Million := Feature_Counts / sum(Feature_Counts) * 10^6]

    # Create a log ratio column in the main counts file using the complementary data set.
    binnedCountsTable[, Log_Ratio := log(Counts_Per_Million/backgroundBinCountsTable$Counts_Per_Million,2)]
  }

  return(binnedCountsTable)

}


# Given a path to a file with information on gene bins, return a binned counts data.table
parseGeneBinData = function(geneBinsCountsFilePath, backgroundFilePath = NA) {

  # Read in the data
  geneBinsCountsTable = fread(geneBinsCountsFilePath)

  # Create columns for counts normalized by total counts
  geneBinsCountsTable[Coding_Strand_Counts == 0, Coding_Strand_Counts := 1]# Pseudo-counts!
  geneBinsCountsTable[Noncoding_Strand_Counts == 0, Noncoding_Strand_Counts := 1]# More Pseudo-counts!
  totalCounts = sum(geneBinsCountsTable$Coding_Strand_Counts) +
                sum(geneBinsCountsTable$Noncoding_Strand_Counts)
  geneBinsCountsTable[,Noncoding_Counts_Per_Million := Noncoding_Strand_Counts / totalCounts * 10^6]
  geneBinsCountsTable[,Coding_Counts_Per_Million := Coding_Strand_Counts / totalCounts * 10^6]

  # Add in a complementary (background) data set (If we have the relevant file).
  if (!is.na(backgroundFilePath)) {
    backgroundCountsTable = fread(backgroundFilePath)

    backgroundCountsTable[Coding_Strand_Counts == 0, Coding_Strand_Counts := 1]# Pseudo-counts!
    backgroundCountsTable[Noncoding_Strand_Counts == 0, Noncoding_Strand_Counts := 1]# More Pseudo-counts!
    totalCounts = sum(backgroundCountsTable$Coding_Strand_Counts) +
      sum(backgroundCountsTable$Noncoding_Strand_Counts)
    backgroundCountsTable[,Coding_Counts_Per_Million := Coding_Strand_Counts / totalCounts * 10^6]
    backgroundCountsTable[,Noncoding_Counts_Per_Million := Noncoding_Strand_Counts / totalCounts * 10^6]

    # Create a log ratio column in the main counts file using the complementary data set.
    geneBinsCountsTable[, Coding_Log_Ratio := log(Coding_Counts_Per_Million/
                                                  backgroundCountsTable$Coding_Counts_Per_Million,2)]
    geneBinsCountsTable[, Noncoding_Log_Ratio := log(Noncoding_Counts_Per_Million/
                                                     backgroundCountsTable$Noncoding_Counts_Per_Million,2)]

    geneBinsCountsTable[, TS_Vs_NTS_Log_Ratio := Noncoding_Log_Ratio - Coding_Log_Ratio]
  }

  return(geneBinsCountsTable)

}


plotGeneBinData = function(geneBinsCountsTable, title = "", xAxisLabel = "Gene Fraction Bin",
                           yAxisLabel = "Log Ratio", ylim = NULL, yData1 = "Coding_Log_Ratio",
                           yData2 = "Noncoding_Log_Ratio", yData3 = "TS_Vs_NTS_Log_Ratio",
                           plotYData3Only = TRUE) {

  if ("Color_Domain" %in% colnames(geneBinsCountsTable)) {
    geneBinPlot = ggplot(geneBinsCountsTable[Color_Domain != "GRAY"], aes(Gene_Fraction, color = Color_Domain))
  } else {
    geneBinPlot = ggplot(geneBinsCountsTable, aes(Gene_Fraction))
  }

  geneBinPlot = geneBinPlot +
    scale_color_identity()

  if (plotYData3Only) {

    geneBinPlot = geneBinPlot +
      geom_line(aes_string(y = yData3)) +
      geom_point(aes_string(y = yData3))

  } else {

    geneBinPlot = geneBinPlot +
      geom_line(aes_string(y = yData1, linetype = shQuote("dashed"))) +
      geom_point(aes_string(y = yData1)) +
      geom_line(aes_string(y = yData2, linetype = shQuote("solid"))) +
      geom_point(aes_string(y = yData2)) +
      scale_linetype_identity(guide = "legend", name = "", breaks = c("solid", "dashed"),
                              labels = c("Transcribed Strand", "Non-Transcribed Strand"))

  }

  geneBinPlot = geneBinPlot +
    labs(title = title, x = xAxisLabel, y = yAxisLabel) +
    coord_cartesian(ylim = ylim) +
    scale_x_continuous(breaks = -2:8 + 0.5, minor_breaks = NULL,
                       labels = c('','',"TSS",'','','','','',"TES",'','')) +
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 15), axis.text = element_text(size = 12),
          legend.text = element_text(size = 12),)

  print(geneBinPlot)

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
createGgplotBinPlots = function(binnedCountsTable, chromosomeSets, yData = "Log_Ratio",
                                yAxisLabel = "Log Ratio", title = "", ylim = NULL, colorPlot = TRUE) {

  for (chromosomeSet in chromosomeSets) {

    plot = ggplot(binnedCountsTable[Chromosome %in% chromosomeSet],
                  aes_string("Bin_Start", yData)) +
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
plotLogRatioDistribution = function(binnedCountsTable, plotType, title = "", yData = "Log_Ratio",
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
  thisPlot = ggplot(binnedCountsTable[Majority_Domain_Color != "GRAY" & Majority_Domain_Color != "WHITE"],
                aes_string("Majority_Domain_Color", yData, color = "Majority_Domain_Color")) +
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
          axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 12), legend.title = element_text(size = 14))

  print(thisPlot)

}


# Graphs median log_ratio counts for different domains across different time points.
# Requires two parallel lists of binned counts tables and time points (as numeric values).
plotLogRatioMediansOverTime = function(binnedCountsTables, timePoints, title = "",
                                       xAxisLabel = "Timepoint", yAxisLabel = "Log Ratio", ylim = NULL) {

  masterMedianTable = rbindlist(mapply(getMedianTable, binnedCountsTables, timePoints, SIMPLIFY = FALSE))

  # Graph a line for each color.
  ggplot(masterMedianTable, aes(Time, Median, color = Majority_Domain_Color)) +
    scale_color_identity() +
    labs(title = title, x = xAxisLabel, y = yAxisLabel) +
    geom_line() + geom_point() +
    coord_cartesian(ylim = ylim) +
    theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title = element_text(size = 15),
          axis.text = element_text(size = 12), legend.text = element_text(size = 12),)

}


# A function for producing a table of median log_ratio values with the domain color and
# time point (as a numeric value) included as columns.
getMedianTable = function(binnedCountsTable, timePoint) {

  return(binnedCountsTable[Majority_Domain_Color != "GRAY" & Majority_Domain_Color != "WHITE",
                           .(Median = median(Log_Ratio), Time = timePoint),
                           Majority_Domain_Color])

}
