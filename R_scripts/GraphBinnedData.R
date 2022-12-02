library(data.table)
library(ggplot2)


# Default text scaling
defaultTextScaling = theme(plot.title = element_text(size = 26, hjust = 0.5),
                           axis.title = element_text(size = 22), axis.text = element_text(size = 18),
                           legend.title = element_text(size = 22), legend.text = element_text(size = 18),
                           strip.text = element_text(size = 22))

# Blank background theme
blankBackground = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Actual domain colors
domainColors = c("BLACK" = "black", "black" = "black", "BLUE" = "blue", "blue" = "blue",
                 "GREEN" = "forestgreen", "green" = "forestgreen",
                 "RED" = "red", "red" = "red", "YELLOW" = "gold2", "yellow" = "gold2",
                 "GRAY" = "gray", "gray" = "gray")


# computes a scaling factor from given raw and background file paths.
# Right now, this function assumes the input file paths have a "Feature_Counts" column to determine counts from.
getScalingFactor = function(rawCountsFilePath, backgroundCountsFilePath) {

  rawCountsTable = fread(rawCountsFilePath)
  backgroundCountsTable = fread(backgroundCountsFilePath)

  return(sum(backgroundCountsTable$Feature_Counts)/sum(rawCountsTable$Feature_Counts))

}


# Given the path to a binned counts file (and potentially paths to files containing info on the domains
# or background binned counts) return a binned counts data.table.
parseBinData = function(binnedCountsFilePath, binnedColorDomainsFilePath = NA,
                        backgroundBinCountsFilePath = NA, scalingFactor = NULL) {

  # Read in the data and set key columns
  binnedCountsTable = fread(binnedCountsFilePath)
  setkey(binnedCountsTable, Chromosome, `Bin_Start-End`)

  # Create a bin start column from the bin range column.
  binnedCountsTable[,Bin_Start := vapply(strsplit(`Bin_Start-End`, '-'),
                                         function(x) as.integer(x[1]), FUN.VALUE = integer(1))]

  # Add a column for the majority color domain if the binned color domains file path was given.
  if (!is.na(binnedColorDomainsFilePath)) {
    binnedColorDomainsTable = fread(binnedColorDomainsFilePath)
    binnedCountsTable = binnedCountsTable[binnedColorDomainsTable]
    setnames(binnedCountsTable, "Domain_Color", "Majority_Domain_Color")
    setkey(binnedCountsTable, Chromosome, `Bin_Start-End`)
  }

  # Add in a complementary (background) data set (If we have the relevant file) and normalize by it.
  # Normalization occurs by taking the base-2 log of the ratio of the raw data to the background data for each bin.
  # Prior to this process, any bins with a zero in either data set are dropped.
  # (Otherwise, these rows cause problems with division or the logarithmic function.)
  # The scaling factor adjusts the normalized values prior to taking the log of the ratio.  A NULL value
  # computes the scaling factor from the given data to center the normalized values on 0.
  if (!is.na(backgroundBinCountsFilePath)) {

    # Read in the background counts and merge them with the raw counts.
    backgroundBinCountsTable = fread(backgroundBinCountsFilePath)
    setnames(backgroundBinCountsTable, "Feature_Counts", "Background_Feature_Counts")
    binnedCountsTable = binnedCountsTable[backgroundBinCountsTable]

    # If no scaling factor was given, compute it from the given data.
    if (is.null(scalingFactor)) {
      scalingFactor = sum(binnedCountsTable$Background_Feature_Counts)/sum(binnedCountsTable$Feature_Counts)
    }

    # Remove rows with 0 counts.
    binnedCountsTable = binnedCountsTable[Feature_Counts > 0 & Background_Feature_Counts > 0]

    # Normalize and scale!
    rawToBackgroundRatio = binnedCountsTable$Feature_Counts / binnedCountsTable$Background_Feature_Counts
    binnedCountsTable[, Scaled_Raw_To_Background_Ratio := rawToBackgroundRatio * scalingFactor]
    binnedCountsTable[, Log_Ratio := log(Scaled_Raw_To_Background_Ratio,2)]
  }

  return(binnedCountsTable)

}


# Given a path to a file with information on gene bins, return a binned counts data.table
parseGeneBinData = function(geneBinsCountsFilePath, backgroundFilePath = NA,
                            scalingFactor = NULL) {

  # Read in the data
  geneBinsCountsTable = fread(geneBinsCountsFilePath)
  if ("Color_Domain" %in% colnames(geneBinsCountsTable)) {
    setkey(geneBinsCountsTable, Color_Domain, Gene_Fraction)
  } else setkey(geneBinsCountsTable, Gene_Fraction)

  # Create columns for normalized counts
  totalCounts = sum(geneBinsCountsTable$Coding_Strand_Counts) +
                sum(geneBinsCountsTable$Noncoding_Strand_Counts)

  # Add in a complementary (background) data set (If we have the relevant file).
  if (!is.na(backgroundFilePath)) {
    backgroundCountsTable = fread(backgroundFilePath)
    setnames(backgroundCountsTable, c("Coding_Strand_Counts","Noncoding_Strand_Counts"),
             c("Background_Coding_Strand_Counts", "Background_Noncoding_Strand_Counts"))
    geneBinsCountsTable = geneBinsCountsTable[backgroundCountsTable]

    # If no scaling factor was given, compute it from the given data.
    if (is.null(scalingFactor)) {
      scalingFactor = sum(geneBinsCountsTable$Background_Coding_Strand_Counts) +
                      sum(geneBinsCountsTable$Background_Noncoding_Strand_Counts) /
                      ( sum(geneBinsCountsTable$Coding_Strand_Counts) +
                        sum(geneBinsCountsTable$Noncoding_Strand_Counts) )
    }

    # Remove rows with 0 counts.
    geneBinsCountsTable = geneBinsCountsTable[Coding_Strand_Counts > 0 & Noncoding_Strand_Counts > 0 &
                                              Background_Coding_Strand_Counts > 0 &
                                              Background_Noncoding_Strand_Counts > 0]

    # Normalize and scale!
    codingRawToBackgroundRatio = geneBinsCountsTable$Coding_Strand_Counts /
                                 geneBinsCountsTable$Background_Coding_Strand_Counts
    noncodingRawToBackgroundRatio = geneBinsCountsTable$Noncoding_Strand_Counts /
                                    geneBinsCountsTable$Background_Noncoding_Strand_Counts

    geneBinsCountsTable[, Scaled_Coding_Ratio := codingRawToBackgroundRatio * scalingFactor]
    geneBinsCountsTable[, Coding_Log_Ratio := log(Scaled_Coding_Ratio,2)]

    geneBinsCountsTable[, Scaled_Noncoding_Ratio := noncodingRawToBackgroundRatio * scalingFactor]
    geneBinsCountsTable[, Noncoding_Log_Ratio := log(Scaled_Noncoding_Ratio,2)]

    geneBinsCountsTable[, TS_Vs_NTS_Log_Ratio := Noncoding_Log_Ratio - Coding_Log_Ratio]
  }

  return(geneBinsCountsTable)

}


# A function for plotting gene fraction data (e.g. TS vs NTS repair rates across genes)
plotGeneBinData = function(geneBinsCountsTable, title = "", xAxisLabel = "Gene Fraction Bin",
                           yAxisLabel = "Log Ratio", ylim = NULL, yData1 = "Coding_Log_Ratio",
                           yData2 = "Noncoding_Log_Ratio", yData3 = "TS_Vs_NTS_Log_Ratio",
                           plotYData3Only = TRUE, flankingBinSize = 356) {

  if ("Color_Domain" %in% colnames(geneBinsCountsTable)) {
    geneBinPlot = ggplot(geneBinsCountsTable[Color_Domain != "GRAY"], aes(Gene_Fraction, color = Color_Domain))
  } else {
    geneBinPlot = ggplot(geneBinsCountsTable, aes(Gene_Fraction))
  }

  geneBinPlot = geneBinPlot +
    scale_color_manual(values = domainColors, guide = "none")

  if (plotYData3Only) {

    geneBinPlot = geneBinPlot +
      geom_line(aes_string(y = yData3), linewidth = 1.25) +
      geom_point(aes_string(y = yData3), size = 2)

  } else {

    geneBinPlot = geneBinPlot +
      geom_line(aes_string(y = yData1, linetype = shQuote("dashed")), linewidth = 1.25) +
      geom_point(aes_string(y = yData1), size = 2) +
      geom_line(aes_string(y = yData2, linetype = shQuote("solid")), linewidth = 1.25) +
      geom_point(aes_string(y = yData2), size = 2) +
      scale_linetype_identity(guide = "legend", name = "", breaks = c("solid", "dashed"),
                              labels = c("Transcribed Strand", "Non-Transcribed Strand"))

  }

  geneBinPlot = geneBinPlot +
    labs(title = title, x = xAxisLabel, y = yAxisLabel) +
    coord_cartesian(ylim = ylim) +
    scale_x_continuous(breaks = -2:8 + 0.5, minor_breaks = NULL,
                       labels = c(-2*flankingBinSize,'',"TSS",'','','','','',"TES",'',2*flankingBinSize)) +
    geom_vline(xintercept = 0.5, linetype = "dashed") + geom_vline(xintercept = 6.5, linetype = "dashed") +
    defaultTextScaling + blankBackground

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

    barplot(relevantCountsTable[, Normalized_Counts],
            names.arg = relevantCountsTable[, Bin_Start], col = "#ca0020", border = "#ca0020",
            xlab = "Chromosome Position", ylab = "Normalized Counts", main = paste(title, chromosome))

    # Then overlay the second data set!
    par(new = T)

    compRelevantCountsTable = compBinCountsTable[Chromosome == chromosome]

    barplot(compRelevantCountsTable[, Normalized_Counts], yaxt = 'n',
            names.arg = compRelevantCountsTable[, Bin_Start], col = "#0571b0", border = "#0571b0")

    # What if we keep track of where they intersect?
    par(new = T)
    barplot(pmin(relevantCountsTable$Normalized_Counts, compRelevantCountsTable$Normalized_Counts), yaxt = 'n',
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

  binnedCountsTable = copy(binnedCountsTable)
  binnedCountsTable[,Bin_Start := Bin_Start / 1000000]

  for (chromosomeSet in chromosomeSets) {

    plot = ggplot(binnedCountsTable[Chromosome %in% chromosomeSet],
                  aes_string("Bin_Start", yData)) +
      geom_bar(stat = "identity", color = "Black", fill = "Black") +
      labs(title = title, x = "Chromosome Position (Mb)", y = yAxisLabel) +
      facet_grid(~Chromosome, space = "free", scales = "free") +
      coord_cartesian(ylim = ylim) +
      defaultTextScaling + blankBackground +
      theme(panel.spacing = unit(0.1, "lines"), legend.position = "none")

    if (colorPlot) {
      plot = plot + geom_bar(aes(color = Majority_Domain_Color, fill = Majority_Domain_Color), stat = "identity") +
        scale_fill_manual(values = domainColors, guide = "none") +
        scale_color_manual(values = domainColors, guide = "none")
    }

    print(plot)

  }
}


# Plot the distribution of log ratio values for each individual color domain throughout the genome.
# Plotted as either a scatter plot or box plot based on given "pseudo constant" values.
SCATTER = 1
BOXPLOT = 2
plotLogRatioDistribution = function(binnedCountsTable, plotType, title = "", yData = "Log_Ratio",
                                    yAxisLabel = "Log Ratio", ylim = NULL, facetColumn = NULL) {

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
    scale_color_manual(values = domainColors, guide = "none")

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

  # Create the facet plot if requested.
  if (!is.null(facetColumn)) {
    thisPlot = thisPlot + facet_grid(facetColumn)
  }

  # Finishing touches
  thisPlot = thisPlot + labs(title = title, y = yAxisLabel) +
    coord_cartesian(ylim = ylim) +
    defaultTextScaling + blankBackground +
    theme(axis.title.x = element_blank())

  print(thisPlot)

}


# Graphs median log_ratio counts for different domains across different time points.
# Requires two parallel lists of binned counts tables and time points (as numeric values).
plotLogRatioMediansOverTime = function(binnedCountsTables, timePoints, title = "",
                                       xAxisLabel = "Timepoint", yAxisLabel = "Log Ratio", ylim = NULL,
                                       lineTypeColumn = NULL, lineTypeMapping = NULL) {

  masterMedianTable = rbindlist(mapply(getMedianTable, binnedCountsTables,
                                       timePoints, MoreArgs = list(lineTypeColumn), SIMPLIFY = FALSE))

  # Graph a line for each color.
  thisPlot = ggplot(masterMedianTable, aes(Time, Median, color = Majority_Domain_Color))

  # If requested, change the linetype based on the values in the given column
  if (!is.null(lineTypeColumn)) {

    thisPlot = thisPlot + aes_string(linetype = lineTypeColumn)

    if (!is.null(lineTypeMapping)) {
      thisPlot = thisPlot + scale_linetype_manual(values = lineTypeMapping, name = NULL)
    }

  }

  thisPlot = thisPlot +
    scale_color_manual(values = domainColors, guide = "none") +
    labs(title = title, x = xAxisLabel, y = yAxisLabel) +
    geom_line(linewidth = 1.25) + geom_point(size = 2) +
    coord_cartesian(ylim = ylim) +
    defaultTextScaling + blankBackground

  print(thisPlot)

}


# A function for producing a table of median log_ratio values with the domain color and
# time point (as a numeric value) included as columns.
getMedianTable = function(binnedCountsTable, timePoint, lineTypeColumn) {

  if (is.null(lineTypeColumn)) {
    return(binnedCountsTable[Majority_Domain_Color != "GRAY" & Majority_Domain_Color != "WHITE",
                             .(Median = median(Log_Ratio), Time = timePoint),
                             Majority_Domain_Color])
  } else {
    return(binnedCountsTable[Majority_Domain_Color != "GRAY" & Majority_Domain_Color != "WHITE",
                             .(Median = median(Log_Ratio), Time = timePoint),
                             by = setNames(list(get("Majority_Domain_Color"), get(lineTypeColumn)),
                                           c("Majority_Domain_Color",lineTypeColumn))])
  }

}
