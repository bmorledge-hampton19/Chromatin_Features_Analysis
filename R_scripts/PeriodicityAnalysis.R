# A series of functions for quantifying and plotting periodicities.
# Trimmed down from mutperiod to accept simple counts data as input.
library(lomb)
library(data.table)
library(ggplot2)

ROTATIONAL = 1
TRANSLATIONAL = 2

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

# Returns a lomb object for the given data.  (Usually nucleosome self-count data)
getLombResult = function(countsTable, rotOrTrans, nucleosomeExclusionBoundary = NA, rotationalPosCutoff = 60,
                        countsColumnName = NULL, plotResult = FALSE, showLspWarnings = FALSE) {

  if (rotOrTrans == ROTATIONAL) {
    lombFrom = 5
    lombTo = 25
  } else if (rotOrTrans == TRANSLATIONAL) {
    lombFrom = 50
    lombTo = 250
  } else stop("Invalid \"rotOrTrans\" argument given.")

  # If no column name for counts was given, try to satisfy one of three default conditions.
  if (is.null(countsColumnName)) {
    if ("Normalized_Both_Strands" %in% colnames(countsTable)) {
      colnames(countsTable)[which(colnames(countsTable) == "Normalized_Both_Strands")] = "Counts"
    } else if ("Both_Strands_Counts" %in% colnames(countsTable)) {
      colnames(countsTable)[which(colnames(countsTable) == "Both_Strands_Counts")] = "Counts"
    } else if (length(colnames(countsTable) == 2 && colnames(CountsTable)[1] == "Dyad_Position")) {
      colnames(countsTable)[2] = "Counts"
    } else stop("No counts column name given and no default conditions satisfied.")

  # If a column name for counts was given, make sure that it exists and then rename it.
  } else {
    if (countsColumnName %in% colnames(countsTable)) {
      colnames(countsTable)[which(colnames(countsTable) == countsColumnName)] = "Counts"
    } else stop("Given counts column name is not present in the column names for the given counts table.")
  }

  if (rotOrTrans == ROTATIONAL) {
    if (!is.na(nucleosomeExclusionBoundary)) stop("nucleosome exclusion boundary should be NA for rotational data.")
    counts = countsTable[Dyad_Position <= rotationalPosCutoff &
                           Dyad_Position >= -rotationalPosCutoff, Counts]
    times = countsTable[Dyad_Position <= rotationalPosCutoff &
                          Dyad_Position >= -rotationalPosCutoff, Dyad_Position]
  } else if (is.na(nucleosomeExclusionBoundary)) {
    counts = countsTable$Counts
    times = countsTable$Dyad_Position
  } else {
    counts = countsTable[Dyad_Position > nucleosomeExclusionBoundary |
                           Dyad_Position < -nucleosomeExclusionBoundary, Counts]
    times = countsTable[Dyad_Position > nucleosomeExclusionBoundary |
                          Dyad_Position < -nucleosomeExclusionBoundary, Dyad_Position]
  }

  if (showLspWarnings) return(lsp(counts, times, lombFrom, lombTo, "period", 100, plot = plotResult))
  else suppressWarnings(return(lsp(counts, times, lombFrom, lombTo, "period", 100, plot = plotResult)))

}


# A helper function for when I want to quickly get lomb results with SNR.
getPeakPeriodicityAndSNR = function(lombResult) {

  # Get the peak periodicity and its associated SNR
  peakPeriodicity = lombResult$peak.at[1]
  noiseBooleanVector = (lombResult$scanned < peakPeriodicity - 0.5
                        | lombResult$scanned > peakPeriodicity + 0.5)
  SNR = lombResult$peak / median(lombResult$power[noiseBooleanVector])

  return(c(peakPeriodicity,SNR))

}


plotLombResult = function(lombData, title = "", xAxisLabel = "Period (bp)",
                          yAxisLabel = "Power", ylim = NULL,
                          colorLabels = NULL, lineTypeLabels = NULL,
                          manualColorOverride = NULL) {

  lombResultPlot = ggplot(lombData, aes(x = Periods, y = Power))

  if ("Color_Domain" %in% colnames(lombData)) {
    lombResultPlot = lombResultPlot + aes(color = Color_Domain)
  }
  if ("Line_Type" %in% colnames(lombData)) {
    lombResultPlot = lombResultPlot + aes(linetype = Line_Type)
  }

  lombResultPlot = lombResultPlot +
    geom_line(linewidth = 2)

  if (is.null(manualColorOverride)) colors = domainColors
  else colors = manualColorOverride

  if (is.null(colorLabels)) {
    lombResultPlot = lombResultPlot + scale_color_manual(values = colors, guide = "none")
  } else {
    lombResultPlot = lombResultPlot +
      scale_color_manual(values = colors, name = '', guide = "legend",
                         breaks = names(colorLabels), labels = colorLabels)
  }

  if (is.null(lineTypeLabels)) {
    lombResultPlot = lombResultPlot + scale_linetype_identity()
  } else {
    lombResultPlot = lombResultPlot +
      scale_linetype_identity(name = '', guide = "legend", breaks = names(lineTypeLabels), labels = lineTypeLabels)
  }

  lombResultPlot = lombResultPlot +
    coord_cartesian(ylim = ylim) +
    labs(title = title, x = xAxisLabel, y = yAxisLabel) +
    defaultTextScaling + blankBackground +
    theme(legend.key.width = unit(2,"cm"))

  print(lombResultPlot)

}


### Modified from plotting tests in run_MutperiodR (might backport later, idk.)

# Defines both individual positions (for coloring) and regions (for grouping with geom_path),
# so this code looks like a hot mess.

# From Cui and Zhurkin, 2010 with 1 added to each side and half-base positions included.
minorInRegions = list(c(4:9-74,5:9-74.5), c(14:19-74,15:19-74.5), c(25:30-74,26:30-74.5), c(36:41-74,37:41-74.5),
                      c(46:51-74,47:51-74.5), c(56:61-74,57:61-74.5), c(66:71-74,67:71-74.5))
minorInRegions = append(minorInRegions, lapply(rev(minorInRegions), function(x) -rev(x)))
names(minorInRegions) = lapply(seq_along(minorInRegions), function(i) paste0("Minor_In_Region_",i))
minorInPositions = unlist(minorInRegions, use.names = FALSE)

minorOutRegions = list(c(9:14-74,10:14-74.5), c(20:24-74,21:24-74.5), c(31:35-74,32:35-74.5),
                       c(41:46-74,42:46-74.5), c(51:56-74,52:56-74.5), c(61:66-74,62:66-74.5))
minorOutRegions = append(minorOutRegions, lapply(rev(minorOutRegions), function(x) -rev(x)))
names(minorOutRegions) = lapply(seq_along(minorOutRegions), function(i) paste0("Minor_Out_Region_",i))
minorOutPositions = unlist(minorOutRegions, use.names = FALSE)

intermediatePositions =
  seq(-73,73,0.5)[!(seq(-73,73,0.5) %in% minorInPositions | seq(-73,73,0.5) %in% minorOutPositions)]

intermediateRegions = list()
currentIntermediateRegion = c()
lastIntermediatePosition = NULL

for (intermediatePosition in intermediatePositions) {

  if (!is.null(lastIntermediatePosition) && intermediatePosition - lastIntermediatePosition != 0.5) {
    intermediateRegions[[length(intermediateRegions) + 1]] = currentIntermediateRegion
    currentIntermediateRegion = c()
  }

  currentIntermediateRegion = append(currentIntermediateRegion,intermediatePosition)
  lastIntermediatePosition = intermediatePosition

}
intermediateRegions = append(intermediateRegions, currentIntermediateRegion)

names(intermediateRegions) = lapply(seq_along(intermediateRegions), function(i) paste0("Intermediate_Region_",i))

linkerRegions = list(Linker_Region_1 = seq(-999,-73,0.5), linker_Region_2 = seq(73,999,0.5))

# Make a data.table with dyad positions and group ID's
rotationalRegionIDs = rbindlist(lapply(
  list(minorInRegions, minorOutRegions, intermediateRegions, linkerRegions), FUN = function(x) {
    rbindlist(
      mapply(function(pos,id) data.table(Dyad_Position = pos, Group_ID = id), x, names(x), SIMPLIFY = FALSE)
    )
  }
))
setkey(rotationalRegionIDs, Dyad_Position)


# Average across 11 base pairs centered on the given position.
smoothValues = function(middlePos, data, dataCol, averagingRadius = 5) {

  positionsToAverage = (middlePos-averagingRadius):(middlePos+averagingRadius)
  valuesToAverage = data[Dyad_Position %in% positionsToAverage][[dataCol]]
  return(mean(valuesToAverage))

}


parsePeriodicityData = function(dataSet, rotationalOnlyCutoff = 60, dataCol = "Normalized_Both_Strands",
                                smoothTranslational = TRUE, fixedNRL = NULL) {

  # If dataSet is a string, get the relevant counts and periodicity data for the given data set name.
  if (is.character(dataSet)) {
    if (dataSet %in% names(mutperiodData$normalizedNucleosomeCountsTables)) {
      countsData = mutperiodData$normalizedNucleosomeCountsTables[[dataSet]]
    } else if (dataSet %in% names(mutperiodData$rawNucleosomeCountsTables)) {
      countsData = mutperiodData$rawNucleosomeCountsTables[[dataSet]]
    } else stop("Unknown data set name.")
    periodicityData = mutperiodData$periodicityResults[Data_Set == dataSet]
    if (dim(periodicityData)[1] > 1) stop("Multiple periodicity data sets found for given name.")
  }
  # Otherwise, the data set that was passed in should just be the counts data.
  else {
    countsData = dataSet
    periodicityData = NULL
  }
  setkey(countsData, Dyad_Position)

  # Determine whether the data is rotational, rotational+linker, or translational.
  rotational = FALSE
  rotationalPlus = FALSE
  translational = FALSE
  if (min(countsData$Dyad_Position) >= -73) {
    rotational = TRUE
  } else if (min(countsData$Dyad_Position) > -999) {
    rotational = TRUE
    rotationalPlus = TRUE
  } else translational = TRUE

  # If only rotational, trim to the cutoff value
  if (rotational && !rotationalPlus) {
    countsData = countsData[Dyad_Position >= -rotationalOnlyCutoff & Dyad_Position <= rotationalOnlyCutoff]
  }

  # Smooth if translational and requested
  if (translational && smoothTranslational) {
    countsData = copy(countsData)
    countsData[, (dataCol) := sapply(countsData$Dyad_Position, smoothValues, data = countsData, dataCol = dataCol)]
  }

  if (rotational) {
    # Color rotational positioning
    countsData[Dyad_Position %in% minorInPositions,
               Periodic_Position_Color := "cyan3"]
               #Periodic_Position_Color := "#1bcc44"]
    countsData[Dyad_Position %in% minorOutPositions,
               Periodic_Position_Color := "darkorange"]
               #Periodic_Position_Color := "#993299"]
    countsData[is.na(Periodic_Position_Color), Periodic_Position_Color := "Black"]

    # Add in group IDs
    countsData = countsData[rotationalRegionIDs, nomatch = NULL]

  }

  if (rotationalPlus) {
    # Color linker DNA in linker+ plots.
    countsData[Dyad_Position <= -73 | Dyad_Position >= 73, Periodic_Position_Color := "deeppink2"]
    #countsData[Dyad_Position <= -73 | Dyad_Position >= 73, Periodic_Position_Color := "gold"]

  }

  if (translational) {

    # Derive linker and nucleosome positions from the expected period of the data.
    if (is.null(fixedNRL)) {
      if (is.null(periodicityData)) stop("No NRL given and no periodicity data available for this data set.")
      nucRepLen = round(periodicityData$Expected_Peak_Periodicity)
    } else {
      nucRepLen = fixedNRL
    }

    nucleosomePositions = sapply(1:10, function(x) return(append((-73+x*nucRepLen):(73+x*nucRepLen),
                                                                 (-72.5+x*nucRepLen):(72.5+x*nucRepLen))))
    nucleosomePositions = append(nucleosomePositions, c(0:73, 0.5:72.5))
    linkerPositions = sapply(0:8, function(x) return( append((74+x*nucRepLen):(-74+(x+1)*nucRepLen),
                                                             (73.5+x*nucRepLen):(-73.5+(x+1)*nucRepLen))))

    # Color translational positioning
    countsData[Dyad_Position %in% nucleosomePositions | -Dyad_Position %in% nucleosomePositions,
               Periodic_Position_Color := "darkblue"]
               #Periodic_Position_Color := "#0571b0"]
    countsData[Dyad_Position %in% linkerPositions | -Dyad_Position %in% linkerPositions,
               Periodic_Position_Color := "deeppink2"]
               #Periodic_Position_Color := "#ca0020"]

    # Group each individual region together

  }

  return(countsData)

}


# Plot the given periodicity data.
# If singleDataSet is true, color the data based on the "color" variable to highlight
# features of the (expected) periodicity.  Otherwise, assume that multiple periodicity
# datasets are present and color them based on a required "Color_Domain" column.
plotPeriodicity = function(countsData, singleDataSet = TRUE,
                           dataCol = "Normalized_Both_Strands", title = "", ylim = NULL,
                           yAxisLabel = "Repair/Damage",
                           xAxisLabel = "Position Relative to Dyad (bp)") {

  if (singleDataSet) {

    periodicityPlot = ggplot(countsData, aes_string("Dyad_Position", dataCol,
                                                         color = "Periodic_Position_Color")) +
      geom_path(linewidth = 1.25, lineend = "round", aes(group = 1))

    # Determine whether the data is rotational, rotational+linker, or translational.
    rotational = FALSE
    rotationalPlus = FALSE
    translational = FALSE
    if (min(countsData$Dyad_Position) >= -73) {
      rotational = TRUE
    } else if (min(countsData$Dyad_Position) > -999) {
      rotational = TRUE
      rotationalPlus = TRUE
    } else translational = TRUE

    if (rotationalPlus) {
      periodicityPlot = periodicityPlot +
        scale_color_identity(name = '', guide = "legend",
                             breaks = c("cyan3", "darkorange", "deeppink2"),
                             #breaks = c("#1bcc44", "#993299", "Gold"),
                             labels = c("Minor-in", "Minor-out", "Linker"))
    } else if (rotational) {
      periodicityPlot = periodicityPlot +
        scale_color_identity(name = '', guide = "legend",
                             breaks = c("cyan3", "darkorange"),
                             #breaks = c("#1bcc44", "#993299"),
                             labels = c("Minor-in", "Minor-out"))
    } else if (translational) {
      periodicityPlot = periodicityPlot +
        scale_color_identity(name = '', guide = "legend",
                             breaks = c("darkblue", "deeppink2"),
                             #breaks = c("#0571b0", "#ca0020"),
                             labels = c("Nucleosome","Linker"))
    }

  } else {

    periodicityPlot = ggplot(countsData, aes_string("Dyad_Position", dataCol,
                                                         color = "Color_Domain")) +
      geom_line(linewidth = 1.25) +
      scale_color_manual(values = domainColors, guide = "none")

  }

  periodicityPlot = periodicityPlot +
    coord_cartesian(ylim = ylim) +
    labs(title = title, x = xAxisLabel, y = yAxisLabel) +
    defaultTextScaling + blankBackground

  print(periodicityPlot)

}


# Combines the parsePeriodicityData and plotPeriodicity functions into one convenience function
parseAndPlotPeriodicity = function(dataSet, rotationalOnlyCutoff = 60, dataCol = "Normalized_Both_Strands",
                                   smoothTranslational = TRUE, fixedNRL = NULL,
                                   title = "", ylim = NULL,
                                   yAxisLabel = "Repair/Damage",
                                   xAxisLabel = "Position Relative to Dyad (bp)") {

  parsedData = parsePeriodicityData(dataSet, rotationalOnlyCutoff, dataCol, smoothTranslational, fixedNRL)
  plotPeriodicity(parsedData, TRUE, dataCol, title, ylim, yAxisLabel, xAxisLabel)

}


# Plot the plus and minus strands (aligned) as two lines on the same graph.
plotPlusAndMinus = function(countsData, title = "", ylim = NULL,
                            yAxisLabel = "Repair/Damage",
                            xAxisLabel = "Position Relative to Dyad (bp)") {

  plusStrandCounts = countsData[[which(grepl("Plus_Strand", colnames(countsData)))[1]]]
  minusStrandCounts = countsData[[which(grepl("Minus_Strand", colnames(countsData)))[1]]]

  strandCountsData = data.table(Dyad_Position = countsData[,Dyad_Position],
                                Plus_Strand_Counts = plusStrandCounts,
                                Minus_Strand_Counts = minusStrandCounts)

  ggplot(strandCountsData, aes(x = Dyad_Position)) +
    geom_line(aes(y = Plus_Strand_Counts, color = "forestgreen"), linewidth = 1.25) +
    geom_line(aes(y = Minus_Strand_Counts, color = "red"), linewidth = 1.25) +
    scale_color_identity(name = '', guide = "legend",
                         breaks = c("forestgreen", "red"),
                         labels = c("Plus Strand", "Minus Strand")) +
    coord_cartesian(ylim = ylim) +
    labs(title = title, x = xAxisLabel, y = yAxisLabel) +
    defaultTextScaling + blankBackground

}


# Add information on timepoint and domain to a given data table.  If none of the expected timepoints or
# none of the expected domains are present, return an empty data.table.
# Also, smooths translational data.
addTimepointAndDomainInfo = function(dataSetName, smoothCols = list("Normalized_Both_Strands"),
                                     expectedTimepoints = c("10m", "30m", "8h", "16h", "24h"),
                                     expectedDomains = c("BLACK", "BLUE", "GREEN", "RED", "YELLOW")) {

  timepoint = names(which(sapply(expectedTimepoints, function(x) grepl(x,dataSetName))))
  domain = names(which(sapply(expectedDomains, function(x) grepl(x,dataSetName))))

  if (length(timepoint) == 0 || length(domain) == 0) return(data.table())
  if (length(timepoint) > 1) stop(paste("Multiple timepoints found in "),dataSetName)
  if (length(domain) > 1) stop(paste("Multiple domains found in "),dataSetName)

  if (dataSetName %in% names(mutperiodData$normalizedNucleosomeCountsTables)) {
    countsData = mutperiodData$normalizedNucleosomeCountsTables[[dataSetName]]
  } else if (dataSetName %in% names(mutperiodData$rawNucleosomeCountsTables)) {
    countsData = mutperiodData$rawNucleosomeCountsTables[[dataSetName]]
  } else stop("Unknown data set name.")

  if (grepl("nuc-group", dataSetName, fixed = TRUE)) {
    countsData = copy(countsData)
    captureOutput = lapply(smoothCols, function(dataCol)
      countsData[, (dataCol) := sapply(countsData$Dyad_Position, smoothValues, data = countsData, dataCol = dataCol)]
    )
  }

  return(countsData[, c("Timepoint", "Domain") := list(rep(timepoint, .N), rep(domain, .N))])

}

# Plot a bunch of figures together using facets, stratified by timepoint on one axis and domains on the other.
plotBulkCountsData = function(bulkCountsData, dataCol = "Normalized_Both_Strands",
                              expectedTimepoints = c("10m", "30m", "8h", "16h", "24h"),
                              title = "", xBreaks = NULL, yBreaks = NULL, ylim = NULL,
                              yAxisLabel = "Repair/Damage",
                              xAxisLabel = "Position Relative to Dyad (bp)") {
  bulkCountsPlot = ggplot(bulkCountsData, aes_string("Dyad_Position", dataCol, color = "Domain")) +
    scale_color_manual(values = domainColors, guide = "none") +
    geom_line() +
    labs(title = title, x = "Position Relative to Dyad (bp)", y = yAxisLabel) +
    facet_grid(factor(Timepoint, levels = expectedTimepoints)~Domain) +
    coord_cartesian(ylim = ylim) + blankBackground +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          strip.background = element_rect(color = "black", linewidth = 1))

  if (is.null(xBreaks)) {
    bulkCountsPlot = bulkCountsPlot +
      scale_x_continuous(breaks = c(round(min(bulkCountsData$Dyad_Position)/2),0,
                                    round(max(bulkCountsData$Dyad_Position)/2)))
  } else {
    bulkCountsPlot = bulkCountsPlot + scale_x_continuous(breaks = xBreaks)
  }

  if (is.null(yBreaks)) {
    bulkCountsPlot = bulkCountsPlot +
      scale_y_continuous(n.breaks = 3)
  } else {
    bulkCountsPlot = bulkCountsPlot + scale_y_continuous(breaks = yBreaks)
  }

  bulkCountsPlot = bulkCountsPlot +
    defaultTextScaling + theme(axis.text.y = element_text(size = 16))

  print(bulkCountsPlot)

}
