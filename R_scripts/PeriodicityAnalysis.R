# A series of functions for quantifying and plotting periodicities.
# Trimmed down from mutperiod to accept simple counts data as input.
library(lomb)
library(ggplot2)

ROTATIONAL = 1
TRANSLATIONAL = 2
# Returns a lomb object for the given data.  (Usually nucleosome self-count data)
getLombResult = function(countsTable, rotOrTrans, nucleosomeExclusionBoundary = NA,
                         showLspWarnings = FALSE) {

  if (rotOrTrans == ROTATIONAL) {
    lombFrom = 5
    lombTo = 25
  } else if (rotOrTrans == TRANSLATIONAL) {
    lombFrom = 50
    lombTo = 250
  } else stop("Invalid \"rotOrTrans\" argument given.")

  if ("Normalized_Both_Strands" %in% colnames(countsTable)) {
    colnames(countsTable)[which(colnames(countsTable) == "Normalized_Both_Strands")] = "Both_Strands_Counts"
  }

  if (is.na(nucleosomeExclusionBoundary)) {
    counts = countsTable$Both_Strands_Counts
  } else {
    counts = countsTable[Dyad_Position > nucleosomeExclusionBoundary |
                           Dyad_Position < -nucleosomeExclusionBoundary, Both_Strands_Counts]
  }

  if (is.na(nucleosomeExclusionBoundary)) {
    times = countsTable$Dyad_Position
  } else {
    times = countsTable[Dyad_Position > nucleosomeExclusionBoundary |
                          Dyad_Position < -nucleosomeExclusionBoundary, Dyad_Position]
  }

  if (showLspWarnings) return(lsp(counts, times, lombFrom, lombTo, "period", 100, plot = FALSE))
  else suppressWarnings(return(lsp(counts, times, lombFrom, lombTo, "period", 100, plot = FALSE)))

}


plotLombResult = function(lombData, title = "", xAxisLabel = "Periods",
                          yAxisLabel = "Power", ylim = NULL,
                          colorLabels = NULL, lineTypeLabels = NULL) {

  lombResultPlot = ggplot(lombData, aes(x = Periods, y = Power))

  if ("Color_Domain" %in% colnames(lombData)) {
    lombResultPlot = lombResultPlot + aes(color = Color_Domain)
  }
  if ("Line_Type" %in% colnames(lombData)) {
    lombResultPlot = lombResultPlot + aes(linetype = Line_Type)
  }

  lombResultPlot = lombResultPlot +
    geom_line(size = 2)

  if (is.null(colorLabels)) {
    lombResultPlot = lombResultPlot + scale_color_identity()
  } else {
    lombResultPlot = lombResultPlot +
      scale_color_identity(name = '', guide = "legend", breaks = names(colorLabels), labels = colorLabels)
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
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 15), axis.text = element_text(size = 12),
          legend.text = element_text(size = 12), legend.key.width = unit(2,"cm"))

  print(lombResultPlot)

}


### Modified from plotting tests in run_MutperiodR (might backport later, idk.)

# From Cui and Zhurkin, 2010 with 1 added to each side and half-base positions included.
minorInPositions = c(4:9-74, 14:19-74, 25:30-74, 36:41-74, 46:51-74, 56:61-74, 66:71-74,
                     5:9-74.5, 15:19-74.5, 26:30-74.5, 37:41-74.5, 47:51-74.5, 57:61-74.5, 67:71-74.5)
minorOutPositions = c(9:14-74, 20:24-74, 31:35-74, 41:46-74 ,51:56-74, 61:66-74,
                      10:14-74.5, 21:24-74.5, 32:35-74.5, 42:46-74.5, 52:56-74.5, 62:66-74.5)

# Average across 11 base pairs centered on the given position.
smoothValues = function(middlePos, data, dataCol, averagingRadius = 5) {

  positionsToAverage = (middlePos-averagingRadius):(middlePos+averagingRadius)
  valuesToAverage = data[Dyad_Position %in% positionsToAverage][[dataCol]]
  return(mean(valuesToAverage))

}


plotPeriodicity = function(dataSet, rotationalOnlyCutoff = 60,
                           smoothTranslational = TRUE, fixedNRL = NULL,
                           dataCol = "Normalized_Both_Strands", title = "", ylim = NULL,
                           yAxisLabel = "Normalized Repair Reads",
                           xAxisLabel = "Position Relative to Dyad (bp)") {

  # If dataSet is a string, get the relevant counts and periodicity data for the given data set name.
  if (is.character(dataSet)) {
    if (dataSet %in% names(mutperiodData$normalizedNucleosomeCountsTables)) {
      countsData = mutperiodData$normalizedNucleosomeCountsTables[[dataSet]]
    } else if (dataSet %in% names(mutperiodData$rawNucleosomeCountsTables)) {
      countsData = mutperiodData$rawNucleosomeCountsTables[[dataSet]]
    } else stop("Unknown data set name.")
    periodicityData = as.list(mutperiodData$periodicityResults[Data_Set == dataSet])
  }
  # Otherwise, the data set that was passed in should just be the counts data.
  else {
    countsData = dataSet
    periodicityData = NULL
  }

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
    countsData[Dyad_Position %in% minorInPositions | -Dyad_Position %in% minorInPositions, Color := "#1bcc44"]
    countsData[Dyad_Position %in% minorOutPositions | -Dyad_Position %in% minorOutPositions, Color := "#993299"]
    countsData[is.na(Color), Color := "Black"]
  }

  if (rotationalPlus) {
    # Color linker DNA in linker+ plots.
    countsData[Dyad_Position <= -73 | Dyad_Position >= 73, Color := "Gold"]
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
    countsData[Dyad_Position %in% nucleosomePositions | -Dyad_Position %in% nucleosomePositions, Color := "#0571b0"]
    countsData[Dyad_Position %in% linkerPositions | -Dyad_Position %in% linkerPositions, Color := "#ca0020"]
  }

  # Plot it!
  periodicityPlot = ggplot(countsData, aes_string("Dyad_Position", dataCol, color = "Color")) +
                    geom_path(size = 1.25, aes(group = 1))
  if (rotationalPlus) {
    periodicityPlot = periodicityPlot +
      scale_color_identity(name = '', guide = "legend",
                           breaks = c("#1bcc44", "#993299", "Gold"),
                           labels = c("Minor-in", "Minor-out", "Linker"))
  } else if (rotational) {
    periodicityPlot = periodicityPlot +
      scale_color_identity(name = '', guide = "legend",
                           breaks = c("#1bcc44", "#993299"),
                           labels = c("Minor-in", "Minor-out"))
  } else if (translational) {
    periodicityPlot = periodicityPlot +
      scale_color_identity(name = '', guide = "legend",
                           breaks = c("#0571b0", "#ca0020"),
                           labels = c("Nucleosome","Linker"))
  }
  periodicityPlot = periodicityPlot +
    coord_cartesian(ylim = ylim) +
    labs(title = title, x = xAxisLabel, y = yAxisLabel) +
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 15), axis.text = element_text(size = 12),
          legend.text = element_text(size = 12))

  print(periodicityPlot)

}


# Add information on timepoint and domain to a given data table.  If none of the expected timepoints or
# none of the expected domains are present, return an empty data.table.
# Also, smooths translational data.
addTimepointAndDomainInfo = function(dataSetName, dataCol = "Normalized_Both_Strands",
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
    countsData[, (dataCol) := sapply(countsData$Dyad_Position, smoothValues, data = countsData, dataCol = dataCol)]
  }

  return(countsData[, c("Timepoint", "Domain") := list(rep(timepoint, .N), rep(domain, .N))])

}

# Plot a bunch of figures together using facets, stratified by timepoint on one axis and domains on the other.
plotBulkCountsData = function(bulkCountsData, dataCol = "Normalized_Both_Strands",
                              expectedTimepoints = c("10m", "30m", "8h", "16h", "24h"),
                              title = "", xBreaks = NULL, ylim = NULL,
                              yAxisLabel = "Normalized Repair Reads",
                              xAxisLabel = "Position Relative to Dyad (bp)") {
  bulkCountsPlot = ggplot(bulkCountsData, aes_string("Dyad_Position", dataCol, color = "Domain")) +
    scale_color_manual(values = c("BLACK" = "black", "BLUE" = "blue", "GREEN" = "forestgreen",
                                  "RED" = "red", "YELLOW" = "gold2"), guide = "none") +
    geom_line() +
    labs(title = title, x = "Position Relative to Dyad (bp)", y = yAxisLabel) +
    facet_grid(factor(Timepoint, levels = expectedTimepoints)~Domain) +
    coord_cartesian(ylim = ylim)
  if (is.null(xBreaks)) {
    bulkCountsPlot = bulkCountsPlot +
      scale_x_continuous(breaks = c(round(min(bulkCountsData$Dyad_Position)/2),0,
                                    round(max(bulkCountsData$Dyad_Position)/2)))
  } else {
    bulkCountsPlot = bulkCountsPlot + scale_x_continuous(breaks = xBreaks)
  }
  bulkCountsPlot = bulkCountsPlot +
    scale_y_continuous(n.breaks = 3) +
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 15), axis.text = element_text(size = 12),
          strip.text = element_text(size = 15))

  print(bulkCountsPlot)

  }
