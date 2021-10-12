# A series of functions for quantifying and plotting periodicities.
# Trimmed down from mutperiod to accept simple counts data as input.
library(lomb)
library(ggplot2)

ROTATIONAL = 1
TRANSLATIONAL = 2
# Returns a lomb object for the given data.  (Usually nucleosome self-count data)
getLombResult = function(countsTable, rotOrTrans, nucleosomeExclusionBoundary = NA) {

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

  return(lsp(counts, times, lombFrom, lombTo, "period", 100, plot = FALSE))

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
