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
                          yAxisLabel = "Power", ylim = NULL) {

  if ("Color_Domain" %in% colnames(lombData)) {
    lombResultPlot = ggplot(lombData, aes(x = Periods, y = Power, color = Color_Domain))
  } else {
    lombResultPlot = ggplot(lombData, aes(x = Periods, y = Power))
  }

  lombResultPlot = lombResultPlot +
    geom_line(size = 2) +
    scale_color_identity() +
    coord_cartesian(ylim = ylim) +
    labs(title = title, x = xAxisLabel, y = yAxisLabel) +
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 15), axis.text = element_text(size = 12),
          legend.text = element_text(size = 12),)

  print(lombResultPlot)

}
