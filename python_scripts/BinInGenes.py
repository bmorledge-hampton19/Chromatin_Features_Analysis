# This script takes a gene designations file and one or more files of features (like mutations) to bin in fractions of that gene.
# It also includes a function for plotting the results.
import os, pandas, math
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.CountThisInThat.Counter import ThisInThatCounter
from benbiohelpers.CountThisInThat.OutputDataStratifiers import AmbiguityHandling
from benbiohelpers.Plotting.PlotnineHelpers import *
from plotnine import *


def binInGenes(featureFilePaths: List[str], geneDesignationsFilePath, flankingBinSize = 0, flankingBinNum = 0, 
               filePathSuffix = "", colorColIndex = None, geneFractionNum = 6):
    """
    Count features (e.g., mutations) on the transcribed and nontranscribed strands of genes and bin them across 6 gene fractions.
    The flankingBinSize and flankingBinNum parameters add additional bins of a constant length on the regions flanking genes (on each side). Importantly,
    these regions must already be a part of the regions given in the gene designations file.
    """

    class BinInGenesCounter(ThisInThatCounter):

        def setUpOutputDataHandler(self):
            super().setUpOutputDataHandler()
            if colorColIndex is not None:
                self.outputDataHandler.addSimpleEncompassingColStratifier(outputName = "Color_Domain", colIndex = colorColIndex)
            self.outputDataHandler.addFeatureFractionStratifier(outputName = "Gene_Fraction", fractionNum = geneFractionNum,
                                                                flankingBinSize = flankingBinSize, flankingBinNum = flankingBinNum)
            self.outputDataHandler.addStrandComparisonStratifier(strandAmbiguityHandling = AmbiguityHandling.tolerate)
            if colorColIndex is None: customStratifyingNames=(None, {True:"Coding_Strand_Counts", False:"Noncoding_Strand_Counts"})
            else: customStratifyingNames=(None, None, {True:"Coding_Strand_Counts", False:"Noncoding_Strand_Counts"})
            self.outputDataHandler.createOutputDataWriter(self.outputFilePath, customStratifyingNames = customStratifyingNames)

    for featureFilePath in featureFilePaths:

        print("\nWorking in", os.path.basename(featureFilePath))

        outputFilePath = featureFilePath.rsplit('.', 1)[0] + "_gene_bins"
        if filePathSuffix: outputFilePath += '_' + filePathSuffix
        outputFilePath += ".tsv"
        metadataFilePath = outputFilePath.rsplit('.',1)[0] + ".metadata"

        counter = BinInGenesCounter(featureFilePath, geneDesignationsFilePath, outputFilePath)
        counter.count()

        # Write metadata to preserve information that is not immediately apparent from the output.
        with open(metadataFilePath, 'w') as metadataFile:
            metadataFile.write(f"Flanking_Bin_Size:\t{flankingBinSize}\n")
            metadataFile.write(f"Flanking_Bins_Each_Side:\t{flankingBinNum}\n")
            metadataFile.write(f"Gene_Designations_File_Path:\t{geneDesignationsFilePath}\n")


domainColors = {"BLACK":"black", "black":"black", "BLUE":"blue", "blue":"blue",
                 "GREEN":"forestgreen", "green":"forestgreen",
                 "RED":"red", "red":"red", "YELLOW":"gold2", "yellow":"gold2",
                 "GRAY":"gray", "gray":"gray"}

# Given a path to a file with information on gene bins, return a binned counts data.table
def parseGeneBinData(geneBinsCountsFilePath, backgroundFilePath = None, scalingFactor = None):

    # Read in the data
    geneBinsCountsTable = pandas.read_table(geneBinsCountsFilePath)

    # Create columns for normalized counts
    totalCounts = sum(geneBinsCountsTable["Coding_Strand_Counts"]) + sum(geneBinsCountsTable["Noncoding_Strand_Counts"])

    # Add in a complementary (background) data set if it was given.
    if backgroundFilePath is not None:
        backgroundCountsTable = pandas.read_table(backgroundFilePath)
        backgroundCountsTable = backgroundCountsTable.rename(columns = {"Coding_Strand_Counts":"Background_Coding_Strand_Counts",
                                                                        "Noncoding_Strand_Counts":"Background_Noncoding_Strand_Counts"})
        geneBinsCountsTable = geneBinsCountsTable.merge(backgroundCountsTable)

        # If no scaling factor was given, compute it from the given data.
        if scalingFactor is None:
            scalingFactor = ((sum(geneBinsCountsTable["Background_Coding_Strand_Counts"]) +
                              sum(geneBinsCountsTable["Background_Noncoding_Strand_Counts"])) /
                             (sum(geneBinsCountsTable["Coding_Strand_Counts"]) +
                              sum(geneBinsCountsTable["Noncoding_Strand_Counts"])))

        # Remove rows with 0 counts.
        geneBinsCountsTable = geneBinsCountsTable.loc[(geneBinsCountsTable["Coding_Strand_Counts"] > 0) &
                                                      (geneBinsCountsTable["Noncoding_Strand_Counts"] > 0) &
                                                      (geneBinsCountsTable["Background_Coding_Strand_Counts"] > 0) &
                                                      (geneBinsCountsTable["Background_Noncoding_Strand_Counts"] > 0)].copy()

        # Normalize, scale, and compute log ratios.
        codingRawToBackgroundRatio = geneBinsCountsTable.Coding_Strand_Counts / geneBinsCountsTable.Background_Coding_Strand_Counts
        noncodingRawToBackgroundRatio = geneBinsCountsTable.Noncoding_Strand_Counts / geneBinsCountsTable.Background_Noncoding_Strand_Counts

        geneBinsCountsTable["Scaled_Coding_Ratio"] = codingRawToBackgroundRatio * scalingFactor
        geneBinsCountsTable["Coding_Log_Ratio"] = geneBinsCountsTable.Scaled_Coding_Ratio.apply(lambda x:math.log(x,2))

        geneBinsCountsTable["Scaled_Noncoding_Ratio"] = noncodingRawToBackgroundRatio * scalingFactor
        geneBinsCountsTable["Noncoding_Log_Ratio"] = geneBinsCountsTable.Scaled_Noncoding_Ratio.apply(lambda x:math.log(x,2))

        geneBinsCountsTable["TS_Vs_NTS_Log_Ratio"] = geneBinsCountsTable.Noncoding_Log_Ratio - geneBinsCountsTable.Coding_Log_Ratio

    else:

        # Remove rows with 0 counts.
        geneBinsCountsTable = geneBinsCountsTable.loc[(geneBinsCountsTable["Coding_Strand_Counts"] > 0) &
                                                      (geneBinsCountsTable["Noncoding_Strand_Counts"] > 0)].copy()

        # Compute log ratio.
        geneBinsCountsTable["TS_Vs_NTS_Log_Ratio"] = (geneBinsCountsTable.Noncoding_Strand_Counts / geneBinsCountsTable.Coding_Strand_Counts).apply(lambda x:math.log(x,2))

    return(geneBinsCountsTable)


def plotGeneBins(geneBinsCountsTable: pandas.DataFrame, title = "", xAxisLabel = "Gene Fraction Bin", yAxisLabel = "Log Ratio", ylim = None, xlim = None,
                 yData1 = "Coding_Log_Ratio", yData2 = "Noncoding_Log_Ratio", yData3 = "TS_Vs_NTS_Log_Ratio",
                 plotYData3Only = True, flankingBinSize = None, flankingBinNum = 0, geneFractionNum = 6):

    if ("Color_Domain" in geneBinsCountsTable.columns):
        geneBinPlot = ggplot(geneBinsCountsTable.loc[geneBinsCountsTable.Color_Domain != "GRAY"], aes("Gene_Fraction", color = "Color_Domain"))
        geneBinPlot = geneBinPlot + scale_color_manual(values = domainColors)
    else:
        geneBinPlot = ggplot(geneBinsCountsTable, aes("Gene_Fraction"))

    if plotYData3Only:

        geneBinPlot = (
            geneBinPlot +
            geom_line(aes(y = yData3), size = 1.25) +
            geom_point(aes(y = yData3), size = 2)
        )

    else:

        geneBinPlot = (geneBinPlot +
            geom_line(aes(y = yData1, linetype = '"dashed"'), size = 1.25) +
            geom_point(aes(y = yData1), size = 2) +
            geom_line(aes(y = yData2, linetype = '"solid"'), size = 1.25) +
            geom_point(aes(y = yData2), size = 2) +
            scale_linetype_identity(guide = "legend", name = " ", breaks = ("dashed", "solid"),
                                    labels = ("Non-Transcribed Strand", "Transcribed Strand"))
        )

    if flankingBinNum > 0:
        geneBinPlot = geneBinPlot + scale_x_continuous(breaks = [i+0.5 for i in range(-flankingBinNum, geneFractionNum + flankingBinNum + 1)], minor_breaks = None,
                                                       labels = [-flankingBinNum*flankingBinSize] + ['']*(flankingBinNum-1) +
                                                                ["TSS"] + ['']*(geneFractionNum-1) + ["TES"] + 
                                                                ['']*(flankingBinNum-1) + [flankingBinNum*flankingBinSize])
    else:
        geneBinPlot = geneBinPlot + scale_x_continuous(breaks = [i+0.5 for i in range(geneFractionNum + 1)], minor_breaks = None,
                                                       labels = ["TSS"] + ['']*(geneFractionNum-1) + ["TES"])

    if xlim is None: xlim = (min(geneBinsCountsTable.Gene_Fraction), max(geneBinsCountsTable.Gene_Fraction)+0.5)

    geneBinPlot = (geneBinPlot +
        labs(title = title, x = xAxisLabel, y = yAxisLabel) +
        coord_cartesian(ylim = ylim, xlim = xlim) +
        geom_vline(xintercept = 0.5, linetype = "dashed") + geom_vline(xintercept = 0.5+geneFractionNum, linetype = "dashed") +
        defaultTextScaling + blankBackground + theme(figure_size = (12,6))
    )

    return geneBinPlot


def main():

    try:
        from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
        workingDirectory = getDataDirectory()
    except ImportError:
        workingDirectory = os.path.dirname(__file__)

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=workingDirectory)
    dialog.createMultipleFileSelector("Feature Files (e.g. mutations):",0,"context_mutations.bed",("Bed Files",".bed"))    
    dialog.createFileSelector("Gene Designations:",1,("Bed Files",".bed"))
    dialog.createCheckbox("Color Domain is present in 7th (index=6) column",2, 0)

    flankDialog = dialog.createDynamicSelector(3, 0)
    flankDialog.initCheckboxController("Gene designations include flanking regions")
    flankSizeDialog = flankDialog.initDisplay(True, "FlankSize")
    flankSizeDialog.createTextField("Flanking bin size:", 0, 0, defaultText="0")
    flankSizeDialog.createTextField("Flanking bin number:", 1, 0, defaultText = "0")
    flankDialog.initDisplayState()

    fileSuffixDialog = dialog.createDynamicSelector(4, 0)
    fileSuffixDialog.initCheckboxController("Custom file suffix")
    suffixDialog = fileSuffixDialog.initDisplay(True, "Suffix")
    suffixDialog.createTextField("File Suffix:", 0, 0, defaultText="_flanked_colored")
    fileSuffixDialog.initDisplayState()

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    if dialog.selections.getToggleStates()[0]: colorColIndex = 6
    else: colorColIndex = None

    if flankDialog.getControllerVar(): 
        flankBinSize = int(dialog.selections.getTextEntries("FlankSize")[0])
        flankBinNum = int(dialog.selections.getTextEntries("FlankSize")[1])
    else: 
        flankBinSize = 0
        flankBinNum = 0

    if fileSuffixDialog.getControllerVar(): fileSuffix = dialog.selections.getTextEntries("Suffix")[0]
    else: fileSuffix = ""

    binInGenes(dialog.selections.getFilePathGroups()[0], dialog.selections.getIndividualFilePaths()[0],
               flankBinSize, flankBinNum, fileSuffix, colorColIndex)


if __name__ == "__main__": main()