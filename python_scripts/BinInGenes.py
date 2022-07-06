# This script takes a gene designations file and one or more files of features (like mutations) to bin in fractions of that gene.
import os
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.CountThisInThat.Counter import ThisInThatCounter
from benbiohelpers.CountThisInThat.OutputDataStratifiers import AmbiguityHandling


def binInGenes(featureFilePaths: List[str], geneDesignationsFilePath, flankingBinSize = 0, flankingBinNum = 0, 
               filePathSuffix = "", colorColIndex = None):

    class BinInGenesCounter(ThisInThatCounter):

        def setUpOutputDataHandler(self):
            super().setUpOutputDataHandler()
            if colorColIndex is not None:
                self.outputDataHandler.addSimpleEncompassingColStratifier(outputName = "Color_Domain", colIndex = colorColIndex)
            self.outputDataHandler.addFeatureFractionStratifier(outputName = "Gene_Fraction", flankingBinSize = flankingBinSize, 
                                                                flankingBinNum = flankingBinNum)
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


def main():

    try:
        from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
        workingDirectory = getDataDirectory()
    except ImportError:
        workingDirectory = os.path.dirname(__file__)

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=workingDirectory)
    dialog.createMultipleFileSelector("Feature Files (e.g. mutations):",0,"bin_me.bed",("Bed Files",".bed"))    
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
    suffixDialog.createTextField("File Suffix:", 0, 0, defaultText="all_genes")
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