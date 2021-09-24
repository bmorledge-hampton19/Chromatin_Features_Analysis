# This script takes a gene designations file and one or more files of features (like mutations) to bin in fractions of that gene.
import os
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
from benbiohelpers.CountThisInThat.Counter import ThisInThatCounter
from benbiohelpers.CountThisInThat.OutputDataStratifiers import AmbiguityHandling

class BinInGenesCounter(ThisInThatCounter):

    def setUpOutputDataHandler(self):
        super().setUpOutputDataHandler()
        self.outputDataHandler.addFeatureFractionStratifier(outputName = "Gene_Fraction")
        self.outputDataHandler.addStrandComparisonStratifier(strandAmbiguityHandling = AmbiguityHandling.tolerate)
        self.outputDataHandler.createOutputDataWriter(self.outputFilePath, customStratifyingNames=(None, {True:"Coding_Strand_Counts", 
                                                                                                          False:"Noncoding_Strand_Counts"}))


def binInGenes(featureFilePaths: List[str], geneDesignationsFilePath):

    for featureFilePath in featureFilePaths:

        print("Working in", os.path.basename(featureFilePath))

        outputFilePath = featureFilePath.rsplit('.', 1)[0] + "_gene_bins.tsv"

        counter = BinInGenesCounter(featureFilePath, geneDesignationsFilePath, outputFilePath)
        counter.count()


def main():

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Feature Files (e.g. mutations):",0,"bin_me.bed",("Bed Files",".bed"))    
    dialog.createFileSelector("Gene Designations:",1,("Bed Files",".bed"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    binInGenes(dialog.selections.getFilePathGroups()[0], dialog.selections.getIndividualFilePaths()[0])


if __name__ == "__main__": main()