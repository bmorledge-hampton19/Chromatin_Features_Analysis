# This script bins RNA reads by chromatin domain, assuming they fall within a gene (determined by given gene designations)
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
from benbiohelpers.CountThisInThat.Counter import ThisInThatCounter
from benbiohelpers.CountThisInThat.OutputDataStratifiers import AmbiguityHandling


def binRNASeqByChromatinDomainInGenes(rNASeqFilePath: str, geneDesignationsFilePath, colorColIndex):

    class BinByCDsInGenesCounter(ThisInThatCounter):

        def setUpOutputDataHandler(self):
            super().setUpOutputDataHandler()
            self.outputDataHandler.addSimpleEncompassingColStratifier(outputName = "Color_Domain", colIndex = colorColIndex)
            self.outputDataHandler.addPlaceholderStratifier(outputName = "RNAseq_Reads")
            self.outputDataHandler.createOutputDataWriter(self.outputFilePath)

    outputFilePath = rNASeqFilePath.rsplit('.', 1)[0] + "_chromatin_domain_counts.tsv"

    counter = BinByCDsInGenesCounter(rNASeqFilePath, geneDesignationsFilePath, outputFilePath)
    counter.count()


def main():

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createFileSelector("RNAseq File:",0,("Bed Files",".bed"))    
    dialog.createFileSelector("Gene Designations (color in 7th column):",1,("Bed Files",".bed"))

    # Run the UI
    dialog.mainloop()

    binRNASeqByChromatinDomainInGenes(dialog.selections.getIndividualFilePaths()[0],
                                      dialog.selections.getIndividualFilePaths()[1], 6)


if __name__ == "__main__": main()