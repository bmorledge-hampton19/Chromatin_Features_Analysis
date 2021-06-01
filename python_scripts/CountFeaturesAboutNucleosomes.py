# This script takes a file of some genome feature and nucleosome dyad centers and determines the 
# density of those features within a 100 bp radius of the dyad centers.
import os
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
from mypyhelper.CountThisInThat.InputDataStructures import EncompassingData
from mypyhelper.CountThisInThat.Counter import ThisInThatCounter
from mypyhelper.CountThisInThat.CounterOutputDataHandler import CounterOutputDataHandler
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog
from typing import List


class NucleosomeData(EncompassingData):

    def setLocationData(self, acceptableChromosomes):
        super().setLocationData(acceptableChromosomes)
        self.strand = '+'


class H1DensityCounter(ThisInThatCounter):

    def setUpOutputDataHandler(self):
        self.outputDataHandler = CounterOutputDataHandler(trackAllEncompassing=True)
        self.outputDataHandler.addEncompassingFeatureStratifier(outputName = "Nucleosome")
        self.outputDataHandler.addPlaceholderStratifier("Feature_Counts")

    def constructEncompassingFeature(self, line) -> NucleosomeData:
        return NucleosomeData(line, self.acceptableChromosomes)


def getNucleosomeH1Density(genomeFeaturesFilePaths: List[str], nucleosomePosFilePath, searchRadius = 100):

    for genomeFeaturesFilePath in genomeFeaturesFilePaths:

        print('\n' + "Working in",os.path.basename(genomeFeaturesFilePath))

        outputFilePath = genomeFeaturesFilePath.rsplit('.',1)[0] + "_nucleosome_stratification.tsv"

        counter = H1DensityCounter(genomeFeaturesFilePath, nucleosomePosFilePath, outputFilePath, 
                                encompassingFeatureExtraRadius = searchRadius)
        counter.count()
        counter.writeResults((None,{None:"Feature_Counts"}))


def main():

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Genome Feature Positions Files:",0,"context_mutations.bed",("Bed Files",".bed"))    
    dialog.createFileSelector("Nucleosome Dyad Center Positions:",1,("Bed Files",".bed"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    getNucleosomeH1Density(dialog.selections.getFilePathGroups()[0], dialog.selections.getIndividualFilePaths()[0])


if __name__ == "__main__": main()