# This script takes a file of some genome feature and nucleosome dyad centers and determines the 
# density of those features within a 100 bp radius of the dyad centers.
import os
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
from benbiohelpers.CountThisInThat.InputDataStructures import EncompassingDataDefaultStrand
from benbiohelpers.CountThisInThat.Counter import ThisInThatCounter
from benbiohelpers.CountThisInThat.CounterOutputDataHandler import CounterOutputDataHandler
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from typing import List


class NucleosomeFeatureCounter(ThisInThatCounter):

    def __init__(self, encompassedFeaturesFilePath, encompassingFeaturesFilePath, 
                 outputFilePath, acceptableChromosomes = None, checkForSortedFiles = True,
                 headersInEncompassedFeatures = False, headersInEncompassingFeatures = False,
                 encompassingFeatureExtraRadius = 0, minEncompassedDistance = 0):
        super().__init__(encompassedFeaturesFilePath, encompassingFeaturesFilePath, 
                 outputFilePath, acceptableChromosomes, checkForSortedFiles,
                 headersInEncompassedFeatures, headersInEncompassingFeatures,
                 encompassingFeatureExtraRadius)

        self.minEncompassedDistance = minEncompassedDistance

    def setUpOutputDataHandler(self):
        self.outputDataHandler = CounterOutputDataHandler()
        self.outputDataHandler.addEncompassingFeatureStratifier(outputName = "Nucleosome")
        self.outputDataHandler.addPlaceholderStratifier("Feature_Counts")

    def constructEncompassingFeature(self, line) -> EncompassingDataDefaultStrand:
        return EncompassingDataDefaultStrand(line, self.acceptableChromosomes)

    def isEncompassedFeatureWithinEncompassingFeature(self, encompassedFeature = None, encompassingFeature = None):

        if encompassedFeature is None: encompassedFeature = self.currentEncompassedFeature
        if encompassingFeature is None: encompassingFeature = self.currentEncompassingFeature

        return (super().isEncompassedFeatureWithinEncompassingFeature(encompassedFeature, encompassingFeature) and
                abs(encompassingFeature.center - encompassedFeature.position) >= self.minEncompassedDistance)


def countFeaturesAboutNucleosomes(genomeFeaturesFilePaths: List[str], nucleosomePosFilePath, onlyCountLinker, searchRadius = 100):

    if onlyCountLinker: minEncompassedDistance = 74
    else: minEncompassedDistance = 0

    for genomeFeaturesFilePath in genomeFeaturesFilePaths:

        print('\n' + "Working in",os.path.basename(genomeFeaturesFilePath))

        outputFilePath = genomeFeaturesFilePath.rsplit('.',1)[0]
        if onlyCountLinker: outputFilePath += "_nucleosome_stratification_linker_only.tsv"
        else: outputFilePath += "_nucleosome_stratification.tsv"

        counter = NucleosomeFeatureCounter(genomeFeaturesFilePath, nucleosomePosFilePath, outputFilePath, 
                                           encompassingFeatureExtraRadius = searchRadius, 
                                           minEncompassedDistance = minEncompassedDistance)
        counter.count()
        counter.writeResults((None,{None:"Feature_Counts"}))


def main():

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Genome Feature Positions Files:",0,"context_mutations.bed",("Bed Files",".bed"))    
    dialog.createFileSelector("Nucleosome Dyad Center Positions:",1,("Bed Files",".bed"))
    dialog.createCheckbox("Only Count Linker", 2, 0)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    countFeaturesAboutNucleosomes(dialog.selections.getFilePathGroups()[0], dialog.selections.getIndividualFilePaths()[0],
                                  dialog.selections.getToggleStates()[0])


if __name__ == "__main__": main()