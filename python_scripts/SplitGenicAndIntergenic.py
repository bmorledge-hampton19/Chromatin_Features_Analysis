# This script takes a file of some genome feature and gene ranges and outputs two files:
# One for genic features and one for intergenic.
import os
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
from benbiohelpers.CountThisInThat.InputDataStructures import ENCOMPASSED_DATA
from benbiohelpers.CountThisInThat.Counter import ThisInThatCounter
from benbiohelpers.CountThisInThat.CounterOutputDataHandler import CounterOutputDataHandler
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.FileSystemHandling.DirectoryHandling import checkDirs
from typing import List


class GenicVsIntergenicCounter(ThisInThatCounter):

    def initOutputDataHandler(self):
        self.outputDataHandler = CounterOutputDataHandler(self.writeIncrementally, trackAllEncompassed = True,
                                                          countNonCountedEncompassedAsNegative = True)

    def setupOutputDataStratifiers(self):
        self.outputDataHandler.addEncompassedFeatureStratifier()
        self.outputDataHandler.addPlaceholderStratifier()


def splitGenicAndIntergenic(genomeFeaturesFilePaths: List[str], geneRegionsFilePath):

    for genomeFeaturesFilePath in genomeFeaturesFilePaths:

        print('\n' + "Working in",os.path.basename(genomeFeaturesFilePath))

        # First, count the number of times that each feature is found within a gene.
        print("Counting features in genic regions...")
        intermediateDirectory = os.path.join(os.path.dirname(genomeFeaturesFilePath),"intermediate_files")
        checkDirs(intermediateDirectory)

        genicCountsOutputFilePath = os.path.join(intermediateDirectory,
                                                 os.path.basename(genomeFeaturesFilePath).rsplit('.',1)[0] + "_genic_counts.bed")
        genicOutputFilePath = genomeFeaturesFilePath.rsplit('.',1)[0] + "_genic.bed"
        intergenicOutputFilePath = genomeFeaturesFilePath.rsplit('.',1)[0] + "_intergenic.bed"
        
        counter = GenicVsIntergenicCounter(genomeFeaturesFilePath, geneRegionsFilePath, genicCountsOutputFilePath,
                                           writeIncrementally = ENCOMPASSED_DATA)
        counter.count()


        # Next, split up the bed entries into the genic and intergenic files based on whether or not they had any counts.
        print("Splitting results into genic and intergenic files...")
        with open(genicCountsOutputFilePath, 'r') as genicCountsOutputFile:
            with open(genicOutputFilePath, 'w') as genicOutputFile:
                with open(intergenicOutputFilePath, 'w') as intergenicOutputFile:

                    for line in genicCountsOutputFile:

                        choppedUpLine = line.strip().split('\t')
                        counts = int(choppedUpLine[-1])
                        if counts < 0:
                            for _ in range(-counts):
                                intergenicOutputFile.write('\t'.join(choppedUpLine[:-1]) + '\n')
                        elif counts > 0: 
                            for _ in range(counts): 
                                genicOutputFile.write('\t'.join(choppedUpLine[:-1]) + '\n')
                        else: raise ValueError("Counts value of 0 found.  Wat??")


def main():

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Genome Feature Positions Files:",0,"context_mutations.bed",("Bed Files",".bed"))    
    dialog.createFileSelector("Gene Ranges File (merged):",1,("Bed Files",".bed"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    splitGenicAndIntergenic(dialog.selections.getFilePathGroups()[0], dialog.selections.getIndividualFilePaths()[0])


if __name__ == "__main__": main()