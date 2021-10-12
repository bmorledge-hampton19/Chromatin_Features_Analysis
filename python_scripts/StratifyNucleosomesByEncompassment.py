# This script takes a file of some encompassing genome feature and tests nucleosomes for inclusion in at least
# one instance of that feature.  (e.g. nucleosomes in genes.)
import os
from posixpath import split
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
from benbiohelpers.CountThisInThat.InputDataStructures import EncompassedDataDefaultStrand
from benbiohelpers.CountThisInThat.Counter import ThisInThatCounter, ENCOMPASSED_DATA
from benbiohelpers.CountThisInThat.CounterOutputDataHandler import CounterOutputDataHandler
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.FileSystemHandling.DirectoryHandling import checkDirs
from typing import List


class EncompassedNucleosomesCounter(ThisInThatCounter):

    def setUpOutputDataHandler(self):
        self.outputDataHandler = CounterOutputDataHandler(self.writeIncrementally, trackAllEncompassed=True)
        self.outputDataHandler.addEncompassedFeatureStratifier(outputName = "Nucleosome")
        self.outputDataHandler.addPlaceholderStratifier(outputName = "Counts")
        self.outputDataHandler.createOutputDataWriter(self.outputFilePath)

    def constructEncompassedFeature(self, line) -> EncompassedDataDefaultStrand:
        return EncompassedDataDefaultStrand(line, self.acceptableChromosomes)


def stratifyNucleosomesByEncompassment(encompassingFeaturesFilePath, nucleosomeFilePath: str):
    
    print("\nWorking in",os.path.basename(nucleosomeFilePath))

    # Create the output file paths.
    baseNucFilePath = os.path.basename(nucleosomeFilePath)
    nucleosomeDir = os.path.dirname(nucleosomeFilePath)
    intermediateDir = os.path.join(nucleosomeDir,"intermediate_files")
    encompassedAndNotNucleosomesFilePath = os.path.join(intermediateDir, baseNucFilePath.rsplit('.',1)[0] +
                                                                         "_encompassed_and_non_encompassed.bed")
    
    encompassedNucleosomeDir = os.path.join(os.path.dirname(nucleosomeDir),nucleosomeDir+"_encompassed")
    encompassedNucleosomesFilePath = os.path.join(encompassedNucleosomeDir, baseNucFilePath.rsplit('.',1)[0] +
                                                                            "_encompassed.bed")

    nonEncompassedNucleosomeDir = os.path.join(os.path.dirname(nucleosomeDir),nucleosomeDir+"_non_encompassed")
    nonEncompassedNucleosomesFilePath = os.path.join(nonEncompassedNucleosomeDir, baseNucFilePath.rsplit('.',1)[0] +
                                                                                  "_non_encompassed.bed")

    checkDirs(intermediateDir, encompassedNucleosomeDir, nonEncompassedNucleosomeDir)

    # Count!
    print("Checking for encompassment...")
    encompassedNucleosomesCounter = EncompassedNucleosomesCounter(nucleosomeFilePath, encompassingFeaturesFilePath,
                                                                  encompassedAndNotNucleosomesFilePath,
                                                                  writeIncrementally = ENCOMPASSED_DATA)                                            
    encompassedNucleosomesCounter.count()

    # Split the results into the encompassed and non-encompassed files.
    print("Splitting results to files for encompassement and non-encompassment.")
    with open(encompassedAndNotNucleosomesFilePath, 'r') as encompassedAndNotNucleosomesFile:
        with open(encompassedNucleosomesFilePath, 'w') as encompassedNucleosomesfile:
            with open(nonEncompassedNucleosomesFilePath, 'w') as nonEncompassedNucleosomesFile:

                for line in encompassedAndNotNucleosomesFile:
                    splitLine = line.split()

                    if int(splitLine[-1]) > 0: encompassedNucleosomesfile.write('\t'.join(splitLine[:-1]) + '\n')
                    else: nonEncompassedNucleosomesFile.write('\t'.join(splitLine[:-1]) + '\n')


def main():

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createFileSelector("Encompassing Feature File:",0,("Bed Files",".bed"))    
    dialog.createFileSelector("Nucleosome Dyad Center Positions:",1,("Bed Files",".bed"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    stratifyNucleosomesByEncompassment(dialog.selections.getIndividualFilePaths()[0], 
                                       dialog.selections.getIndividualFilePaths()[1])


if __name__ == "__main__": main()