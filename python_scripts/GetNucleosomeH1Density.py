# This script takes a file of H1 midpoints and nucleosome dyad centers and determines the 
# density of H1 proteins within a 100 bp radius of the dyad centers.
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
from mypyhelper.CountThisInThat.InputDataStructures import EncompassingData
from mypyhelper.CountThisInThat.Counter import ThisInThatCounter
from mypyhelper.CountThisInThat.CounterOutputDataHandler import CounterOutputDataHandler
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog


class NucleosomeData(EncompassingData):

    def setLocationData(self, acceptableChromosomes):
        super().setLocationData(acceptableChromosomes)
        self.strand = '+'


class H1DensityCounter(ThisInThatCounter):

    def setUpOutputDataHandler(self):
        self.outputDataHandler = CounterOutputDataHandler()
        self.outputDataHandler.addEncompassingFeatureStratifier(outputName = "Nucleosome")
        self.outputDataHandler.addPlaceholderStratifier("H1_Counts")

    def constructEncompassingFeature(self, line) -> NucleosomeData:
        return NucleosomeData(line, self.acceptableChromosomes)


def getNucleosomeH1Density(h1CentersFilePath, nucleosomePosFilePath, outputFilePath, searchRadius = 100):

    counter = H1DensityCounter(h1CentersFilePath, nucleosomePosFilePath, outputFilePath, encompassingFeatureExtraRadius = searchRadius)
    counter.count()
    counter.writeResults((None,{None:"Proximal_H1_Count"}))


def main():

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createFileSelector("H1 centers:",0,("Bed Files",".bed"))    
    dialog.createFileSelector("Nucleosome Dyad Center Positions:",1,("Bed Files",".bed"))
    dialog.createFileSelector("Output File:",2,("TSV Files",".tsv"), newFile = True)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    getNucleosomeH1Density(dialog.selections.getIndividualFilePaths()[0], dialog.selections.getIndividualFilePaths()[1],
                           dialog.selections.getIndividualFilePaths()[2])


if __name__ == "__main__": main()