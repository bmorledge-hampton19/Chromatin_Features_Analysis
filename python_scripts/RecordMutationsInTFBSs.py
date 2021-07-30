# This script counts mutations that occur in transcription factor binding sites (TFBSs),
# keeping track of those counts for every genome position with > 1 mutation and the TFBSs encompassing that location.

from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
from benbiohelpers.CountThisInThat.Counter import ThisInThatCounter
from benbiohelpers.CountThisInThat.InputDataStructures import TfbsData
from benbiohelpers.CountThisInThat.CounterOutputDataHandler import AmbiguityHandling, CounterOutputDataHandler
from benbiohelpers.CountThisInThat.SupplementalInformation import TfbsSupInfoHandler

class MutationsInTfbsCounter(ThisInThatCounter):

    def setUpOutputDataHandler(self):
        self.outputDataHandler = CounterOutputDataHandler()
        self.outputDataHandler.addEncompassedFeatureStratifier("Mutation Position", False)
        #self.outputDataHandler.addStrandComparisonStratifier(AmbiguityHandling.record)
        self.outputDataHandler.addPlaceholderStratifier(AmbiguityHandling.record)
        self.outputDataHandler.addSupplementalInformationHandler(TfbsSupInfoHandler, 0)

    def constructEncompassingFeature(self, line) -> TfbsData:
        return TfbsData(line, self.acceptableChromosomes)


def recordMutationsInTFBSs(mutationPosFilePath, tFBSPosFilePath, outputFilePath):

    counter = MutationsInTfbsCounter(mutationPosFilePath, tFBSPosFilePath, outputFilePath)
    counter.count()
    counter.writeResults( (None,{None:"Counts"}) )
    #counter.writeResults( (None,{True:"Same_Strand",False:"Opposite_Strand",None:"Ambiguous"}) )


def main():

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createFileSelector("Bed Mutation Data:",0,("Bed Files",".bed"))    
    dialog.createFileSelector("TFBS Positions:",1,("Bed Files",".bed"))
    dialog.createFileSelector("Output File:",2,("TSV Files",".tsv"), newFile = True)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    recordMutationsInTFBSs(dialog.selections.getIndividualFilePaths()[0], dialog.selections.getIndividualFilePaths()[1],
                           dialog.selections.getIndividualFilePaths()[2])


if __name__ == "__main__": main()