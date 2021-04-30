# This script takes a mutation file and the coordinates for binding motifs
# and calculates how many mutations occurred at positions within those motifs.
# NOTE:  Both input files must be sorted for this script to run properly. 
#        (Sorted first by chromosome (string) and then by nucleotide position (numeric))

import os, warnings
from typing import List
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (Metadata, generateFilePath, getDataDirectory,
                                                                  DataTypeStr, getAcceptableChromosomes)

class MutationData:

    def __init__(self, line, acceptableChromosomes):

        # Read in the next line.
        choppedUpLine = line.strip().split()

        # Assign variables
        self.chromosome = choppedUpLine[0] # The chromosome that houses the mutation.
        self.position = float(choppedUpLine[1]) # The position of the mutation in its chromosome. (0 base)
        self.strand = choppedUpLine[5] # Either '+' or '-' depending on which strand houses the mutation.

        # Make sure the mutation is in a valid chromosome.
        if not self.chromosome in acceptableChromosomes:
            raise ValueError(choppedUpLine[0] + " is not a valid chromosome for the mutation trinuc file.")


# Contains data on a single binding motif position obtained by reading the next available line in a given file.
class BindingMotifData:

    def __init__(self, line):

        # Read in the next line.
        choppedUpLine = line.strip().split()

        self.chromosome = choppedUpLine[0]
        self.startPos = int(choppedUpLine[1]) # 0 base
        self.endPos = int(choppedUpLine[2]) - 1 # Still 0 base 
        self.strand = choppedUpLine[5] # The strand with the binding motif


# Uses the given binding motif positions file and mutation file to count the number of mutations 
# in those regions for each position.  
# Generates a new file to store these results.
# NOTE:  It is VITAL that both files are sorted, first by chromosome number and then by starting coordinate.
#        (Sorted first by chromosome (string) and then by nucleotide position (numeric))
#        This code is pretty slick, but it will crash and burn and give you a heap of garbage as output if the inputs aren't sorted.
class CountsFileGenerator:

    def __init__(self, mutationFilePath, bindingMotifsFilePath, 
                 bindingMotifsMutationCountsFilePath, acceptableChromosomes):

        # Open the mutation and binding motif positions files to compare against one another.
        self.mutationFile = open(mutationFilePath, 'r')
        self.bindingMotifsFile = open(bindingMotifsFilePath,'r')

        # Store the other arguments passed to the constructor
        self.acceptableChromosomes = acceptableChromosomes
        self.bindingMotifsMutationCountsFilePath = bindingMotifsMutationCountsFilePath

        # Dictionaries holding the number of mutations found at each position in the binding motif.
        # Key is an integer giving a 1-based position of the mutation relative to the center (rounded down) of the binding motif.
        # Each position in the reverse dictionary is the complementary base at the same position in the binding motif strand.
        self.bindingMotifStrandMutationCounts = dict()
        self.reverseMotifStrandMutationCounts = dict()

        # Half of the binding motif length, rounded down.
        self.halfBindingMotifLength = None

        # Keeps track of mutations that matched to a binding motif to check for overlap.
        self.mutationsInPotentialOverlap: List[MutationData] = list()

        # The mutation and binding motif currently being investigated.
        self.currentMutation: MutationData = None
        self.bindingMotif: BindingMotifData = None


    # Reads in the next mutation from the mutation data into currentMutation
    def readNextMutation(self) -> MutationData:

        # Read in the next line.
        nextLine = self.mutationFile.readline()

        # Check if EOF has been reached.
        if len(nextLine) == 0: self.currentMutation = None
        # Otherwise, read in the next mutation.
        else: self.currentMutation = MutationData(nextLine, self.acceptableChromosomes)

    
    # Reads in the next bindingMotif from the file of binding motif positions
    def readNextBindingMotif(self) -> BindingMotifData:

        # Read in the next line.
        nextLine = self.bindingMotifsFile.readline()

        # Check if EOF has been reached.
        if len(nextLine) == 0: 
            self.bindingMotif = None
        # Otherwise, read in the next mutation.
        else:
            self.bindingMotif = BindingMotifData(nextLine)

        # Check for mutations in overlapping regions between this binding motif and the last one.
        if self.bindingMotif is not None: self.checkMutationsInOverlap() 


    # Takes a mutation object and binding motif object which have unequal chromosomes and read through data until they are equal.
    def reconcileChromosomes(self):
        
        chromosomeChanged = False # A simple flag to determine when to inform the user that a new chromosome is being accessed

        # Until the chromosomes are the same for both mutations and binding motifs, read through the one with the eariler chromosome.
        while (self.currentMutation is not None and self.bindingMotif is not None and 
               self.currentMutation.chromosome != self.bindingMotif.chromosome):
            chromosomeChanged = True
            if self.currentMutation.chromosome < self.bindingMotif.chromosome: self.readNextMutation()
            else: self.readNextBindingMotif()

        if chromosomeChanged and self.bindingMotif is not None and self.currentMutation is not None: 
            print("Counting in",self.bindingMotif.chromosome)


    # Determines whether or not the current mutation is past the range of the current binding motif.
    def isMutationPastBindingMotif(self):

        if self.currentMutation is None:
            return True
        elif self.currentMutation.position > self.bindingMotif.endPos:
            return True
        elif not self.currentMutation.chromosome == self.bindingMotif.chromosome:
            return True
        else: 
            return False


    # A function which checks to see if the current mutation falls within the current binding motif and if it does, adds it to the dictionary.
    # (Assumes that the current mutation is not beyond the binding motif because this is checked first by isMutationPastBindingMotif.)
    def addMutationIfInBindingMotif(self, mutation: MutationData, bindingMotif: BindingMotifData, checkingOverlap):

        if mutation.position >= bindingMotif.startPos:

            # Assign this mutation to its position in the binding motif
            relativeMutPos = mutation.position - bindingMotif.startPos - self.halfBindingMotifLength
            if mutation.strand == bindingMotif.strand: 
                self.bindingMotifStrandMutationCounts[relativeMutPos] += 1
            else: 
                self.reverseMotifStrandMutationCounts[relativeMutPos] += 1
            
            # Add the mutation to the list of mutations in the current binding motif (if not already checking for overlap)
            if not checkingOverlap: self.mutationsInPotentialOverlap.append(mutation)


    # Check to see if any previous mutations called for previous binding motifs are present in the current binding motif due to overlap.
    def checkMutationsInOverlap(self):    

        # First, get rid of any mutations that fall before the start position of the new binding motif.
        self.mutationsInPotentialOverlap = [mutation for mutation in self.mutationsInPotentialOverlap 
                                            if mutation.position >= self.bindingMotif.startPos and 
                                            mutation.chromosome == self.bindingMotif.chromosome]

        # Next, add all remaining mutations to the binding motif counts for the current binding motif.
        for mutation in self.mutationsInPotentialOverlap:           
            self.addMutationIfInBindingMotif(mutation, self.bindingMotif, True)


    # Count all mutations within binding motifs and assign a relative position to them.
    def count(self):
        # Get data on the first mutation and binding motif and reconcile their chromosomes if necessary to start things off.
        # If either the mutation file or binding motif file is empty, make sure to bypass the check.
        self.readNextMutation()
        self.readNextBindingMotif()
        if self.currentMutation is None or self.bindingMotif is None:
            warnings.warn("Empty Mutation or Binding Motif Positions file.  Output will most likely be unhelpful.")
        elif self.bindingMotif.chromosome == self.currentMutation.chromosome: 
            print("Counting in",self.bindingMotif.chromosome)
        else: self.reconcileChromosomes()

        # Set up the mutation count dictionaries using the length of the current motif (assume it's constant across all motifs)
        bindingMotifLength = self.bindingMotif.endPos - self.bindingMotif.startPos + 1
        self.halfBindingMotifLength = bindingMotifRangeEnd = int(bindingMotifLength / 2)
        bindingMotifRangeStart = -bindingMotifRangeEnd
        if bindingMotifLength % 2 == 1: bindingMotifRangeEnd += 1

        self.intPositions = list()
        self.halfPositions = list()
        for i in range(bindingMotifRangeStart, bindingMotifRangeEnd):
            self.bindingMotifStrandMutationCounts[i] = 0
            self.reverseMotifStrandMutationCounts[i] = 0
            self.intPositions.append(i)
            if i < bindingMotifRangeEnd:
                self.bindingMotifStrandMutationCounts[i+0.5] = 0
                self.reverseMotifStrandMutationCounts[i+0.5] = 0
                self.halfPositions.append(i+0.5)

        # The core loop goes through each binding motif one at a time and checks mutation positions against it until 
        # one exceeds its rightmost position or is on a different chromosome (or mutations are exhausted).  
        # Then, the next binding motif is checked, then the next, etc. until no binding motifs remain.
        while self.bindingMotif is not None:

            # Read mutations until the mutation is past the range of the current binding motif.
            while not self.isMutationPastBindingMotif():

                # Check and see if we need to add the mutation to our lists.
                self.addMutationIfInBindingMotif(self.currentMutation, self.bindingMotif, False)
                #Get data on the next mutation.
                self.readNextMutation()

            # Read in a new binding motif.
            self.readNextBindingMotif()

            # Reconcile the mutation data and binding motif data to be sure
            # that they are looking at the same chromosome for the next iteration
            self.reconcileChromosomes()

        # Close the input files.
        self.mutationFile.close()
        self.bindingMotifsFile.close()


    def writeResults(self):

        # Write the results to the output file.
        with open(self.bindingMotifsMutationCountsFilePath,'w') as BMMutationCountsFile:
            
            # Write the headers to the file.
            BMMutationCountsFile.write('\t'.join(("Motif_Position", "Motif_Strand_Mutation_Counts", "Reverse_Strand_Mutation_Counts")) + '\n')

            # Determine whether or not to include just half positions, just int positions, or all positions in the final output.
            motifPosRange = self.intPositions
            hasHalfPositions = False

            # Are there half position counts?
            for i in self.halfPositions:
                if self.bindingMotifStrandMutationCounts[i] > 0 or self.reverseMotifStrandMutationCounts[i] > 0:
                    motifPosRange = self.halfPositions
                    hasHalfPositions = True
                    break

            # If there were half position counts, are there also int position counts?
            if hasHalfPositions:
                for i in self.intPositions:
                    if self.bindingMotifStrandMutationCounts[i] > 0 or self.reverseMotifStrandMutationCounts[i] > 0:
                        motifPosRange = self.bindingMotifStrandMutationCounts.keys()
                        break

            # Write the counts.
            for pos in motifPosRange:
                BMMutationCountsFile.write('\t'.join((str(pos), str(self.bindingMotifStrandMutationCounts[pos]), 
                                                        str(self.reverseMotifStrandMutationCounts[pos]))) + '\n')


# Main functionality starts here.
def countInBindingMotifs(mutationFilePaths, bindingMotifsFilePaths):

    bindingMotifsMutationCountsFilePaths = list() # A list of paths to the output files generated by the function

    # Loop through each given mutation file path, creating a corresponding binding motifs mutation count file for each.
    for mutationFilePath in mutationFilePaths:
        for bindingMotifsFilePath in bindingMotifsFilePaths:

            print("\nWorking with",os.path.basename(mutationFilePath),"and",os.path.basename(bindingMotifsFilePath))

            # Make sure we have the expected file type.
            if not DataTypeStr.mutations in os.path.basename(mutationFilePath): 
                raise ValueError("Mutation file should have \"" + DataTypeStr.mutations + "\" in the name.")
            
            # Get metadata and use it to generate a path to the nucleosome positions file.
            metadata = Metadata(mutationFilePath)

            # Generate the output file path for mutation counts.
            binder = os.path.basename(bindingMotifsFilePath).rsplit("binding_motifs",1)[0]
            if "binding_motifs" not in os.path.basename(bindingMotifsFilePath):
                warnings.warn("\"binding_motifs\" not found in basename of binding motifs file.  The output file's name is probably a garbled mess.")

            bindingMotifsMutationCountsFilePath = generateFilePath(directory = metadata.directory,
                                                                    dataGroup = metadata.dataGroupName, 
                                                                    fileExtension = ".tsv", dataType = binder + "binding_motif_mutation_counts")
            bindingMotifsMutationCountsFilePaths.append(bindingMotifsMutationCountsFilePath)

            # Ready, set, go!
            counter = CountsFileGenerator(mutationFilePath, bindingMotifsFilePath, bindingMotifsMutationCountsFilePath, 
                                        getAcceptableChromosomes(metadata.genomeFilePath))
            counter.count()
            counter.writeResults()

    return bindingMotifsMutationCountsFilePaths


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Mutation Files:",0,DataTypeStr.mutations + ".bed",("Bed Files",".bed"))
    dialog.createMultipleFileSelector("Binding Motifs Files:", 1, "binding_motifs.bed", ("Bed Files",".bed"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    countInBindingMotifs(dialog.selections.getFilePathGroups()[0], dialog.selections.getFilePathGroups()[1])

if __name__ == "__main__": main()