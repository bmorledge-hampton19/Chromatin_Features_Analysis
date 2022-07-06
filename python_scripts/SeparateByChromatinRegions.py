# This script takes a bed file and the coordinates for chromatin regions
# and splits the rows in the bed file into new files for each domain.
# NOTE:  Both input files must be sorted for this script to run properly. 
#        (Sorted first by chromosome (string) and then by nucleotide position (numeric))

import os, warnings
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog, Selections
from benbiohelpers.FileSystemHandling.DirectoryHandling import checkDirs


class MutationData:

    def __init__(self, line: str):

        # Read in the next line.
        choppedUpLine = line.split()

        # Assign variables
        self.line = line
        self.chromosome = choppedUpLine[0] # The chromosome that houses the mutation.
        self.position = float(choppedUpLine[1]) # The position of the mutation in its chromosome. (0 base)
        self.domainName = None # The chromatin domain this mutation is a part of.


# Contains data on a single domain range obtained by reading the next available line in a given file.
class DomainData:

    def __init__(self, line: str):

        # Read in the next line.
        choppedUpLine = line.strip().split()

        self.chromosome = choppedUpLine[0]
        self.startPos = int(choppedUpLine[1]) # 0 base
        self.endPos = int(choppedUpLine[2]) - 1 # Still 0 base (I think?  I don't actually know whether or not the domain file is bed-formatted)
        self.domainName = choppedUpLine[3].replace('/',"_or_")


# Uses the given domain ranges file and mutation file to split mutations by domain. 
# Generates a folder of new files to store these results.
# NOTE:  It is VITAL that both files are sorted, first by chromosome number and then by starting coordinate.
#        (Sorted first by chromosome (string) and then by nucleotide position (numeric))
#        This code is pretty slick, but it will crash and burn and give you a heap of garbage as output if the inputs aren't sorted.
class DomainSplitter:

    def __init__(self, mutationFilePath, domainRangesFilePath):

        # Open the mutation and gene positions files to compare against one another.
        self.mutationFile = open(mutationFilePath, 'r')
        self.domainRangesFile = open(domainRangesFilePath,'r')

        # Set up the file system for outputting files for different domains..
        self.domainOutputFiles = dict()
        self.domainOutputFolder = os.path.join(os.path.dirname(mutationFilePath),os.path.basename(domainRangesFilePath).rsplit('.',1)[0])
        checkDirs(self.domainOutputFolder)
        self.domainOutputFilePathBasename = os.path.basename(mutationFilePath).rsplit('.',1)[0]

        # Keeps track of mutations that matched to a domain to check for overlap.
        self.mutationsInPotentialOverlap: List[MutationData] = list()

        # The mutation and gene currently being investigated.
        self.currentMutation: MutationData = None
        self.currentDomain: DomainData = None


    # Reads in the next mutation from the mutation data into currentMutation
    def readNextMutation(self) -> MutationData:

        # Read in the next line.
        nextLine = self.mutationFile.readline()

        # Check if EOF has been reached.
        if len(nextLine) == 0: 
            self.currentMutation = None
        # Otherwise, read in the next mutation.
        else:
            self.currentMutation = MutationData(nextLine)

    
    # Reads in the next domain from the domain ranges file into current domain
    def readNextDomain(self) -> DomainData:

        # Read in the next line.
        nextLine = self.domainRangesFile.readline()

        # Check if EOF has been reached.
        if len(nextLine) == 0: 
            self.currentDomain = None
        # Otherwise, read in the next mutation.
        else:
            self.currentDomain = DomainData(nextLine)

        # Check for mutations in overlapping regions between this gene and the last one.
        if self.currentDomain is not None: self.checkMutationsInOverlap() 


    # Takes a mutation object and domain object which have unequal chromosomes and read through data until they are equal.
    def reconcileChromosomes(self):
        
        chromosomeChanged = False # A simple flag to determine when to inform the user that a new chromosome is being accessed

        # Until the chromosomes are the same for both mutations and domains, read through the one with the eariler chromosome.
        while (self.currentMutation is not None and self.currentDomain is not None and 
               self.currentMutation.chromosome != self.currentDomain.chromosome):
            chromosomeChanged = True
            if self.currentMutation.chromosome < self.currentDomain.chromosome: self.readNextMutation()
            else: self.readNextDomain()

        if chromosomeChanged and self.currentDomain is not None and self.currentMutation is not None: 
            print("Binning by domain in",self.currentDomain.chromosome)


    # Determines whether or not the current mutation is past the range of the current domain.
    def isMutationPastDomain(self):

        if self.currentMutation is None:
            return True
        elif self.currentMutation.position > self.currentDomain.endPos:
            return True
        elif not self.currentMutation.chromosome == self.currentDomain.chromosome:
            return True
        else: 
            return False


    # A function which checks to see if the current mutation falls within the current domain and if it does, adds it to the list.
    # (Assumes that the current mutation is not beyond the domain because this is checked first by isMutationPastDomain.)
    def addMutationIfInDomain(self):

        if self.currentMutation.position >= self.currentDomain.startPos:
            
            # Set the mutation's domain name to this domain's name.
            self.currentMutation.domainName = self.currentDomain.domainName

            # Add the mutation to the list of mutations in the current domain
            self.mutationsInPotentialOverlap.append(self.currentMutation)


    # Check to see if any previous mutations called for previous domains are present in the current domain due to overlap.
    # If not, also write them to the relevant domain file.
    def checkMutationsInOverlap(self):    

        # First, flag any mutations that fall before the start position of the new domain to be written.
        mutationsToWrite = [mutation for mutation in self.mutationsInPotentialOverlap 
                            if mutation.position < self.currentDomain.startPos or 
                            mutation.chromosome != self.currentDomain.chromosome]

        self.mutationsInPotentialOverlap = list(set(self.mutationsInPotentialOverlap) - set(mutationsToWrite))

        # Next, check all remaining mutations to see if their previous domain assignment matches with the new domain.
        for mutation in self.mutationsInPotentialOverlap:

            # Is the mutation within the upper bound of the domain (endPos)?
            # Does the name for the previous domain(s) match this domain?
            if (mutation.position <= self.currentDomain.endPos and 
                mutation.domainName != self.currentDomain.domainName):

                # If not, switch the domain to ambiguous (None).
                mutation.domainName = None

        # Remove any mutations that were found to be ambiguous.
        self.mutationsInPotentialOverlap = [mutation for mutation in self.mutationsInPotentialOverlap 
                                            if mutation.domainName is not None]

        # Write the flagged mutations.
        for mutation in mutationsToWrite: self.writeMutationToDomainFile(mutation)


    # Write the mutation to its assigned domain.
    def writeMutationToDomainFile(self, mutation: MutationData):

        # First, determine if we have a new chromatin domain.
        # If we do, we need to set up a new output file for it.
        if not mutation.domainName in self.domainOutputFiles:

            domainOutputFilePath = os.path.join(self.domainOutputFolder, self.domainOutputFilePathBasename + '_' + mutation.domainName + "_domain.bed")
            self.domainOutputFiles[mutation.domainName] = open(domainOutputFilePath, 'w')

        # Now, write the mutation's line to the relevant file.
        self.domainOutputFiles[mutation.domainName].write(mutation.line)


    # Split given mutations into domain ranges present in the given file.  (Or, drop them if they don't belong to exactly one domain.)
    def splitByDomains(self):

        # Get data on the first mutation and domain and reconcile their chromosomes if necessary to start things off.
        # If either the mutation file or domain ranges file is empty, make sure to bypass the check.
        self.readNextMutation()
        self.readNextDomain()
        if self.currentMutation is None or self.currentDomain is None:
            warnings.warn("Empty mutation or domain ranges file.  Output will most likely be unhelpful.")
        elif self.currentDomain.chromosome == self.currentMutation.chromosome: 
            print("Binning by domain in",self.currentDomain.chromosome)
        else: self.reconcileChromosomes()

        # The core loop goes through each domain range one at a time and checks mutation positions against it until 
        # one exceeds its rightmost position or is on a different chromosome (or mutations are exhausted).  
        # Then, the next domain range is checked, then the next, etc. until no ranges remain.
        while self.currentDomain is not None:

            # Read mutations until the mutation is past the range of the current domain.
            while not self.isMutationPastDomain():

                # Check and see if we need to add the mutation to our lists.
                self.addMutationIfInDomain()
                #Get data on the next mutation.
                self.readNextMutation()

            # Read in a new domain.
            self.readNextDomain()

            # Reconcile the mutation data and domain data to be sure that they are looking at the same chromosome for the next iteration
            self.reconcileChromosomes()

        # Close the input files.
        self.mutationFile.close()
        self.domainRangesFile.close()

        for domain in self.domainOutputFiles:
            self.domainOutputFiles[domain].close()


# Main functionality starts here.
def separateByChromatinRegions(mutationFilePaths, domainRangesFilePath: str):

    # Loop through each given mutation file path, splitting it up based on the domain ranges given in the relevant file path.
    for mutationFilePath in mutationFilePaths:

        print("\nWorking with",os.path.basename(mutationFilePath))

        # Make sure we have the expected file type.
        if not "context_mutations" in os.path.basename(mutationFilePath): 
            warnings.warn("Mutation file is expected to have \"" + "context_mutations" + "\" in the name.  Are you sure this is the right file type?")

        # Ready, set, go!
        counter = DomainSplitter(mutationFilePath, domainRangesFilePath)
        counter.splitByDomains()


def main():

    try:
        from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
        workingDirectory = getDataDirectory()
    except ImportError:
        workingDirectory = os.path.dirname(__file__)

    # Create the Tkinter UI
    with TkinterDialog(workingDirectory=workingDirectory) as dialog:
        dialog.createMultipleFileSelector("File(s) to separate:",0, "context_mutations.bed",("Bed Files",".bed"))
        dialog.createFileSelector("Domain Range File:", 1, ("Bed File",".bed"))

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    mutationFilePaths = selections.getFilePathGroups()[0] # A list of mutation file paths
    domainRangesFilePath = selections.getIndividualFilePaths()[0] # The gene positions file path

    separateByChromatinRegions(mutationFilePaths, domainRangesFilePath)

if __name__ == "__main__": main()