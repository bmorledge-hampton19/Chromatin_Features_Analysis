# This script takes a file of gene designations (sorted) and combines any overlapping regions.
# Can either preserve strand information and discard ambiguous regions or 
# discard strand information to preserve ambiguous regions.
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory


def mergeGeneRanges(geneDesignationsFilePaths: List[str], preserveAmbiguousStrandRegions):

    for geneDesignationsFilePath in geneDesignationsFilePaths:

        print("Merging gene ranges for",geneDesignationsFilePath)

        # First, condense all overlapping gene regions and remove any ambiguous regions.
        mergedGeneRangesFilePath = geneDesignationsFilePath.rsplit('.', 1)[0] + "_merged.bed"

        currentGeneRangeChromosome = None
        currentGeneRangeStart = None
        currentGeneRangeEnd = None
        currentGeneRangeStrand = None

        with open(geneDesignationsFilePath, 'r') as geneDesignationsFile:
            with open(mergedGeneRangesFilePath, 'w') as mergedGeneRangesFile:
                for line in geneDesignationsFile:

                    # Parse out the gene range info from the current line.
                    choppedUpLine = line.strip().split('\t')
                    lineChromosome = choppedUpLine[0]
                    lineGeneStart = int(choppedUpLine[1])
                    lineGeneEnd = int(choppedUpLine[2])
                    lineStrand = choppedUpLine[5]

                    # Unless we are starting a new gene range, check to see if the gene region on this line overlaps with the current one.
                    if currentGeneRangeChromosome is not None:

                        # If they overlap, expand the current range and check to see if strands match.
                        if currentGeneRangeChromosome == lineChromosome and lineGeneStart < currentGeneRangeEnd:
                            
                            currentGeneRangeEnd = lineGeneEnd
                            if currentGeneRangeStrand is not None and currentGeneRangeStrand != lineStrand: currentGeneRangeStrand = '.'

                        # If they don't overlap, check to see if the strand designation for the current region is unambiguous, then write it.
                        # Also, keep in mind to expand the ranges by one bp on either side for trinucleotide context at the borders.
                        else:

                            if currentGeneRangeStrand != '.' or preserveAmbiguousStrandRegions:
                                mergedGeneRangesFile.write('\t'.join((currentGeneRangeChromosome, str(currentGeneRangeStart),
                                                                      str(currentGeneRangeEnd), '.', '.', currentGeneRangeStrand)) + '\n')
                            
                            # Don't forget to reset the chromosome variable to flag the rest for reassignment!
                            currentGeneRangeChromosome = None


                    # If we are starting to look at a new gene range, assign all the values from this line.
                    if currentGeneRangeChromosome is None:
                        currentGeneRangeChromosome = lineChromosome
                        currentGeneRangeStart = lineGeneStart
                        currentGeneRangeEnd = lineGeneEnd
                        currentGeneRangeStrand = lineStrand

                # Do one last check so we don't miss the last gene range.
                if currentGeneRangeStrand != '.' or preserveAmbiguousStrandRegions:
                    mergedGeneRangesFile.write('\t'.join((currentGeneRangeChromosome, str(currentGeneRangeStart), 
                                                          str(currentGeneRangeEnd), '.', '.', currentGeneRangeStrand)) + '\n')


def main():

    # Create a simple dialog for selecting the gene expression file.
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Gene Designations Files", 0, "gene_designations.bed", ("Bed File", ".bed"))
    dialog.createCheckbox("Preserve Ambiguous Regions", 1, 0)
    dialog.mainloop()

    if dialog.selections is None: quit()

    # Retrieve the selections and pass the relevant arguments to the primary function.
    mergeGeneRanges(dialog.selections.getFilePathGroups()[0], dialog.selections.getToggleStates()[0])


if __name__ == "__main__": main()