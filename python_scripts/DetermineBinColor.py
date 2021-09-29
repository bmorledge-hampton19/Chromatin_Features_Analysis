from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.CountThisInThat.InputDataStructures import EncompassingData, EncompassingDataDefaultStrand, ColorDomainData
from typing import List
import os


# This function takes a bed file of chromatin color domains and assigns a color to bins regularly spaced to cover the whole genome.
# Bins are colored based on the majority domain coverage present in that region.
# If no domain achieves minimum coverage, it defaults to gray.
# NOTE: input files must be sorted by chromosome ID (alphabetically) and feature start position (numerically).
def determineRegularBinColors(colorDomainsFilePath, chromSizesFilePath, binSize, minimumCoverage = 0.5):

    # Retrieve information on the sizes of the chromosomes being used.
    chromSizes = dict()
    with open(chromSizesFilePath, 'r') as chromSizesFile:
        for line in chromSizesFile:
            chromID, chromSize = line.split()
            chromSizes[chromID] = int(chromSize)

    print("\nWorking in:", os.path.basename(colorDomainsFilePath))

    # Generate an output file path
    binnedFeaturesFilePath = colorDomainsFilePath.rsplit('.', 1)[0] + '_' + str(binSize) + "bp_binned.tsv"

    # Prepare for binning!
    bins = dict() # A nested dictionary.  The first key is for chromosome number, the second is for the start of the bin.
    with open(colorDomainsFilePath, 'r') as colorDomainsFile:

        # Read in the first line of the input file.
        choppedUpLine = colorDomainsFile.readline().split()
        if not choppedUpLine: domainChrom = None
        else: 
            domainChrom = choppedUpLine[0]
            assert domainChrom in chromSizes, "Unrecognized chromosome: " + domainChrom
            domainStartPos = int(choppedUpLine[1])
            domainEndPos = int(choppedUpLine[2]) - 1
            domainColor = choppedUpLine[3]

        # Iterate through the chromosome, determining the color of each bin.
        for binChrom in chromSizes:

            print("Binning in",binChrom)

            bins[binChrom] = dict()
            binStart = 0

            # Bin until there is no more domain data for the current chromosome AND all bins have been initialized for the current chromosome.
            while (domainChrom is not None and domainChrom == binChrom) or binStart < chromSizes[binChrom]:
                
                # Initialize the dictionary and set the minimum coverage level as GRAY.
                encompassedBasesByColor = dict()
                encompassedBasesByColor["GRAY"] = binSize*minimumCoverage

                # Are we within the current bin?  If so, check to see how much of the bin is encompassed.
                while domainChrom is not None and domainChrom == binChrom and domainStartPos < binStart + binSize:

                    if domainEndPos < binStart + binSize and domainStartPos >= binStart:
                        encompassedBases = (domainEndPos - domainStartPos + 1)
                    elif domainEndPos < binStart + binSize:
                        encompassedBases = (domainEndPos - binStart + 1)
                    elif domainStartPos >= binStart:
                        encompassedBases = (binStart + binSize - domainStartPos)
                    else: encompassedBases = binSize

                    encompassedBasesByColor[domainColor] = encompassedBasesByColor.setdefault(domainColor, 0) + encompassedBases

                    # Read in the next domain, UNLESS it extends into the next bin.  In that case, continue onto the next bin.
                    # NOTE: This means that this code doesn't handle overlapping domains very well, so I'm just assuming they don't occur very frequently.
                    if domainEndPos < binStart + binSize:
                        choppedUpLine = colorDomainsFile.readline().split()
                        if not choppedUpLine: domainChrom = None
                        else: 
                            domainChrom = choppedUpLine[0]
                            domainStartPos = float(choppedUpLine[1])
                            domainEndPos = int(choppedUpLine[2])
                            domainColor = choppedUpLine[3]
                    else: break

                # Assign the majority color to the bin.
                maxCoverage = max(encompassedBasesByColor.values())
                maxColors = [key for key, value in encompassedBasesByColor.items() if value == maxCoverage]
                if len(maxColors) > 1: bins[binChrom][binStart] = "GRAY"
                else: bins[binChrom][binStart] = maxColors[0]

                # Increment the binStart.
                binStart += binSize

            # We should never exit a bin chromosome while we still have a feature of that chromosome... Right?  *Sigh* Better double check...
            assert domainChrom != binChrom, ("Chromosome " + domainChrom + " bin exited before assigning feature starting at " + 
                                                str(domainStartPos) + ".  Are the chrom.sizes incorrect?")
            # Also, make sure we recognize the new chromosome from the chrom.sizes file.
            assert domainChrom is None or domainChrom in chromSizes, "Unrecognized chromosome: " + domainChrom


    # Write the results of the binning.
    with open(binnedFeaturesFilePath, 'w') as binnedFeaturesFile:

        print("Writing results...")

        # Write headers
        binnedFeaturesFile.write('\t'.join(("Chromosome","Bin_Start-End","Domain_Color")) + '\n')

        # Write the bins and feature counts.
        for chromosome in bins:
            for binStart in bins[chromosome]:

                binnedFeaturesFile.write('\t'.join((chromosome, str(binStart)+'-'+str(binStart+binSize-1), bins[chromosome][binStart])) + '\n')


# This function takes a bed file of chromatin color domains and a bed file of specified regions
# and assigns a color to each region based on majority coverage, defaulting to gray if no domain achieves minimum coverage.
# NOTE: input files must be sorted by chromosome ID (alphabetically) and feature start position (numerically).
def determineSpecifiedBinColors(colorDomainsFilePath, featureFilePath: str, minimumCoverage = 0.5):

    print("\nWorking in:", os.path.basename(colorDomainsFilePath))

    # Generate an output file path
    coloredFeaturesFilePath = featureFilePath.rsplit('.', 1)[0] + "_color_domain_designations.tsv"

    # Prepare for binning!
    with open(colorDomainsFilePath, 'r') as colorDomainsFile:
        with open(featureFilePath, 'r') as featureFile:
            with open(coloredFeaturesFilePath, 'w') as coloredFeaturesFile:

                # Read in the first line of the color domains file.
                colorDomainsFileLine = colorDomainsFile.readline()
                if not colorDomainsFileLine: currentDomainData = None
                else:  currentDomainData = ColorDomainData(colorDomainsFileLine, None)

                # Iterate through the features, determining the color for each.
                currentChrom = "Not a chromosome name yet."
                validDomainsForFeature: List[ColorDomainData] = list()
                for featureFileLine in featureFile:

                    featureData = EncompassingDataDefaultStrand(featureFileLine, None)

                    if featureData.chromosome != currentChrom: 
                        currentChrom = featureData.chromosome
                        print("Assigning domain colors in", currentChrom)
                    
                    # Until the next domain is fully beyond the current feature, add it to the valid domains list.
                    while not isACompletelyPastB(currentDomainData, featureData):
                        validDomainsForFeature.append(currentDomainData)
                        colorDomainsFileLine = colorDomainsFile.readline()
                        if not colorDomainsFileLine: currentDomainData = None
                        else:  currentDomainData = ColorDomainData(colorDomainsFileLine, None)

                    # Iterate through valid Domains and discard any that have been completely passed.
                    validDomainsForFeature = [validDomainData for validDomainData in validDomainsForFeature if not isACompletelyPastB(featureData, validDomainData)]

                    # Get ready to assign a color using the valid domains!
                    # Initialize the dictionary and set the minimum coverage level as GRAY.
                    encompassedBasesByColor = dict()
                    encompassedBasesByColor["GRAY"] = (featureData.endPos - featureData.startPos + 1)*minimumCoverage

                    for validDomainData in validDomainsForFeature:

                        # Check for non-overlap.
                        # The first assertion shouldn't be possible, but the second check IS technically possible if the previous feature 
                        # had a greater end pos, causing previously overlapping domains to now be ahead of the feature entirely.
                        assert not isACompletelyPastB(featureData, validDomainData), "Passed domain encountered.  This should not be possible..."
                        if isACompletelyPastB(validDomainData, featureData): continue
                        
                        # Determine range of overlap by finding the start and end of the overlap.
                        if validDomainData.startPos < featureData.startPos: overlapStartPos = featureData.startPos
                        else: overlapStartPos = validDomainData.startPos

                        if validDomainData.endPos > featureData.endPos: overlapEndPos = featureData.endPos
                        else: overlapEndPos = validDomainData.endPos

                        # Update the dictionary with the number of encompassed bases for the given color.
                        overlap=overlapEndPos-overlapStartPos+1
                        encompassedBasesByColor[validDomainData.color] = encompassedBasesByColor.setdefault(validDomainData.color, 0) + overlap

                    # Assign the majority color to the bin.
                    maxCoverage = max(encompassedBasesByColor.values())
                    maxColors = [key for key, value in encompassedBasesByColor.items() if value == maxCoverage]
                    if len(maxColors) > 1: featureColor = "GRAY"
                    else: featureColor = maxColors[0]

                    # Write the result to the output file.
                    coloredFeaturesFile.write(featureFileLine[:-1] + '\t' + featureColor + '\n')


def isACompletelyPastB(A: EncompassingData, B: EncompassingData):
    if A is None or B is None: return True
    else: return A.chromosome > A.chromosome or A.startPos > B.endPos


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=os.path.dirname(__file__))
    dialog.createFileSelector("Chromatin Domains File:", 0, ("Bed Files", ".bed"))

    binnerTypeDS = dialog.createDynamicSelector(1, 0)
    binnerTypeDS.initDropdownController("Bins are...", ("Regular", "Specific Ranges"))
    regularBinsDialog = binnerTypeDS.initDisplay("Regular", "Regular")
    specificRangeBinsDialog = binnerTypeDS.initDisplay("Specific Ranges", "Specific Ranges")

    regularBinsDialog.createFileSelector("Chromosome Sizes File:", 0, ("Text File",".txt"))
    regularBinsDialog.createDropdown("Bin Size (bp):", 1, 0, ("1000","10000","100000","1000000"))

    specificRangeBinsDialog.createFileSelector("Ranges to bin:", 0, ("Bed File", ".bed"))

    binnerTypeDS.initDisplayState()

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    if binnerTypeDS.getControllerVar() == "Regular":
        determineRegularBinColors(dialog.selections.getIndividualFilePaths()[0], dialog.selections.getIndividualFilePaths("Regular")[0],
                                  int(dialog.selections.getDropdownSelections("Regular")[0]))
    elif binnerTypeDS.getControllerVar() == "Specific Ranges":
        determineSpecifiedBinColors(dialog.selections.getIndividualFilePaths()[0],
                                    dialog.selections.getIndividualFilePaths("Specific Ranges")[0])

if __name__ == "__main__": main()