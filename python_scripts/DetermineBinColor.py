from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from typing import List
import os


# This function takes a bed file of chromatin color domains and assigns a color to each bin, defaulting to gray if no domain covers that area.
# Bins are colored based on the majority color present in that region.
# NOTE: input files must be sorted by chromosome ID (alphabetically) and feature start position (numerically).
def binAcrossGenome(colorDomainsFilePath, chromSizesFilePath, binSize, minimumCoverage = 0.5):

    # Retrieve information on the sizes of the chromosomes being used.
    chromSizes = dict()
    with open(chromSizesFilePath, 'r') as chromSizesFile:
        for line in chromSizesFile:
            chromID, chromSize = line.split()
            chromSizes[chromID] = int(chromSize)

    print("Working in:", os.path.basename(colorDomainsFilePath))

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

        # Iterate through the chromosome, tracking how many features start within each bin.
        for binChrom in chromSizes:

            print("Binning in",binChrom)

            bins[binChrom] = dict()
            binStart = 0

            # Remain in the bin until the chromosomes don't match AND all bins have been initialized.
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

                    # Read in the next chromosome, UNLESS this feature extends into the next bin.  In that case, continue onto the next bin.
                    # NOTE: This means that this code doesn't handle overlapping regions very well, so I'm just assuming they don't occur very frequently.
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


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=os.path.dirname(__file__))
    dialog.createFileSelector("Chromatin Domains File:", 0, ("Bed Files", ".bed"))
    dialog.createFileSelector("Chromosome Sizes File:", 1, ("Text File",".txt"))
    dialog.createDropdown("Bin Size (bp):", 2, 0, ("1000","10000","100000","1000000"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    binAcrossGenome(dialog.selections.getIndividualFilePaths()[0], dialog.selections.getIndividualFilePaths()[1],
                    int(dialog.selections.getDropdownSelections()[0]))

if __name__ == "__main__": main()