from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import DataTypeStr
from typing import List
import os


# This function takes a bed file of genome coordinates and bins them across each chromosome using the specified bin size.
# NOTE: input files must be sorted by chromosome ID (alphabetically) and feature start position (numerically).  Only the start position is used when binning.
def binAcrossGenome(genomeFeatureFilePaths: List[str], chromSizesFilePath, binSize):

    # Retrieve information on the sizes of the chromosomes being used.
    chromSizes = dict()
    with open(chromSizesFilePath, 'r') as chromSizesFile:
        for line in chromSizesFile:
            chromID, chromSize = line.split()
            chromSizes[chromID] = int(chromSize)

    for genomeFeatureFilePath in genomeFeatureFilePaths:

        print("Working in:", os.path.basename(genomeFeatureFilePath))

        # Generate an output file path
        binnedFeaturesFilePath = genomeFeatureFilePath.rsplit('.', 1)[0] + '_' + str(binSize) + "bp_binned.tsv"

        # Prepare for binning!
        bins = dict() # A nested dictionary.  The first key is for chromosome number, the second is for the start of the bin.
        with open(genomeFeatureFilePath, 'r') as genomeFeatureFile:

            # Read in the first line of the input file.
            choppedUpLine = genomeFeatureFile.readline().split()
            if not choppedUpLine: featureChrom = None
            else: 
                featureChrom = choppedUpLine[0]
                assert featureChrom in chromSizes, "Unrecognized chromosome: " + featureChrom
                featureStartPos = float(choppedUpLine[1])

            # Iterate through the chromosome, tracking how many features start within each bin.
            for binChrom in chromSizes:

                print("Binning in",binChrom)

                bins[binChrom] = dict()
                binStart = 0

                # Remain in the bin until the chromosomes don't match AND all bins have been initialized.
                while (featureChrom is not None and featureChrom == binChrom) or binStart < chromSizes[binChrom]:
                    bins[binChrom][binStart] = 0

                    # Are we in the current bin?  If so, increment this bin and read in a new feature.
                    while featureChrom is not None and featureChrom == binChrom and featureStartPos < binStart + binSize:

                        bins[binChrom][binStart] += 1

                        choppedUpLine = genomeFeatureFile.readline().split()
                        if not choppedUpLine: featureChrom = None
                        else: 
                            featureChrom = choppedUpLine[0]
                            featureStartPos = float(choppedUpLine[1])

                    # Increment the binStart.
                    binStart += binSize

                # We should never exit a bin chromosome while we still have a feature of that chromosome... Right?  *Sigh* Better double check...
                assert featureChrom != binChrom, ("Chromosome " + featureChrom + " bin exited before assigning feature starting at " + 
                                                  str(featureStartPos) + ".  Are the chrom.sizes incorrect?")
                # Also, make sure we recognize the new chromosome from the chrom.sizes file.
                assert featureChrom is None or featureChrom in chromSizes, "Unrecognized chromosome: " + featureChrom


        # Write the results of the binning.
        with open(binnedFeaturesFilePath, 'w') as binnedFeaturesFile:

            print("Writing results...")

            # Write headers
            binnedFeaturesFile.write('\t'.join(("Chromosome","Bin_Start-End","Feature_Counts")) + '\n')

            # Write the bins and feature counts.
            for chromosome in bins:
                for binStart in bins[chromosome]:

                    binnedFeaturesFile.write('\t'.join((chromosome, str(binStart)+'-'+str(binStart+binSize-1), str(bins[chromosome][binStart]))) + '\n')


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=os.path.dirname(__file__))
    dialog.createMultipleFileSelector("Genome Feature Files:", 0, DataTypeStr.mutations + ".bed", ("Bed Files", ".bed"))
    dialog.createFileSelector("Chromosome Sizes File:", 1, ("Text File",".txt"))
    dialog.createDropdown("Bin Size (bp):", 2, 0, ("1000","10000","100000","1000000"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    binAcrossGenome(dialog.selections.getFilePathGroups()[0], dialog.selections.getIndividualFilePaths()[0],
                    int(dialog.selections.getDropdownSelections()[0]))

if __name__ == "__main__": main()