# This script takes a nucleosome positions file and the upper and lower quartile nucleosomes from H1 density data 
# (or other data, I suppose), and creates new nucleosome positioning files containing only those lower/upper quartile nucleosomes.
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulBioinformaticsFunctions import parseFastaDescription
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import checkDirs, getDataDirectory
import os, subprocess
from typing import List


# NOTE: If sloppyCopy is true, this function doesn't actually check to see if the nucleosome positions in the quartile files are present in the
#       base nucleosome file.  It just converts them to bed format.  This has the potential to cause pRoBlEmS.
def getQuartileNucleosomePositions(quartileFilePaths: List[str], nucPosDir: str, stratificationType, sloppyCopy):
    
    # If this isn't just a sloppy copy, create a dictionary containing each line in the root nucPos file for the corresponding location ID
    # Maybe this could cause memory issues, but I think it should be fine since the nucleosome maps are usually not too big.
    if not sloppyCopy:
        nucPosLines = dict()
        with open(os.path.join(nucPosDir,os.path.basename(nucPosDir)+".bed")) as nucPosFile:

            for line in nucPosFile:
                chromosome, startPos, endPos = line.split()[:3]
                key = (chromosome + ':' + str(float(startPos)) + '-' + str(float(endPos) - 1) + "(+)")
                nucPosLines[key] = line

    for quartileFilePath in quartileFilePaths:

        print("Working in",os.path.basename(quartileFilePath))

        if "lower_quartile" in os.path.basename(quartileFilePath):
            quartile = "lower_quartile"
        elif "upper_quartile" in os.path.basename(quartileFilePath):
            quartile = "upper_quartile"
        else:
            raise ValueError("Quartile Designation not found in file name: " + os.path.basename(quartileFilePath))

        # Create a name for the new nucleosome data based on the name of the quartile file and the parent nucleosome data.
        nucleosomeDataName = '_'.join((os.path.basename(nucPosDir),stratificationType,quartile))
        outputNucPosFilePath = os.path.join(os.path.dirname(nucPosDir), nucleosomeDataName, nucleosomeDataName + ".bed")
        checkDirs(os.path.dirname(outputNucPosFilePath))
        
        with open(quartileFilePath, 'r') as quartileFile:
            quartileFile.readline() # Get rid of headers
            with open(outputNucPosFilePath, 'w') as outputNucPosFile:

                # Parse out the location information and convert it to bed format.
                for line in quartileFile:

                    if sloppyCopy:
                        chromosome, startPos, endPos, strand = parseFastaDescription(line.split()[0])
                        outputNucPosFile.write('\t'.join((chromosome, startPos, str(float(endPos) + 1), '.', '.', strand)) + '\n')
                    else:
                        outputNucPosFile.write(nucPosLines[line.split()[0]])

        # Sort the output
        subprocess.run(("sort","-k1,1","-k2,3n",outputNucPosFilePath,"-o",outputNucPosFilePath), check = True)


def main():

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Quartile Files:",0, "quartile.tsv", ("Tab Separated Files",".tsv"))    
    dialog.createFileSelector("Nucleosome Directory:",1, directory=True)
    dialog.createDropdown("Stratification Type", 2, 0, ["h1 density", "other"])
    dialog.createCheckbox("Sloppy Copy?", 3, 0)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    getQuartileNucleosomePositions(dialog.selections.getFilePathGroups()[0], dialog.selections.getIndividualFilePaths()[0],
                                   dialog.selections.getDropdownSelections()[0].replace(' ', '_'), dialog.selections.getToggleStates()[0])


if __name__ == "__main__": main()