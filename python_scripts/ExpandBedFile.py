from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog
import os


# This function takes a bed file expands its coordinates to encompass extra bases on each side.
# Returns the file path to the expanded bed file.
def expandBedFile(baseBedFilePath: str, expansionRadius = 50):

    # Generate paths to the files needed to expand the binding motif.
    expandedBedFilePath = baseBedFilePath.rsplit('.',1)[0] + "_expanded.bed"

    # Expand the bed coordinates.
    print("Expanding nucleosome coordinates...")
    with open(baseBedFilePath,'r') as baseBedFile:
        with open(expandedBedFilePath, 'w') as expandedBedFile:

            # Write the expanded positions to the new file, one line at a time.
            for line in baseBedFile:
                choppedUpLine = line.strip().split('\t')
                choppedUpLine[1] = str(int(choppedUpLine[1]) - 50)
                choppedUpLine[2] = str(int(choppedUpLine[2]) + 50)

                # Write the results to the expansion file as long as it is not before the start of the chromosome.
                if int(choppedUpLine[1]) > -1: expandedBedFile.write('\t'.join(choppedUpLine) + '\n')
                else: print("Nucleosome at chromosome", choppedUpLine[0], "with expanded start pos", choppedUpLine[1],
                            "extends into invalid positions.  Skipping.")                                      

    return expandedBedFilePath

def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=os.path.dirname(__file__))
    dialog.createFileSelector("Bed File:", 0, ("Bed Files",".bed"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    expandBedFile(dialog.selections.getIndividualFilePaths()[0])

if __name__ == "__main__": main()