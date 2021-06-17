from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
import os


# Takes the txt file from the spivakov paper looks for any and all entries associated with the given transcription factors.
# Valid lines are converted to bed format.
def parseSpivakovToBed(spivakovFilePath: str, acceptableTFs = ("CTCF",)):

    # Create an output file path
    bedOutputFilePath = spivakovFilePath.rsplit('.',1)[0] + ".bed"

    with open(spivakovFilePath, 'r') as spivakovFile:
        with open(bedOutputFilePath, 'w') as bedOutputFile:

            # Read through each line, searching for valid transcription factors.
            # If any are found, write the line in bed format to the output file.
            # NOTE: The start and stop sites of the Spivakov file are 1-based
            for line in spivakovFile:
                choppedUpLine = line.split()

                for acceptableTF in acceptableTFs:
                    if acceptableTF in choppedUpLine[5]:
                        if choppedUpLine[3] == '1': strand = '+'
                        elif choppedUpLine[3] == "-1": strand = '-'
                        else: raise ValueError("Unexpected strand designation found: " + choppedUpLine[3])
                        bedOutputFile.write('\t'.join(("chr"+choppedUpLine[0], str(int(choppedUpLine[1])-1), choppedUpLine[2], '.', '.', strand)) + '\n')
                        continue


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=os.path.dirname(__file__))
    dialog.createFileSelector("Spivakov File:", 0, ("Text File",".txt"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    parseSpivakovToBed(dialog.selections.getIndividualFilePaths()[0])

if __name__ == "__main__": main()