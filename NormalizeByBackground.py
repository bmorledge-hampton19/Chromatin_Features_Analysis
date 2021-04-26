from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
import os


# This function normalizes one or more raw counts files.
# If "headers" is true, the first line in the files is not used during normalization, and 
# the first line from the raw file is preserved in the normalized output file.
# "columnsToNormalize" describes which columns are normalized across the two files.
# All other columns are preserved in the state present in the raw counts file.
def normalizeByBackground(rawCountsFilePaths: str, backgroundCountsFilePath, headers = True, columnsToNormalize = (1,2)):

    # Iterate through the raw counts file paths, normalizing for each one.
    for rawCountsFilePath in rawCountsFilePaths:

        # Create the output file.
        normalizedFilePath = rawCountsFilePath.rsplit('.',1)[0] + "_normalized." + rawCountsFilePath.rsplit('.',1)[1]

        with open(rawCountsFilePath, 'r') as rawCountsFile:
            with open(backgroundCountsFilePath, 'r') as backgroundCountsFile:
                with open(normalizedFilePath, 'w') as normalizedFile:

                    # If present, trim the headers.
                    if headers:
                        normalizedFile.write(rawCountsFile.readline())
                        backgroundCountsFile.readline()

                    for rawLine in rawCountsFile:
                        backgroundLine = backgroundCountsFile.readline()

                        # Double check to make sure the background file hasn't ended prematurely.
                        assert backgroundLine, "Background file ended before raw counts file."

                        rawChoppedUpLine = rawLine.split()
                        backgroundChoppedUpLine = backgroundLine.split()

                        normalizedChoppedUpLine = list()
                        for i in range(len(rawChoppedUpLine)):
                            if i in columnsToNormalize:
                                # Add a pseudocount if necessary to prevent dividing by zero.
                                if backgroundChoppedUpLine[i] == "0": backgroundChoppedUpLine[i] = "1"
                                normalizedChoppedUpLine.append(str(int(rawChoppedUpLine[i])/int(backgroundChoppedUpLine[i])))
                            else: normalizedChoppedUpLine.append(rawChoppedUpLine[i])

                        normalizedFile.write('\t'.join(normalizedChoppedUpLine) + '\n')
                    
                    assert not backgroundCountsFile.readline(), "Raw file ended before background counts file."


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Raw Counts Files:",0, "binding_motif_mutation_counts.bed", ("Bed Files",".bed"), ("TSV files", ".tsv"))
    dialog.createFileSelector("Background Counts File:", 1, ("Bed Files",".bed"), ("TSV files", ".tsv"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    normalizeByBackground(dialog.selections.getFilePathGroups()[0], dialog.selections.getIndividualFilePaths()[0])

if __name__ == "__main__": main()