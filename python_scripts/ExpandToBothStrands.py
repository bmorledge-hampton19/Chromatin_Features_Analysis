from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
import os


# Given a bed file, write a new bed file where each entry is duplicated, with the new entries being assigned to each strand.
def expandToBothStrands(bedFilePath: str):

    newBedFilePath = bedFilePath.rsplit('.',1)[0] + "_stranded.bed"

    with open(bedFilePath, 'r') as bedFile:
        with open(newBedFilePath, 'w') as newBedFile:

            for line in bedFile:
                chromosome, startPos, endPos = line.split()[:3]
                
                for strand in ('+','-'): 
                    newBedFile.write('\t'.join((chromosome, startPos, endPos, '.', '.', strand)) + '\n')



def main():

    # Create a simple dialog for selecting the gene designation files.
    dialog = TkinterDialog(workingDirectory=os.path.dirname(__file__))
    dialog.createFileSelector("Bed File:", 0, ("bed file", ".bed"))

    dialog.mainloop()

    if dialog.selections is None: quit()

    expandToBothStrands(dialog.selections.getIndividualFilePaths()[0],)


if __name__ == "__main__": main()