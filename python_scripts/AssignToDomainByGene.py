# Uses gene names to assign rows of data (e.g. RPKM) to color domains.
import os
from pickle import TRUE
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog


def assignToDomainByGene(coloredGeneDesignationsFilePath: str, colorlessGeneDataFilePath: str,
                         geneIDindex = 0, omitGrayDomain = TRUE, addSecondaryID = True):

    # Create an output file path
    coloredGeneDataFilePath = colorlessGeneDataFilePath.rsplit('.',1)[0] + "_colored.tsv"

    # Read all the gene designations (and potentially secondary names) into a dictionary.
    # (This shouldn't cost too much memory, right?)
    print("Retrieving information from gene designations file...")
    geneColorsByID = dict()
    secondaryNameByID = dict()
    with open(coloredGeneDesignationsFilePath, 'r') as coloredGeneDesignationsFile:
        for line in coloredGeneDesignationsFile:
            splitLine = line.strip().split('\t')
            if splitLine[6] != "GRAY" or not omitGrayDomain:
                geneColorsByID[splitLine[3]] = splitLine[6]
                if addSecondaryID:
                    secondaryNameByID[splitLine[3]] = splitLine[4]

    
    # Use the created dictionary(ies) to generate the output file.
    print("Writing new information to colored gene data file.")
    with open(colorlessGeneDataFilePath, 'r') as colorlessGeneDataFile:
        with open(coloredGeneDataFilePath, 'w') as coloredGeneDataFile:

            # First, create the headers for the output file.
            headers = colorlessGeneDataFile.readline().strip().split('\t')
            headers.append("Color_Domain")
            if addSecondaryID: headers.append("Secondary_ID")
            coloredGeneDataFile.write('\t'.join(headers)+'\n')

            # Next, get additional data from the previously created dictionaries and add them
            # to all relevant rows.
            for line in colorlessGeneDataFile:

                splitLine = line.strip().split('\t')
                geneID = splitLine[geneIDindex]

                geneColor = geneColorsByID.get(geneID,None)
                if geneColor is None: continue
                else: splitLine.append(geneColor)

                if addSecondaryID: splitLine.append(secondaryNameByID[geneID])

                coloredGeneDataFile.write('\t'.join(splitLine)+'\n')
                

def main():

    dialog = TkinterDialog(workingDirectory=os.path.dirname(__file__))
    dialog.createFileSelector("Colored Gene Designations:", 0, ("Bed File",".bed"))
    dialog.createFileSelector("Colorless Gene Data (e.g. RPKM)", 1, ("Tab separated files",".tsv"))

    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    assignToDomainByGene(dialog.selections.getIndividualFilePaths()[0],
                         dialog.selections.getIndividualFilePaths()[1])

if __name__ == "__main__": main()