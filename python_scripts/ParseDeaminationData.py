# This script takes the CPD and deamination data files, validates them, and converts them to a more standard bed format
import os
from typing import List

from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.FileSystemHandling.AddSequenceToBed import addSequenceToBed
from benbiohelpers.DNA_SequenceHandling import isPurine
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory, getAcceptableChromosomes


def parseDeaminationData(cPDFilePaths: List[str], deaminationFilePaths: List[str], genomeFastaFilePath):
    """
    See script header.
    """

    acceptableChromosomes = getAcceptableChromosomes(genomeFastaFilePath)

    for cPDFilePath in cPDFilePaths:

        print("\nWorking in",os.path.basename(cPDFilePath))

        # Create a name for the parsed output file.
        cPDParsedFilePath = cPDFilePath.rsplit('.',1)[0] + "_parsed.bed"

        # Create a path to the output file with only cytosine positions.
        cPDCytosinePositionsFilePath = cPDFilePath.rsplit('.',1)[0] + "_cytosines.bed"

        with open(cPDFilePath, 'r') as cPDFile:
            with open (cPDParsedFilePath, 'w') as cPDParsedFile:

                print("Parsing original file...")

                cPDFile.readline() # Skip the header line.

                for line in cPDFile:

                    choppedUpLine = line.split()
                    
                    # Record position information, inferring strand from the given gap sequencing value 
                    # and extending positions to encompass the full CPD sequence.
                    chromosome = choppedUpLine[0]

                    # Check for invalid chromosomes
                    if chromosome not in acceptableChromosomes: 
                        print("Skipping invalid chromosome:", chromosome)
                        continue 

                    if isPurine(choppedUpLine[3]):
                        strand = '-'
                        position0 = str( int(choppedUpLine[1]) - 1 )
                        position1 = choppedUpLine[2]
                    else:
                        strand = '+'
                        position0 = choppedUpLine[1]
                        position1 = str( int(choppedUpLine[2]) + 1 )

                    # Make sure there is aggreement among the weird sequence columns as to the CPD sequence
                    assert choppedUpLine[5][-2:] == choppedUpLine[6][-2:] == choppedUpLine[7][:2] == choppedUpLine[8][:2], line
                    cPD = choppedUpLine[5][-2:]

                    cPDParsedFile.write('\t'.join((chromosome, position0, position1, cPD, '.', strand)) + '\n')

        # Derive the sequences directly from the positions and make sure it matches the cPD value obtained previously.
        # At the same time, create the file with only the cytosine positions in CPDs.
        print("Deriving sequence context from given fasta file...")
        addSequenceToBed(cPDParsedFilePath, genomeFastaFilePath, 4)
        with open(cPDParsedFilePath, 'r') as cPDParsedFile:
            with open(cPDCytosinePositionsFilePath, 'w') as cPDCytosinePositionsFile:

                print("Validating original file and trimming to single base cytosine positions...")

                for line in cPDParsedFile:

                    choppedUpLine = line.split()
                    assert choppedUpLine[3] == choppedUpLine[4], line

                    for i, base in enumerate(choppedUpLine[3]):
                        if base == 'C':

                            if choppedUpLine[5] == '+':
                                position0 = str( int(choppedUpLine[1]) + i )
                            else:
                                position0 = str( int(choppedUpLine[1]) + 1 - i )
                            position1 = str( int(position0) + 1)

                            cPDCytosinePositionsFile.write('\t'.join((choppedUpLine[0], position0, position1, '.', '.', choppedUpLine[5])) + '\n')


    for deaminationFilePath in deaminationFilePaths:

        print("\nWorking in",os.path.basename(deaminationFilePath))

        # Create a name for the parsed output file.
        deaminationParsedFilePath = deaminationFilePath.rsplit('.',1)[0] + "_parsed.bed"

        # Create a path to the output file with only cytosine positions in dipy contexts.
        dipyDeaminationPositionsFilePath = deaminationFilePath.rsplit('.',1)[0] + "_dipy_cytosines.bed"

        with open(deaminationFilePath, 'r') as deaminationFile:
            with open(deaminationParsedFilePath, 'w') as deaminationParsedFile:

                print("Parsing original file...")

                deaminationFile.readline() # Skip the header line.

                for line in deaminationFile:

                    choppedUpLine = line.split()

                    # Check for non-cytosine positions, and ignore them if found.
                    if choppedUpLine[3] in ('A','T'): continue

                    # Record position information for all other rows.  Expand to trinucleotide context.
                    chromosome = choppedUpLine[0]

                    # Check for invalid chromosomes
                    if chromosome not in acceptableChromosomes: 
                        print("Skipping invalid chromosome:", chromosome)
                        continue 

                    position0 = str( int(choppedUpLine[1]) - 1 )
                    position1 = str( int(choppedUpLine[2]) + 1 )
                    if choppedUpLine[3] == 'C': strand = '+'
                    else: strand = '-'
                    trinuc = choppedUpLine[5]

                    deaminationParsedFile.write('\t'.join((chromosome, position0, position1, trinuc, '.', strand)) + '\n')

        # Derive the trinucleotide context sequences directly from the positions and make sure it matches the sequence obtained from the file.
        # At the same time, create the file with only the cytosine positions with an adjacent pyrimidine.
        print("Deriving sequence context from given fasta file...")
        addSequenceToBed(deaminationParsedFilePath, genomeFastaFilePath, substitutionPosition = 4)
        with open(deaminationParsedFilePath, 'r') as deaminationParsedFile:
            with open(dipyDeaminationPositionsFilePath, 'w') as dipyDeaminationPositionsFile:

                print("Validating original file and trimming cytosines without adjacent pyrimidines...")

                for line in deaminationParsedFile:

                    choppedUpLine = line.split()
                    
                    assert choppedUpLine[3] == choppedUpLine[4], line

                    if isPurine(choppedUpLine[3][0]) and isPurine(choppedUpLine[3][2]): continue
                    else:

                        position0 = str( int(choppedUpLine[1]) + 1 )
                        position1 = str( int(choppedUpLine[2]) - 1 )

                        dipyDeaminationPositionsFile.write('\t'.join((choppedUpLine[0], position0, position1, '.', '.', choppedUpLine[5])) + '\n')


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("CPD Files:",0,"CPD_data.bed",("Bed Files",".bed"))
    dialog.createMultipleFileSelector("Deamination Files:", 1,"deamination_data.bed", ("Bed Files",".bed"))
    dialog.createFileSelector("Genome Fasta File:", 2, ("Fasta File",".fa"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    parseDeaminationData(dialog.selections.getFilePathGroups()[0], dialog.selections.getFilePathGroups()[1],
                         dialog.selections.getIndividualFilePaths()[0])


if __name__ == "__main__": main()