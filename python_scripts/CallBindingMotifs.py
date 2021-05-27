import MOODS.parsers
import MOODS.tools
import MOODS.scan

from nucperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog
from nucperiodpy.helper_scripts.UsefulBioinformaticsFunctions import bedToFasta, FastaFileIterator

import os
from typing import List


# This function uses the MOODS package to identify the centers of binding motifs among ChIp-seq peak calls.
def callBindingMotifs(peakRegionFilePaths: List[str], genomeFilePath, bindingMotifFilePath):

    # Pre-process the binding motif.
    bg = MOODS.tools.flat_bg(4)
    pseudocount = 0.0001
    motifMatrix = MOODS.parsers.pfm_to_log_odds(bindingMotifFilePath, bg, pseudocount)
    threshold = MOODS.tools.threshold_from_p(motifMatrix, bg, 0.0001)
    
    # Get the size of the binding motif.
    with open(bindingMotifFilePath, 'r') as bindingMotifFile:
        bindingMotifFile.readline() # Skip the header line.
        bindingMotifLength = len(bindingMotifFile.readline().split())

    for peakRegionFilePath in peakRegionFilePaths:

        print()
        print("Working in", os.path.basename(peakRegionFilePath))

        # Generate a semi-intelligent name for the resulting motif centers file.
        motifCentersFilePath = peakRegionFilePath.rsplit('.',1)[0]
        if motifCentersFilePath.endswith("peak_regions"):
            motifCentersFilePath = motifCentersFilePath.rsplit("_peak_regions",1)[0] 
        motifCentersFilePath += "_" + os.path.basename(bindingMotifFilePath).rsplit('.',1)[0] + "_binding_motifs.bed"

        # Generate a file path for a corresponding fasta file to convert the bed file to.
        peakRegionSequencesFilePath = peakRegionFilePath.rsplit('.',1)[0] + ".fa"

        # Write some metadata
        metadataFilePath = os.path.join(os.path.dirname(peakRegionFilePath), ".metadata")
        with open(metadataFilePath, 'w') as metadataFile:
            metadataFile.write("Binding_Motif_File_Path: " + bindingMotifFilePath + '\n')
            metadataFile.write("Genome_File_Path: " + genomeFilePath + '\n')

        # Generate the fasta file if it doesn't already exist.
        if os.path.exists(peakRegionSequencesFilePath):
            print("Fasta file already exists.")
        else:
            print("Fasta file not found.  Generating...")
            bedToFasta(peakRegionFilePath, genomeFilePath, peakRegionSequencesFilePath)

        # Scan for the motif in all the given DNA sequences.
        with open(peakRegionSequencesFilePath, 'r') as peakRegionSequencesFile:
            with open(motifCentersFilePath, 'w') as motifCentersFile:
                for fastaEntry in FastaFileIterator(peakRegionSequencesFile):
                    results = MOODS.scan.scan_dna(fastaEntry.sequence, (motifMatrix,), bg, (threshold,), 7)

                    # Output the results to the motif centers file.
                    for rs in results:
                        for r in rs: 

                            # Determine the start and end positions depending on the strand the motif was found on.
                            if fastaEntry.strand == '+':
                                startPos = str(r.pos + int(fastaEntry.startPos))
                                endPos = str(r.pos + int(fastaEntry.startPos) + bindingMotifLength)
                            else:
                                startPos = str(int(fastaEntry.endPos) - r.pos - bindingMotifLength)
                                endPos = str(int(fastaEntry.endPos) - r.pos)


                            motifCentersFile.write('\t'.join((fastaEntry.chromosome, startPos, endPos, '.', str(r.score), fastaEntry.strand)) + '\n')


def main():

    # Create a simple dialog for selecting the gene designation files.
    dialog = TkinterDialog(workingDirectory=os.path.dirname(__file__))
    dialog.createMultipleFileSelector("Peak Region Bed Files:", 0, "peak_regions.bed", 
                                      ("Bed Files", ".bed"))
    dialog.createFileSelector("Genome Fasta File:", 1, ("fasta File", ".fa"))
    dialog.createFileSelector("Binding Motif File:", 2, ("pfm File", ".pfm"))

    dialog.mainloop()

    if dialog.selections is None: quit()

    callBindingMotifs(dialog.selections.getFilePathGroups()[0], dialog.selections.getIndividualFilePaths()[0],
                      dialog.selections.getIndividualFilePaths()[1])


if __name__ == "__main__": main()