# This script takes one or more bed files of transcription factor binding sites and a file of motif offsets and generates
# standardized bed files of single-nucleotide motif midpoints.
import os, subprocess
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.FileSystemHandling.DirectoryHandling import getTempDir
from benbiohelpers.FileSystemHandling.RemoveDuplicates import removeDuplicates
from benbiohelpers.FileSystemHandling.AddSequenceToBed import addSequenceToBed

def getTFBS_MidpointsFromOffsets(TFBS_FilePaths: List[str], offsetsFilePath: str, genomeFastaFilePath = None,
                                 retainSequence = True, removeDups = True):
    """
    Takes one or more bed files of transcription factor binding sites and a file of motif offsets and
    generates standaradized bed files of single-nucleotide motif midpoints.

    Transcription factor name is expected to be in the 7th column and will replace the 5th
    column so that the 7th column can (optionally) be used to store the original sequence.

    By default, duplicate entries (those with the same midpoint and transcription factor) will be reduced to a single entry.
    """

    # Generate a dictionary of offsets for the different motifs.
    offsets = dict()
    with open(offsetsFilePath, 'r') as offsetsFile:
        offsetsFile.readline() # Skip headers
        for line in offsetsFile:
            motif, offset = line.strip().split('\t')
            offsets[motif] = int(offset)

    # Use the offsets to find the midpoints of the given TFBSs.
    for TFBS_FilePath in TFBS_FilePaths:

        print(f"\nWorking with {os.path.basename(TFBS_FilePath)}...")

        # Create the output file path, as well as the necessary temporary intermediate paths.
        baseName = os.path.basename(TFBS_FilePath).rsplit(".bed",1)[0]
        TFBS_MidpointFilePath = os.path.join(os.path.dirname(TFBS_FilePath),baseName + "_midpoints.bed")

        reformattedTFBS_FilePath = os.path.join(getTempDir(TFBS_FilePath), baseName+"_reformatted.bed")
        preDedupTFBS_MidpointFilePath = os.path.join(getTempDir(TFBS_FilePath), baseName+"_midpoints_pre_dedup.bed")

        # First, reformat the original TFBS file by writing the TF name to the 5th column and then putting the sequence of the
        # TFBS in the 7th column (if requested).
        print("Moving TF name to 5th column...")
        with open(TFBS_FilePath, 'r') as TFBS_File, open(reformattedTFBS_FilePath, 'w') as reformattedTFBS_File:
            for line in TFBS_File:
                splitLine = line.strip().split('\t')
                splitLine[4] = splitLine[6]
                reformattedTFBS_File.write('\t'.join(splitLine[:6])+'\n')

        if retainSequence:
            print("Adding original motif sequence to file...")
            addSequenceToBed(reformattedTFBS_FilePath, genomeFastaFilePath)

        # Determine which output file path to use, based on whether or not we are deduplicating afterwards.
        if removeDups: outputFilePath = preDedupTFBS_MidpointFilePath
        else: outputFilePath = TFBS_MidpointFilePath

        # Calculate binding site midpoints using the given offsets.
        print("Calculating and writing midpoints...")
        missingOffsets = set()
        with open(reformattedTFBS_FilePath, 'r') as reformattedTFBS_File, open(outputFilePath, 'w') as outputFile:
            for line in reformattedTFBS_File:
                splitLine = line.strip().split('\t')

                # Make sure we actually have an offset for this motif. If not, skip it.
                if splitLine[3] not in offsets:
                    if splitLine[3] not in missingOffsets:
                        missingOffsets.add(splitLine[3])
                        print(f"No offset found for motif {splitLine[3]}. Skipping all related binding sites.")
                    continue

                # Calculate the preliminary midpoint.
                midpoint = (int(splitLine[1]) + int(splitLine[2])-1)/2 # 0-based

                # For binding sites with a half-base midpoint, round up if it is on the + strand and down if on the -.
                if int(midpoint) != midpoint:
                    if splitLine[5] == '+': midpoint += 0.5
                    else: midpoint -= 0.5

                # Add the offset (or the negative offset if on the minus strand)
                if splitLine[5] == '+': midpoint += offsets[splitLine[3]]
                else: midpoint -= offsets[splitLine[3]]

                # substitute the midpoint in the bed file.
                splitLine[1] = str(int(midpoint))
                splitLine[2] = str(int(midpoint + 1))

                outputFile.write('\t'.join(splitLine) + '\n')

        # Sort the result
        subprocess.check_call(("sort", "-k1,1", "-k2,2n", "-k3,3n", "-k5,5", "-k6,6",
                               "-s", "-o", outputFilePath, outputFilePath))

        # Remove duplicates, if requested.
        if removeDups:
            print("Removing duplicate entries based on TF name...")
            deduppedTFBS_MidpointFilePath = removeDuplicates([preDedupTFBS_MidpointFilePath], [0, 1, 2, 4, 5], verbose = False)[0]
            os.replace(deduppedTFBS_MidpointFilePath, TFBS_MidpointFilePath)


def main():
    with TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__), "..","data"), title = "Get TFBS Midpoints") as dialog:
        dialog.createMultipleFileSelector("TFBS files:", 0, "TFBS.bed", ("Bed files", "*.bed"))
        dialog.createFileSelector("Midpoints file:", 1, ("Tab-delimited file", ".tsv"))
        with dialog.createDynamicSelector(2, 0) as retainSequenceDS:
            retainSequenceDS.initCheckboxController("Retain binding site sequence")
            retainSequenceDS.initDisplay(True, "retainSequence").createGenomeSelector(0, 0)
        dialog.createCheckbox("Remove duplicates", 3, 0)

    if retainSequenceDS.getControllerVar(): genomeFastaFilePath = dialog.selections.getGenomes("retainSequence", "fasta")[0]
    else: genomeFastaFilePath = None

    getTFBS_MidpointsFromOffsets(dialog.selections.getFilePathGroups()[0], dialog.selections.getIndividualFilePaths()[0],
                                 genomeFastaFilePath, retainSequenceDS.getControllerVar(),
                                 dialog.selections.getToggleStates()[0])

if __name__ == "__main__": main()