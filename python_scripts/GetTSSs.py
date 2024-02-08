# Takes a file of gene designations and converts it to a file of single-nucleotide transcription start sites (TSSs)
import os
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog


def getTSSs(geneDesignationsFilePaths: List[str]) -> List[str]:
    """
    For each gene designations file path, creates a new file containing the single-base TSS for each gene.
    """

    tssFilePaths: List[str] = list()

    for geneDesignationsFilePath in geneDesignationsFilePaths:

        print(f"\nWorking in {os.path.basename(geneDesignationsFilePath)}...")

        # Create an output file path for the 
        if geneDesignationsFilePath.endswith("gene_designations.bed"):
            tssFilePath = geneDesignationsFilePath.rsplit("gene_designations.bed", 1)[0] + "TSSs.bed"
        else: tssFilePath = geneDesignationsFilePath.rsplit(".bed", 1)[0] + "_TSSs.bed"
        tssFilePaths.append(tssFilePath)

        # Create a set of TSS sites to make sure no duplicate entries slip through.
        TSS_Sites = set()

        # Iterate through the file, writing each TSS based on whether the gene is on the plus or minus strand.
        with open(geneDesignationsFilePath, 'r') as geneDesignationsFile, open(tssFilePath, 'w') as tssFile:

            for line in geneDesignationsFile:

                splitLine = line.strip().split('\t')

                if splitLine[5] == '+': splitLine[2] = str(int(splitLine[1])+1)
                elif splitLine[5] == '-': splitLine[1] = str(int(splitLine[2])-1)
                else: print("Warning: Found line without + or - strand designation. Skipping."); continue

                TSS_Site = (splitLine[0],splitLine[1],splitLine[2],splitLine[5])
                if TSS_Site in TSS_Sites: print(f"Warning: Found duplicate TSS Site: {TSS_Site}\nSkipping."); continue
                else: TSS_Sites.add(TSS_Site)

                tssFile.write('\t'.join(splitLine)+'\n')

    return tssFilePaths


def main():
    
    with TkinterDialog(workingDirectory = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..", "data"),
                       title = "Get TSSs") as dialog:
        dialog.createMultipleFileSelector("Gene designations Files", 1, "gene_designations.bed")

    getTSSs(dialog.selections.getFilePathGroups()[0])


if __name__ == "__main__": main()