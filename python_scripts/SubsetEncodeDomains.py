# This script subsets ENCODE chromatin domain files (e.g., to isolate euchromatin domains).
import os
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog

def subsetEncodeDomains(encodeDomainsFilePaths: List[str], subsetDomainsFilePaths: List[str]):

    for encodeDomainsFilePath in encodeDomainsFilePaths:

        print(f"\nWorking with {os.path.basename(encodeDomainsFilePath)}...")

        # Derive the basename for the given file path.
        outputDir = os.path.dirname(encodeDomainsFilePath)
        inputBasename = os.path.basename(encodeDomainsFilePath).rsplit('.', 1)[0]
        if inputBasename.endswith("chromatin_domains"):
            outputBasename = os.path.basename(encodeDomainsFilePath).rsplit("_chromatin_domains",1)[0]
        else: outputBasename = inputBasename

        for subsetDomainsFilePath in subsetDomainsFilePaths:

            print(f"Subsetting using {os.path.basename(subsetDomainsFilePath)}")

            # Create the output file path.
            outputFilePath = os.path.join(outputDir,outputBasename + '_' + os.path.basename(subsetDomainsFilePath).rsplit('.',1)[0] + ".bed")

            with open(encodeDomainsFilePath, 'r') as encodeDomainsFile,\
                 open(subsetDomainsFilePath, 'r') as subsetDomainsFile,\
                 open(outputFilePath, 'w') as outputFile:
                
                subsetDomains = [line.strip() for line in subsetDomainsFile]

                for line in encodeDomainsFile:
                    if line.split('\t')[3] in subsetDomains: outputFile.write(line)


def main():

    with TkinterDialog(workingDirectory=os.path.dirname(os.path.dirname(__file__)), title = "Subset Encode Domains") as dialog:
        dialog.createMultipleFileSelector("Encode Domains Files", 0, "chromatin_domains.bed", ("Bed files", ".bed"))
        dialog.createMultipleFileSelector("Domain Subsets Files", 1, "domains.txt", ("Text files", ".txt"))

    subsetEncodeDomains(dialog.selections.getFilePathGroups()[0], dialog.selections.getFilePathGroups()[1])

if __name__ == "__main__": main()