# This script takes a list of genes and counts the occurrence of features within the first 500 nucleotides of those genes,
# using the "count this in that" module.
import os
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
from benbiohelpers.CountThisInThat.InputDataStructures import EncompassingData
from benbiohelpers.CountThisInThat.Counter import ThisInThatCounter
from benbiohelpers.CountThisInThat.CounterOutputDataHandler import CounterOutputDataHandler
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from typing import List

class GeneStart(EncompassingData):

    def setLocationData(self, acceptableChromosomes):
        super().setLocationData(acceptableChromosomes)
        

class NucleosomeFeatureCounter(ThisInThatCounter):

    def __init__(self, encompassedFeaturesFilePath, encompassingFeaturesFilePath, 
                 outputFilePath, acceptableChromosomes = None, checkForSortedFiles = True,
                 headersInEncompassedFeatures = False, headersInEncompassingFeatures = False,
                 encompassingFeatureExtraRadius = 0, minEncompassedDistance = 0):
        super().__init__(encompassedFeaturesFilePath, encompassingFeaturesFilePath, 
                 outputFilePath, acceptableChromosomes, checkForSortedFiles,
                 headersInEncompassedFeatures, headersInEncompassingFeatures,
                 encompassingFeatureExtraRadius)

        self.minEncompassedDistance = minEncompassedDistance

    def setUpOutputDataHandler(self):
        self.outputDataHandler = CounterOutputDataHandler()
        self.outputDataHandler.addEncompassingFeatureStratifier(outputName = "Nucleosome")
        self.outputDataHandler.addPlaceholderStratifier("Feature_Counts")

    def constructEncompassingFeature(self, line) -> EncompassingDataDefaultStrand:
        return EncompassingDataDefaultStrand(line, self.acceptableChromosomes)

    def isEncompassedFeatureWithinEncompassingFeature(self, encompassedFeature = None, encompassingFeature = None):

        if encompassedFeature is None: encompassedFeature = self.currentEncompassedFeature
        if encompassingFeature is None: encompassingFeature = self.currentEncompassingFeature

        return (super().isEncompassedFeatureWithinEncompassingFeature(encompassedFeature, encompassingFeature) and
                abs(encompassingFeature.center - encompassedFeature.position) >= self.minEncompassedDistance)