# Chromatin_Features_Analysis
 A suite of scripts for analyzing DNA damage/repair/mutagenesis with respect to chromatin features like transcription factors and histone modifications
***

## Table of Contents
1. [Project Overview](#project-overview)
2. [Dependencies and File System Naming Conventions](#dependencies-and-file-system-naming-conventions)
3. [Preparing Data for Analysis](#preparing-data-for-analysis)
4. [Analysis and Figure Generation](#analysis-and-figure-generation)
5. [Data Availability](#data-availability)
***

## Project Overview

#### Accompanying Literature
The code in this repository was primarily used to produce the findings in [NO LINK YET](). This paper contains additional insight into the analysis using this code.

#### Background
f
***

## Dependencies and File System Naming Conventions

#### A Disclaimer...
Admittedly, the software in this repository was largely not designed with other users in mind. Despite this, the core analysis is still very much reproducible after taking a few steps to set up a specific environment and following a few conventions pertaining to the project's file system structure. In the event that the following steps are insufficient to reproduce the desired results, you are more than welcome to [post an issue](../../issues) or email me directly at b.morledge-hampton@wsu.edu.

#### Dependencies
Besides the code in this repository, the analysis requires that the following packages are installed (The most up-to-date version of each package is recommended):
- Python:
  - [benbiohelpers](https://github.com/bmorledge-hampton19/benbiohelpers)
  - [mutperiod](https://github.com/bmorledge-hampton19/mutperiod)
- R:
  - [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
  - [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
- Must be available through command line (Recommend installing through apt, where possible):
  - [bedtools](https://bedtools.readthedocs.io/en/latest/)
  - [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  - [samtools](http://www.htslib.org/)
  - [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
  
Although not required, the scripts in [this repository](https://github.com/bmorledge-hampton19/XR-seq_Analysis) may be helpful for [trimming](https://github.com/bmorledge-hampton19/XR-seq_Analysis/blob/main/TrimAdaptorSequences.py) and [aligning](https://github.com/bmorledge-hampton19/XR-seq_Analysis/blob/main/TrimmedFastqToSam.py) sequencing reads and [parsing them to bed](https://github.com/bmorledge-hampton19/XR-seq_Analysis/blob/main/SamToBed.py). If desired, these three operations can be chained together using the [AlignXRSeqReads.py](https://github.com/bmorledge-hampton19/XR-seq_Analysis/blob/main/AlignXRSeqReads.py) script.

#### Directory structure naming conventions
After cloning the deamination_determination repository, you will need to create a directory at the top level called "data". Within this data directory, 

This data can be found in the [Data Availability](#data-availability) section

When running the R notebooks contained in this repository, they will ask you to provide a "bioinformatics directory". This is merely the directory containing the cloned repository (i.e. one directory up from "deamination_determination/"). Once this directory is selected, the underlying file system is inferred using the naming conventions in this section.

#### Data file naming conventions
When obtaining fastq sequencing files (or sam alignment files), a strict file naming system needs to be maintained in order for the rest of the pipeline to run without errors. The file names should be made up of the following identifying information, separated by underscores:
- Organism or cell type ID
- Lesion type ("CPD", "6-4", or "cisplatin")
- Timepoint (in shorthand format, such as "10min" for minutes or "1h" for hours)
- Repitition number (In the format "rep#" where '#' is the repitition number or "all_reps" when combined)

Here is an example of a valid sam file name: "Dm_CPD_1h_rep1.sam"

#### Periodicity Analysis
f
***

## Preparing Data for Analysis

#### ???
f
***

## Analysis and Figure Generation
If the data has been prepared correctly, running the analyses and generating figures is easy: Simply run/knit the [R notebook](R/) in the repository. Also, when running the notebook, keep in mind the following details:

#### Defining the "Bioinformatics_Projects" Directory
The notebook requires a path to the "Bioinformatics_Projects" directory. (More information on this directory is available at the end of the [directory structure naming conventions section](#directory-structure-naming-conventions)). In brief, this is the parent directory to the cloned chromatin_features_analysis repository. If you are running the R notebooks in interactive mode, you will be able to select this directory through a file browser dialog. If you are knitting the notebook, you will select this directory from a list of choices in the Shiny UI for the notebook's parameters. If the desired directory is not listed, you will need to add it manually to the R notebook (generally on line ????).

#### Figure Output
When knitting figure generation notebooks, the figures are automatically output to the R/????? directory, which will be created if it doesn't already exist.
***

## Data Availability
The sequencing data used for the analyses supported by this repository can be found [here](NO-LINK-YET)
