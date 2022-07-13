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
Much of this analysis builds on [previous work](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8349767/) examining mutation rates relative to nucleosome positions. This project utilizes the software from that previous work in tandem with a suite of new scripts to examine DNA damage and repair data relative to nucleosomes in Drosophila as well as stratifying that data across [five distinct chromatin domains](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3119929/). As a result, the python scripts contained in this repository are directed at organizing genomic features such as nucleosomes with respect to these chromatin domains in preparation for comparing them against DNA damage and repair rates. The R scripts and notebook in this repository are then used to format and plot the data.
***

## Dependencies and File System Naming Conventions

#### A Disclaimer...
Admittedly, the software in this repository was largely not designed with other users in mind. This documentation is really meant more for the author of the software (me!) than anyone. Despite this, the core analysis SHOULD still be reproducible following the instructions below. In the event that the following steps are insufficient to reproduce the desired results, you are more than welcome to [post an issue](../../issues) or email me directly at b.morledge-hampton@wsu.edu.

#### Dependencies
Besides the code in this repository, the analysis requires that the following packages are installed (The most up-to-date version of each package is recommended):
- Python:
  - [mutperiod](https://github.com/bmorledge-hampton19/mutperiod)
  - [benbiohelpers](https://github.com/bmorledge-hampton19/benbiohelpers) (Note that this should install automatically with mutperiod).
- R:
  - [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
  - [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
- Must be available through command line (Recommend installing through apt, where possible):
  - [bedtools](https://bedtools.readthedocs.io/en/latest/)
  - [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  - [samtools](http://www.htslib.org/)
  - [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
 
##### Optional Dependencies
Although not required, the scripts in [this repository](https://github.com/bmorledge-hampton19/xrlesionfinder/tree/main/xrlesionfinder) may be helpful for [trimming](https://github.com/bmorledge-hampton19/xrlesionfinder/blob/main/xrlesionfinder/AlignmentAndFormatting/TrimAdaptorSequences.py) and [aligning](https://github.com/bmorledge-hampton19/xrlesionfinder/blob/main/xrlesionfinder/AlignmentAndFormatting/TrimmedFastqToSam.py) sequencing reads and [parsing them to bed](https://github.com/bmorledge-hampton19/xrlesionfinder/blob/main/xrlesionfinder/AlignmentAndFormatting/SamToBed.py). If desired, these three operations can be chained together using the [AlignXRSeqReads.py](https://github.com/bmorledge-hampton19/xrlesionfinder/blob/main/xrlesionfinder/AlignmentAndFormatting/AlignXRSeqReads.py) script.

#### File System Naming Conventions
This repository's analysis depends on strict directory and file path naming conventions. These are outlined below and are later expanded upon in the [Preparing Data for Analysis](#preparing-data-for-analysis) section:

##### This Repository's Data Directory
After cloning the chromatin_features_analysis repository, you will need to create a directory at the top level called "data". Within this data directory, two directories must be present:
- __RNAseq/__ This directory will contain RPKM values for later analysis.
- __chromatin_domains/drosophila/__  This directory will contain bed-formatted files of regions labeled as one of the five colored chromatin domains.

##### The bioinformatics directory
When running the R notebook contained in this repository, it expects data to be contained within a specific bioinformatics directory. This is merely the parent directory containing this cloned repository as well as data from mutperiod (i.e. one directory up from chromatin_features_analysis/). Typically, this directory is called "Bioinformatics_Projects/". Once this directory is selected, the underlying file system is inferred using the naming conventions in this section.

##### The mutperiod_data directory
Data generated by mutperiod is stored in its own directory independent of the chromatin_features_analysis project directory. This directory is created the first time mutperiod requests file input. If possible, this directory should be created within Bioinformatics_Projects/mutperiod/. (The directory that will be created is Bioinformatics_Projects/mutperiod/mutperiod_data/.) Next, a drosophila_data/ directory should be created under mutperiod_data/ and the following four directories should be created under drosophila_data/:
- __damage/__ This directory will contain nucleosome periodicity results for cellular DNA CPD-seq data.
- __naked/__ This directory will contain nucleosome periodicity results for naked DNA CPD-seq data.
- __nucleosomes/__ This directory will contain nucleosome periodicity (repeat length) data.
- __repair/__ This directory will contain nucleosome periodicity results for XR-seq data.

Finally, under the mutperiod_data/\_\_external_data\ directory (which should have been created automatically), create a directory named "dm6/". This is where you will store the dm6 genome assembly and various other data files relating to it, such as nucleosome maps.

#### Data file naming conventions
When obtaining fastq sequencing files (or sam alignment files), a strict file naming system needs to be maintained in order for the rest of the pipeline to run without errors. I'll apologize in advance for the inconsistencies between these conventions. At this point, it just is what it is.


For XR-seq files, the file names should be made up of the following identifying information, separated by underscores:
- Organism ID. (For this analysis, this will always be "Dm").
- Lesion type (For this analysis, this will always be "CPD")
- Timepoint (in shorthand format, such as "10m" for minutes or "1h" for hours)
- Repetition number (In the format "rep#" where '#' is the repetition number or "all" when combined)

Here is an example of a valid XR-seq file name: "Dm_CPD_1h_all.fastq"


For CPD-seq files, the cellular damage files should be constructed in the following manner: "WT_0hr_UV_\[REPETITION\]" where \[REPETITION\] is replaced by "rep1", "rep2", or "aggregate". Meanwhile, naked DNA damage files should be constructed as such: "WT_naked_DNA" with no indicator of repetition.

#### Deviating from naming conventions
It is still possible to run the analysis without following these instructions to the letter. This structure is imposed so that the R notebook within this repository can be run with minimal user input. However, different naming conventions can be used as long as the corresponding names within the notebook are changed.

***

## Preparing Data for Analysis
This section details the necessary steps to prepare the data for analysis and figure generation. It is admittedly an arduous, specific, and unforgiving process, as no streamlined pipeline was ever created to accommodate it. Despite this, that process is documented below in the event that it needs to be replicated.

#### Acquiring Data
The following data are necessary to run this analysis and often need to be incorporated into the file system in a specific manner. The data files themselves can be found at the links in the [data availability](#data-availability) section.
- __The dm6 Genome__: This should be placed in the mutperiod_data/\_\_external_data/dm6 directory and should be called "dm6.fa".
- __dm6 chromosome sizes__: File location and name are irrelevant.
- __CPD-seq data sets__: These files should be accessible from the [mutperiod_data](#the-mutperiod_data-directory) directory, and their names should follow the [proper naming conventions](#data-file-naming-conventions).
- __XR-seq sequencing data__: These files should be downloaded in fastq form and be accessible from the [mutperiod_data](#the-mutperiod_data-directory) directory. Their names should follow the [proper naming conventions](#data-file-naming-conventions). Optionally, [this script](https://github.com/bmorledge-hampton19/xrlesionfinder/blob/main/xrlesionfinder/AlignmentAndFormatting/SRA_ToFastq.py) can be used to download the sequencing data more efficiently by providing the relevant SRA accession IDs in a newline-separated text format and optionally providing custom names (in the same format) for the resulting fastq files. Note that running this script requires that the [SRA Toolkit](https://github.com/ncbi/sra-tools) is installed and accessible through your PATH environment variable.
- __Nucleosome Map__: A directory for this file should be created under mutperiod_data/\_\_external_data/dm6 with the name "S2_nucleosome_map". The map itself should be named "S2_nucleosome_map.bed". If desired, the map can be recreated from the MNase data linked below.
- __Chromatin Domain Designations__: This file should be placed within [this repository's data directory](#this-repository's-data-directory) under chromatin_domains/drosophila/ and should be called drosophila_chromatin_color_domains.bed.
- __Gene Region Designations__: A directory for this file should be created under mutperiod_data/\_\_external_data/dm6 with the name "gene_designations". The file itself should be named "dm6_BDGP6_89_all_genes.bed".
- __RPKM Values__: This file should be placed within [this repository's data directory](#this-repository's-data-directory) under RNAseq/ and should be called drosophila_S2_gene_RPKM.tsv (it will have to be converted from the xls format it is hosted as).


#### Coloring Genome-Wide Bins.
The current Drosophila chromatin color domains file contains irreguluarly spaced regions assigned to different chromatin domains and doesn't fully cover the Drosophila genome. However, later steps in analysis will expect chromatin domain designations in regularly spaced, genome-wide bins. The transformation to this format can be accomplished using the [DetermineBinColor](https://github.com/bmorledge-hampton19/Chromatin_Features_Analysis/blob/main/python_scripts/DetermineBinColor.py) script. When running the script, make sure that the "Bins are..." option is set to "regular" and that the relevant chromatin domain designations file and chromosome sizes file are given. Furthermore, the chromatin domain designations input file must be sorted, first by chromosome ID (first column) and then by start position (second column), and the chromosome sizes file must be sorted by chromosome ID. Later analysis expects that colored bins of 10,000 and 100,000 bp are present, so the script needs to run twice, once for each bin size. All output files will automatically be written in the necessary location and named appropriately.

#### Filtering, "Coloring", and Expanding Gene Designations
For this study, our analysis was limited to better characterized genes, as defined by the presence of a common name or annotation symbol other than the FlyBase ID in the gene designations file. This can be found in the fifth column (index 4) of the file if it is present (The column will be blank if no common name is present). This repository does not contain a script for filtering out these rows, but this can be easily done using any matrix-like data structure that allows for stratification of rows, such as awk for command line or pandas for python. The resulting file should be given the "\_named_genes" suffix (i.e. dm6_BDGP6_89_named_genes.bed).

Once the file has been filtered, the genes can be "colored" using the color domain designations. Similarly to how genome-wide bins were colored in the previous step, the [DetermineBinColor](https://github.com/bmorledge-hampton19/Chromatin_Features_Analysis/blob/main/python_scripts/DetermineBinColor.py) script can be used to color genes, provided that both input files are sorted by chromosome ID first and start position second. Ensure that the "Bins are..." option is set to "Specific Ranges" and the relevant files are provided.

Lastly, the genomic positions in the colored gene designations file need to be expanded for the sake of a later step which examines flanking regions. The [ExpandBedFile](https://github.com/bmorledge-hampton19/Chromatin_Features_Analysis/blob/main/python_scripts/ExpandBedFile.py) script can achieve this easily. When run, set the "Expansion Radius" option to "1068". (This is half the median gene length from the filtered genes file.)

#### Stratifying and Subsetting Nucleosome Maps
The analysis supported by this repository examines nucleosomes with respect to different chromatin domains. In order to achieve this, the original nucleosome map needs to be stratified by those chromatin domains. This can be achieved using the [SeparateByChromatinRegions](https://github.com/bmorledge-hampton19/Chromatin_Features_Analysis/blob/main/python_scripts/SeparateByChromatinRegions.py) script along with some manual file system management. This script should be run with the S2 nucleosome map and chromatin domain files. (A warning may be printed to the terminal about using unexpected input. This can be ignored.) After running, the script will produce a folder named after the chromatin domain file, containing the stratified nucleosome maps. This folder will be present in the directory containing the original nucleosome map, but its contents will need to be moved to directories of their own, as mutperiod expects one nucleosome map per directory. These new directories will also need to be named after their respective files. (E.g. the S2_nucleosome_map_BLACK_domain.bed file should be moved to mutperiod_data/\_\_external_data/dm6/S2_nucleosome_map_BLACK_domain)

Later steps in the analysis expect random subsets of the larger stratified nucleosome maps (BLACK, BLUE, and YELLOW) that are equal in size to the smaller stratified nucleosome maps (GREEN and RED). This repository does not contain scripts to perform this subsetting, but any program which performs random sampling without replacement should be fine. For this analysis, the "shuf" command line tool was used. Each of the three nucleosome maps listed above needs to be subsetted into a 47843 row sample and a 21857 row sample. The former should be named with the suffix "\_subset" and the latter with the suffix "\_mini_subset". (E.g. S2_nucleosome_map_BLACK_domain_mini_subset.bed") As with the previous step, each of these new nucleosome maps needs its own appropriately named directory within the dm6/ directory.

For the sake of keeping file names short but distinct, mutperiod needs to be told specifically how to name files associated with each nucleosome map. This is done by including an "append_to_data_name.txt" file in the relevant nucleosome map directory. For this analysis, the original S2 nucleosome map represents the default map, and does not need this naming file, but all of its derivative maps do. Every other nucleosome map directory needs to define the text it will append to data names. This text should be the suffix that sets that nucleosome map apart from the default map and should be given as a single line in the append_to_data_name.txt file. (E.g. the "S2_nucleosome_map_BLUE_domain_subset" directory will need to contain an "append_to_data_name.txt" file with a single line of text that reads "\_BLUE_domain_subset".)

#### Combining Data Repetitions
Later steps in the analysis expect data with the "all" or "aggregate" identifiers, representing data which has combined multiple repetitions for the same timepoint. It is a good idea to combine these repitions now, before setting up the directories for mutperiod. The [CombineReps.py](https://github.com/bmorledge-hampton19/benbiohelpers/blob/main/python/benbiohelpers/FileSystemHandling/CombineReps.py) script in benbiohelpers can be used to accomplish this. Just make sure to change the "Combined Repetition String" parameter to just "all" for XR-seq data and "aggregate" for cellular CPD-seq data. Also note that this script assumes that the 2 repetitions are present in the same directory, so you may need to reorganize how the files are structured at this point. It may be easiest to do this in tandem with the next step.

#### Preparing Directories For mutperiod
Later steps in the analysis use mutperiod, which expects exactly one directory per data set. So for each data set, a directory needs to be created beneath the mutperiod_data/ directory, following the same [naming conventions](#data-file-naming-conventions) as the accompanying data set. (Remember, only the "all" files containing both repetitions for a given timepoint are required for the following steps in this analysis). These directories can be present within a subdirectory of the mutperiod_data/ directory. It is recommended that these directories are localized to the drosophila_data/ directory, or some other subdirectory to isolate them from other projects in mutperiod. Keep in mind that any data stratification not directly supported by mutperiod will have to be accompanied by the creation of additional directories. This will become relevant when [splitting across genic and intergenic positions](#split-across-genic-and-intergenic-positions).

#### Aligning XR-seq Data
For this analysis, reads were aligned using [bowtie2.4.5](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) with default parameters after trimming adaptor sequences with [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). Samtools and bedtools were then used to convert to bam and bed format respectively. See the [Optional Dependencies](#optional-dependencies) section for some handy-dandy scripts to help with this process!

#### Calling Lesions From XR-seq Data
Mutperiod has a specialized function for reducing XR-seq reads to single-nucleotide resolution, although this function is hidden from the command line interface. To access this function, run the [ParseXRSeq.py](https://github.com/bmorledge-hampton19/mutperiod/blob/master/python_packages/mutperiod/mutperiodpy/input_parsing/ParseXRSeq.py) script found within the input_parsing module in mutperiod. In the resulting UI, give the XR-seq data in bed format, generated during the above [alignment step](#Aligning-XR-seq-Data). For the lesion call parameter file, give [this file](https://github.com/bmorledge-hampton19/mutperiod/blob/master/python_packages/mutperiod/mutperiodpy/input_parsing/drosophila_CPD_call_params.tsv) which is also found in the mutperiod input_parsing module. Finally, select the dm6.fa genome file where prompted. When run, the program should output files containing the presumed lesion positions at single-nucleotide resolution as well as create the expected mutperiod directory structure for each file.

#### Parsing CPD-seq Data For mutperiod.
Mutperiod also has a function for parsing the CPD-seq bed files with dinucleotide lesion positions. Run the [ParseStandardBed](https://github.com/bmorledge-hampton19/mutperiod/blob/master/python_packages/mutperiod/mutperiodpy/input_parsing/ParseStandardBed.py) script, which is part of the mutperiod input_parsing module. Provide the relevant bed files and genome fasta. The script should format the data for mutperiod from there.

#### Expanding Lesion Context
Later analysis expects that the called XR-seq lesions have trinucleotide context surrounding them. This can be achieved using the [ExpandContext.py](https://github.com/bmorledge-hampton19/mutperiod/blob/master/python_packages/mutperiod/mutperiodpy/ExpandContext.py) script, which is part of mutperiod. Simply run the script, select the "drosophila_data" directory containing all the directories with called XR-seq (but NOT CPD-seq) data formatted for mutperiod, and select "trinuc/quadrunuc" as the expansion context. The data should be converted to trinucleotide context and be suitable for the rest of the pipeline.

#### Split Across Genic and Intergenic Positions
In order to split the data accurately across genic and intergenic regions, overlapping gene designations need to be merged. This can be done using the [MergeGeneRanges.py](https://github.com/bmorledge-hampton19/Chromatin_Features_Analysis/blob/main/python_scripts/MergeGeneRanges.py) script in this repository. When running the script, select the gene ranges file, and make sure to check the box labeled "Preserve Ambiguous Regions". The resulting merged genes file can then be used in tandem with any bed files to split the rows into genic and intergenic output files via the [SplitGenicAndIntergenic.py](https://github.com/bmorledge-hampton19/Chromatin_Features_Analysis/blob/main/python_scripts/SplitGenicAndIntergenic.py) script in this repository. Later steps in this analysis expect these split file types, so this script should be used to generate them by running it and providing it with the directory containing the XR-seq and CPD-seq data files (The relevant files will be acquired automatically.)

As mentioned [previously](#preparing-directories-for-mutperiod), data stratification not directly supported by mutperiod has to be accompanied by manual stratification of the data files into their own individual directories. New directories need to be created for the new genic and intergenic variants of the newly created files. These directories should be named the same as their parent directories but with a "\_genic" or "\_intergenic" suffix. Lastly, these directories have to be prepared for the mutperiod pipeline. This can be achieved using the [ParsePreparedInput.py](https://github.com/bmorledge-hampton19/mutperiod/blob/master/python_packages/mutperiod/mutperiodpy/input_parsing/ParsePreparedInput.py) script from the mutperiod input_parsing module.

#### Binning Across the Genome
The [BinAcrossGenome](https://github.com/bmorledge-hampton19/Chromatin_Features_Analysis/blob/main/python_scripts/BinAcrossGenome.py) script is used to bin damage/repair counts across the genome. This script needs to be run on all damage (naked and cellular) and repair data for 10,000 and 100,00 bp bins. The relevant input files end with "context_mutations.bed" (despite actually being repair or damage data) in the mutperiod_data sub directories for each data set. These files will be auto-acquired if any directory above them is selected in the BinAcrossGenome UI. The output files will automatically be named appropriately but will still need to be manually placed in the relevant "damage/", "naked/", or "repair/", subdirectories of the mutperiod_data directory.

#### Binning Within Genes Across Strands
In order to assess the transcriptional asymmetry of repair and damage, these data need to be binned across the transcribed and non-transcribed strands of genes. This is accomplished by running the [BinInGenes](https://github.com/bmorledge-hampton19/Chromatin_Features_Analysis/blob/main/python_scripts/BinInGenes.py) script across all the relevant repair and damage files. The script will have to be run twice for each file, once with a color designation (check the "Color domain is present..." box and give "\_flanked_colored" as the custom suffix) and once without a color designation (give "\_flanked" as the custom suffix). In either case, the filtered, colored, and expanded gene designations file generated [previously](#filtering-coloring-and-expanding-gene-designations) should be given and the "gene designations include flanking regions" box should be checked, with bin size set to 356 and bin number set to 3. Similar to the [previous step](#binning-across-the-genome), input files should be those with the "context_mutations.bed" suffix and the output files will need to be manually relocated to the relevant folders.

#### Periodicity Analysis With mutperiod
Much of the following analysis depends on periodicity data generated by mutperiod. Thankfully, this pipeline is fairly streamlined. Full documentation on this pipeline can be found at the [mutperiod repository](https://github.com/bmorledge-hampton19/mutperiod), but the process will also be summarized here as well. For starters, the main pipeline needs to be run, either through [this script](https://github.com/bmorledge-hampton19/mutperiod/blob/master/python_packages/mutperiod/mutperiodpy/RunAnalysisSuite.py) or by running `mutperiod mainPipeline` on the command line. In the UI that follows, all the nucleosome maps created [previously](#stratifying-and-subsetting-nucleosome-maps) need to be selected (this must be done one file at a time), and the normalization method should be set to "Custom Background". All the available toggles need to be checked. (There should be 5 once they are all checked.) For the "Bed Mutation Files" field, select all the aggregate cellular damage directories (those containing data with the WT_0hr_UV prefix) and for the "Custom Background Directory" select the naked damage directory, WT_naked_DNA/. Then, repeat the pipeline with the combined repetitions of the repair data as the "Bed Mutation Files", and the cellular damage data as the "Custom Background Directory".

The main pipeline detailed above will produce nucleosome-relative counts, but these still need to be processed into periodicity data. This can be accomplished either through [this script](https://github.com/bmorledge-hampton19/mutperiod/blob/master/python_packages/mutperiod/mutperiodpy/RunNucleosomeMutationAnalysis.py) or by running `mutperiod periodicityAnalysis` on the command line. In the resulting UI, ALL the repair and damage nucleosome counts files resulting from the main mutperiod pipeline need to be selected. If all the relevant directories are stored within a single larger directory, such as the drosophila_data/ or even the mutperiod_data/ directory, only that directory needs to selected. The output file should be named "full_color_domain_analysis_combined_5.rda" and should be placed directly in the mutperiod_data/drosophila_data/ directory.

#### Retrieving Nucleosome Repeat Length Data
As a byproduct of producing periodicity data from repair and damage datasets, mutperiod creates files to calculate and record nucleosome repeat length using nucleosome positions relative to one another. The relevant files are found in the individual nucleosome map directories within mutperiod_data/\_\_external_data/dm6 and can be recognized by the file suffix: "nuc-group_self_raw_nucleosome_mutation_counts.tsv" These files (one for each nucleosome map) need to be manually moved into the mutperiod_data/drosophila_data/nucleosomes/ directory for later analysis. 

And that should be it! The data files are prepared for the final stages of analysis and figure generation.
***

## Analysis and Figure Generation
If the data has been prepared correctly, running the analyses and generating figures is easy: Simply run/knit the [R notebook](R_scripts/GraphBinnedData.Rmd) in the repository. Also, when running the notebook, keep in mind the following details:

#### Defining the "Bioinformatics_Projects" Directory
The notebook requires a path to the "Bioinformatics_Projects" directory. (More information on this directory is available in the [bioinformatics directory](#the-bioinformatics-directory) section.) In brief, this is the parent directory to the cloned chromatin_features_analysis repository. If you are running the R notebooks in interactive mode, you will be able to select this directory through a file browser dialog. If you are knitting the notebook, you will select this directory from a list of choices in the Shiny UI for the notebook's parameters. If the desired directory is not listed, you will need to add it manually to the R notebook (on line 9).

#### Figure Output
When knitting figure generation notebooks, the figures are automatically output to the R/Drosophila_chromatin_domains_analysis_output directory, which will be created if it doesn't already exist.
***

## Data Availability
The data used for the analyses supported by this repository can be found at the following locations (All data relates to Drosophila Melanogaster):
- dm6 genome: [https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz)
- dm6 chromosome sizes: [https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes](https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes)
- CPD-seq data: [NO-LINK-YET](NO-LINK-YET)
- XR-seq data: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138846](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138846)
- MNase map: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1200479](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1200479)
- Nucleosome map: [NO-LINK-YET](NO-LINK-YET)
- Chromatin domain designations: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22069](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22069)
- Gene region designations: [??????](??????)
- RPKM values: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM480160](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM480160)
