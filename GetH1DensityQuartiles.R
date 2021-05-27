library(data.table)

# Get the h1 density file and read it into a data.table
h1DensityFilePath = choose.files(multi = FALSE, caption = "H1 Density File")
h1DensityTable = fread(h1DensityFilePath)
h1DensityTable = h1DensityTable[order(-h1DensityTable[[2]])]

# Construct file names for the upper and lower density quartiles.
h1DensityUpperQuartileFilePath = paste0(strsplit(h1DensityFilePath, '.', fixed = T)[[1]][1],
                                        "_upper_quartile.tsv")
h1DensityLowerQuartileFilePath = paste0(strsplit(h1DensityFilePath, '.', fixed = T)[[1]][1],
                                        "_lower_quartile.tsv")

# Get the upper and lower quartiles for the data.
h1DensityUpperQuartileTable = h1DensityTable[1:(nrow(h1DensityTable)/4)]
h1DensityLowerQuartileTable = h1DensityTable[(nrow(h1DensityTable)*3/4 + 1):nrow(h1DensityTable)]

fwrite(h1DensityUpperQuartileTable, h1DensityUpperQuartileFilePath, sep = '\t')
fwrite(h1DensityLowerQuartileTable, h1DensityLowerQuartileFilePath, sep = '\t')
