library(data.table)

# Read in the data and derive a name for it
binnedCountsFilePath = choose.files(multi = FALSE)
binnedCountsTable = fread(binnedCountsFilePath)
binnedCountsDataName = strsplit(basename(binnedCountsFilePath),'.', fixed = T)[[1]][1]
#binnedCountsDataName = strsplit(binnedCountsDataName, "_1000", fixed = T)[[1]][1]

# Create a bin start column from the bin range column.
binnedCountsTable[,Bin_Start := vapply(strsplit(`Bin_Start-End`, '-'),
                                       function(x) as.integer(x[1]), FUN.VALUE = integer(1))]

# Create a column for counts normalized by total counts
binnedCountsTable[,Counts_Per_Million := Feature_Counts / sum(Feature_Counts) * 10^6]

# Add a column for the majority color domain.
binnedChromatinDomainColorsFilePath = choose.files(multi = FALSE)
binnedChromatinDomainColorsTable = fread(binnedChromatinDomainColorsFilePath)
binnedCountsTable[,Majority_Domain_Color := binnedChromatinDomainColorsTable$Domain_Color]

# Add in a complementary data set.
compBinCountsFilePath = choose.files(multi = FALSE)
compBinCountsTable = fread(compBinCountsFilePath)
compBinCountsDataName = strsplit(basename(compBinCountsFilePath),'.', fixed = T)[[1]][1]

compBinCountsTable[,Bin_Start := vapply(strsplit(`Bin_Start-End`, '-'),
                                       function(x) as.integer(x[1]), FUN.VALUE = integer(1))]

compBinCountsTable[,Counts_Per_Million := Feature_Counts / sum(Feature_Counts) * 10^6]

# Create a log ratio column in the main counts file using the complementary data set.
binnedCountsTable[, Log_Ratio := log(Counts_Per_Million/compBinCountsTable$Counts_Per_Million,2)]

# Create a graph for each chromosome (One data set only)
for (chromosome in unique(binnedCountsTable$Chromosome)) {

  # Skip those funky chromosome fragments...
  if (grepl('_',chromosome)) next

  relevantCountsTable = binnedCountsTable[Chromosome == chromosome]

  barplot(relevantCountsTable[, Feature_Counts],
          names.arg = relevantCountsTable[, Bin_Start],
          xlab = "Chromosome Position", ylab = "Counts", main = paste0(binnedCountsDataName,"_",chromosome))


}

# Default comparison is damage and repair.
firstDataSet = "Damage"
compDataSet = "Repair"

# Create a graph for each chromosome (Repair and Damage both)
for (chromosome in unique(binnedCountsTable$Chromosome)) {

  # Skip those funky chromosome fragments...
  if (grepl('_',chromosome)) next

  # Plot the first data set...
  relevantCountsTable = binnedCountsTable[Chromosome == chromosome]

  barplot(relevantCountsTable[, Counts_Per_Million],
          names.arg = relevantCountsTable[, Bin_Start], col = "#ca0020", border = "#ca0020",
          xlab = "Chromosome Position", ylab = "Counts Per Million Total", main = paste(comparisonTitle, chromosome))

  # Then overlay the second data set!
  par(new = T)

  compRelevantCountsTable = compBinCountsTable[Chromosome == chromosome]

  barplot(compRelevantCountsTable[, Counts_Per_Million], yaxt = 'n',
          names.arg = compRelevantCountsTable[, Bin_Start], col = "#0571b0", border = "#0571b0")

  # What if we keep track of where they intersect?
  par(new = T)
  barplot(pmin(relevantCountsTable$Counts_Per_Million, compRelevantCountsTable$Counts_Per_Million), yaxt = 'n',
          names.arg = compRelevantCountsTable[, Bin_Start], col = "#683968", border = "#683968")


  # Finally, add a legend
  legend("topleft", c(firstDataSet, compDataSet), col=c("#ca0020", "#0571b0"),
         lwd=c(3,3), cex = 0.8, bg = "white")

}


# Create a graph for each chromosome that uses the log ratio between the two data sets.
for (chromosome in unique(binnedCountsTable$Chromosome)) {

  # Skip those funky chromosome fragments...
  if (grepl('_',chromosome)) next

  # Plot the first data set...
  relevantCountsTable = binnedCountsTable[Chromosome == chromosome]

  barplot(relevantCountsTable[, Log_Ratio], names.arg = relevantCountsTable[, Bin_Start],
          xlab = "Chromosome Position", ylab = "Counts Per Million Log Ratio (Damage/Repair)",
          col = relevantCountsTable$Majority_Domain_Color, border = relevantCountsTable$Majority_Domain_Color,
          main = paste(comparisonTitle, chromosome))

}
