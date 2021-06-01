# This script compares feature counts across nucleosomes in pre-defined quartile ranges.
library(data.table)
library(ggplot2)

# Get the quartile range files and convert them to data.tables
upperQuartileFilePath = choose.files(multi = FALSE, caption = "Upper Quartile File")
upperQuartileTable = fread(upperQuartileFilePath)

lowerQuartileFilePath = choose.files(multi = FALSE, caption = "Lower Quartile File")
lowerQuartileTable = fread(lowerQuartileFilePath)

# Get the feature counts, read them into a table, and split up that table based on the quartile ranges.
featureCountsFilePath = choose.files(multi = FALSE, caption = "Feature Counts File")
featureCountsTable = fread(featureCountsFilePath)

# If desired, get a background counts to normalize the feature counts.
backgroundCountsFilePath = choose.files(multi = FALSE, caption = "Background Counts File")
if (length(backgroundCountsFilePath) > 0) {
  backgroundCountsTable = fread(backgroundCountsFilePath)
  backgroundCountsTable[Feature_Counts == 0, Feature_Counts := 1]# Pseudo-counts!
  featureCountsTable[,Feature_Counts := Feature_Counts / backgroundCountsTable$Feature_Counts]
}

# Split the data into lower and upper quartile ranges.
featureCountsTable[Nucleosome %in% upperQuartileTable$Nucleosome, Quartile := "Upper Quartile"]
featureCountsTable[Nucleosome %in% lowerQuartileTable$Nucleosome, Quartile := "Lower Quartile"]

title = "PLACEHOLDER TITLE"
xAxislabel = "Quartile_Range"
yAxisLabel = "Feature_Counts"
# Graph the data
ggplot(featureCountsTable[!is.na(Quartile)], aes(Quartile, Feature_Counts)) +
  geom_jitter(width = 0.2, shape = 1, size = 2) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 2, colour = "red") +
  labs(title = title,
       x = xAxislabel, y = yAxisLabel) +
  theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 15), axis.title.x = element_blank())
