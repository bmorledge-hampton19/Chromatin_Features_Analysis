library(data.table)
library(ggplot2)

# Create a bar plot of dipyrimidine nucleotide frequencies from the given one-dimensional data.
frequencyPlot = function(sequences, title = "Sequence Frequency") {

  frequencies = data.table(table(sequences))
  colnames(frequencies) = c("Sequences", "Count")
  ggplot(data = frequencies, aes(x = Sequences, y = Count)) +
    geom_bar(stat = "identity") +
    labs(title = title)

}

# plot the CPD frequencies and deamination trinucleotide context frequencies
dataFilePath = file.path("..","data", "deamination_data", "CPDs", "control_CPD_data_parsed.bed")
data = fread(dataFilePath)
frequencyPlot(data[[4]], title = "CPD Frequencies")

dataFilePath = file.path("..","data", "deamination_data", "DA_0h", "0h_deamination_data_parsed.bed")
data = fread(dataFilePath)
frequencyPlot(data[[4]], title = "Deamination 0h Trinuc Frequencies")

dataFilePath = file.path("..","data", "deamination_data", "DA_24h", "24h_deamination_data_parsed.bed")
data = fread(dataFilePath)
frequencyPlot(data[[4]], title = "Deamination 24h Trinuc Frequencies")

dataFilePath = file.path("..","data", "deamination_data", "DA_48h", "48h_deamination_data_parsed.bed")
data = fread(dataFilePath)
frequencyPlot(data[[4]], title = "Deamination 48h Trinuc Frequencies")
