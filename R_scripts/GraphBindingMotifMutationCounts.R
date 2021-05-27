library(data.table)

countsFilePath = choose.files(multi = FALSE)
countsTable = fread(countsFilePath)
countsDataName = strsplit(basename(countsFilePath),'.', fixed = T)[[1]][1]
countsDataName = strsplit(countsDataName, "aggregate_", fixed = T)[[1]][2]

barplot(countsTable[[2]], names.arg = countsTable[[1]], xlab = "Motif Position",
        ylab = "Normalized Counts", main = paste(countsDataName, "(Motif Strand)"))

barplot(countsTable[[3]], names.arg = countsTable[[1]], xlab = "Motif Position",
        ylab = "Normalized Counts", main = paste(countsDataName, "(Opposite Strand)"))
