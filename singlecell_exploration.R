library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(Matrix)
library(ggplot2)


load('/Users/ameyateli/LAs BeST/Single Cell/sce.Rdata')

# use R to plot the distribution of a single gene across all the cells 
# pick diff genes, play around with plots etc. 
# assay, functions in the document 
# try a few transformations to see something 

# exploring the data, 
class(sce)

Assays(sce)

dim(sce) # dimensions of count matrix

sce[1:10, 1:3] # examining matrix

sce[7, 1:3]


counts <- counts(sce) # accessing the counts matrix


# looking at a specific cell/cells and its totals 
counts_per_cell <- Matrix::colSums(counts)

cat("counts per cell: ", counts_per_cell[1:5], "\n") ## counts for first 5 cells

# cell number 7
cat("counts per cell: ", counts_per_cell[7], "\n") ## counts for cell 7



# looking at a specific gene and its totals 
counts_per_gene <- Matrix::rowSums(counts)

# gene number 7
cat("counts per gene: ", counts_per_gene[7], "\n") ## counts for cell 7


# barplot of row7: gene7
gene7 <- counts[7, ]

barplot(gene7, main = "Times Gene 7 Was Expressed Per Cell")

hist(gene7, main = "Frequency of Cells expressing gene 7", breaks=11, right = FALSE)
hist(log10(gene7), main = "Frequency of Cells expressing gene 7")

# barplot of row 47: gene 47
gene47 <- counts[47, ]

barplot(gene47, main = "Times Gene 47 Was Expressed Per Cell")

hist(gene47, main = "Frequency of Cells expressing gene 47", breaks=11, right = FALSE)
hist(log10(gene47), main = "Frequency of Cells expressing gene 47")

boxplot(gene47, main = "Frequency of Cells expressing gene 47")


# both gene 7 and gene 47
boxplot(gene7, gene47, main = "Frequency of Cells expressing gene 7 and gene 47")


# testing different genes 
hist(counts[6000,])


# histogram of MT-ATP6
curr_row <- counts[c("MT-ATP6"),]
hist(curr_row, right = FALSE, breaks = max(curr_row)+1, main = "Frequency of Cells expressing MT-ATP6")


# exploring metadata
cell_info <- colData(sce)
head(cell_info)


# 
hist(log(gene7 + 1), main = "Frequency of Cells expressing gene 7", breaks=11, right = FALSE)
hist(log(gene7), main = "Frequency of Cells expressing gene 7", breaks=11, right = FALSE)




# Load the package
library(scater)


plotColData(sce, y = "nGene")
plotColData(sce, y = "stage")
plotColData(sce, y = "type")
plotColData(sce, y = "nUMI")
plotColData(sce, y = "doublet_score")


colnames(colData(sce))


plotColData(sce, x = "stage", y = "nGene", 
            color_by = "stage") + 
  scale_color_manual(values = c(
    "early" = "maroon", 
    "locAdv" = "gold", 
    "metast" = "tomato"
  ))


# normalization 

sce <- logNormCounts(sce)
assayNames(sce)

log_counts <- logcounts(sce) # accessing the counts matrix

# after normalization 
log_gene47 <- log_counts[47,]
hist(log_gene47, main = "Frequency of Cells expressing gene 47", breaks=11, right = FALSE)

# before normalization to compare 
gene47 <- counts[47, ]
hist(gene47, main = "Frequency of Cells expressing gene 47", breaks=11, right = FALSE)





