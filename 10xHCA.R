rm(list = ls())

suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(scran))

file = "ImmuneCells"

## read expression matrix (gene x cell)
matrix_dir = paste0("/Users/izayaorihara/Documents/Genomics/PeakEnreachment/data/hca/",file,"/") # local
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
feature.path <- paste0(matrix_dir, "features.tsv.gz")
cell.path <- paste0(matrix_dir, "cells.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")

mat <- readMM(file = matrix.path)

feature.names <- read.delim(feature.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names <- read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
cells.names <- read.delim(cell.path,
                         header = TRUE,
                         stringsAsFactors = FALSE) %>% select(file_uuid)

colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1


## seperate by samples
ids <- unique(cells.names$file_uuid)
i=2
id <- which(cells.names$file_uuid == ids[i])
# donor 1
mat.id <- mat[,id]
dim(mat.id)
# length(unique(colnames(mat.donor1))) # no duplicated cells




##################
## Filter Cells ##
##################
## exploratory
# Visualize QC metrics as a violin plot
# data <- CreateSeuratObject(counts = mat.id, project = file)
# VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


## remove cells with too few or too many features detected
nFeature <- colSums(mat.id!=0)
kid <- which(nFeature > 200 & nFeature < 6000) # adjusted based on vlnplot
mat.id <- mat.id[,kid]


## filter cells with high perc of mt read cnts (cells that have >= 20% mitochondrial counts were removed)
# unique(feature.names$V5)
mt.index <- which(feature.names$V5 == "chrM")
mt.perct <- colSums(mat.id[mt.index,]) / colSums(mat.id)
# hist(mt.perct)
kid <- which(mt.perct < .2) # adjusted based on hist
mat.id <- mat.id[,kid]


## protein_coding genes were retained
mat.id <- mat.id[which(feature.names$V4 == "protein_coding"),]


## cells with less than 500 protein coding genes detected were removed
nFeature <- colSums(mat.id != 0)
# hist(nFeature)
kid <- which(nFeature > 500)
mat.id <- mat.id[,kid]




##################
## Filter Genes ##
##################
## remove Mitochondrial and KI/GL genes
chr <- c(paste0("chr",1:22), "chrX", "chrY")
kid <- which(feature.names$V5[feature.names$V4 == "protein_coding"] %in% chr)
mat.id <- mat.id[kid,]

## remove not expressed genes
mat.id <- mat.id[rowSums(mat.id) > 0,]

## remove duplicated genes
if(length(unique(rownames(mat.id))) != dim(mat.id)[1]) {
  gn <- rownames(mat.id)  
  rs <- rowSums(mat.id)
  kid <- sapply(unique(gn),function(sid) {
    tmp <- which(gn==sid)
    if (length(tmp)==1) {
      tmp
    } else {
      tmp[which.max(rs[tmp])]
    }
  })
  mat.id <- mat.id[kid,]
  row.names(mat.id) <- gn[kid]
  mat.id <- mat.id[!grepl('^MT-',row.names(mat.id)),]
  mat.id <- round(mat.id)
}

## remove low expr genes
kid <- which(rowMeans(mat.id > 0) >= 0.01)
mat.id <- mat.id[kid,]
dim(mat.id)



##################
## Save Files   ##
##################
mat.id <- round(mat.id)
saveRDS(mat.id, paste0("/Users/izayaorihara/Documents/Genomics/PeakEnreachment/data/hca/",file,"/processed/genebycell_",i,".rds"))


## normalization
sce <- SingleCellExperiment(list(counts=mat.id))
if (ncol(mat.id) < 21){
  sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5), sizes=c(5,10,15,20)) # Window sliding is repeated with different window sizes to construct the linear system, as specified by sizes. By default, the number of cells in each window ranges from 21 to 101.
} else {
  sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5))  
} # calculate size factor
sf <- sizeFactors(sce) # extract size factor
normmat <- sweep(mat.id, 2, sf, '/') # devided by size factor (normalization)


## log2 transformation
normmat <- log2(normmat + 1)
saveRDS(normmat, paste0("/Users/izayaorihara/Documents/Genomics/PeakEnreachment/data/hca/", file, "/processed/norm_genebycell_",i,".rds"))



## save file for BIRD prediction
gene_id <- rownames(normmat)
rownames(normmat) <- NULL
data <- cbind(gene_id, as.data.frame(normmat))
write.table(data, paste0("/Users/izayaorihara/Documents/Genomics/PeakEnreachment/data/hca/",file,"/processed/norm_genebycell_",i,".txt"),
            append = FALSE, sep = "\t", dec = ".", quote = FALSE, row.names = FALSE, col.names = TRUE)

# take random 500 cells
sample <- sample(dim(mat.id)[2], 500)
data <- cbind(gene_id, as.data.frame(normmat)[,sample])
write.table(data, paste0("/Users/izayaorihara/Documents/Genomics/PeakEnreachment/data/hca/",file,"/processed/norm_genebycell_500_",i,".txt"),
            append = FALSE, sep = "\t", dec = ".", quote = FALSE, row.names = FALSE, col.names = TRUE)



## BIRD command
cat("/Users/izayaorihara/Documents/Genomics/BIRD-master/BIRD_predict -b /Users/izayaorihara/Documents/Genomics/BIRD-master/human_hg19_model.bin -i /Users/izayaorihara/Documents/Genomics/PeakEnreachment/data/hca/",file,"/processed/norm_genebycell_500_1.txt -o /Users/izayaorihara/Documents/Genomics/PeakEnreachment/data/hca/",file,"/processed/norm_genebycell_500_1_out.txt",sep = "")




