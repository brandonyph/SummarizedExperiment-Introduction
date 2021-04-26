
#https://biodatascience.github.io/compbio/bioc/SE.html

library(SummarizedExperiment)
colData <- data.frame(sample=factor(1:6),
                      condition=factor(c("A","A","B","B","C","C")),
                      treated=factor(rep(0:1,3)))
colData


#BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
txdb <- EnsDb.Hsapiens.v86
g <- genes(txdb)
g <- keepStandardChromosomes(g, pruning.mode="coarse")
rowRanges <- g[1:10]

exprs <- matrix(rnorm(6 * 10), ncol=6, nrow=10)
se <- SummarizedExperiment(assay=list("exprs"=exprs),
                           colData=colData,
                           rowRanges=rowRanges)
se


assayNames(se)

mcols(se)$score <- rnorm(10)
mcols(se)

#Using the ranges of a SummarizedExperiment

query <- GRanges("1", IRanges(25000,40000))
se.sub <- se[overlapsAny(se, query), ]
se.sub <- se[se %over% query,]
rowRanges(se.sub)
assay(se.sub)
seqinfo(se)

#Downloading SummarizedExperiment data
url <- "http://duffel.rail.bio/recount/SRP046226/rse_gene.Rdata"
download.file(url, "rse_gene.Rdata")
load("rse_gene.Rda")

rse <- rse_gene

colData(rse)[,1:6]

class(rse$characteristics)

rse$characteristics[1:3]

rse$characteristics[[1]]

rse$condition <- sapply(rse$characteristics, `[`, 3)
rse$treatment <- sapply(rse$characteristics, `[`, 4)

rowRanges(rse)

seqinfo(rse)

##########################################################
#Visualizing count matrix data in a SummarizedExperiment
library(DESeq2)
dds <- DESeqDataSet(rse, ~condition + treatment)

vsd <- vst(dds, blind=FALSE)

library(matrixStats)
rv <- rowVars(assay(vsd))

anno.col <- as.data.frame(colData(vsd)[,c("condition","treatment")])
anno.col

library(pheatmap)
pheatmap(assay(vsd)[head(order(rv, decreasing=TRUE),100),],
         annotation_col=anno.col,
         show_rownames=FALSE, show_colnames=FALSE)

sessionInfo()






