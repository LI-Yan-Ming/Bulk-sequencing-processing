rm(list=ls())

setwd("Z:/CT Surgery Lab/Yanming Li/CellLine_RIPseq/antiTFAM_293T/DiffBind/")

####### Annotate TFAM-RNA peaks to genes, using ChIPseeker ###########
setwd("Z:/CT Surgery Lab/Yanming Li/CellLine_RIPseq/antiTFAM_293T/macs2/")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

files <- list.files(path = "./", pattern = "*.narrowPeak$")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

for (i in 1:length(files)){
  peak <- readPeakFile(files[[i]])
  peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Hs.eg.db")
  write.table(as.data.frame(peakAnno), paste0(files[i], "_peaks_Anno.txt"), quote=F, sep="\t")
}

########## check the TFAM-RNA peaks' genomic location ############################
peaks = lapply(files, readPeakFile)
plotPeakProf2(peaks, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body",
              TxDb = txdb, facet = "row", nbin = 400)


################### DiffBind Analysis of 7 samples ##########################################
## There must be at least two sample groups with at least three replicates ##
setwd("Z:/CT Surgery Lab/Yanming Li/CellLine_RIPseq/antiTFAM_293T/DiffBind/")

library(DiffBind)

files <- list.files(path = "./", pattern = "*.narrowPeak$")
SampleID <- stringr::str_split_fixed(files, "[.]", 2)[,1]
SampleID <- substr(SampleID, 1, nchar(SampleID)-6)
bam <- list.files(path = "./", pattern = "*_rmBlackList.bam$")

# Sample sheet
samples <- data.frame(
  SampleID = SampleID,
  Condition = stringr::str_split_fixed(files, "_", 4)[,1],
  Treatment = stringr::str_split_fixed(files, "_", 4)[,2],
  Replicates = stringr::str_split_fixed(files, "_", 4)[,3],
  Peaks = files,
  bamReads = c(bam[8], bam[9], bam[10], bam[11], bam[1], bam[2], bam[3]),
  bamControl = c(bam[6], bam[6], bam[7], bam[7], bam[4], bam[4], bam[5]),
  stringsAsFactors = FALSE
)

# Load data
dba_obj <- dba(sampleSheet = samples, scoreCol=5)
dba_obj <- dba.count(dba_obj)
plot(dba_obj)
dba_obj <- dba.normalize(dba_obj)

# diff analysis
dba_obj <- dba.contrast(dba_obj,design="~Condition + Treatment")
dba_obj <- dba.analyze(dba_obj)
dba.show(dba_obj, bContrasts=TRUE)

saveRDS(dba_obj, "TFAM_RIPseq_7samples_DiffBind_dba.rds")

# View results
results <- dba.report(dba_obj)
print(results)

# Visualization
dba.plotHeatmap(dba_obj)
dba.plotPCA(dba_obj, label=DBA_ID)
dba.plotVolcano(dba_obj) #


#################### DiffBind, CCCP vs CTRL ###########################
dba_obj <- readRDS("TFAM_RIPseq_7samples_DiffBind_dba.rds")

results <- dba.report(dba_obj, contrast = 2, bFlip=TRUE)
print(results)
saveRDS(results, "TFAM_RIPseq_CCCPvsCTRL_grange.rds")
sum(results$Fold>0)
sum(results$Fold<0)

dba.plotHeatmap(dba_obj)
dba.plotPCA(dba_obj, label=DBA_ID) #contrast=2, 

plot(dba_obj, contrast=2)
dba.plotVolcano(dba_obj, contrast=2, bFlip=TRUE) 
dba.plotBox(dba_obj, contrast=2)

## profile plot
rep <- dba.report(dba_obj, contrast=2, bFlip=TRUE)
repList <- GRangesList(Gain=rep[rep$Fold>0,],Loss=rep[rep$Fold<0,])
profiles <- dba.plotProfile(dba_obj, sites=repList, 
                            samples = list(CTRL=dba_obj$masks$Ctrl, CCCP=dba_obj$masks$CCCP),
                            merge=NULL)
dba.plotProfile(profiles)

############### Annotate CCCPvsCTRL peaks to genes, using ChIPseeker and clusterProfile ####################
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(stringr)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
data_background <- toTable(org.Hs.egSYMBOL)

rep <- dba.report(dba_obj, contrast=2, bFlip=TRUE)

peakAnno <- annotatePeak(rep[rep$Fold>0,], tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")
write.table(as.data.frame(peakAnno), "peakAnno_CCCPvsCTRL_up.txt", quote=F, sep="\t")
test1 = bitr(as.data.frame(peakAnno)$SYMBOL, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
each_up_GO <- enrichGO(test1$ENTREZID,OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", 
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = data_background$gene_id, 
                       qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE, pool = FALSE)
each_up_GO@result$Description <- str_remove(each_up_GO@result$Description, ",")
write.csv(each_up_GO, "GO_peakAnno_CCCPvsCTRL_up.csv", quote=F)

peakAnno <- annotatePeak(rep[rep$Fold<0,], tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")
write.table(as.data.frame(peakAnno), "peakAnno_CCCPvsCTRL_down.txt", quote=F, sep="\t")
test1 = bitr(as.data.frame(peakAnno)$SYMBOL, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
each_up_GO <- enrichGO(test1$ENTREZID,OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", 
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = data_background$gene_id, 
                       qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE, pool = FALSE)
each_up_GO@result$Description <- str_remove(each_up_GO@result$Description, ",")
write.csv(each_up_GO, "GO_peakAnno_CCCPvsCTRL_down.csv", quote=F)





