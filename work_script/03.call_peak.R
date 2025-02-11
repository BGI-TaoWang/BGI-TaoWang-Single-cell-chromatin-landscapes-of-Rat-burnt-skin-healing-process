library(ArchR)
library(dplyr)
library(BSgenome)
library(AnnotationDbi)
library(parallel)

print(Sys.time())

obj_path = "/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result/proj_recluster.rds"
cluster_name = "subcelltype"
outpath = "/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result"

print(getwd())

############create genome annotation file##############
library("BSgenome.Rnorvegicus.UCSC.rn7")
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Rnorvegicus.UCSC.rn7)
org <- loadDb("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/01.ATAC_rn7/ArchR1/00.data/annotation/orgDB_Rattus.sqlite")
txdb <- loadDb("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/01.ATAC_rn7/ArchR1/00.data/annotation/TxDb.Rnorvegicus.UCSC.rn7.sqlite")
geneAnnotation <- createGeneAnnotation(TxDb = txdb , OrgDb = org)
##############read ArchR proj ###################
setwd(outpath)
projHeme2 <- readRDS(obj_path)
addArchRThreads(threads = 6) 

projHeme3 <- addGroupCoverages(ArchRProj = projHeme2, groupBy = cluster_name)

pathToMacs2 <- "/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/USER/daixi/miniconda3/envs/scanpy/bin/macs2"
projHeme3 <- addReproduciblePeakSet(
    ArchRProj = projHeme3, 
    groupBy = cluster_name,  
    pathToMacs2 = pathToMacs2,  
    genomeAnnotation = genomeAnnotation,
    geneAnnotation = geneAnnotation,
    genomeSize = 2647899415 
)
projHeme4 <- addPeakMatrix(projHeme3)
saveRDS(projHeme4,"proj_peak_call.rds")


