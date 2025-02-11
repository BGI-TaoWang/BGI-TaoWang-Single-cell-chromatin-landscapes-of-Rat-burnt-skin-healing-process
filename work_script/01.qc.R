library(dplyr)
library(ArchR)
library(Seurat)
library(BSgenome)
library(AnnotationDbi)
library(stringr)
set.seed(1234)


sample_path = "/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/00.data/all_sample.txt"
out_path = "/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result"
threads = 6
minTSS = 2


sample_list = read.table(sample_path,header = F)
sample_id = unlist(sample_list[,1]) %>% gsub("hot","burn",.)  %>% gsub("SS","NC",.)
fragment_name = unlist(sample_list[,2])
fragments = paste0("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/00.data/fragments/",fragment_name,".fragments.tsv.gz")

if (!dir.exists(out_path)) {dir.create(out_path,recursive = T)}
setwd(out_path)
addArchRThreads(threads = threads) 

#create genome file
library("BSgenome.Rnorvegicus.UCSC.rn7")
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Rnorvegicus.UCSC.rn7)

org <- loadDb("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/01.ATAC_rn7/00.data/annotation/orgDB_Rattus.sqlite")
txdb <- loadDb("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/01.ATAC_rn7/00.data/annotation/TxDb.Rnorvegicus.UCSC.rn7.sqlite")
geneAnnotation <- createGeneAnnotation(TxDb = txdb , OrgDb = org)

# read ArchR fragment
ArrowFiles <- createArrowFiles(
  inputFiles = fragments,
  sampleNames = fragment_name,
  minTSS = minTSS, #Dont set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1,
   force = TRUE
)
# create ArchR project
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  copyArrows = F,  #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)

# filter cell 
idxPass <- which((proj$TSSEnrichment >= 4) & (proj$nFrags >= 1000))
cellsPass <- proj$cellNames[idxPass]
proj  = proj[cellsPass, ]

# add meta
names(sample_id) = fragment_name
proj@cellColData$AnalysisName = proj@cellColData$Sample
proj@cellColData$condition = sample_id[as.vector(proj@cellColData$Sample)]
proj@cellColData$Time = sapply(str_split(proj@cellColData$condition,"-"),function(x) x[1])
proj@cellColData$term = sapply(str_split(proj@cellColData$condition,"-"),function(x) x[2])
proj@cellColData$Position = sapply(str_split(proj@cellColData$condition,"-"),function(x) x[3])

#filter Doublets
proj <- filterDoublets(proj) #filterRatio = 1
p1 <- plotGroups(
  ArchRProj = proj,
  groupBy = "AnalysisName",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "ridges"
)

p2 <- plotGroups(
  ArchRProj = proj,
  groupBy = "AnalysisName",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p3 <- plotGroups(
  ArchRProj = proj,
  groupBy = "AnalysisName",
  colorBy = "cellColData",
  name = "nFrags",
  plotAs = "ridges"
)

p4 <- plotGroups(
  ArchRProj = proj,
  groupBy = "AnalysisName",
  colorBy = "cellColData",
  name = "nFrags",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 20)

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list( 
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30
)

proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)


proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "LSI_Clusters",
  resolution = 0.8
)

proj <- addClusters(
  input = proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Harmony_Clusters",
  resolution = 0.8
)

proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)


##################plot LSI-cluster######################
p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "AnalysisName", embedding = "UMAP")
p6 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "LSI_Clusters", embedding = "UMAP")
ggAlignPlots(p5, p6, type = "h")
plotPDF(p5,p6, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 10)

##################plot Harmony cluster######################
p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "AnalysisName", embedding = "UMAPHarmony")
p6 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Harmony_Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p5, p6, type = "h")
plotPDF(p5,p6, name = "Plot-UMAPHarmony-Sample-Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 10)

##################plot sample distribution######################
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Time", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "condition", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Time-Condition.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 10)
                        
saveArchRProject(ArchRProj = proj, outputDirectory = out_path, load = FALSE)
saveRDS(proj, file = 'proj_filter.rds')
