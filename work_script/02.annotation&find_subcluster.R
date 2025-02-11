###细分###
library(ArchR)
library(Seurat)

proj = readRDS("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result/proj_anno.rds")
outpath = "/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result"
setwd(outpath)

#############################Find marker###################################
markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "LSI_Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
for (i in 1:length(markerList)) {
  markerList[[i]]$cluster = names(markerList)[i]
}
marker_df = do.call(rbind, markerList)
write.csv(marker_df, "anno/LSI_Clusters_markers.csv")

# #############################画marker基因的heatmap########################
heatmapGS <- markerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  labelMarkers = markerGenes,
  transpose = TRUE
)

plotPDF(heatmapGS, name = "anno/GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE)


#############Find KC subcluster########################
cellsPass <- proj$cellNames[which((proj$LSI_Clusters == "C11"))]
proj_sub  = proj[cellsPass, ]
proj_sub <- addClusters(
  input = proj_sub,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "LSI_Clusters_sub",
  resolution = 0.2,
  force=T
)

proj_sub$LSI_Clusters_sub = paste0("C11_sub@",proj_sub$LSI_Clusters_sub)
proj$LSI_Clusters[match(proj_sub$cellNames,proj$cellNames)] = proj_sub$LSI_Clusters_sub
##########add annotation################################
anno = read.csv("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result/anno/anno0801.csv")
proj = readRDS("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result/proj_recluster0801.rds")
proj$celltype = anno[as.vector(match(proj@cellColData$LSI_Clusters,anno$LSI_Clusters)),"subcelltype"]

#########Find FIB subcluster##########
cells = c("FB_reti","FB_reti_profibrosis")
cellsPass <- proj_filter$cellNames[which((proj_filter$celltype %in% cells ))]
proj_fb  = proj_filter[cellsPass, ]

proj_fb <- addClusters(
  input = proj_fb,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "LSI_Clusters_sub",
  resolution = 0.15,
  force=T
)

fb_anno = c("C4" = "FB_papi","C2" = "FB_Tnc","C1"="FB_pro_inf","C3"="FB_COL")
proj_fb$celltype = fb_anno[proj_fb$LSI_Clusters_sub]
proj_filter$celltype[match(proj_fb$cellNames,proj_filter$cellNames)] = proj_fb$celltype
saveRDS(proj_filter,"proj_recluster0802end.rds")

########filter unknow cell ####################
cells = setdiff(unique(proj$celltype),c("unknow1","unknow2"))
cellsPass <- proj$cellNames[which((proj$celltype %in% cells ))]
proj_filter  = proj[cellsPass, ]

proj_filter <- addUMAP(
    ArchRProj = proj_filter, 
    reducedDims = "IterativeLSI", 
    name = "UMAP_filter_unknow", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force=T
)
saveRDS(proj_filter,"proj_recluster.rds")
