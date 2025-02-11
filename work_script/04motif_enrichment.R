library(ArchR)
library(ggplot2)
library(patchwork)
library("BSgenome.Rnorvegicus.UCSC.rn7") 

proj = readRDS("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result/proj_peak_call.rds")
path="/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result/"
outpath='/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result/motif_enrich/'

##add  motif annotation
library(chromVARmotifs)
data("mouse_pwms_v1")
proj = addMotifAnnotations(ArchRProj = proj, motifPWMs = mouse_pwms_v1,species = "mouse",annoName="motif_mouse")
dir.create(outpath)
setwd(path)
dir.create(outpath)
setwd(path)


#filter cell lower than 1000
lower_cell = names(table(proj$celltype))[table(proj$celltype) < 100]
cellsPass = proj$cellNames[!(proj@cellColData$celltype %in% lower_cell)]
proj_filter <- proj[cellsPass, ]

#get celltype marker peak
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

#TF motif enrichment
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "motif_mouse",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
saveRDS(enrichMotifs,paste0(outpath,"enrichMotifs_mouse.rds"))
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 10,cutOff = 15, transpose = TRUE)
pdf(paste0(outpath,"Motifs-Enriched-Marker-Heatmap_rat.pdf"), width = 20, height = 8)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()
saveRDS(projHeme4,"proj_motif_enrichment.rds")



