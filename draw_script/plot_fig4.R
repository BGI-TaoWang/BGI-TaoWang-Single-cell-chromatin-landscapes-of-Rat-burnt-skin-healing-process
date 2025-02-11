#Fig2
library(dplyr)
library(ArchR)
library(stringr)
library(ggplot2)
set.seed(1234)

setwd("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result")
figure_path = "Figure4"
dir.create(figure_path)
##fig4a
proj = readRDS("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result/motif_enrich/proj_motif_enrichment.rds")
enrichMotifs = readRDS("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result/motif_enrich/enrichMotifs_mouse.rds")
mt <- plotEnrichHeatmap(enrichMotifs, n = 10,cutOff = 15, transpose = TRUE,returnMatrix=TRUE)

mt_col <- str_remove(colnames(mt), "ENSMUSG\\d+_LINE\\d+_")
mt_col <- str_remove(mt_col, "XP_\\d+_LINE\\d+_")
mt_col <- str_remove(mt_col, "BAE\\d+_LINE\\d+_")
colnames(mt) = mt_col

library(pheatmap)
co <- colorRampPalette(c("#E6E7E8", "#3B94FC","#6B44C6","#71138B" ,"#000000"))(100)
a = pheatmap(mt,cluster_cols = F,cluster_rows = F,color =co,main="TF enrichment",angle_col=180)
pdf(paste0(figure_path,'/TF_enrichment.pdf'),width = 16,height = 8)
a 
dev.off()


##fig4b footprint
mt <- plotEnrichHeatmap(enrichMotifs, n = 10,cutOff = 15, transpose = TRUE,returnMatrix=TRUE)
markerMotifs <- str_remove(colnames(mt), " \\(\\d+\\)")
motifPositions <- getPositions(proj,name="motif_mouse")
BSgenome.Rnorvegicus.UCSC.rn7 = readRDS("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/00.data/BSgenome.Rnorvegicus.UCSC.rn7.rds")
seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = "celltype"
)
saveRDS(seFoot,paste0(figure_path,"/seFoot0820.rds"))
source("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/01.script/colorlist.R")

plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj, 
  normMethod = "Subtract",
  plotName = "all_Footprints0819",
  addDOC = FALSE,
  smoothWindow = 5,
  pal=colorlist
)