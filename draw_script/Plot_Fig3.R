#Fig3
library(dplyr)
library(ArchR)
library(stringr)
library(ggplot2)
set.seed(1234)

setwd("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result")
figure_path = "Figure3"
dir.create(figure_path)
##fig3a\c
proj = readRDS("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result/proj_recluster.rds")
source("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/01.script/colorlist.R")


Muscle = c("pericytes","smooth_muscle_cells")
Endothelial = c("venular_endothelial-2","venular_endothelial-3","lymphatic_endothelia","venular_endothelial-1")
Fibroblast = c("FB_COL","FB_pro_inf","FB_papi","FB_Tnc","FB-myo","Dermal_Papilla")
Immnue = c("Innate lymphoid cell","Neutrophil","monocyte","DC")
Epidermal_Keratinocyte = c("IFE_basel","IFE_basel2","IFE_Spinous")
Hair_follicle_keratinocyte = c("HFSC","IRS","Cortex","ORS","Melanocytes","TAC","HFSC_active","Sebaceous_glands" )
all_celltype = c(Muscle,Fibroblast,"Neuronal_cell",Endothelial,Epidermal_Keratinocyte,Hair_follicle_keratinocyte,Immnue)

proj$bulk_labels = proj$celltype
proj$bulk_labels[proj$celltype %in% Muscle] = "Muscle"
proj$bulk_labels[proj$celltype %in% Endothelial] = "Endothelial"
proj$bulk_labels[proj$celltype %in% Fibroblast] = "Fibroblast"
proj$bulk_labels[proj$celltype %in% Immnue] = "Immnue"
proj$bulk_labels[proj$celltype %in% Epidermal_Keratinocyte] = "Epidermal_Keratinocyte"
proj$bulk_labels[proj$celltype %in% Hair_follicle_keratinocyte] = "Hair_follicle_keratinocyte"

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "celltype", embedding = "UMAP_filter_unknow",pal = colorlist) 
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "bulk_labels", embedding = "UMAP_filter_unknow",pal = colorlist_bulk) 

figure_name = paste0(figure_path, "/Fig3_umap.pdf")
ggsave(figure_name , p1|p2 , width = 15 , height = 8,limitsize=F) 


##fig3b\d
df = as.data.frame(proj@cellColData)
df$sample = paste(df$term,df$Time,sep= "-")
df$time <- str_extract(df$condition, "\\d+[hD]")
df$time_num <- as.numeric(str_remove(df$time, "[hD]"))
df$time_unit <- str_extract(df$time, "[hD]")
df$time_num[df$time_unit == "h"] <- df$time_num[df$time_unit == "h"] / 24
df$sample <- factor(df$sample, levels = rev(unique(df$sample[order(df$time_num)])))
df$celltype <- factor(df$celltype,levels = names(colorlist))
df$bulk_labels <- factor(df$bulk_labels,levels = rev(names(colorlist_bulk)))

#比例图
p_bar1 <- df %>% ggplot(aes(x= sample,fill=celltype))+
  geom_bar(stat = "count",col="black",position = "fill")+
  scale_fill_manual(values = colorlist)+
  theme_classic()+
  labs(y= "Cell Percent")+
  theme(text = element_text(size = 10,face="bold"),axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1), plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 15))+
  theme(panel.border = element_rect(fill = NA,color = "black",size=1.5,linetype = "solid"))+
  coord_flip()

p_bar2 <- df %>% ggplot(aes(x= sample,fill=bulk_labels))+
  geom_bar(stat = "count",col="black",position = "fill")+
  scale_fill_manual(values = colorlist_bulk)+
  theme_classic()+
  labs(y= "Cell Percent")+
  theme(text = element_text(size = 10,face="bold"),axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1), plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 15))+
  theme(panel.border = element_rect(fill = NA,color = "black",size=1.5,linetype = "solid"))+
  coord_flip()

ggsave(paste0(figure_path,"/Fig3bd_cell_percent_barplot.pdf"),(p_bar1|p_bar2),width = 20,height = 5,limitsize=F)

##fig3e track plot
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
anno = read.csv("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result/anno/anno0801.csv")
marker = anno$marker %>% paste(.,sep = ",") %>% strsplit(.,split = ",") %>% unlist() %>%  unique()
marker = marker[marker != ""]

marker_add = c("Enpp2","Crabp1","Cspg4","Colec12","Robo2","Col23a1","Tnc","Postn","Sfrp2","Runx1",
              "Ccl2","Cxcl1","Tgfb2","Cxcl3","Il6","Igf1","Col1a1","Col1a2","Plac8","Pi16","Dpt",
              "Procr","Stat3","Ptgs2","Irf4","Pcsk5","Clec2g")
marker = c(marker,marker_add)

p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "celltype", 
    useGroups= all_celltype,
    pal= colorlist,
    geneSymbol = marker,
    upstream = 20000,
    downstream = 20000
)

plotPDF(plotList = p, 
    name = paste0(figure_path,"_Plot_Tracks_all_Marker_Genes_short_range.pdf"), 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)



#Fig3d
library(Seurat)
library(ComplexHeatmap)
GeneScoreMatrix <- getMatrixFromProject(
                ArchRProj = proj,
                useMatrix = "GeneScoreMatrix"
        )
gene_score <- assays(GeneScoreMatrix)$GeneScoreMatrix
rownames(gene_score) <- rowData(GeneScoreMatrix)$name
mat <- log(gene_score + 1)
obj <- CreateSeuratObject(
                counts = mat,
                assay = 'GeneScore',
                project = 'ATAC',
                min.cells = 1,
                meta.data = as.data.frame(proj@cellColData))
obj <- ScaleData(obj, verbose = FALSE)
saveRDS(obj,paste0(figure_path,"/Seurat_obj0812.rds"))

#caclulate average genescore
obj$celltype = factor(obj$celltype ,levels = all_celltype)
Idents(obj) = 'celltype'
ave_condition = AverageExpression(obj)
anno = read.csv("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result/Figure3/dot_plot_marker.csv")
anno_list <- c()
for (i in 1:nrow(anno)){
        celltype <- anno[i, 'celltype']
        markers <- unlist(strsplit(anno[i, 'marker'], split = ","))
        if (!all(markers %in% rownames(obj@assays$GeneScore))){
                print(paste0("Warnning ", markers[!markers %in% rownames(obj@assays$GeneScore)], " do not exist "))
                markers <- markers[markers %in% rownames(obj@assays$GeneScore)]
        }

        if (length(markers)>0){
                names(markers) <- rep(celltype, length(markers))
                anno_list <- c(anno_list, markers)
        }
}

data = ave_condition$GeneScore[unique(anno_list),]


marker_vec = lapply(rownames(data),function(x){
  i  = data[x,][data[x,] == max(data[x,])]
  names(x) = names(i)
  return(x)
}) %>% unlist()
marker_vec = lapply(colnames(data),function(x){
  tmp = marker_vec[names(marker_vec) == x]
  return(tmp)
}) %>% unlist()
data = data[marker_vec,]
co <- colorRampPalette(c("#4760AD", "white", "#B22837"))(100)
a = pheatmap(data,cluster_cols = F,cluster_rows = F,scale = "row",color =co,main="genescore")
pdf(paste0(figure_path,'/celltype_marker0821.pdf'),width = 6,height = 15)
a 
dev.off()


##Fig3c FeaturePlot
obj = readRDS("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result/Figure3/Seurat_obj0821.rds")
obj <- obj %>% 
  NormalizeData %>%  
  FindVariableFeatures %>%
  ScaleData %>%
  RunPCA %>% 
  RunUMAP(1:30)

UMAP <- getEmbedding(proj, embedding = "UMAP_filter_unknow", returnDF = TRUE)
UMAP <- as.matrix(UMAP)
UMAP = as.data.frame(UMAP)
obj@reductions$umap@cell.embeddings[,1] = UMAP$`IterativeLSI#UMAP_Dimension_1`[match(rownames(obj@reductions$umap@cell.embeddings),rownames(UMAP))]
obj@reductions$umap@cell.embeddings[,2] = UMAP$`IterativeLSI#UMAP_Dimension_2`[match(rownames(obj@reductions$umap@cell.embeddings),rownames(UMAP))]

color_palette <- colorRampPalette(c("#4760AD", "white", "#B22837"))(100)

draw_marker <-  function(obj,marker_list,out_path){
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(stringr)
    library(viridis)
    if (!dir.exists(out_path)) {dir.create(out_path)}
    DefaultAssay(obj) = "GeneScore"
    cells <- unique(names(marker_list))
    marker_umap <- lapply(cells,function(cell){
        cell_marker = marker_list[names(marker_list) == cell]
        cell_marker <- cell_marker[cell_marker %in% rownames(obj)]
        if(length(cell_marker) ==0) {return(NULL)}
        nrow <- ceiling(length(cell_marker)/4) #
        plot <- FeaturePlot(object= obj,features = cell_marker,raster=T,combine = FALSE)
        fix_col = scale_color_viridis()
        p2 <- lapply(plot, function (x) x + fix_col)
        plot = CombinePlots(p2,ncol=4)
        file = file.path(out_path,paste0(cell,"_featureplot.pdf"))
        ggsave(filename = file,plot,width = 20,height = 5*nrow,limitsize=T,dpi = 100)    
    })
}

outpath= "/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result/Figure3/marker_FP"
draw_marker(obj,anno_list,outpath)