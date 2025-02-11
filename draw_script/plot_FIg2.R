#Fig2
library(dplyr)
library(ArchR)
library(stringr)
library(ggplot2)
set.seed(1234)

outpath = "/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result"
figure_path = "Figure2"
dir.create(figure_path)
proj <- readRDS("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result/proj_filter.rds")

#对样本进行排序
df = as.data.frame(proj@cellColData)
df$time <- str_extract(df$condition, "\\d+[hD]")
df$time_num <- as.numeric(str_remove(df$time, "[hD]"))
df$time_unit <- str_extract(df$time, "[hD]")
df$time_num[df$time_unit == "h"] <- df$time_num[df$time_unit == "h"] / 24
df$condition <- factor(df$condition, levels = unique(df$condition[order(df$time_num)]))


#Fig.2a
source("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/01.script/colorlist.R")
colorlist_sample = colorlist_sample[1:length(unique(proj@cellColData$condition))]
p1 = ggGroup(x = proj@cellColData$condition,y = proj@cellColData$TSSEnrichment,
            groupOrder=levels(df$condition),
            ylabel = "TSSEnrichment",
            xlabel = "library",
            plotAs="violin",
            baseSize = 18,
            size =0.5,
            #ratioYX = 0.2,
            pal = colorlist_sample
            )
p1 = p1 + scale_y_continuous(n.breaks=10)
ggsave(filename= paste0(figure_path,"/TSS_violin.pdf"),p1,width = 10)

p2 = ggGroup(x = proj@cellColData$condition,y = log(proj@cellColData$nFrags,base = 10),
            groupOrder=levels(df$condition),
            ylabel = "log10(nFragments)",
            xlabel = "library",
            plotAs="violin",
            baseSize = 18,
            size =0.5,
            #ratioYX = 0.2,
            pal = colorlist_sample
            )
p2 = p2 #+ scale_y_continuous(n.breaks=10)
ggsave(filename= paste0(figure_path,"/log10(nFragments)_violin.pdf"),p2,width = 10)


##Fig2b:TSS 分布图
p = plotTSSEnrichment(ArchRProj = proj,pal= colorlist_sample)
figure_name = paste0(figure_path, "/Fig2c_TSS_distribution.pdf")
ggsave(figure_name , p , width = 5 , height = 5,limitsize=F) 


#Fig.2c
xlim <- c(2.5,max(log(proj@cellColData$nFrags,base = 10))*1.1  )
ylim <- c(0,max(proj@cellColData$TSSEnrichment)*1.1  )
# 定义一个颜色调色板
color_palette <- colorRampPalette(c("#4760AD", "white", "#B22837"))(100)
p3 = ggPoint( y= proj@cellColData$TSSEnrichment,x = log(proj@cellColData$nFrags,base = 10),
              colorDensity =TRUE,
              xlim = xlim,
              ylim = ylim,
              baseSize = 18,
              legendSize = 15,
              xlabel = "log10(nFragments)",
              ylabel = "TSSEnrichment",
              pal = color_palette
)+
  geom_vline(xintercept = 3, linetype = "dotted", color = "black") + 
  geom_hline(yintercept = 4, linetype = "dotted", color = "black") 

ggsave(filename= paste0(figure_path,"/point_filtering_plots.pdf"),p3)



##fig2d
proj = readRDS("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/Skin_standby/03.ATAC_burn_nc/02.result/proj_recluster0802end.rds")
proj$sample = paste(proj$term,proj$Time,sep="-")
names(colorlist_samples2) = unique(proj$sample)
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "sample", embedding = "UMAP_filter_unknow",pal = colorlist_samples2,baseSize = 0,size=1.2) 
figure_name = paste0(figure_path, "/Fig2_umap.pdf")
ggsave(figure_name , p1 , width = 10 , height = 10,limitsize=F) 

##fig2e
p_bar3 = df %>% ggplot(aes(x= sample,fill=sample))+
  geom_bar(stat = "count",col="black")+
  scale_fill_manual(values = colorlist_samples2)+
  theme_classic()+
  labs(y= "Cell Number")+
  theme(text = element_text(size = 10,face="bold"),axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1), plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 15))+
  theme(panel.border = element_rect(fill = NA,color = "black",size=1.5,linetype = "solid"))+
  coord_flip()
ggsave(paste0(figure_path,"/Fig2_cell_num_barplot.pdf"),p_bar3,width = 8,height = 5,limitsize=F)


