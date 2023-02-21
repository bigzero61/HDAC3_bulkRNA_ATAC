options(stringsAsFactors = FALSE)
rm(list=ls())
library(ggplot2)
library(DESeq2)
library(DiffBind)
library(ggplot2)
#BiocManager::install("org.Hs.eg.db")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Mm.eg.db)
#library(org.Hs.eg.db)
library(ggraph)
library(ChIPseeker)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(ggVennDiagram)
library(tidyverse)
library(EnhancedVolcano)
library(ggplotify)
library(patchwork)
library(pheatmap)
library(VennDiagram)
library(regioneR)
setwd("H:/ATAC_seq")
color <- c("#808080","#8babd3","#ff8080","#8bc164") #four colors for four groups


# Functions ---------------------------------------------------------------

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


# Processing --------------------------------------------------------------


dbObj <- dba(sampleSheet="sample_sheet.csv")
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
#dbObj <- readRDS("dbObj.rds")
dbObj <- dba.normalize(dbObj)
#batch factor stored in [conditon], using block to remove batch effect.
dbObj <- dba.contrast(dbObj,design = ~ Condition + Treatment)
dbObj <- dba.contrast(dbObj,design = ~ Condition + Treatment,contrast=c("Treatment","WT_MCAO","WT_sham"))
dbObj <- dba.contrast(dbObj,design = ~ Condition + Treatment,contrast=c("Treatment","KO_MCAO","WT_MCAO"))
dbObj <- dba.contrast(dbObj,design = ~ Condition + Treatment,contrast=c("Treatment","KO_sham","WT_sham"))
dbObj <- dba.contrast(dbObj,design = ~ Condition + Treatment,contrast=c("Treatment","KO_MCAO","KO_sham"))



#dbObj <- dba.contrast(dbObj,design = ~ Treatment,minMembers = 2) single factor
dbObj <- dba.analyze(dbObj,method=DBA_DESEQ2)
#To use the results from the blocked analysis, 
#you need to specify method=DBA_DESEQ2_BLOCK (or method=DBA_EDGER_BLOCK) 
#when calling dba.report() or any of the plotting functions.!!!★
## ABOVE is an out-of-date method, just including batch in design formula as we did in RNA-seq analysis.

saveRDS(dbObj, file = "dbObj.rds")
#dbObj <- readRDS("dbObj.rds")

dba.show(dbObj,bContrasts = T)
comp <- dba.report(dbObj, contrast = 2,th = 0.01
                   )
data <- as.data.frame(comp)
#comp_pos <- dba.report(dbObj, contrast = 1,fold = 0.5,bGain = T) 


#Visualization
p6 <- dba.plotHeatmap(dbObj)
pcaplot <- dba.plotPCA(dbObj,DBA_TREATMENT,vColors = add.alpha(color,0.8) )
update(pcaplot,lwd=2,cex = 2)
grid.draw()
p7 <- as.ggplot(p7)

# annotation --------------------------------------------------------------
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

peak.comp  <-  annotatePeak(
  peak = comp,
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  annoDb = "org.Mm.eg.db"
)
View(as.data.frame(peak.comp))
data_anno_wtvssham <- as.data.frame(peak.comp) %>% dplyr::filter(grepl("Promoter", annotation)) %>% 
  dplyr::filter(Fold > 0.5) %>% 
  dplyr::select(geneId,SYMBOL,Fold)

datamerge <- base::merge(data_anno_wtvssham,data_anno_kovswt,by = "geneId")
shammcao <- enrichPathway(gene=datamerge$geneId,pvalueCutoff=1, readable=T,organism = "mouse")
p10 <- EnhancedVolcano(data_anno,
                lab = NA,
                x = 'Fold',
                #selectLab = label,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                #drawConnectors = TRUE,
                #boxedLabels = TRUE,
                colAlpha = 4/5,
                title="KO sham vs WT sham",subtitle = NULL,
                #widthConnectors = 1,
                y = 'FDR')
volcano_atac <- (plot_spacer() + p8) / 
  (p9 + p10)
ggsave(plot=volcano_atac,filename = "volcano_atac.pdf",device = "pdf",bg=NULL,width=12,height=12,units = "in")
lfc.cutoff <- 0.5
FDR.cutoff <- 0.01
data_anno$change = ifelse(data_anno$FDR < FDR.cutoff & abs(data_anno$Fold) >= lfc.cutoff, 
                    ifelse(data_anno$Fold> lfc.cutoff ,'Up','Down'),
                    'Stable')
keyvals <- ifelse(
  data_anno$Fold < -lfc.cutoff & data_anno$FDR < FDR.cutoff, 'royalblue',
  ifelse(data_anno$Fold > lfc.cutoff & data_anno$FDR < FDR.cutoff, '#FF622F',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#FF622F'] <- 'Up'
names(keyvals)[keyvals == 'black'] <- 'Stable'
names(keyvals)[keyvals == 'royalblue'] <- 'Down'
volcano_plot_sham <- EnhancedVolcano(data_anno,
                                     lab = NA,
                                     x = 'Fold',
                                     #selectLab = label,
                                     pCutoff = FDR.cutoff,
                                     FCcutoff = lfc.cutoff,
                                     #drawConnectors = TRUE,
                                     #boxedLabels = TRUE,
                                     caption = paste0("total = ", nrow(data_anno), " genes"),
                                     colAlpha = 0.4,
                                     colCustom = keyvals,
                                     ylab=bquote(~-Log[10] ~ italic(p.adj)),
                                     legendLabels=c('NS',
                                                    paste0('|log2FC| > ',lfc.cutoff),
                                                    paste0("p.adj < ",FDR.cutoff),
                                                    paste0('|log2FC| > ',lfc.cutoff," & FDR < ",
                                                           FDR.cutoff)),
                                     title="KO sham vs WT sham", #note
                                     pointSize = 3,
                                     subtitle = paste0("Up: ",nrow(data_anno[grep(pattern="Up",data_anno$change),]),"  Down: ",nrow(data_anno[grep(pattern="Down",data_anno$change),])),
                                     #widthConnectors = 1,
                                     xlim = c(-2,2),
                                     y = 'FDR')






pie <- plotAnnoPie(peak.comp)
ggsave(plot=pie,filename = "pie.svg",device = "svg",bg=NULL,width=5,height=5,units = "in")

p9 <- ChIPseeker::upsetplot(peak.comp,vennpie=TRUE)
gene <- seq2gene(comp, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene,organism = "mouse")
edox <- pairwise_termsim(pathway2)
treeplot(edox)
emapplot(edox)


geneList <- as.data.frame(peak.comp) %>% 
  dplyr::filter(grepl("Promoter", annotation)) %>% 
  .$Fold
names(geneList) = as.data.frame(peak.comp) %>% 
  dplyr::filter(grepl("Promoter", annotation)) %>% 
  .$geneId
geneList = sort(geneList, decreasing = TRUE)
rea_kovs <- gsePathway(geneList, 
                       pvalueCutoff = 0.2,
                       organism = "mouse",
                       pAdjustMethod = "BH", 
                       verbose = FALSE)
edox2 <- pairwise_termsim(rea_kovs)
p5 <- emapplot(edox2)
treeplot(edox2)
gseaplot2(rea_kovs, geneSetID = c(7,64))

# integrate with RNA seq --------------------------------------------------

comp_pos <- peak.comp %>% as.data.frame(.) %>% filter(Fold>0.5 & FDR < 0.01) %>% 
  filter(grepl("Promoter", annotation)) %>% dplyr::select(geneId,Fold)

comp_neg <- peak.comp %>% as.data.frame(.) %>% filter(Fold < -0.5 & FDR < 0.01) %>% 
  filter(grepl("Promoter", annotation)) %>% dplyr::select(geneId,Fold)
deg_pos <- deg %>% filter(change=="Up") %>% dplyr::select(ENTREZID,log2FoldChange)
deg_neg <- deg %>% filter(change=="Down") %>% dplyr::select(ENTREZID,log2FoldChange)
pos <- list(`ATAC seq` = comp_pos$geneId, `RNA seq`= deg_pos$ENTREZID)
neg <- list(`ATAC seq` = comp_neg$geneId, `RNA seq`= deg_neg$ENTREZID)
p1 <- ggvenn::ggvenn(pos,stroke_size = 0.5,set_name_size = 4) +scale_fill_jco()
p2 <- ggvenn::ggvenn(neg,stroke_size = 0.5,set_name_size = 4) +scale_fill_jco()
merge <- peak.comp %>% as.data.frame(.) %>% filter(abs(Fold)>0.5 & FDR < 0.01) %>% 
  #filter(grepl("Promoter", annotation))%>% 
  dplyr::select(geneId,Fold)
degg <- filter(deg,change=="Up" | change=="Down") %>% dplyr::select(ENTREZID,log2FoldChange)
merge <- merge(x=degg,y=merge,by.x="ENTREZID",by.y="geneId")
merge$vs <- merge$log2FoldChange/merge$Fold
merge$change = ifelse(merge$vs >0 , 
                    ifelse(merge$log2FoldChange >0 ,'Up','Down'),
                    'Inconsistent')
merge$change <- factor(merge$change,levels = c("Up","Down","Inconsistent"))

merge$x <- merge$log2FoldChange*merge$Fold

cor <- cor.test(merge$log2FoldChange,merge$Fold,method = "pearson")
p3 <- ggplot(data = merge,aes(x=log2FoldChange,y=Fold))+
  geom_point(color="black",aes(fill=factor(change),alpha=abs(x)),size=4,shape=21)+
  geom_smooth(method = lm,color="#FA9191",fill="#FFDEDE")+scale_fill_npg()+theme_bw()+
  geom_hline(yintercept = 0,linetype="dashed")+
  #geom_rect(aes(xmin = -Inf,xmax = -0.58,ymin = -Inf,ymax = -0.58),fill = "red",alpha = .01)+
  geom_vline(xintercept = 0,linetype="dashed")+
  xlab("RNA-seq log2FC")+ylab("ATAC-seq log2FC")+labs(alpha="RNA*ATAC",fill="change")+
  stat_cor(method = "pearson",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))
  #facet_wrap(~change)+
  #annotate("text",label=paste0("r2 = ",cor[["estimate"]][["cor"]][1],"\np-value = ",cor[["p.value"]][1]),x=-4,y=1)

olap_rea <- enrichPathway(gene=dplyr::filter(merge,change=="Down")$ENTREZID,pvalueCutoff=1, readable=T,organism = "mouse")
olap_rea <- enrichGO(gene=dplyr::filter(merge,change=="Down")$ENTREZID,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = 'ENTREZID',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2)

p4 <- dotplot(olap_rea)

plot_spacer + p6 + plot_layout(widths=c(2,1))

(p8 + p9 + plot_layout(widths = c(1,2))) /
  (p5 + (p1 / p2) + plot_layout(widths = c(2,1))) /
  p3 + p4 + plot_annotation(tag_levels = "A")
  
comp_neg <- peak.comp %>% as.data.frame(.) %>% filter(Fold< -0.5 & FDR < 0.01)
merge <- merge(x=comp_neg,y= dplyr::select(deg_neg,ENTREZID),by.x="geneId",by.y="ENTREZID")
merge <- dplyr::select(merge,1:12)
merge <- dplyr::select(merge,-geneId)
merge <- merge[order(merge$seqnames, merge$start),]
head(merge)
rea_overlap_neg <- enrichPathway(unique(merge$geneId),organism = "mouse",pvalueCutoff = 1,pAdjustMethod = "fdr")
View(as.data.frame(rea_overlap_neg))
write.table(merge, file= "H:/ATAC_seq/kovswtmcao_overlap_neg_p0.01.bed", sep="\t", quote=F, row.names = F,col.names = F)
kk <- enrichKEGG(gene    = unique(merge$geneId),
                 keyType = 'kegg',
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)
View(as.data.frame(kk))
ego <- enrichGO(gene         = unique(merge$ENSEMBL),
                OrgDb         = org.Mm.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2)
cnetplot(ego)
dotplot(ego, showCategory=20)


# Custom GSEA analysis for overlap genes ----------------------------------
rna <- deg %>% filter(padj < 0.05)
atac <- as.data.frame(peak.comp) %>% filter(FDR < 0.05 & abs(Fold) > 1) %>% dplyr::select(geneId)
overlap <- merge(x=rna,y=atac,by.x="ENTREZID",by.y="geneId") %>% distinct(ENTREZID,.keep_all = T)
geneList <- overlap$log2FoldChange
names(geneList) = overlap$ENTREZID
geneList = sort(geneList, decreasing = TRUE)
library(msigdbr)
msigdbr_inflam_dataset <- msigdbr(species = "Mus musculus",category = "C2",subcategory ="CGP" ) %>% 
  dplyr::select(gs_name, entrez_gene, gene_symbol)
gsea_inflam <- GSEA(geneList = geneList,TERM2GENE = msigdbr_inflam_dataset[,c(1,2)],pvalueCutoff = 1)
gseaplot2(gsea_inflam,c("COATES_MACROPHAGE_M1_VS_M2_UP",
                        "COATES_MACROPHAGE_M1_VS_M2_DN"))


rna_mcao <- deg %>% filter(padj < 0.05)
atac_mcao <- as.data.frame(peak.comp) %>% filter(FDR < 0.01) %>% dplyr::select(geneId)
overlap <- merge(x=rna_mcao,y=atac_mcao,by.x="ENTREZID",by.y="geneId") %>% distinct(ENTREZID,.keep_all = T)
geneList <- rna_mcao$log2FoldChange
names(geneList) = rna_mcao$ENTREZID
geneList = sort(geneList, decreasing = TRUE)
gsea_inflam_mcao <- GSEA(geneList = geneList,TERM2GENE = msigdbr_inflam_dataset[,c(1,2)],pvalueCutoff = 1)
gseaplot2(gsea_inflam_mcao,c("GSE25123_CTRL_VS_IL4_STIM_MACROPHAGE_DN",
                        "GSE9988_LPS_VS_CTRL_TREATED_MONOCYTE_UP"))
test <- merge(dplyr::select(filter(as.data.frame(gsea_inflam),p.adjust < 0.05),Description),
              dplyr::select(filter(as.data.frame(gsea_inflam_mcao),p.adjust < 0.05),Description),by = "Description")


custom_lps <- read_tsv("H:/CustomGeneSet/LPS.tsv")
custom_lps <- custom_lps %>% filter(adj.P.Val < 0.05) %>% dplyr::select(ID,Gene.ID,Gene.symbol,logFC) %>% 
  dplyr::arrange(-logFC) %>% na.omit(.) %>% head(.,n=200) %>% rename(.,entrez_gene = Gene.ID, gs_name = ID) %>% 
  mutate(.,gs_name = "GSE49329 LPS_vs_untreated_UP") %>% dplyr::arrange(Gene.symbol)

custom_il4 <- read_tsv("H:/CustomGeneSet/IL4.tsv")
custom_il4 <- custom_il4 %>% filter(adj.P.Val < 0.05) %>% dplyr::select(ID,Gene.ID,Gene.symbol,logFC) %>% 
  dplyr::arrange(-logFC) %>% na.omit(.) %>% head(.,n=200) %>% rename(.,entrez_gene = Gene.ID, gs_name = ID) %>% 
  mutate(.,gs_name = "GSE49329 IL4_vs_untreated_UP") %>% dplyr::arrange(Gene.symbol)
custom <- rbind(custom_lps,custom_il4)
gsea_mcao_custom <- GSEA(geneList = geneList,TERM2GENE = custom[,c(1,2)],pvalueCutoff = 1)
gseaplot2(gsea_mcao_custom,c("GSE49329 LPS_vs_untreated_UP",
                             "GSE49329 IL4_vs_untreated_UP"))
gsea_kovswt_custom <- GSEA(geneList = geneList,TERM2GENE = custom[,c(1,2)],pvalueCutoff = 1)
gseaplot2(gsea_kovswt_custom,c("GSE49329 LPS_vs_untreated_UP",
                             "GSE49329 IL4_vs_untreated_UP"))
rld <- rlog(dds)
countmatrix <- assay(rld)
countmatrix <- as.data.frame(countmatrix)
countmatrix$gene_id <- rownames(countmatrix)
countmatrix <- merge(countmatrix,filter(an, gene_type == "protein_coding")[,c(1,2)],by="gene_id")
custom_20 <- rbind(head(arrange(custom_lps,-logFC),n=20),head(arrange(custom_il4,-logFC),n=20))
custom_core <- as.data.frame(gsea_kovswt_custom)["GSE49329 LPS_vs_untreated_UP","core_enrichment"] %>% 
  strsplit(.,split = "/") %>% .[[1]] %>% data.frame(entrez_gene = .) %>% 
  merge(.,custom[,c(2,3)],by = "entrez_gene")
countmatrix_custom_m1m2 <- merge(countmatrix,as.data.frame(unique(custom_core$Gene.symbol)),
                                 by.x="gene_name",by.y="unique(custom_core$Gene.symbol)")

rownames(countmatrix_custom_m1m2) <- countmatrix_custom_m1m2$gene_name
loc <- match(custom_20$Gene.symbol,rownames(countmatrix_custom_m1m2))
countmatrix_custom_m1m2 <- countmatrix_custom_m1m2[loc,]
countmatrix_custom_m1m2 <- countmatrix_custom_m1m2 %>% 
  dplyr::select(WT_S1,WT_S2,WT_S3,KO_S1,KO_S2,KO_S3,WT_M1,WT_M2,WT_M3,KO_M1,KO_M2,KO_M3) %>% 
  as.matrix(.)
countmatrix_custom_m1m2 <- countmatrix_custom_m1m2 - rowMeans(countmatrix_custom_m1m2)
pheatmap(countmatrix_custom_m1m2,cluster_rows = F,cluster_cols = F)





# regioneR ----------------------------------------------------------------
atac <- peak.comp %>% as.data.frame(.) %>% 
  filter(Fold < 0) %>% 
  filter(grepl("Promoter", annotation)) %>% 
  filter(annotation == "Promoter (1-2kb)" | annotation == "Promoter (<=1kb)") %>% 
  .$geneId

pathway2 <- enrichPathway(atac,organism = "mouse")
edox <- setReadable(go, 'org.Mm.eg.db', 'ENTREZID')
go <- enrichGO(gene=atac,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = 'ENTREZID',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2)
cnetplot(edox)



cnetplot(edox)

rna <- deg %>% filter(change == "Down") %>% .$ENTREZID


peak.comp %>% as.data.frame(.) %>% filter(Fold < 0) %>% 
  filter(grepl("Promoter", annotation)) %>% 
  #dplyr::select(seqnames,start,end) %>% 
  arrange(seqnames,start) %>% 
  write.table(.,file = "H:/ATAC_seq/mcao.kovswt_promoterdown.bed",sep="\t", quote=F, row.names = F,col.names = F)

peak.comp %>% as.data.frame(.) %>% filter(Fold < -0.5) %>% 
  filter(grepl("Promoter", annotation)) %>% 
  dplyr::select(seqnames,start,end) %>% 
  arrange(seqnames,start) %>% 
  write.table(.,file = "H:/ATAC_seq/mcao.kovswt_promoterdown1101.bed",sep="\t", quote=F, row.names = F,col.names = F)

A <- toGRanges(A = "H:/ATAC_seq/mcao.kovswt_promoterdown1101.bed",genome = "mm10")
B <- toGRanges(A = "H:/ATAC_seq/mcaokovswt0422rnadown.bed",genome = "mm10")
#pt <- overlapPermTest(A=A, B=B, ntimes=200,genome = "mm10", verbose = TRUE)
pt2 <- permTest(A=A, B=B, ntimes=100, alternative="auto", 
                verbose=TRUE, genome="mm10", 
                evaluate.function=meanDistance, randomize.function=randomizeRegions, non.overlapping=F)
plot(pt2)
lz <- localZScore(pt=pt, A=A, B=B)
plot(lz)




# Homer Motif ----------------------------------------------------------------

library("RCurl")
library("XML")
library("magrittr")
library("rvest")
homer <- readHTMLTable("H:/ATAC_seq/result/kovswt_promoter_down_nobg_motif/homerResults.html")
homer <- as.data.frame(homer[["NULL"]])
homer <- homer[1:18,]
homer$`Best Match/Details` <- c("Sfpi1","NFY","Sp2","BORIS","Fos","GFY","CREB3L4","PKNOX1","NFkB-p65","RUNX","GFX","TBP","JGL","OSR2","HLF","CDX2","NRF1","KHDRBS1")
homer$`-log P-pvalue` <- -as.numeric(homer$`log P-pvalue`)
homer <- homer %>% arrange(-`-log P-pvalue`)
homer$`Best Match/Details` <- factor(homer$`Best Match/Details`,levels = rev(as.character(homer$`Best Match/Details`)))
ggplot(homer,aes(x = `Best Match/Details`,y = `-log P-pvalue`))+
  geom_bar(stat = "identity",fill = "blue")+
  coord_flip()+
  ggpubr::theme_pubr()
