options(stringsAsFactors = FALSE)

# RNA seq universal script ------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(DESeq2)
library(org.Mm.eg.db)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(ReactomePA)
library(do)
library(ggpubr)
library(enrichplot)
library(export)
library(ggVennDiagram)
library(ggvenn)
library(ggsci)
library(pheatmap)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(ReactomePA)
library(circlize)
color <- c("#000000","#0000ff", "#ff4040","#00c000") #four colors for four groups
setwd("H:/RNASEQ")
# Functions ---------------------------------------------------------------
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
geneBarplot <- function(genepoint){
  if(!("an" %in% ls())){
    load("gtf_gene.Rdata")
    an <-  gtf_gene[,c("gene_name","gene_id","gene_type")]
  }
  if(!("dds" %in% ls())){
    load("dds.RData")
  }
  test <- plotCounts(dds, gene=an[which(an$gene_name==genepoint),"gene_id"], intgroup="group", returnData=TRUE)
  p <- ggplot(test,aes(group, count))+
    geom_bar(stat = "summary",fun=mean,aes(fill=group),show.legend = FALSE)+
    #scale_y_log10()+
    geom_point(show.legend=FALSE,size=6,fill="white",shape=21)+
    scale_fill_manual(values = color)+theme_bw()+
    ggtitle(genepoint)
  return(p)
}
geneHeatmap <- function(genelist,matrix = rld,type = "symbol",group_col = color, heatmap_col,heatmap_name = "Heatmap"){
  matrix <- as.data.frame(matrix)
  matrix$gene <- rownames(matrix)
  if(!("an" %in% ls())){
    load("gtf_gene.Rdata")
    an <-  gtf_gene[,c("gene_name","gene_id","gene_type")]
  }
  matrix <- merge(matrix,an[,c(1,2)],by.x = "gene", by.y="gene_id")
  gene.df <- bitr(matrix$gene_name, fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db)
  matrix <- merge(x=matrix,y=gene.df,by.x="gene_name",by.y="SYMBOL")
  if(type == "symbol"){
    matrix <- matrix %>% filter(gene_name %in% genelist)
  }
  if(type == "ID"){
    matrix <- matrix %>% filter(ENTREZID %in% genelist)
  }
  rownames(matrix) <- matrix$gene_name
  matrix <- matrix %>% 
    dplyr::select(WT_S1,WT_S2,WT_S3,KO_S1,KO_S2,KO_S3,WT_M1,WT_M2,WT_M3,KO_M1,KO_M2,KO_M3) %>% 
    as.matrix(.)
  matrix_scale <- as.data.frame(t(apply(matrix,1,scale)))
  colnames(matrix_scale) <- colnames(matrix)
  matrix_scale <- as.matrix(matrix_scale)
  names(group_col) <- c("WT sham","KO sham","WT MCAO","KO MCAO")
  anno <- HeatmapAnnotation(group = factor(c(rep("WT sham",3),rep("KO sham",3),
                                             rep("WT MCAO",3),rep("KO MCAO",3)),
                                           levels = c("WT sham","KO sham","WT MCAO","KO MCAO")),
                            col = list(group = group_col),
                            border = T) 
  p <- Heatmap(matrix_scale, name = "z-score", col = heatmap_col,
                          #rect_gp = gpar(col = "white", lwd = 2),
                          column_title = heatmap_name,cluster_columns = F,
                          top_annotation = anno)
  return(p)
}
acquire_core_geneid <- function(term_name,enrich_object,top = "ALL",deg_data = deg){
  coregene <- as.data.frame(enrich_object)[term_name,"core_enrichment"] %>% strsplit(split = "/") %>% .[[1]]
  if(top == "ALL"){
    return(coregene)
  } else{
    if(as.data.frame(enrich_object)[term_name,"NES"] > 0){
      gene <- filter(deg_data,ENTREZID %in% coregene) %>% slice_max(log2FoldChange, n = top) %>% .$ENTREZID
      return(gene)
    }
    if(as.data.frame(enrich_object)[term_name,"NES"] < 0){
      gene <- filter(deg_data,ENTREZID %in% coregene) %>% slice_min(log2FoldChange, n = top) %>% .$ENTREZID
      return(gene)
    }
  }
}
genelistify <- function(deg,p.adj = 1){
  geneList <- deg$log2FoldChange
  names(geneList) = deg %>% filter(padj < p.adj) %>% .$ENTREZID
  geneList = sort(geneList, decreasing = TRUE)
  return(geneList)
}
geneBarplot_removebatch <- function(genepoint,matrix=mat){
  result <- t(matrix[which(matrix$gene_name==genepoint),2:13]) %>% as.data.frame(.)
  colnames(result) <- "Counts"
  result$group <- rownames(result) %>% gsub(".{1}$","",.)
  result$group <- factor(result$group,levels = c("WT_S","KO_S","WT_M","KO_M"))
  result <- group_by(result,group)
  plot <- ggplot(result,aes(x=group,y=Counts)) + 
    geom_bar(stat = "summary",fun=mean,aes(fill=group),show.legend = FALSE) + 
    geom_point(show.legend=FALSE,size=6,fill="white",shape=21)+
    scale_fill_manual(values = color)+theme_bw()+
    ggtitle(genepoint)
  return(plot)
}

# Establish DESEQ2 object -------------------------------------------------
#One way ANOVA
{
directory="count" #directory storing .count files
samplesheet <- read.table("H:/RNASEQ/samplesheet.csv",sep=",",header = T)
samplesheet$group <- paste0(samplesheet$condition,".",samplesheet$treatment)
samplesheet[,3:7] <- lapply(samplesheet[,3:7],as.factor)#3:7 refers to all the columns for grouping
samplesheet$treatment <- relevel(samplesheet$treatment,"WT")
samplesheet$condition <- relevel(samplesheet$condition,"sham")
samplesheet$group <- factor(samplesheet$group,levels = c("sham.WT","sham.KO","MCAO.WT","MCAO.KO"))
design= ~ batch +group
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = samplesheet,directory = directory,
                                       design= design)
rm(list=c("directory","samplesheet"))
dds <- DESeq(ddsHTSeq)
}
#save(dds,file = "dds.RData")
load("dds.RData")
vsd <- vst(dds, blind=F)
mat <- assay(vsd)
mm <- model.matrix(~group, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
assay(vsd) <- mat
pcaData <- plotPCA(vsd,intgroup="group",returnData = T)
pcaData$group <- factor(pcaData$group,levels = c("sham.WT","sham.KO","MCAO.WT","MCAO.KO"))
percentVar <- round(100 * attr(pcaData, "percentVar"))
p6 <- ggplot(pcaData, aes(PC1, PC2)) +
  geom_point(size=6,aes(fill=group),alpha=0.4,stroke=0.1,shape=21) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  ylim(c(-20,15))+
  theme_bw()+
  theme(legend.position = "top",legend.box = "vertical",legend.margin=margin())+
  scale_fill_manual(values = color)
rld <- rlog(dds,blind = F)
rld <- assay(rld)
rld <- limma::removeBatchEffect(rld, batch=vsd$batch, design=mm)

#saveRDS(rld,"rld_count_matrix")
#save(rld,file = "rld.RData")

mat_rmbatch <- limma::removeBatchEffect(log2(counts(dds, normalized=TRUE)+1), batch=vsd$batch, design=mm)
mat_rmbatch <- 2^mat_rmbatch
mat_rmbatch <- as.data.frame(mat_rmbatch)
mat_rmbatch <- merge(mat_rmbatch,an[,1:2],by.x = 0,by.y = "gene_id")
geneBarplot_removebatch("Lgals3",mat_rmbatch)

# Correlation -------------------------------------------------------------

matrix <- mat_rmbatch[,c(5:7,11:14)]
matrix <- matrix[which(matrix$gene_name %in% c("Spi1","Mki67","Aif1","Il1b","Arg1","Mrc1","Tnf","Nfkb1","Trem2")),]
rownames(matrix) <- matrix$gene_name
matrix <- log2(as.matrix(matrix[,-7]))
matrix <- t(matrix)
corp <- rcorr(matrix)
corrplot(corp$r,p.mat = corp$P,
         #type = "upper",
         tl.col = "black",
         sig.level = c(0.001,0.01,0.05),
         insig = "label_sig",
         pch.col = "white",
         #method = "color",
         pch.cex = 1.2,is.corr = T)
# Comparisons results -----------------------------------------------------
{
res <- results(dds,contrast = c("group","sham.KO","sham.WT"),alpha=0.05)
res <- lfcShrink(dds=dds,contrast = c("group","sham.KO","sham.WT"),type = "ashr",res=res)
#saveRDS(res,"final.res_MCAO.KO_vs_MCAO.WT")
#saveRDS(res,"final.res_MCAO.WT_vs_sham.WT")
#saveRDS(res,"final.res_sham.KO_vs_sham.WT")

res <- readRDS("final.res_MCAO.KO_vs_MCAO.WT")
res <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
}
deg <- res
# Prepare DEG datasets ----------------------------------------------------
{
#load("gtf_gene.Rdata")
an <-  gtf_gene[,c("gene_name","gene_id","gene_type")] #取出gtf文件中这三列
deg <- merge(res,an,by.x = "gene", by.y="gene_id")
deg <- deg[deg$gene_type=="protein_coding",]
deg <- na.omit(deg)
gene.df <- bitr(deg$gene_name, fromType = "SYMBOL",
                toType = c("ENTREZID", "ENSEMBL"),
                OrgDb = org.Mm.eg.db)
deg <- merge(x=deg,y=gene.df,by.x="gene_name",by.y="SYMBOL")
deg <- dplyr::filter(deg, ENTREZID != "NA")
deg <- deg[which(duplicated(deg$ENTREZID) == F), ]
#saveRDS(deg,"final.deg_unfiltered_MCAO.WT_vs_sham.WT")
}
deg <- readRDS("final.deg_unfiltered_MCAO.KO_vs_MCAO.WT")

{
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

deg$change = ifelse(deg$padj < padj.cutoff & abs(deg$log2FoldChange) >= lfc.cutoff, 
                    ifelse(deg$log2FoldChange> lfc.cutoff ,'Up','Down'),
                    'Stable')
}
label <- deg %>% arrange(-abs(log2FoldChange),padj) %>% head(.,20) %>% .$gene_name
keyvals <- ifelse(
  deg$log2FoldChange < -lfc.cutoff & deg$padj < padj.cutoff, 'royalblue',
  ifelse(deg$log2FoldChange > lfc.cutoff & deg$padj < padj.cutoff, '#FF622F',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#FF622F'] <- 'Up'
names(keyvals)[keyvals == 'black'] <- 'Stable'
names(keyvals)[keyvals == 'royalblue'] <- 'Down'
volcano_plot_sham <- EnhancedVolcano(deg,
                      lab = deg$gene_name,
                      x = 'log2FoldChange',
                      selectLab = label,
                      pCutoff = padj.cutoff,
                      FCcutoff = lfc.cutoff,
                      #drawConnectors = TRUE,
                      #boxedLabels = TRUE,
                      caption = paste0("total = ", nrow(deg), " genes"),
                      colAlpha = 0.4,
                      colCustom = keyvals,
                      ylab=bquote(~-Log[10] ~ italic(p.adj)),
                      legendLabels=c('NS',
                                     paste0('|FC| > ',round(2^lfc.cutoff,1)),
                                     paste0("p.adj < ",padj.cutoff),
                                     paste0('|FC| > ',round(2^lfc.cutoff,1)," & p.adj < ",
                                            padj.cutoff)),
                      title="WT MCAO vs WT sham", #note
                      pointSize = 3,
                      subtitle = paste0("Up: ",nrow(deg[grep(pattern="Up",deg$change),]),"  Down: ",nrow(deg[grep(pattern="Down",deg$change),])),
                      #widthConnectors = 1,
                      xlim = c(-7.5,12.5),
                      y = 'padj')
ggsave(plot=volcano_plot,filename = "volcano_rna_wtvswtsham.svg",device = "svg",bg=NULL,width=5,height=5,units = "in")
geneBarplot("Plk1")
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("#0f769e", "#FFFFFF", "#c72e40"))
col_fun(seq(-3, 3))
geneHeatmap(genelist = c("Spi1","Nfya","Sp2","Ctcfl","Fos","Creb3l4","Pknox1","Rela","Runx1","Tbp","Nub1","Osr2","Hlf","Nrf1","Khdrbs1"),group_col = add.alpha(color,alpha = 0.8),heatmap_col = col_fun,heatmap_name = "M2-like Genes")
wm_genes <- readxl::read_excel("H:/RNASEQ/mmc4.xlsx",sheet = "Venn Analysis",col_names = TRUE)
geneHeatmap(genelist = wm_genes$`Common Genes (39)`[1:39])
# Downstream analysis (enrichment) ----------------------------------------

geneList <- genelistify(deg)

rea <- gsePathway(geneList, 
           pvalueCutoff = 0.05,
           organism = "mouse",
           pAdjustMethod = "BH", 
           verbose = FALSE)
View(as.data.frame(rea))
rea <- as.data.frame(rea) %>% dplyr::select(ID,NES,p.adjust) %>% 
  rename(NES_WTMCAOvsSHAM = NES, p.adjust_WTMCAOvsSHAM = p.adjust)
merge_reactome <- merge(dplyr::select(as.data.frame(rea_kovswt),ID,Description,NES,p.adjust),rea,by = "ID")
merge_reactome <- merge_reactome %>% filter(NES*NES_WTMCAOvsSHAM < 0)
merge_reactome <- merge_reactome %>% 
  rename(`KO MCAO vs WT MCAO` = p.adjust, `WT MCAO vs WT sham` = p.adjust_WTMCAOvsSHAM) %>% 
  pivot_longer(.,cols = c("KO MCAO vs WT MCAO","WT MCAO vs WT sham"),names_to = "group",values_to = "p.adjust") %>% 
  rename(`KO MCAO vs WT MCAO` = NES, `WT MCAO vs WT sham` = NES_WTMCAOvsSHAM) %>% 
  pivot_longer(.,cols = c("KO MCAO vs WT MCAO","WT MCAO vs WT sham"),names_to = "group",values_to = "NES",names_repair = "unique") %>% 
  filter(group...3 == group...5) %>% rename(group = group...3) %>% dplyr::select(-group...5)
merge_reactome$Description <- factor(merge_reactome$Description,levels = arrange(filter(merge_reactome,NES < 0),NES)$Description)

ggplot(data = merge_reactome,aes(x=Description,y=NES))+
  geom_bar(stat = "identity",aes(fill=factor(group),alpha=p.adjust))+
  scale_alpha(range = c(0.7,0.1))+theme_pubr()+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())+
  geom_text(aes(label=Description),stat = "summary")+scale_fill_npg()+
  coord_flip()


rea_kovswt <- rea

coregene <- acquire_core_geneid("R-MMU-68877",rea,top = 20,deg_data = deg) #note deg_data
geneHeatmap(coregene,type = "ID")

# ATAC seq universal script -----------------------------------------------

library(DiffBind)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ggraph)
library(ChIPseeker)
library(enrichplot)
library(ggVennDiagram)
library(tidyverse)
library(ggplotify)
library(patchwork)
library(pheatmap)
library(VennDiagram)

setwd("H:/ATAC_seq")
dbObj <- dba(sampleSheet="sample_sheet.csv")
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
saveRDS(dbObj, file = "dbObj.rds")
#dbObj <- readRDS("dbObj.rds")
dbObj <- dba.normalize(dbObj)
#batch factor stored in [conditon], using block to remove batch effect.
dbObj <- dba.contrast(dbObj,design = ~ Condition + Treatment,minMembers = 2)
dbObj <- dba.contrast(dbObj,design = ~ Condition + Treatment,contrast=c("Treatment","WT_MCAO","WT_sham"),minMembers = 2)
dbObj <- dba.contrast(dbObj,design = ~ Condition + Treatment,contrast=c("Treatment","KO_MCAO","WT_MCAO"),minMembers = 2)
dbObj <- dba.contrast(dbObj,design = ~ Condition + Treatment,contrast=c("Treatment","KO_sham","WT_sham"),minMembers = 2)



#dbObj <- dba.contrast(dbObj,design = ~ Treatment,minMembers = 2) single factor
dbObj <- dba.analyze(dbObj,method=DBA_DESEQ2)

