options(stringsAsFactors = F)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
# 自定义函数 -------------------------------------------------------------------
CoxUniVar <- function(dat, status, times, var, digit = 3){
    dat[[status]] <- as.factor(dat[[status]])
    subgroup <- levels(as.factor(dat[[status]]))
    subgroup1 <- paste0(subgroup[2], " vs ", subgroup[1])
    dat[[status]] <- as.numeric(dat[[status]])
    formu <- as.formula(paste0("Surv(",times,",",status,") ~", var))
    fit <- coxph(formu,data= dat)
    unisum <- summary(fit)
    HR <- exp(coef(fit))[1]
    HR <- round(HR, digit)
    ci <- exp(confint(fit))[1,]
    ci <- round(ci, digit)
    cito <- paste0(ci[1], " - ", ci[2])
    p <- unisum$coefficients[1, "Pr(>|z|)"]
    p.val <- p
    p <- ifelse(p < 0.001, "< 0.001", round(p, 3))
    var1 <- names(exp(coef(fit)))[1]
    result <- c(var1, status,subgroup1, HR, cito, p.val, p)
    names(result) <- c("var", "group", "subgroup","HR", "95%CI", "p","p.val")
    return(result)
}
head5 <- function(dat){dat[1:5,1:5]}

# 加载必须数据 ------------------------------------------------------------------

idmap <- read.delim("~/project/gencode.v22.annotation.gene.probeMap", header = T)
rownames(idmap) <- idmap$id
head(idmap)
geneFamily1 <- c("HIF1A", "EPAS1", "HIF3A", "ARNT", "ARNT2", "EGLN2",
                 "EGLN1", "EGLN3", "VHL", "CREBBP", "EP300", "HIF1AN"
)
geneFamily2 <- c("HNRNPL","HNRNPA2B1","HNRNPH2","HNRNPUL1","HNRNPH1",
                 "HNRNPCL1","HNRNPLL","HNRNPR","HNRNPD","HNRNPU",
                 "HNRNPA1L2","HNRNPA1P33","HNRNPUL2","HNRNPF","HNRNPH3",
                 "HNRNPC","HNRNPA0","HNRNPM","HNRNPUL2-BSCL2","HNRNPAB",
                 "HNRNPA1","HNRNPDL")
geneinfo <- idmap[idmap$gene %in% geneFamily2,]
write.csv(geneinfo, "HNRGene.csv", row.names = F, quote = F)
geneFamily3 <- paste0("SRSF", 1:12)
geneFamily1 <- intersect(geneFamily1,idmap$gene)
geneFamily2 <- intersect(geneFamily2,idmap$gene)
geneFamily3 <- intersect(geneFamily3, idmap$gene)
targetGene <- list(oxy = geneFamily1, HNRNPL = geneFamily2,
                   SRSF = geneFamily3)

# figure1A 寻找pan癌当中的差异分析结果 ------------------------------------------------

### 前期的数据整理
difDatIdenx <- Sys.glob("~/project/TCGADat/DeseqDiff/*")
diffExprCancerIndex <- c("HNSC", "LUSC", "LIHC", "COAD", "READ", "LUAD", 
                         "STAD", "CHOL", "UCEC", "BLCA", "BRCA", "ESCA", 
                         "THCA", "PRAD", "KICH", "KIRC", "KIRP")
difDatIdenx <- str_subset(difDatIdenx, paste(diffExprCancerIndex, collapse = "|"))

## 获得差异表达的logFC
geneFamilyDiffDat <- list()
num <- 0
for(i in targetGene){
    geneFamilyDat1 <- data.frame()
    num <- num + 1
    for(j in difDatIdenx){
        dat <- vroom::vroom(gzfile(j))
        geneFamilyDat <- dat %>% filter(gsm %in% i) %>% 
            select(gsm, log2FoldChange) %>% mutate(cancertype = str_extract(j, "(?<=Diff\\/).+?(?=_)"))
        geneFamilyDat1 <- rbind(geneFamilyDat1, geneFamilyDat)
    }
    geneFamily1Dat2 <- spread(geneFamilyDat1, key = "cancertype", value = "log2FoldChange")
    geneFamilyDiffDat[[num]] <- geneFamily1Dat2
}
geneFamilyDiffDat1 <- map(geneFamilyDiffDat, function(x){
    x <- as.data.frame(x)
    rownames(x) <- x[,1]
    x[-1]
})
### 基因的logfc值进行转换
logfcCat <- lapply(geneFamilyDiffDat1, function(dat){
    dat1 <- apply(dat, 2, function(x){
        cut(x, breaks = c(-Inf, -1, -0.5, 0.5, 1, Inf),
            labels = c("< -1", "-1 - -0.5", "-0.5 - 0.5", "0.5 - 1", "> 1"))
    })
    rownames(dat1) <- rownames(dat)
    return(dat1)
})
# ### 把基因的P值提取为数据数值
geneFamilyDiffP <- list()
num <- 0
for(i in targetGene){
    geneFamilyDat1 <- data.frame()
    num <- num + 1
    for(j in difDatIdenx){
        dat <- vroom::vroom(j)
        geneFamilyDat <- dat %>% filter(gsm %in% i) %>%
            select(gsm, padj) %>% mutate(cancertype = str_extract(j, "(?<=Diff\\/).+?(?=_)"))
        geneFamilyDat1 <- rbind(geneFamilyDat1, geneFamilyDat)
    }
    geneFamily1Dat2 <- spread(geneFamilyDat1, key = "cancertype", value = "padj")
    geneFamilyDiffP[[num]] <- geneFamily1Dat2
}
geneFamilyDiffP1 <- map(geneFamilyDiffP, function(x){
    x <- as.data.frame(x)
    rownames(x) <- x[,1]
    x[-1]
})
# ### 差异表达结果当中P值>0.05的进行标记
logfcCat1 <- map2(logfcCat, geneFamilyDiffP1, function(logfc, p){
    logfc[p >= 0.05] <- "P >= 0.05"
    return(logfc)
})
unique(matrix(logfcCat, ncol  = 1))
#
# # ### 绘图
unique(matrix(logfcCat1[[1]], ncol  = 1))
# ## 每个分类定义一个颜色
col_cat <- c("> 1" = "#A80C3A", "0.5 - 1" = "#ED5E57", "-0.5 - 0.5" = "#DDD3D2",
             "-1 - -0.5" = "#6B9AB7", "< -1" = "#2F5B89", "P >= 0.05" = "white")
cell_fun <- function(logfc, dataP, logfcCutoff = 1, PCutoff = 0.05,
                     darkcol = "black", lightcol = "white", digit = 2, fontsize  = 6){
    function(j, i, x, y, width, height, fill){
        if(abs(logfc[i,j]) > logfcCutoff & dataP[i,j] < PCutoff){
            grid.text(round(logfc, digit)[i, j], x, y,
                      gp = gpar(fontsize = fontsize, col  = lightcol))
        }else{
            grid.text(round(logfc, digit)[i, j], x, y,
                      gp = gpar(fontsize = fontsize, col  = darkcol))
        }
    }
}
## 定义注释信息的颜色
an_col <- c("#ADDD8E", "#1D91C0", "#DD3498")
names(an_col) <- unique(genetype$Group)
row_an <-  HeatmapAnnotation(type = genetype$Group, ##注释信息的内容。
                             show_annotation_name = F, ## 是否显示注释的标题
                             col = list(type = an_col), ## 注释信息的颜色
                             show_legend = T,  ## 是否显示注释信息的说明
                             annotation_legend_param = list(title = "Group"
                                                            #nrow = 1
                                                            ), ## 注释信息图例的标题
                             which = "row" #对行或者列进行注释 
)
p1 <- Heatmap(matrix = logfcCat1[[1]],
        name = "logFC", #主要图例的标题
        rect_gp = gpar(col = "grey", lwd = 1), ##边框颜色的变化
        col = col_cat,
        row_names_side = "right",
        cell_fun = cell_fun(geneFamilyDiffDat1[[1]], geneFamilyDiffP1[[1]]),
        row_names_gp = gpar(fontsize = 10),left_annotation = row_an,
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 45, 
        #heatmap_legend_param = list(nrow = 1),
        #left_annotation = row_an
        )
pdf("oxyFigure1A.pdf", width = 8, height = 6)
p1
dev.off()
write.csv(geneFamilyDiffDat1[[1]], "./oxy/figure 1/oxyDiffLogFc.csv", quote = F)
write.csv(geneFamilyDiffP1[[1]], "./oxy/figure 1/oxyDiffP.csv", quote = F)
write.csv(logfcCat1[[1]], "./oxy/figure 1/oxylogfcCat1.csv", quote = F)
# figure 1B某一个基因的差异表达结果 ---------------------------------------------------
library(ggsci)
tpmindex <- Sys.glob("~/project/TCGADat/tpmDat/*gz")
HIF3AIndex <- diffExprCancerIndex[-8]
TpmDatIdenx <- str_subset(tpmindex, paste(HIF3AIndex, collapse = "|"))
HIF3Aplot <- lapply(TpmDatIdenx[-6], function(x){
    example <- vroom::vroom(gzfile(x))
    colnames(example)[1] <- "id"
    cancertype <- str_extract(x, "(?<=tpmDat/).+?(?=_)")
    HIF3A <- example[example$id == "ENSG00000100644.15", ]
    HIF3A$cancertype <- cancertype
    HIF3ALong <- HIF3A %>% select(-id) %>% gather(key = "sampleid", value = tpm) %>% 
        mutate(sampletype = ifelse(str_detect(sampleid, "0..$"), "cancer", "normal"))
    HIF3ALong$tpm <- as.numeric(HIF3ALong$tpm)
    p1 <- ggplot(HIF3ALong, aes(x = sampletype, y = tpm, group = sampletype)) +
        geom_point(alpha = 0.4, size = 3, stroke = 0, na.rm = TRUE,
                   position = position_jitterdodge(dodge.width = 0.4,jitter.height = 0),
                   aes(color = sampletype)) +
        geom_violin(width = 0.5,
                    alpha = 0.2, fill = "white", na.rm = TRUE) +
        stat_boxplot(notch = F,
                     geom = "boxplot", width = 0.3, alpha = 0.2, fill = "white",
                     na.rm = TRUE)  +
        theme_classic() + scale_color_npg() +
        guides(color = F) + labs(y = "HIF1A expression", x = cancertype) + scale_y_log10()
    return(p1)
})
cowplot::plot_grid(plotlist = HIF3Aplot, ncol = 4)
ggsave("./oxy/figure 1/HIF1A.pdf", width = 8, height = 6)



# figure 2A突变统计 -----------------------------------------------------------

### 前期数据处理
mutateIndex <- Sys.glob("~/project/TCGADat/MutateDat/*")
MutateDat <- list()
num <- 0
for(i in targetGene){
    num <- num + 1
    MutaFreq <- c()
    for(j in mutateIndex){
        example <- read.delim(gzfile(j))
        nsample <- length(unique(example$Sample_ID))
        MutaFreq1 <- sapply(i, function(x){
            nrow(example[example$gene == x, ])/nsample
        })
        MutaFreq <- rbind(MutaFreq, MutaFreq1)
    }
    MutateDat[[num]] <- MutaFreq
}
MutateDat <- lapply(MutateDat, function(x){
    rownames(x) <- str_extract(mutateIndex, "(?<=-).+?(?=\\.)")
    return(x)
})
MutateDat1 <- map(MutateDat, t)
names(MutateDat1) <- names(targetGene)
write.csv(MutateDat1[[2]], "HNRMutation.csv", quote = F)
col_funMut <- colorRamp2(c(0, 0.03, 0.06), c("#2185C5", "#FFF6E5","#FF7F66"))
an_col <- c("#DDD5AF", "#F4666D", "#4A81AA")
names(an_col) <- unique(genetype$Group)
row_an <-  HeatmapAnnotation(type = genetype$Group, ##注释信息的内容。
                             show_annotation_name = F, ## 是否显示注释的标题
                             col = list(type = an_col), ## 注释信息的颜色
                             show_legend = T,  ## 是否显示注释信息的说明
                             annotation_legend_param = list(title = "Group", nrow = 1), ## 注释信息图例的标题
                             which = "row" #对行或者列进行注释 
)
row_an
p1 <- Heatmap(MutateDat1[[1]], name = "Mutation Frequency", 
              heatmap_legend_param = list(direction  = "horizontal"),
              col = col_funMut,rect_gp = gpar(col = "grey", lwd = 1), left_annotation = row_an,
              row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
pdf(paste0("./oxy/figure 2/figure2A.pdf"), height = 5, width = 8)
draw(p1, heatmap_legend_side = "bottom")
dev.off()


# figure 2B 某一个肿瘤的突变 ------------------------------------------------------

library(TCGAbiolinks)
mut <- GDCquery_Maf(tumor = "UCEC", pipelines = "mutect2")
library(maftools)
maf <- read.maf(maf = mut)
col = RColorBrewer::brewer.pal(n = 9, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')
pdf("./oxy/figure 2/figure2B.pdf", width = 8, height = 5)
oncoplot(maf = maf, genes = geneFamily1, colors = col,writeMatrix = T)
dev.off()


# figure 2C ccle突变数据结果 ----------------------------------------------------
library(circlize)
CCLEMutation <- data.table::fread("./HNRNPL/CCLE_mutations.csv",
                                  data.table = F)
CCLEInfo <- data.table::fread("./HNRNPL/sample_info.csv",
                              data.table = F)
CCLECount <- table(CCLEInfo$lineage)
OxyMutate <- CCLEMutation %>% 
    filter(Hugo_Symbol %in% targetGene[[1]]) %>% 
    select(Hugo_Symbol, Tumor_Sample_Barcode) %>% 
    left_join(CCLEInfo[c("DepMap_ID", "lineage")], 
              by = c("Tumor_Sample_Barcode" = "DepMap_ID"))
OxyMutate1 <- OxyMutate[OxyMutate$lineage != "",]
OxyCount <- OxyMutate1 %>% group_by(Hugo_Symbol, lineage) %>% 
    dplyr::count()
OxyCount$All <- CCLECount[OxyCount$lineage]
OxyCount$Percent <- OxyCount$n/OxyCount$All*100

TargetTissue <- c("breast", "gastric", "colorectal", "kidney",
                  "lung", "bone", "ovary", "skin", "fibroblast",
                  "liver")
CCLECount[TargetTissue]
OxyHeat <- OxyCount %>% select(1:2,5) %>% 
    filter(lineage %in% TargetTissue) %>% 
    spread(key = "lineage", value = "Percent")
OxyHeat[is.na(OxyHeat)] <- 0
OxyHeat <- as.data.frame(OxyHeat)
OxyHeat <- column_to_rownames(OxyHeat, var = "Hugo_Symbol")
OxyHeat1 <- sapply(OxyHeat, function(x){
    as.vector(x)/100
})
rownames(OxyHeat1) <- rownames(OxyHeat)
write.csv(OxyHeat1, "./oxy/figure 2/OxyHeat1.csv", quote = F)
OxyHeat2 <- t(OxyHeat1)
### 绘图
pdf("./oxy/figure 2/figure 2C.pdf", width = 5,height = 5)
col_fun = colorRamp2(c(0, 0.2), c("#FFFFBF", "red"))
circos.par(canvas.xlim =c(-1.1,1.1),
           canvas.ylim = c(-1.2,1.2),
           cell.padding = c(0, 0, 0, 0),
           track.height = 0.08,
           gap.degree = 5)
circos.initialize(factors = rownames(OxyHeat1), xlim = c(0,1))
circos.track(rownames(OxyHeat1), ylim = c(0,1),
             #bg.border = NA, 
             panel.fun = function(x,y){
                 circos.text(CELL_META$xcenter,
                             CELL_META$cell.ylim[2] + uy(5,"mm"),
                             CELL_META$sector.index,
                             cex = 0.6,
                             niceFacing = T
                 )
                 sector.index = CELL_META$sector.index
                 m = OxyHeat2[[sector.index]]
                 col_mat = col_fun(m)
                 circos.rect(0,0,1,1,border = col_mat,
                             col = col_mat)
             }
)
for(i in 2:10){
    circos.track(rownames(HNRCCLEHeat1), ylim = c(0,1),
                 #bg.border = NA, 
                 panel.fun = function(x,y){
                     # circos.text(CELL_META$xcenter,
                     #             CELL_META$cell.ylim[2] + uy(1,"mm"),
                     #             CELL_META$sector.index,
                     #             cex = 0.5,
                     # )
                     sector.index = CELL_META$sector.index
                     m = HNRCCLEHeat2[[sector.index]][i]
                     col_mat = col_fun(m)
                     circos.rect(0,0,1,1,border = col_mat,
                                 col = col_mat)
                 }
    )
}
#circos.rect()
#circos.rect(borders())
#plot(HNRCCLEHeat1[,1], HNRCCLEHeat1[,2])
circos.clear()
dev.off()



# figure 2D 拷贝数变异 ---------------------------------------------------------
cnvindex <- Sys.glob("~/project/TCGADat/CNVDatCategory/*")
CNVDat <- list()
num <- 0
for(i in targetGene){
    gainRes <- c()
    lostRes <- c()
    num <- num + 1
    for (j in cnvindex) {
        example <- vroom::vroom(gzfile(j))
        example1 <- example %>% filter(`Gene Symbol` %in% i) %>% gather(key = "sampleID", value = "value", -1) %>% 
            spread(key = "Gene Symbol", value = "value")
        sampleSum <- nrow(example1)
        gain <- sapply(example1[-1], function(x){
            sum(x > 0)/sampleSum
        })
        lost <- sapply(example1[-1], function(x){
            sum(x < 0)/sampleSum
        })
        gainRes <- rbind(gainRes, gain)
        lostRes <- rbind(lostRes, lost)
    }
    together <- list(gain = gainRes, lost = lostRes)
    together <- lapply(together, function(x){
        rownames(x) <- str_extract(cnvindex, "(?<=gory\\/).+?(?=_)")
        return(x)
    })
    CNVDat[[num]] <- together
}
CNVDat1 <- lapply(CNVDat, function(x) lapply(x, t))
names(CNVDat1) <- names(targetGene)
CNVDatLong <- lapply(CNVDat1, function(x){
    aaa <- lapply(x, function(y){
        y %>% as.data.frame() %>% rownames_to_column(var = "Gene") %>% 
            gather(key = "Cancer", value = "CNVFreq", -Gene)
    })
    res <- cbind(aaa[[1]], aaa[[2]][,3])
    colnames(res)[3:4] <- c("Gain", "Loss")
    return(res)
})
write.csv(CNVDatLong[[1]], "./oxy/figure 2/HNRCNV.csv", row.names = F, quote = F)
UpColor <- colorRamp2(breaks = c(0, 0.5), colors = c("#FFF6E5","#FF7F66"))
DnColor <- colorRamp2(breaks = c(0, 0.5), colors = c("#FFF6E5","#2185C5"))   
DiagFunc <- function(up, down){
    function(j, i, x, y, width, height, fill){
        grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width), 
                     unit.c(y - 0.5*height, y + 0.5*height, y + 0.5*height),
                     gp = gpar(fill = DnColor(down[i, j]), col = "grey")) 
        grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width), 
                     unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height),
                     gp = gpar(fill = UpColor(up[i, j]), col = "grey"))
    }
}

dat <- CNVDat1[[1]]
an_col <- c("#DDD5AF", "#F4666D", "#4A81AA")
names(an_col) <- unique(genetype$Group)
row_an <-  HeatmapAnnotation(type = genetype$Group, ##注释信息的内容。
                             show_annotation_name = F, ## 是否显示注释的标题
                             col = list(type = an_col), ## 注释信息的颜色
                             show_legend = T,  ## 是否显示注释信息的说明
                             annotation_legend_param = list(title = "Group", nrow = 1), ## 注释信息图例的标题
                             which = "row" #对行或者列进行注释 
)
p1 <- Heatmap(dat$gain, column_title = "Copy number variation across cancer types", 
              rect_gp = gpar(type = "none"), show_heatmap_legend = F,  
              cluster_rows = F, cluster_columns = T, ##绘制空的热图框
              left_annotation = row_an, ##添加左侧注释信息
              cell_fun = DiagFunc(up = dat$gain, down = dat$lost) ## 绘制表格内的内容
)

col_fun = colorRamp2(c(-0.5, 0, 0.5), c("#2185C5", "#FFF6E5","#FF7F66")) ##自定义颜色信息
lgd <- Legend(title = "CNV Frequency", ## 注释的标题
              col_fun = col_fun, ## 注释的颜色
              at = c(-0.5,-0.25,0,0.25,0.5), ## 注释刻度的分组
              labels = c("0.5","Loss","0","Gain","0.5"),  ## 注释刻度的重命名
              direction = "horizontal" ## 注释的方向
) 
pdf(paste0("./oxy/figure 2/figure2D.pdf"), height = 5, width = 8)
draw(p1, annotation_legend_list = lgd, annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom", merge_legend = TRUE)
dev.off()




# panCancer数据提取 -----------------------------------------------------------

options(stringsAsFactors = F)
load("~/project/PanCaner/FigureYa71ssGSEA/easy_input_immunity.rdata")
library(clusterProfiler)
immunitGym <- lapply(immunity, function(x){
    result <- bitr(x, 
                   fromType="ENTREZID", toType=c("SYMBOL","ENSEMBLPROT"), OrgDb="org.Hs.eg.db") 
    return(result)
})
immnuGym <- plyr::ldply(immunitGym, data.frame)
names(immnuGym)[1] <- "type"
genecode23 <- read.delim("~/project/gencode.v23.metadata.HGNC", header = F)
genecode23 <- genecode23[!duplicated(genecode23$V1),]
immnuGym1 <- left_join(immnuGym, genecode23, by = c("SYMBOL" = "V2"))
immnuGym2 <- immnuGym1[!duplicated(immnuGym1$SYMBOL),]
immnuGym3 <- immnuGym2[!is.na(immnuGym2$V1),]
targetgene <- c(geneFamily1, geneFamily2)
targetgeneENSID <- genecode23[genecode23$V2 %in% targetgene,]
targetgeneENSID <- targetgeneENSID[!duplicated(targetgeneENSID$V2),]
targetgeneENSID1 <- targetgeneENSID$V1

targetid <- c(immnuGym3$V1, targetgeneENSID1)


library(data.table)
rownum <- nrow(fread("tcga_Kallisto_tpm", select = 1L, header = F))
rownum


## 确定有多个个15000
index <- seq(0,rownum, by = 15000)
##构建循环
result <- c()
for(i in index){
    dat <- fread("tcga_Kallisto_tpm", nrows = 15000,skip = i, header = F, data.table = F)
    result1 <- dat[dat$V1 %in% targetid, ]
    result <- rbind(result, result1)
}
result[1:5,1:5]
colna <- fread("tcga_Kallisto_tpm", nrows = 1, header = F, data.table = F)
colna1 <- unlist(colna[1,])
colnames(result) <- colna1
colnames(result)[1] <- "ensid"
rownames(result) <- result$ensid
result1 <- result[,-1]
result1 <- apply(result1, 2, as.numeric)
rownames(result1) <- rownames(result)
result1 <- 2^result1 - 0.001
resultImmue <- result[immnuGym3$V1,]
resultImmue <- resultImmue[!is.na(resultImmue$ensid),]
rownames(immnuGym3) <- immnuGym3$V1
resultImmue$entrez <- immnuGym3[resultImmue$ensid,]$ENTREZID
resultImmue <- resultImmue[-1]
rownames(resultImmue) <- NULL
resultImmue <- column_to_rownames(resultImmue, var = "entrez")
library(GSVA)
tcga_gsva <- as.data.frame(t(gsva(as.matrix(resultImmue), immunity, method = "ssgsea")))
write.csv(tcga_gsva, "immune.csv", quote = F)
write.csv(t(ImmnuDat), "ImmuDat.csv", quote = F)

resultTargetGene <- result[targetgeneENSID$V1,]
resultTargetGene <- resultTargetGene[,-1]
resultTargetGene[1:5,1:5]


# figure 3A  --------------------------------------------------------------
### 数据提取
genecode23 <- read.delim("~/project/gencode.v23.metadata.HGNC", header = F)
genecode23 <- genecode23[!duplicated(genecode23$V1),]
targetgeneENSID <- genecode23[genecode23$V2 %in% rownames(genetype),]
targetgeneENSIDNoDuplicated <- targetgeneENSID[!duplicated(targetgeneENSID$V2),]
rownames(targetgeneENSIDNoDuplicated) <- targetgeneENSIDNoDuplicated$V1
library(data.table)
rownum <- nrow(fread("~/project/TCGADat/panCancer/tcga_Kallisto_tpm", select = 1L, header = F))
rownum

## 确定有多个个15000
index <- seq(0,rownum, by = 15000)
##构建循环
OxyDat <- c()
for(i in index){
    dat <- fread("~/project/TCGADat/panCancer/tcga_Kallisto_tpm", nrows = 15000,skip = i, header = F, data.table = F)
    result1 <- dat[dat$V1 %in% targetgeneENSIDNoDuplicated$V1, ]
    OxyDat <- rbind(OxyDat, result1)
}
OxyDat[1:5,1:5]
colna <- fread("~/project/TCGADat/panCancer/tcga_Kallisto_tpm", nrows = 1, header = F, data.table = F)
colna1 <- unlist(colna[1,])
colnames(OxyDat) <- colna1
colnames(OxyDat)[1] <- "ensid"
rownames(OxyDat) <- OxyDat$ensid
OxyDat <- OxyDat[,-1]
OxyDat1 <- OxyDat
OxyDat1 <- apply(OxyDat1, 2, as.numeric)
rownames(OxyDat1) <- rownames(OxyDat)
OxyDat1 <- 2^OxyDat1 - 0.001
OxyDat1 <- as.data.frame(OxyDat1)
OxyDat1$gsm <- targetgeneENSIDNoDuplicated[rownames(OxyDat1),]$V2
OxyDat1 <- OxyDat1[!duplicated(OxyDat1$gsm),]
rownames(OxyDat1) <- OxyDat1$gsm
OxyDat1 <- OxyDat1[,-10536]
OxyDatT <- as.data.frame(t(OxyDat1))
### 导入cacnerhallexpr数据
cancerHall <- read.csv("~/project/TCGADat/panCancer/PanCancerHall.csv", row.names = 1)
library(psych)
CorRes <- corr.test(cancerHall, OxyDatT, method = "spearman")
CorResR <- as.data.frame(CorRes$r)
CorResP <- as.data.frame(CorRes$p)
write.csv(CorResR, "./oxy/figure3/CancerCorGeneR.csv", quote = F)
write.csv(CorResP, "./oxy/figure 3/CancerCorGeneP.csv", quote = F)

### 到处ctytoscape数据
CorResR[CorResP >= 0.05 | abs(CorResR) < 0.3] <- NA
class(CorResR)
CytCor <- rownames_to_column(CorResR, var = "PathWay") %>% 
    gather(key = "gene", value = "R2", -PathWay) %>% filter(!is.na(R2))
CytCor$PathWay <- str_split(CytCor$PathWay, pattern = "_", n = 2,simplify = T)[,2]
write.csv(CytCor, "./oxy/figure 3/CytCor.csv", quote = F, row.names = F)

colnames(CytCor) <- c("from", "to", "R2")
CytCor1 <- CytCor %>% mutate(group = ifelse(R2 > 0, "Pos", "Neg"),
                             ModiPa = ifelse(R2 > 0, paste0("Pos_", from),
                                             paste0("Neg_", from))) %>% 
    select(2:4, from = ModiPa)
head(CytCor1)
### 定义各个点的位置
CytCorInfo <- CytCor1 %>% select(from, to) %>% 
    gather(key = "Direction", value = "Iterm") %>% distinct() %>% 
    mutate(Group = ifelse(str_detect(Iterm, "^Pos"), "Pos",
                          ifelse(str_detect(Iterm, "^Neg"), "Neg", "Gene")))
write.csv(CytCor1, "./oxy/figure 3/CytCor.csv", row.names = F, quote = F)
CytCorInfo$Group1 <- CytCorInfo$Group
CytCorInfo$Group1[28:35] <- genetype[CytCorInfo$Iterm[28:35],]$Group
write.csv(CytCorInfo, "./oxy/figure 3/CytCorInfo.csv", quote = F, row.names = F)

# figure 3B做图 -------------------------------------------------------------
NumCount <- map2(CorResR, CorResP, function(x, y){
    c(sum(x > 0.3 & y < 0.05), sum(x < -0.3 & y < 0.05))
})
NumCount1 <- do.call(rbind, NumCount)
colnames(NumCount1) <- c("PosCount", "NegCount")
NumCount1 <- as.data.frame(NumCount1)
NumCount1$gene <- rownames(NumCount1)
NumCount1$NegCount <- -NumCount1$NegCount

NumCount1Long <- NumCount1 %>% filter(PosCount != 0 | NegCount != 0) %>% 
    gather(key = "Group", value = "sum", -gene)
NumCount1Long$Group <- ifelse(NumCount1Long$Group == "PosCount", 
                              "Positive", "Negative")
NumCount1Long$Group <- factor(NumCount1Long$Group, 
                              levels = c("Positive", "Negative"))
write.csv(NumCount1Long, "./oxy/figureCancerCount.csv", row.names = F, quote = F)
ggplot(NumCount1Long, aes(gene, sum, fill = Group)) + 
    geom_bar(stat = "identity") + theme_classic() + 
    scale_y_continuous(breaks = c(-10,-5,0,5,10), 
                       labels = c(10,5,0,5,10)) + 
    ylab("Number of Pathway")  +  xlab(NULL) + 
    scale_fill_manual(values = c("#FB8072","#80B1D3")) + 
    labs(fill = NULL) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 12), 
          axis.title.y = element_text(size = 13),
          legend.position = "top") 
ggsave("./oxy/figure 3/PathwayCount.pdf", width = 6, height = 4)


# figure 3C 基因内的相关分析图 -----------------------------------------------------
MarkerName <- str_split(CytCor1$from, pattern = "_", n = 2, simplify = T)
sort(table(MarkerName[,2]),decreasing = T)
TargetMarker <- c("PROTEIN_SECRETION", "UV_RESPONSE_DN", "ANDROGEN_RESPONSE")
MarkerR <- CorResR[paste0("HALLMARK_", TargetMarker),] %>% t() %>% 
    as.data.frame() %>% rownames_to_column(var = "Gene")
colnames(MarkerR) <- c( "Gene", TargetMarker)
MarkerRLong <- gather(MarkerR, key = "Pathway", value = "r", -Gene)
MarkerP <- CorResP[paste0("HALLMARK_", TargetMarker),] %>% t() %>% 
    as.data.frame() %>% rownames_to_column(var = "Gene")
colnames(MarkerP) <-  c( "Gene", TargetMarker)
MarkerPLong <- gather(MarkerP, key = "Pathway", value = "p", -Gene)
identical(MarkerPLong$Pathway, MarkerRLong$Pathway)
identical(MarkerPLong$Gene, MarkerRLong$Gene)
MarkCorDf <- cbind(MarkerRLong, MarkerPLong$p)
MarkCorDf <- MarkCorDf[c(2,1,3,4)]
colnames(MarkCorDf) <- c("x", "y", "r", "p")

exaDat <- OxyDatT
colnames(exaDat) <- letters[1:12]
aaa <- letters[1:12]
names(aaa) <- colnames(OxyDatT)
ExadatAn <- MarkCorDf
ExadatAn$y <- aaa[MarkCorDf$y]
TargetMarkerCode <- paste0("Group", 1:3)
names(TargetMarkerCode) <- TargetMarker
ExadatAn$x <- TargetMarkerCode[MarkCorDf$x]
library(ggcor)
corr <- fortify_cor(exaDat, type = "upper", show.diag = TRUE,
                    cor.test = TRUE, cor.test.method = "spearman",
                    cluster.type = "all")
xname <- c("e", "f", "h", "c", "l", "a", "g", 
           "b", "d", "i", "j", "k")
names(xname) <- 1:12
yname <- rev(xname)
names(yname) <- 1:12
corr1 <- corr
corr1$x <- xname[corr$x]
corr1$y <- yname[corr$y]
corr1$x <- sapply(corr1$x, function(x) names(which(aaa == x)))
corr1$y <- sapply(corr1$y, function(x) names(which(aaa == x)))
write.csv(corr1, "GeneCor.csv", row.names = F, quote = F)
collll <- c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", 
            "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")
ggcor(corr, xlim = c(-10, 30), fill.colours = rev(collll))  +  
    add_link(ExadatAn, diag.label = T, mapping = aes(colour  = r)) +
    add_diaglab(angle = 45)  + geom_pie2() +
    scale_color_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#F46D43") + 
    remove_axis("y") 
ggsave("corres.pdf", width = 6, height = 6)



# figure 3D 基因和免疫浸润的关系 ----------------------------------------------------
ImmnueTargetDatT <- read.csv("~/project/TCGADat/panCancer/ImmnueTargetDatT.csv",
                             row.names = 1)
ImmuCor <- corr.test(ImmnueTargetDatT, OxyDatT, method = "spearman")
ImmuCorR <- ImmuCor$r
ImmuCorP <- ImmuCor$p
write.csv(ImmuCorP, "./oxy/figure 3/ImmuCorP.csv", quote = F)
write.csv(ImmuCorR, "./oxy/figure 3/ImmuCorR.csv", quote = F)
ImmuCorR[ImmuCorP >= 0.05 | abs(ImmuCorR) < 0.3] <- NA
ImmuCorR <- as.data.frame(ImmuCorR)
ImmuCorRLong <- ImmuCorR %>% rownames_to_column(var = "ImmuGene") %>% 
    gather(key = "gene", value = "R2", -ImmuGene) %>% filter(!is.na(R2)) %>% 
    mutate(Group = 1)
colnames(ImmuCorRLong) <- c("from", "to", "R2", "Group")
write.csv(ImmuCorRLong, "./oxy/figure 3/ImmuCorRLong.csv", row.names = F, quote = F)
SelfGene <- tibble(from = "", to = c(unique(ImmuCorRLong$to)), 
                   R2 = "", Group = 2)
ImmuEdge <- rbind(ImmuCorRLong, SelfGene)

###制作node的注释信息
ImmuDegree <- table(ImmuCorRLong$from)
ImmueNodeInfo <- tibble(node = names(ImmuDegree),
                        Group = 1,
                        Degree = as.vector(ImmuDegree)
)
immnuGym3 <- read.csv("~/project/TCGADat/panCancer/immnuGym3.csv")
ImmueNodeInfo1 <- left_join(ImmueNodeInfo, immnuGym3[c(1,3)],
                            by = c("node" = "SYMBOL"))
ImmueNodeInfo1$type[51] <- "CD8 T cells"
OxyNodeDegree <- table(ImmuCorRLong$to)
OxyNodeInfo <- tibble(node = c(names(OxyNodeDegree),""),
                      Group = 2, 
                      Degree = c(as.vector(OxyNodeDegree),0),
                      type = c(names(OxyNodeDegree),""),)
OxyNodeGene <- unique(OxyNodeInfo$node)
NodeInfo <- rbind(ImmueNodeInfo1, OxyNodeInfo)
head(NodeInfo)
write.csv(NodeInfo, "./oxy/figure 3/NodeInfo.csv", row.names = F, quote = F)
library(igraph)
library(ggraph)
graph <- graph_from_data_frame(ImmuEdge, vertices = NodeInfo)
layout <- create_layout(graph, layout = 'circle')
ImmuGene <- read.csv("~/project/TCGADat/panCancer/ImmuGene.csv")
ImmuGene <- unlist(ImmuGene$x)
ImmuGeneNUm <- length(ImmuGene)
##设置指定颜色
TargetCol <- read.table("~/project/TCGADat/panCancer/allcol_minus.txt")
TargetCol <- TargetCol[[1]]
ImmueCol <- c(Assembly = "#efb306", C = '#FFAA5D', `C++` = "#e8351e", JavaScript = "#cd023d", 
              Java = "#852f88", R = "#4e54ac", Python = "#0f8096", Ruby = "#7db954", 
              SQL = "#17a769")
# library(rPlotter)
# mycol<-extract_colours("ColorFig.jpeg",num_col = 30)

OxyCol <- c(an_col[genetype[unique(NodeInfo$type)[11:20],]$Group],"white")
NodeCOl <- c(ImmueCol, OxyCol)
names(NodeCOl) <- unique(NodeInfo$type)
EdgeCol <- NodeCOl[10:30]
outer_circle <- layout %>%
    filter(Group == 1) %>%
    mutate(type = factor(type, names(NodeCOl)[1:9])) %>%
    arrange(type, desc(name)) %>%
    mutate(
        x = cos((row_number() - 1) / ImmuGeneNUm * 4.65 * pi)*0.8,
        y = sin((row_number() - 1) / ImmuGeneNUm * 4.65 * pi)*0.8
    )

inner_circle <- layout %>%
    filter(Group ==2) %>%
    #mutate(type = factor(type, names(ImmueCol)[-10])) %>%
    #arrange(type, desc(name)) %>%
    mutate(
        x = cos((row_number() - 1) / 10 * 2 * pi) * 0.5,
        y = sin((row_number() - 1) / 10 * 2 * pi) * 0.5
    )
layout[] <- bind_rows(outer_circle, inner_circle) %>%
    arrange(.ggraph.index)
p1 <- ggraph(layout) +
    geom_edge_diagonal0(
        aes(edge_color = node1.type),
        edge_width = 0.3, show.legend = FALSE,
        edge_alpha = 0.2
    ) + 
    geom_node_point(
        aes(size = Degree, color = type,),
        alpha = 0.6, show.legend = F
    ) + 
    geom_node_text(
        aes(
            x = 1.0175 * x,
            y = 1.0175 * y,
            label = name,
            angle = -((-node_angle(x, y) + 90) %% 180) + 90,
            filter = !(name %in% OxyNodeGene)
        ),
        size = 1.7, hjust = 'outward') + 
    geom_node_text(
        aes(
            x = x,
            y = y,
            label = name,
            filter = name %in% OxyNodeGene
        ),
        size = 2
    ) + 
    scale_edge_color_manual(values = NodeCOl[1:10]) +
    scale_color_manual(values = NodeCOl) + 
    scale_size_area(max_size = 10) + 
    scale_edge_alpha_manual(values = c(0.15, 1)) +
    theme_void() 
p1
p1  + annotate("text", x = 0.7, y = 0.2, label = "CD8 T cells", 
               color = "#efb306", angle = 300#, vjust = 0.4
) + 
    annotate("text", x = 0.2, y = 0.7, label = "Macrophages",
             color = "#eb990c", angle = -10) + 
    annotate("text", x = -0.45, y = 0.55, label = "T helper cells",
             color = "#e8351e", angle =30) + 
    annotate("text", x = -0.65, y = 0.3, label = "NK cells",
             color = "#cd023d", angle =30) + 
    annotate("text", x = -0.7, y = -0.2, label = "Th1 cells",
             color = "#852f88", angle =-60) + 
    annotate("text", x = 0, y = -0.7, label = "T cells",
             color = "#4e54ac", angle =0) + 
    annotate("text", x = 0.35, y = -0.6, label = "B cells",
             color = "#0f8096", angle =10) + 
    annotate("text", x = 0.58, y = -0.4, label = "Neutrophils",
             color = "#7db954", angle =30) + 
    annotate("text", x = 0.7, y = -0.2, label = "Mast cells",
             color = "#17a769", angle =30, size = 3) + 
    xlim(c(-0.9,0.9)) + ylim(c(-0.9,0.9))

ggsave("./oxy/figure 3/ImmuNetWork.pdf", width = 6, height = 6)

# figure 4A 生存分析热图 --------------------------------------------------------

SurDatIndex <- x("~/project/TCGADat/MedianSurv/*")
CoxRes <- list()
num <- 0
for(i in targetGene){
    num <- num + 1
    Surres <- list()
    for(j in SurDatIndex){
        id <- str_extract(j, "(?<=Surv/)\\w+(?=Median)")
        surDat <- read.csv(gzfile(j))
        surDat <- surDat[-1]
        targetSurDat <- surDat[surDat$gsm %in% i, ]
        Surres[[id]] <- targetSurDat
    }
    CoxRes[[num]] <- Surres
}
HNRCox <- plyr::ldply(CoxRes[[2]], data.frame)
write.csv(HNRCox, "HNRCox.csv", row.names = F, quote = F)
CoxRes1 <- lapply(CoxRes, function(x){
    res <- plyr::ldply(x, data.frame)
    names(res)[1] <- "CancerType"
    return(res)
})
CoxResHR <- lapply(CoxRes1, function(x){
    x %>% select(CancerType:HR) %>% mutate(HR1 = round(HR, 2)) %>% 
        select(-HR) %>% spread(key = "CancerType", value = "HR1") %>% 
        column_to_rownames(., var = "gsm")
})
write.csv(CoxResHR[[1]], "./oxy/figure 4/CoxResHR.csv", quote = F)
CoxResP <- lapply(CoxRes1, function(x){
    x %>% select(CancerType:gsm, p.val) %>% 
        spread(key = "CancerType", value = "p.val") %>% 
        column_to_rownames(., var = "gsm")
})
write.csv(CoxResP[[1]], "./oxy/figure 4/CoxResP.csv", quote = F)

CoxResCate <- map2(CoxResHR, CoxResP, function(data1, data2){
    logfcCat <- apply(data1, 2, function(x){
        cut(x, breaks = c(-1,  1, Inf),
            labels = c("Low Risk", "High Risk"))
    })
    logfcCat[data2 >= 0.05] <- "P >= 0.05"
    rownames(logfcCat) <- rownames(data1)
    return(logfcCat)
})
col_cat <- c("High Risk" = "#ED5E57",  "Low Risk" = "#6B9AB7", "P >= 0.05" = "#d9d9d9")

cell_fun <- function(logfc, dataP, logfcCutoff = 1, PCutoff = 0.05, 
                     darkcol = "black", lightcol = "white", digit = 2, fontsize  = 6){
    function(j, i, x, y, width, height, fill){
        if(dataP[i,j] < PCutoff){
            grid.text(format(round(logfc, digit)[i, j], nsmall = digit), x, y, 
                      gp = gpar(fontsize = fontsize, col  = lightcol))
        }else{
            grid.text(format(round(logfc, digit)[i, j], nsmall = digit), x, y, 
                      gp = gpar(fontsize = fontsize, col  = darkcol))
        }
    }
}
names(CoxResCate) <- names(targetGene)

pdf("OxySur.pdf", width = 8, height = 4)
CoxResCate[[1]] <- CoxResCate[[1]][genetype$Gene.Symbol,]
row_an <-  HeatmapAnnotation(type = genetype$Group, ##注释信息的内容。
                             show_annotation_name = F, ## 是否显示注释的标题
                             col = list(type = an_col), ## 注释信息的颜色
                             show_legend = T,  ## 是否显示注释信息的说明
                             annotation_legend_param = list(title = "Group"), ## 注释信息图例的标题
                             which = "row" #对行或者列进行注释 
)
pdf("./oxy/figure 4/OxySur.pdf", width = 8, height = 4)
Heatmap(CoxResCate[[1]], 
        rect_gp = gpar(lwd = 1, col = "grey"), col = col_cat,
        left_annotation = row_an)
dev.off()
write.csv(CoxResCate[[1]], "./oxy/figure 4/CoxResCate.csv", quote = F)
# figure 4B 森林图绘图---------------------------------------------------------------

### 绘制森林图
#### 绘制HNRNPC
### 统计各个癌的样本数
CoxTpmIndex <- str_subset(tpmindex, 
                          paste(COxCancerId, collapse = "|"))
sampleN <- lapply(CoxTpmIndex, function(x){
    dat1 <- read.csv(gzfile(x), nrows = 1)
    Sample <- sum(str_detect(colnames(dat1), "0..$"))
    Cancer <- str_extract(x, "(?<=tpmDat/).+?(?=_)")
    c(Cancer, Sample)
})
SampleDat <- do.call(rbind, sampleN)
SampleDat1 <- as.data.frame(SampleDat)
colnames(SampleDat1) <- c(".id",  "Sample")
### 绘制HNRNPC的森林图
HIF1ACox <- lapply(CoxRes[[1]], function(x){
    x[x$gsm == "HIF1A",]
})
HIF1ACox1 <- plyr::ldply(HIF1ACox, data.frame)
COxCancerId <- HIF1ACox1$`.id`
HIF1AJoin <- left_join(HIF1ACox1, SampleDat1, by  = ".id")
rownames(HIF1AJoin) <- HIF1AJoin$.id

HIF1AJoin$HR <- round(HIF1AJoin$HR,2)
HIF1AJoin$low95 <- round(HIF1AJoin$low95, 2)
HIF1AJoin$up95 <- round(HIF1AJoin$up95, 2)
HIF1AFor <- data.frame(Cancer = rownames(HIF1AJoin),
                       Patients = HIF1AJoin$Sample,
                       HR = paste0(format(HIF1AJoin$HR, nsmall = 2), 
                                   "(", format(HIF1AJoin$low95, nsmall = 2),
                                   "-", format(HIF1AJoin$up95, nsmall = 2),
                                   ")"))

CoxUl2 <- cbind(c("\nCancer", NA,NA, HIF1AFor$Cancer, NA),
                c("\nPatients", NA, NA, HIF1AFor$Patients, NA),
                c("Hazard Ratio\n(95% CI)", NA, NA, HIF1AFor$HR, NA))
HIF1AJoin$Sample <- as.numeric(HIF1AJoin$Sample)
HIF1AJoin$Percent <- (HIF1AJoin$Sample/10165)*100
pdf("./oxy/figure 4/HIF1A.pdf", width=5,height = 8)
forestplot(labeltext=CoxUl2, #图中的文本
           mean=c(NA,NA,1,HIF1AJoin$HR, NA),#HR
           lower=c(NA,NA,1,HIF1AJoin$low95, NA), #95%置信区间下限
           upper=c(NA,NA,1,HIF1AJoin$up95, NA),#95%置信区间上限
           #title="Hazard Ratio of different echocardiography indicator for CVD",
           graph.pos=3,#图在表中的列位置
           clip = c(0,2),
           graphwidth =  "auto",#unit(.4,"npc"),#图在表中的宽度比例
           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
           col=fpColors(box="#fc9272", lines="#969696", zero = "#969696",
                        axes = "#969696"),#box颜色
           boxsize=c(NA,NA,NA,HIF1AJoin$Percent,NA)/10,#box大小根据样本量设置
           lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           #grid = structure(c(HNRNPCJoin$HR[2]), gp = gpar(col = "black", lty=2,lwd=2)),#虚线(可多条)及其横坐标、颜色、线宽
           xticks = c(0, 0.5, 1,1.5, 2),#横坐标刻度根据需要可随意设置
           lwd.xaxis=2,#X轴线宽
           xlab= "HIF1A",#X轴标题
           hrzl_lines=list("3" = gpar(lwd=2, col="#969696"),#第三行顶部加黑线，引号内数字标记行位置
                           #"4" = gpar(lwd=60,lineend="butt", columns=c(1:4), col="#99999922"),#加阴影，弱项不建议使用
                           "36" = gpar(lwd=2, col="#969696")),#最后一行底部加黑线,""中数字为nrow(data)+5
           txt_gp=fpTxtGp(label=gpar(cex=0.7),#各种字体大小设置
                          ticks=gpar(cex=0.6),
                          xlab=gpar(cex = 0.6),
                          title=gpar(cex = 0.6)),
           #is.summary = c(T,rep(F,27)),#首行字体类型设置
           lineheight = unit(.5,"cm"),#固定行高
           align=c("l","c","c"),#每列文字的对齐方式，偶尔会用到
           cex=10,
           colgap = unit(0,"cm"),#列间隙
           mar=unit(rep(1.25, times = 4), "cm"),#图形页边距
           new_page = T#是否新页
)      

dev.off()

### 绘制VHL的森林图
VHLCox <- lapply(CoxRes[[1]], function(x){
    x[x$gsm == "VHL",]
})
HNRNPCox1 <- plyr::ldply(VHLCox, data.frame)
COxCancerId <- HNRNPCox1$`.id`
HNRNPCJoin <- left_join(HNRNPCox1, SampleDat1, by  = ".id")
rownames(HNRNPCJoin) <- HNRNPCJoin$.id

HNRNPCJoin$HR <- round(HNRNPCJoin$HR,2)
HNRNPCJoin$low95 <- round(HNRNPCJoin$low95, 2)
HNRNPCJoin$up95 <- round(HNRNPCJoin$up95, 2)
HNRNPFor <- data.frame(Cancer = rownames(HNRNPCJoin),
                       Patients = HNRNPCJoin$Sample,
                       HR = paste0(format(HNRNPCJoin$HR, nsmall = 2), 
                                   "(", format(HNRNPCJoin$low95, nsmall = 2),
                                   "-", format(HNRNPCJoin$up95, nsmall = 2),
                                   ")"))

CoxUl2 <- cbind(c("\nCancer", NA,NA, HNRNPFor$Cancer, NA),
                c("\nPatients", NA, NA, HNRNPFor$Patients, NA),
                c("Hazard Ratio\n(95% CI)", NA, NA, HNRNPFor$HR, NA))
HNRNPCJoin$Sample <- as.numeric(HNRNPCJoin$Sample)
HNRNPCJoin$Percent <- (HNRNPCJoin$Sample/10165)*100



pdf("./oxy/figure 4/VHL.pdf", width=5,height = 8)
forestplot(labeltext=CoxUl2, #图中的文本
           mean=c(NA,NA,1,HNRNPCJoin$HR, NA),#HR
           lower=c(NA,NA,1,HNRNPCJoin$low95, NA), #95%置信区间下限
           upper=c(NA,NA,1,HNRNPCJoin$up95, NA),#95%置信区间上限
           #title="Hazard Ratio of different echocardiography indicator for CVD",
           graph.pos=3,#图在表中的列位置
           clip = c(0,2),
           graphwidth =  "auto",#unit(.4,"npc"),#图在表中的宽度比例
           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
           col=fpColors(box="#fc9272", lines="#969696", zero = "#969696",
                        axes = "#969696"),#box颜色
           boxsize=c(NA,NA,NA,HNRNPCJoin$Percent,NA)/10,#box大小根据样本量设置
           lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           #grid = structure(c(HNRNPCJoin$HR[2]), gp = gpar(col = "black", lty=2,lwd=2)),#虚线(可多条)及其横坐标、颜色、线宽
           xticks = c(0, 0.5, 1,1.5, 2),#横坐标刻度根据需要可随意设置
           lwd.xaxis=2,#X轴线宽
           xlab= "HNRNPA1",#X轴标题
           hrzl_lines=list("3" = gpar(lwd=2, col="#969696"),#第三行顶部加黑线，引号内数字标记行位置
                           #"4" = gpar(lwd=60,lineend="butt", columns=c(1:4), col="#99999922"),#加阴影，弱项不建议使用
                           "36" = gpar(lwd=2, col="#969696")),#最后一行底部加黑线,""中数字为nrow(data)+5
           txt_gp=fpTxtGp(label=gpar(cex=0.7),#各种字体大小设置
                          ticks=gpar(cex=0.6),
                          xlab=gpar(cex = 0.6),
                          title=gpar(cex = 0.6)),
           #is.summary = c(T,rep(F,27)),#首行字体类型设置
           lineheight = unit(.5,"cm"),#固定行高
           align=c("l","c","c"),#每列文字的对齐方式，偶尔会用到
           cex=10,
           colgap = unit(0,"cm"),#列间隙
           mar=unit(rep(1.25, times = 4), "cm"),#图形页边距
           new_page = T#是否新页
)      
dev.off()



# figure 4D 聚类生存分析 ----------------------------------------------------------
KIRCindex <- str_subset(tpmindex, "KIRC")

KIRCDat <- vroom::vroom(gzfile(KIRCindex))
colnames(KIRCDat)[1] <- "id"
KIRCDat <- as.data.frame(KIRCDat)
rownames(KIRCDat) <- KIRCDat$id
cancerid <- str_detect(colnames(KIRCDat), "0..$")
LIHCCancerDat <- KIRCDat[cancerid]

SurDIffGene <- c("EPAS1", "ARNT", "ARNT2", "EGLN1", 
                 "CREBBP", "EP300", "HIF1AN")
SurDiffTarget <- idmap[idmap$gene %in% SurDIffGene,]

ClusterDat <- LIHCCancerDat[SurDiffTarget$id,]
ClusterDat$gsm <- SurDiffTarget[rownames(ClusterDat),]$gene
rownames(ClusterDat) <- NULL
ClusterDat <- column_to_rownames(ClusterDat, var = "gsm")

LIHCSurInfo <- read.delim(gzfile("~/project/TCGADat/survivalDat/TCGA-KIRC.survival.tsv.gz"))
rownames(LIHCSurInfo) <- LIHCSurInfo$sample
ToId <- intersect(colnames(ClusterDat), LIHCSurInfo$sample)
library(ClassDiscovery)
ClusterDat1 <- ClusterDat[,ToId]
hcs <- hclust(distanceMatrix(as.matrix(ClusterDat1), "pearson"), "ward.D") 
group <- cutree(hcs,k=2)
LIHCSurInfo1 <- LIHCSurInfo[ToId,]
LIHCSurInfo1$Clust2 <- paste0("C",group[rownames(LIHCSurInfo1)])
library(survival)
library(survminer)
my.surv <- Surv(LIHCSurInfo1$X_OS, LIHCSurInfo1$X_OS_IND)
fit <- survfit(my.surv ~ LIHCSurInfo1$Clust2)

data.survdiff <- survdiff(my.surv ~ LIHCSurInfo1$Clust2)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")

p2 <- ggsurvplot(fit, data = LIHCSurInfo1 ,
           #ggtheme = theme_bw(), #想要网格就运行这行
           conf.int = F, #不画置信区间，想画置信区间就把F改成T
           conf.int.style = "step",#置信区间的类型，还可改为ribbon
           censor = T, #不显示观察值所在的位置
           palette = "npg", #线的颜色对应高、低
           risk.table = T,
           risk.table.col="strata",
           #legend.title = i,#基因名写在图例题目的位置
           font.legend = 11,#图例的字体大小
           #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
           
           # #在图例上标出高低分界点的表达量，和组内sample数量
           # legend.labs=c("HNRNP High",
           #               "HNPNP Low"),
           
           #在左下角标出pvalue、HR、95% CI
           #太小的p value标为p < 0.001
           pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))
p2
ggsave("./oxy/figure 4/ClusterSur.pdf", width = 6,height = 6)

# figure 4C 热图绘制 ----------------------------------------------------------

phenoDat <- read.delim(gzfile("~/project/TCGADat/phenotypeData/TCGA-KIRC.GDC_phenotype.tsv.gz"))
rownames(phenoDat) <- phenoDat$submitter_id.samples
phenoDat1 <- phenoDat[ToId,]
phenoDat2 <- phenoDat1[c( "submitter_id.samples",
                          "age_at_initial_pathologic_diagnosis", "gender.demographic"
)]
phenoDat3 <- inner_join(phenoDat2, LIHCSurInfo1, 
                        by = c("submitter_id.samples" = "sample"))
phenoDat4 <- phenoDat3[c(1:4,9)]
colnames(phenoDat4) <- c("sample", "Age", "Gender", "Event", "Group")
phenoDat4$Group <- ifelse(phenoDat4$Group == "C1", "High", "Low")
phenoDat4$Event <- ifelse(phenoDat4$Event == 0, "Alive", "Dead")
phenoDat4$Gender <- str_to_title(phenoDat4$Gender)

phenoDat4 <- column_to_rownames(phenoDat4, var = "sample")
n <- t(scale(t(ClusterDat1)))
n[n>2] = 2
n[n<-2] = -2
annColors <- list()
annColors[["Gender"]] <- c("Male"="#8DD3C7","Female"="#BEBADA")
annColors[["Event"]] <- c("Alive"="#FDB462","Dead"="#D9D9D9")
annColors[["Group"]] <- c("High"="#FB8072","Low"="#80B1D3")
annColors[["Age"]] <- c("#D0D1E6", "#3690C0")
col <- colorRampPalette(c("#3C7DAF", "#EAF4F1","#FFFCBA", "#E83140"))(20)
#Heatmap(n)
library(pheatmap)
pheatmap(n, show_rownames = T, show_colnames = F, annotation_col = phenoDat4,
         cluster_cols = hcs,
         color = col, 
         annotation_colors = annColors)
library(ggplotify)
g <- as.ggplot(Heat)
p2 <- as.ggplot(p2)
cowplot::plot_grid(g, p2, ncol = 1)
