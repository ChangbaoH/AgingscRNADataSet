##############packages########
install.packages('Seurat')
install.packages('dplyr')
install.packages('readr')
install.packages('tidyverse')
install.packages('dplyr')
install.packages('Matrix')
install.packages('cowplot')
install.packages('patchwork')
install.packages('ggplot2')
install.packages('igraph')

library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library("igraph")
library("readr")
library("dplyr")
library(Seurat)
#############胰腺
expr <- read_tsv("/data03/HCB/aging/Pancreas/Adult_Pancreas/exprMatrix.tsv")
expr <- as.data.frame(expr)
meta <- read_tsv("/data03/HCB/aging/Pancreas/Adult_Pancreas/meta.tsv")
meta <- as.data.frame(meta)
meta <- cbind(meta,rep(NA,length(rownames(meta))))
colnames(meta)[length(colnames(meta))] <- "Sample_ID"
for (i in 1:length(rownames(meta))) {
  meta$Sample_ID[i] <- strsplit(x=as.character(meta[i,1]),split='_')[[1]][2]
}
meta1 <- meta
meta1 <- distinct(meta1,patient_ID,age,Sample_ID,.keep_all = T)
for (i in 1:length(rownames(meta1))) {
  meta_temp <- meta[which(meta[,2]==meta1[i,2]&meta[,10]==meta1[i,10]),]
  temp <- expr[,meta_temp[,1]]
  write.csv(temp,paste("/data03/HCB/aging/Pancreas/counts/",meta1[i,2],"_",meta1[i,10],"_",meta_temp[i,4],"years.csv",sep = ""))
}

expr1 <- read_tsv("/data03/HCB/aging/Pancreas/Neonatal_Pancreas/exprMatrix.tsv",row.names=1)
expr1 <- as.data.frame(expr1)
meta <- read_tsv("/data03/HCB/aging/Pancreas/Neonatal_Pancreas/meta.tsv")
meta <- as.data.frame(meta)
meta1 <- meta
meta1 <- distinct(meta1,patient_ID,.keep_all = T)

for (i in 1:length(rownames(meta1))) {

  meta_temp <- meta[which(meta[,3]==meta1[i,3]),]
  temp <- expr1[,meta_temp[,1]]
  colnames(temp) <- paste(colnames(temp),"youth",sep = "_")
  if(i==1){  
	scdata <- temp
  }else{
  
  }
  write.csv(temp,paste("/data03/HCB/aging/Pancreas/counts/",meta1[i,3],"_0_years.csv",sep = ""))
}
#############视网膜
#### 2. 读入原始表达数据 ####
bm1 <- Read10X("/data03/HCB/aging/retina/data3/")
bm2 <- Read10X("~/Project/others/小张聊科研/单细胞课程/Day2课程资料/BoneMarrow/BM2/")
colnames(bm1) <- paste(colnames(bm1),"BM1",sep = "_")
colnames(bm2) <- paste(colnames(bm2),"BM2",sep = "_")

barcodes <- read_tsv("barcodes.tsv")
bar <- barcodes[,1:2]
##############blood
gene <- read_tsv('')
gene <- as.data.frame(gene)
UMI1[1:5,1:5]
UMI1 <- cbind(rep(NA,length(rownames(UMI1))),UMI1)
colnames(UMI1)[1] <- 'gene'
UMI1[,1] <- gene[rownames(UMI1),2] 
barcodes <- read.table("03.Cell.Barcodes.txt")
UMI1 <- as.data.frame(UMI1)


SC1 <- gsub("-",".",c('gene',barcodes[which(barcodes[,2]=="SC1"),1]))  
temp <- UMI1[,SC1]
write.csv(temp,"/data03/HCB/aging/blood/SC1_110s.csv")
SC2 <- gsub("-",".",c('gene',barcodes[which(barcodes[,2]=="SC2"),1]))  
temp <- UMI1[,SC2]
write.csv(temp,"/data03/HCB/aging/blood/SC2_110s.csv")
SC3 <- gsub("-",".",c('gene',barcodes[which(barcodes[,2]=="SC3"),1]))  
temp <- UMI1[,SC3]
write.csv(temp,"/data03/HCB/aging/blood/SC3_110s.csv")
SC4 <- gsub("-",".",c('gene',barcodes[which(barcodes[,2]=="SC4"),1]))  
temp <- UMI1[,SC4]
write.csv(temp,"/data03/HCB/aging/blood/SC4_110s.csv")
SC5 <- gsub("-",".",c('gene',barcodes[which(barcodes[,2]=="SC5"),1]))  
temp <- UMI1[,SC5]
write.csv(temp,"/data03/HCB/aging/blood/SC5_110s.csv")
SC6 <- gsub("-",".",c('gene',barcodes[which(barcodes[,2]=="SC6"),1]))  
temp <- UMI1[,SC6]
write.csv(temp,"/data03/HCB/aging/blood/SC6_110s.csv")
SC7 <- gsub("-",".",c('gene',barcodes[which(barcodes[,2]=="SC7"),1]))  
temp <- UMI1[,SC7]
write.csv(temp,"/data03/HCB/aging/blood/SC7_110s.csv")

CT1 <- gsub("-",".",c('gene',barcodes[which(barcodes[,2]=="CT1"),1]))  
temp <- UMI1[,CT1]
write.csv(temp,"/data03/HCB/aging/blood/CT1_50s.csv")
CT2 <- gsub("-",".",c('gene',barcodes[which(barcodes[,2]=="CT2"),1]))  
temp <- UMI1[,CT2]
write.csv(temp,"/data03/HCB/aging/blood/CT2_70s.csv")
CT3 <- gsub("-",".",c('gene',barcodes[which(barcodes[,2]=="CT3"),1]))  
temp <- UMI1[,CT3]
write.csv(temp,"/data03/HCB/aging/blood/CT3_60s.csv")
CT4 <- gsub("-",".",c('gene',barcodes[which(barcodes[,2]=="CT4"),1]))  
temp <- UMI1[,CT4]
write.csv(temp,"/data03/HCB/aging/blood/CT4_70s.csv")
CT5 <- gsub("-",".",c('gene',barcodes[which(barcodes[,2]=="CT5"),1]))  
temp <- UMI1[,CT5]
write.csv(temp,"/data03/HCB/aging/blood/CT5_80s.csv")

prefetch SRR1482463 -O output
prefetch SRR1482463 -O output
prefetch SRR1482463 -O output



library(Seurat)
library(SeuratData)
library(SeuratDisk)
# Converting from AnnData to Seurat via h5Seurat
Convert('Full_obj_raw_counts_nosoupx.h5ad', "h5seurat",overwrite = TRUE,assay = "RNA")
scRNA <- LoadH5Seurat("Full_obj_raw_counts_nosoupx.h5seurat")
save(scRNA,file="/data/HCB/aging/ileum/raw_SeuratObject.RData")


Convert('pediatric_RAWCOUNTS_cellxgene.h5ad', "h5seurat",overwrite = TRUE,assay = "RNA")
scRNA1 <- LoadH5Seurat("pediatric_RAWCOUNTS_cellxgene.h5seurat")
save(scRNA1,file="/data/HCB/aging/ileum/raw_SeuratObject1.RData")

Convert('global_raw.h5ad', "h5seurat",overwrite = TRUE,assay = "RNA")
scRNA <- LoadH5Seurat("global_raw.h5seurat")
save(scRNA,file="/data/HCB/aging/heart/raw_SeuratObject.RData")


gunzip GSM4504182_LZ003matrix.csv.gz
gunzip GSM4504184_LZ005matrix.csv.gz
gunzip GSM4504185_LZ007matrix.csv.gz
gunzip GSM4504186_LZ008matrix.csv.gz
gunzip GSM4504187_LZ009matrix.csv.gz
gunzip GSM4504189_LZ011matrix.csv.gz
gunzip GSM4504191_LZ013matrix.csv.gz
gunzip GSM4504192_LZ014matrix.csv.gz
gunzip GSM4504193_LZ015matrix.csv.gz
gunzip GSM4504194_LZ016matrix.csv.gz
##########################################################
temp1 <- read.csv("/data03/HCB/aging/Pancreas/counts/AFES365_11_53years.csv")
temp2 <- read.csv("/data03/HCB/aging/Pancreas/counts/st19061908_0_years.csv")
temp3 <- read.csv("/data03/HCB/aging/Pancreas/counts/TUM_C1_24_77years.csv")


temp1[1:5,1:5]
temp2[1:5,1:5]
temp3[1:5,1:5]
















####  Code Description              ####
#---  1. Written by WoLin @ 2019.03.15，last update 19.07.14 ---#
#---  2. Analysis for single sample    ---#
#---  3. Support 10X data & expression matrix ---#
#---  4. Need to change:sample data, dims ---#

#### 1. 加载分析使用的工具包 ####
install.packages('Seurat')
install.packages('ggplot2')
install.packages('cowplot')
install.packages('Matrix')
install.packages('dplyr')

library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)

#### 2. 读入原始表达数据 ####
#以下两种方式二选一
#10X 数据
bm1 <- Read10X("/data03/HCB/aging/test/BoneMarrow/BM1/")
bm2 <- Read10X("/data03/HCB/aging/test/BoneMarrow/BM2/")
colnames(bm1) <- paste(colnames(bm1),"BM1",sep = "_")
colnames(bm2) <- paste(colnames(bm2),"BM2",sep = "_")

#将所有读入的数据合并成一个大的矩阵
#合并时需注意行名一致
#既有10X的数据又有表达矩阵的数据，全部转换为表达矩阵再进行合并
#关于矩阵合并请见单独的矩阵合并脚本“merge_matrix.R”
experiment.data <- cbind(bm1,bm2)

#创建一个文件夹用于写分析结果
sam.name <- "result"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

#### 3. 创建Seurat分析对象 ####
experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project = "scRNA tutorial", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")
#将数据写到文件中一边后续分析使用
save(experiment.aggregate,file=paste0("./",sam.name,"/",sam.name,"_raw_SeuratObject.RData"))

#### 4. 数据概览 & QC ####
#查看SeuratObject中的对象
slotNames(experiment.aggregate)
#assay
experiment.aggregate@assays
experiment.aggregate@meta.data
experiment.aggregate@active.assay
experiment.aggregate@active.ident
experiment.aggregate@neighbors
experiment.aggregate@reductions
experiment.aggregate@images
experiment.aggregate@project.name
experiment.aggregate@misc
experiment.aggregate@version
experiment.aggregate@commands
experiment.aggregate@tools

aging_counts<-experiment.aggregate[["RNA"]]@counts

experiment.aggregate@meta.data[1:5,]

#细胞及细胞中基因与RNA数量
dim(experiment.aggregate@meta.data)
View(experiment.aggregate@meta.data)

##QC：统计线粒体基因在每个细胞中的占比
experiment.aggregate[["percent.mt"]] <- PercentageFeatureSet(experiment.aggregate,pattern = "^MT-")
pdf(paste0("./",sam.name,"/QC-VlnPlot.pdf",width = 8,height = 4.5))
VlnPlot(experiment.aggregate, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)
dev.off()

##QC：统计基因数，RNA，线粒体基因分布
gene.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nFeature_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nCount_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$percent.mt,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq))
colnames(freq.combine) <- c("Count_Gene","Count_RNA","MT_percent")
write.table(freq.combine,file = paste0(sam.name,"/QC-gene_frequency.txt"),quote = F,sep = "\t")
rm(gene.freq,rna.freq,mt.freq)

##QC：基因数与线粒体基因以及RNA数量的分布相关性
plot1 <- FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(paste0("./",sam.name,"/QC-FeatureScatter.pdf"),width = 8,height = 4.5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()
rm(plot1,plot2)

#### 5. 筛选细胞 ####
cat("Before filter :",nrow(experiment.aggregate@meta.data),"cells\n")
experiment.aggregate <- subset(experiment.aggregate, 
                               subset = 
                                 nFeature_RNA > 500 & 
                                 nCount_RNA > 1000 & 
                                 nCount_RNA < 20000 &
                                 percent.mt < 5)
cat("After filter :",nrow(experiment.aggregate@meta.data),"cells\n")

#### 6. 表达量标准化 ####
experiment.aggregate <- NormalizeData(experiment.aggregate, 
                                      normalization.method = "LogNormalize",
                                      scale.factor = 10000)

#计算表达量变化显著的基因FindVariableFeatures
experiment.aggregate <- FindVariableFeatures(experiment.aggregate, 
                                             selection.method = "vst",
                                             nfeatures = 2000)

#展示标准化之后的整体表达水平
top10 <- head(x = VariableFeatures(experiment.aggregate), 10)
plot1 <- VariableFeaturePlot(experiment.aggregate)
plot2 <- LabelPoints(plot = plot1, points = top10)
pdf(file = paste0(sam.name,"/Norm-feature_variable_plot.pdf"),width = 8,height = 5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()

#### 7. 均一化与PCA ####
#均一化（需要一点时间）
experiment.aggregate <- ScaleData(
  object = experiment.aggregate,
  do.scale = FALSE,
  do.center = FALSE,
  vars.to.regress = c("percent.mt"))

#PCA计算
experiment.aggregate <- RunPCA(object = experiment.aggregate, 
                               features = VariableFeatures(experiment.aggregate),
                               verbose = F)

#PCA结果展示-1
pdf(paste0("./",sam.name,"/PCA-VizDimLoadings.pdf"),width = 7,height = 5)
VizDimLoadings(experiment.aggregate, dims = 1:2, reduction = "pca")
dev.off()

#PCA结果展示-2
pdf(paste0("./",sam.name,"/PCA-DimPlot.pdf"),width = 5,height = 4)
DimPlot(experiment.aggregate, reduction = "pca")
dev.off()

#PCA结果展示-3
pdf(paste0("./",sam.name,"/PCA-DimHeatmap.pdf"),width = 5,height = 4)
DimHeatmap(experiment.aggregate, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()

#### 8. 确定细胞类群分析PC ####
#耗时较久
experiment.aggregate <- JackStraw(experiment.aggregate, num.replicate = 100,dims = 40)
experiment.aggregate <- ScoreJackStraw(experiment.aggregate, dims = 1:40)
pdf(paste0("./",sam.name,"/PCA-JackStrawPlot_40.pdf"),width = 6,height = 5)
JackStrawPlot(object = experiment.aggregate, dims = 1:40)
dev.off()
#确定用于细胞分群的PC
dim.use <- 1:9

#### 9. 细胞分群TSNE算法 ####
#TSNE算法
experiment.aggregate <- FindNeighbors(experiment.aggregate, dims = dim.use)
experiment.aggregate <- FindClusters(experiment.aggregate, resolution = 0.5)

experiment.aggregate <- RunTSNE(experiment.aggregate, dims = dim.use, do.fast = TRUE)
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = experiment.aggregate, pt.size=0.5,label = T)
dev.off()
write.table(experiment.aggregate@meta.data,file = paste0("./",sam.name,"/",sam.name,"_cells_details_tsne_",max(dim.use),"PC.txt"))

#按照数据来源分组展示细胞异同
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_SamGroup_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = experiment.aggregate, group.by="orig.ident", pt.size=0.5)
dev.off()

#### 10. 计算marker基因 ####
all.markers <- FindAllMarkers(experiment.aggregate, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/",sam.name,"_marker_genes_tsne_",max(dim.use),"PC.txt"),sep="\t",quote = F)

# 遍历每一个cluster然后展示其中前4个基因
marker.sig <- all.markers %>% filter(p_val_adj <= 0.05)
for(cluster in unique(marker.sig$cluster)){
  cluster.markers <- FindMarkers(experiment.aggregate, ident.1 = cluster, min.pct = 0.25)
  cluster.markers <- as.data.frame(cluster.markers) %>% 
    mutate(Gene = rownames(cluster.markers))
  cl4.genes <- cluster.markers %>% arrange(desc(avg_logFC))
  cl4.genes <- cl4.genes[1:min(nrow(cl4.genes),4),"Gene"]
  
  #VlnPlot
  pvn <- VlnPlot(experiment.aggregate, features = cl4.genes)
  pdf(paste0("./",sam.name,"/MarkerGene-VlnPlot_cluster",cluster,"_tsne_",max(dim.use),"PC.pdf"),width = 7,height = 6)
  print(pvn)
  dev.off()
  
  #Feater plot 
  pvn <- FeaturePlot(experiment.aggregate,features=cl4.genes)
  pdf(paste0("./",sam.name,"/MarkerGene-FeaturePlot_cluster",cluster,"_tsne_",max(dim.use),"PC.pdf"),width = 7,height = 6)
  print(pvn)
  dev.off()
}
rm(cl4.genes,cluster.markers,pvn)

#热图展示Top marker基因
top5 <- marker.sig %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
pdf(paste0("./",sam.name,"/MarkerGene-Heatmap_all_cluster_tsne_",max(dim.use),"PC.pdf"),width = 8,height = 5)
DoHeatmap(experiment.aggregate, features = top5$gene,size = 2)
dev.off()

install.packages("harmony")
library(harmony)

2. Seurat结合Harmony的整合流程
这种整合方法很简单，而且占内存少，速度快

library(harmony)

bm1 <- Read10X("/data03/HCB/aging/test/BoneMarrow/BM1/")
bm2 <- Read10X("/data03/HCB/aging/test/BoneMarrow/BM2/")
colnames(bm1) <- paste(colnames(bm1),"BM1",sep = "_")
colnames(bm2) <- paste(colnames(bm2),"BM2",sep = "_")
testdf <- cbind(bm1,bm2)
test.seu <- CreateSeuratObject(counts = testdf) %>%
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData()
test.seu <- RunPCA(test.seu, npcs = 50, verbose = FALSE)
test.seu@meta.data$patient=str_replace(test.seu$orig.ident,"_.*$","")

先运行Seurat标准流程到PCA这一步，然后就是Harmony整合，可以简单把这一步理解为一种新的降维

test.seu=test.seu %>% RunHarmony("patient", plot_convergence = TRUE)

> test.seu
An object of class Seurat 
33538 features across 6746 samples within 1 assay 
Active assay: RNA (33538 features)
 2 dimensional reductions calculated: pca, harmony
1
2
3
4
5
6
接着就是常规聚类降维，都是基于Harmony的Embeddings矩阵

test.seu <- test.seu %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

test.seu <- test.seu %>% 
  RunTSNE(reduction = "harmony", dims = 1:30)
1
2
3
4
5
6
7
看看效果

p3 <- DimPlot(test.seu, reduction = "tsne", group.by = "patient", pt.size=0.5)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
p4 <- DimPlot(test.seu, reduction = "tsne", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
fig_tsne <- plot_grid(p3, p4, labels = c('patient','ident'),align = "v",ncol = 2)
ggsave(filename = "tsne3.pdf", plot = fig_tsne, device = 'pdf', width = 27, height = 12, units = 'cm')
————————————————
版权声明：本文为CSDN博主「TOP生物信息」的原创文章，遵循CC 4.0 BY-SA版权协议，转载请附上原文出处链接及本声明。
原文链接：https://blog.csdn.net/qq_38774801/article/details/112292947


rm(list = ls())
bm1 <- Read10X("/data03/HCB/aging/test/BoneMarrow/BM1/")
bm2 <- Read10X("/data03/HCB/aging/test/BoneMarrow/BM2/")
colnames(bm1) <- paste(colnames(bm1),"BM1",sep = "_")
colnames(bm2) <- paste(colnames(bm2),"BM2",sep = "_")
experiment.data <- cbind(bm1,bm2)
scRNA <- CreateSeuratObject(
  experiment.data,
  project = "scRNA tutorial", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")

scRNA@meta.data$orig.ident
scRNA@meta.data$percent.mt
scRNA@meta.data$percent.rb
scRNA@meta.data$percent.HB

cellinfo <- subset(scRNA@meta.data, select = c("orig.ident", "percent.mt", "percent.rb", "percent.HB"))
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 



cellinfo <- subset(scRNA@meta.data, select = c("orig.ident", "percent.mt", "percent.rb", "percent.HB"))
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
2.2 数据标准化（和锚点整合不同，不需拆分样本，直接标准化）
### SCT标准化数据
scRNA <- SCTransform(scRNA)
2.3 使用harmony整合数据

### PCA
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)

### 整合方法1：单个样本间进行整合（推荐，效果更好）
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 
# group.by.vars参数是设置按哪个分组来整合
# max.iter.harmony设置迭代次数，默认是10。运行RunHarmony结果会提示在迭代多少次后完成了收敛。
#⚠️RunHarmony函数中有个lambda参数，默认值是1，决定了Harmony整合的力度。lambda值调小，整合力度变大，反之。（只有这个参数影响整合力度，调整范围一般在0.5-2之间）


作者：Hayley笔记
链接：http://events.jianshu.io/p/159df5c4ff21
来源：简书
著作权归作者所有。商业转载请联系作者获得授权，非商业转载请注明出处。
2.1 准备数据
rm(list = ls())

scRNA <- readRDS("scRNA.rds")
cellinfo <- subset(scRNA@meta.data, select = c("orig.ident", "percent.mt", "percent.rb", "percent.HB"))
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
2.2 数据标准化（和锚点整合不同，不需拆分样本，直接标准化）
### SCT标准化数据
scRNA <- SCTransform(scRNA)
2.3 使用harmony整合数据

### PCA
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)

### 整合方法1：单个样本间进行整合（推荐，效果更好）
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 
# group.by.vars参数是设置按哪个分组来整合
# max.iter.harmony设置迭代次数，默认是10。运行RunHarmony结果会提示在迭代多少次后完成了收敛。
#⚠️RunHarmony函数中有个lambda参数，默认值是1，决定了Harmony整合的力度。lambda值调小，整合力度变大，反之。（只有这个参数影响整合力度，调整范围一般在0.5-2之间）


作者：Hayley笔记
链接：http://events.jianshu.io/p/159df5c4ff21
来源：简书
著作权归作者所有。商业转载请联系作者获得授权，非商业转载请注明出处。







########################################################################################################################################Seurat 标准流程
#标准 Seurat 工作流采用原始的单细胞表达数据，旨在数据中查找clusters。此过程包括数据标准化和高变基因选择、数据归一化、高变基因的PCA、共享近邻图形的构建以及使用模块优化进行聚类。最后，我们使用 t-SNE 在二维空间中可视化我们的clusters。
pbmc <- experiment.aggregate


pbmc.counts <- Read10X(data.dir = "~/Downloads/pbmc3k/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.counts)
pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc)
pbmc <- FindNeighbors(object = pbmc)
pbmc <- FindClusters(object = pbmc)
pbmc <- RunTSNE(object = pbmc)
DimPlot(object = pbmc, reduction = "tsne")
#########################################################################################################################################Seurat对象交互
#自 Seurat v4.0 以来，对 Seurat 对象进行了改进，并增加了用户交互的新方法。还为常见任务引入简单的功能，例如取子集和合并。

# Get cell and feature names, and total numbers
colnames(x = pbmc)
Cells(object = pbmc)
rownames(x = pbmc)
ncol(x = pbmc)
nrow(x = pbmc)
# Get cell identity classes
Idents(object = pbmc)
levels(x = pbmc)

# Stash cell identity classes
pbmc[["old.ident"]] <- Idents(object = pbmc)
pbmc <- StashIdent(object = pbmc, save.name = "old.ident")

# Set identity classes
Idents(object = pbmc) <- "CD4 T cells"
Idents(object = pbmc, cells = 1:10) <- "CD4 T cells"

# Set identity classes to an existing column in meta data
Idents(object = pbmc, cells = 1:10) <- "orig.ident"
Idents(object = pbmc) <- "orig.ident"

# Rename identity classes
pbmc <- RenameIdents(object = pbmc, `CD4 T cells` = "T Helper cells")
# Subset Seurat object based on identity class, also see ?SubsetData
subset(x = pbmc, idents = "B cells")
subset(x = pbmc, idents = c("CD4 T cells", "CD8 T cells"), invert = TRUE)

# Subset on the expression level of a gene/feature
subset(x = pbmc, subset = MS4A1 > 3)

# Subset on a combination of criteria
subset(x = pbmc, subset = MS4A1 > 3 & PC1 > 5)
subset(x = pbmc, subset = MS4A1 > 3, idents = "B cells")

# Subset on a value in the object meta data
subset(x = pbmc, subset = orig.ident == "Replicate1")

# Downsample the number of cells per identity class
subset(x = pbmc, downsample = 100)
# Merge two Seurat objects
merge(x = pbmc1, y = pbmc2)
# Merge more than two Seurat objects
merge(x = pbmc1, y = list(pbmc2, pbmc3))
################################################################################################################################################数据访问
##在 Seurat 中访问数据很简单，使用明确定义的取子集代码可以快速查找所需的数据。

# View metadata data frame, stored in object@meta.data
pbmc[[]]

# Retrieve specific values from the metadata
pbmc$nCount_RNA
pbmc[[c("percent.mito", "nFeature_RNA")]]

# Add metadata, see ?AddMetaData
random_group_labels <- sample(x = c("g1", "g2"), size = ncol(x = pbmc), replace = TRUE)
pbmc$groups <- random_group_labels
# Retrieve or set data in an expression matrix ('counts', 'data', and 'scale.data')
GetAssayData(object = pbmc, slot = "counts")
pbmc <- SetAssayData(object = pbmc, slot = "scale.data", new.data = new.data)
# Get cell embeddings and feature loadings
Embeddings(object = pbmc, reduction = "pca")
Loadings(object = pbmc, reduction = "pca")
Loadings(object = pbmc, reduction = "pca", projected = TRUE)
# FetchData can pull anything from expression matrices, cell embeddings, or metadata
FetchData(object = pbmc, vars = c("PC_1", "percent.mito", "MS4A1"))
#############################################################################################################################################Seurat的可视化
#默认情况下，所有绘图功能都将返回 ggplot2 绘图，从而允许使用 ggplot2 轻松定制。

# Dimensional reduction plot for PCA or tSNE
DimPlot(object = pbmc, reduction = "tsne")
DimPlot(object = pbmc, reduction = "pca")

# Dimensional reduction plot, with cells colored by a quantitative feature
FeaturePlot(object = pbmc, features = "MS4A1")

# Scatter plot across single cells, replaces GenePlot
FeatureScatter(object = pbmc, feature1 = "MS4A1", feature2 = "PC_1")
FeatureScatter(object = pbmc, feature1 = "MS4A1", feature2 = "CD3D")

# Scatter plot across individual features, repleaces CellPlot
CellScatter(object = pbmc, cell1 = "AGTCTACTAGGGTG", cell2 = "CACAGATGGTTTCT")

VariableFeaturePlot(object = pbmc)

# Violin and Ridge plots
VlnPlot(object = pbmc, features = c("LYZ", "CCL5", "IL32"))
RidgePlot(object = pbmc, feature = c("LYZ", "CCL5", "IL32"))

# Heatmaps
DoHeatmap(object = pbmc, features = heatmap_markers)
DimHeatmap(object = pbmc, reduction = "pca", cells = 200)

# New things to try!  Note that plotting functions now return ggplot2 objects, so you can add themes, titles, and options
# onto them
VlnPlot(object = pbmc, features = "MS4A1", split.by = "groups")
DotPlot(object = pbmc, features = c("LYZ", "CCL5", "IL32"), split.by = "groups")
FeaturePlot(object = pbmc, features = c("MS4A1", "CD79A"), blend = TRUE)
DimPlot(object = pbmc) + DarkTheme()
DimPlot(object = pbmc) + labs(title = "2,700 PBMCs clustered using Seurat and viewed\non a two-dimensional tSNE")
##Seurat 提供了许多预制的主题，可添加到 ggplot2 中，以便快速定制

#主题	功能
DarkTheme
设置带有白色文本的黑色背景
FontSize
为绘图的各个元素设置字体大小
NoAxes
删除轴和轴文本
NoLegend
删除所有图例元素
RestoreLegend
删除后恢复图例
RotatedAxis
旋转 x 轴标签
# Plotting helper functions work with ggplot2-based scatter plots, such as DimPlot, FeaturePlot, CellScatter, and
# FeatureScatter
plot <- DimPlot(object = pbmc) + NoLegend()

# HoverLocator replaces the former `do.hover` argument It can also show extra data throught the `information` argument,
# designed to work smoothly with FetchData
HoverLocator(plot = plot, information = FetchData(object = pbmc, vars = c("ident", "PC_1", "nFeature_RNA")))

# FeatureLocator replaces the former `do.identify`
select.cells <- FeatureLocator(plot = plot)

# Label points on a ggplot object
LabelPoints(plot = plot, points = TopCells(object = pbmc[["pca"]]), repel = TRUE)
######################################################################################################################################不同assay之间切换
#使用 Seurat，您可以轻松地在单细胞水平的不同assays 之间切换（例如来自 CITE-seq 的 ADT 计数，或整合/批次校正的数据）。大多数功能现在需要assays 参数，但可以设置默认assay以避免重复语句。

cbmc <- CreateSeuratObject(counts = cbmc.rna)
# Add ADT data
cbmc[["ADT"]] <- CreateAssayObject(counts = cbmc.adt)
# Run analyses by specifying the assay to use
NormalizeData(object = cbmc, assay = "RNA")
NormalizeData(object = cbmc, assay = "ADT", method = "CLR")

# Retrieve and set the default assay
DefaultAssay(object = cbmc)
DefaultAssay(object = cbmc) <- "ADT"
DefaultAssay(object = cbmc)

# Pull feature expression from both assays by using keys
FetchData(object = cbmc, vars = c("rna_CD3E", "adt_CD3"))

# Plot data from multiple assays using keys
FeatureScatter(object = cbmc, feature1 = "rna_CD3E", feature2 = "adt_CD3")
###########################################################################################################################Seurat v2.X与v4.X命令的区别
Seurat v2.X	Seurat v4.X
object@data	GetAssayData(object = object)[1]
object@raw.data	GetAssayData(object = object, slot = "counts")
object@scale.data	GetAssayData(object = object, slot = "scale.data")
object@cell.names	colnames(x = object)[2]
rownames(x = object@data)[3]	rownames(x = object)[4]
object@var.genes	VariableFeatures(object = object)[5]
object@hvg.info	HVFInfo(object = object)[6]
object@assays$assay.name	object[["assay.name"]]
object@dr$pca	object[["pca"]]
GetCellEmbeddings(object = object, reduction.type = "pca")	Embeddings(object = object, reduction = "pca")
GetGeneLoadings(object = object, reduction.type = "pca")	Loadings(object = object, reduction = "pca")
AddMetaData(object = object, metadata = vector, col.name = "name")	object$name <- vector
object@meta.data$name	object$name
object@idents	Idents(object = object)[7]
SetIdent(object = object, ident.use = "new.idents")	Idents(object = object) <- "new.idents"
SetIdent(object = object, cells.use = 1:10, ident.use = "new.idents")	Idents(object = object, cells = 1:10) <- "new.idents"
StashIdent(object = object, save.name = "saved.idents")	object$saved.idents <- Idents(object = object)
levels(x = object@idents)[8]	levels(x = object)[9]
RenameIdent(object = object, old.ident.name = "old.ident", new.ident.name = "new.ident")	RenameIdents(object = object, "old.ident" = "new.ident")
WhichCells(object = object, ident = "ident.keep")	WhichCells(object = object, idents = "ident.keep")
WhichCells(object = object, ident.remove = "ident.remove")	WhichCells(object = object, idents = "ident.remove", invert = TRUE)
WhichCells(object = object, max.cells.per.ident = 500)[10]	WhichCells(object = object, downsample = 500)[11]
WhichCells(object = object, subset.name = "name", low.threshold = low, high.threshold = high)	WhichCells(object = object, expression = name > low & name < high)[12]
FilterCells(object = object, subset.names = "name", low.threshold = low, high.threshold = high)	subset(x = object, subset = name > low & name < high)[13]
SubsetData(object = object, subset.name = "name", low.threshold = low, high.threshold = high)	subset(x = object, subset = name > low & name < high)[14]
MergeSeurat(object1 = object1, object2 = object2)	merge(x = object1, y = object2)[15]






pdf(paste0("./",result.name,"/subcluster-DotPlot.pdf"),width = 28,height = 12)
DotPlot(exp.seurat, features = c("CSF1R", "CD74", "NRGN", "VCAN","FLT1","AMBP", "SPI1", "MRC1", "TMEM119", "CX3CR1", "CLDN5","VTN","GLUL","SOX9","AQP4","GJA1","NDRG2","GFAP","ALDH1A1","ALDH1L1", "VIM", "PTGDS", "PDGFRA", "PCDH15", "OLIG2", "OLIG1", "PLP1", "MAG", "MOG", "MOBP", "MBP", "SATB2", "SLC17A7", "SLC17A6", "GAD1", "GAD2", "SLC32A1", "SNAP25", "STMN2", "RBFOX3"),cols = c("blue", "red"))
dev.off()


exp.seurat <- RenameIdents(exp.seurat, `0` = "excitatory neurons", `1` = "excitatory neurons", `2` = "excitatory neurons", `3` = "excitatory neurons", `4` = "inhibitory neurons", `5` = "excitatory neurons", `6` = "inhibitory neurons", `7` = "excitatory neurons", `8` = "inhibitory neurons", `9` = "astrocytes", `10` = "excitatory neurons", `11` = "inhibitory neurons", `12` = "oligodendrocytes", `13` = "excitatory neurons", `14` = "astrocytes", `15` = "OPCs", `16` = "excitatory neurons", `17` = "endothelial cells", `18` = "inhibitory neurons", `19` = "astrocytes")
save(exp.seurat,file=paste0("./",result.name,"/",result.name,"_Celltype_SeuratObject_final.RData"))