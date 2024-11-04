rm(list = ls())   
options(stringsAsFactors = F)


################	Differentially Methylated Probe (DMP) Analysis############
###P_adj < 0.05 and log|FC| > 0.2###
library("FactoMineR")
library("factoextra")
library(clusterProfiler)
library(ggplot2)
library(limma)
library(org.Hs.eg.db)
library(pheatmap)
library(ChAMP)

#load raw data#
testDir="X/x/x"
myLoad <- champ.load(testDir,arraytype="EPIC")
save(myLoad,file='myLoad.Rdata')
pd<-myLoad$pd[,3:5]
rownames(pd)<-pd$Sample_Name
group_list<-as.data.frame(pd$Sample_Group)
rownames(group_list)<-pd$Sample_Name
group_list$Sample_Group<-pd$Sample_Group
names(group_list)<-c("Sample_Group")
save(pd,group_list,file = 'pd_grouplist.Rdata')

#QC#
champ.QC(as.matrix(myLoad$beta),pheno = group_list$Sample_Group,Rplot = T)
myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC")
save(myNorm,file = 'myNorm.Rdata')

#beta value#
beta.m<-myNorm
plot(density(beta.m[,c(1,2,3,4,5,6,7,8)]),
     col="#FD8D62",lwd=2.5,bty="l" ,ylim=c(0,4.5),
     xlab="β",ylab="Density",main="Density Plot")+
  lines(density(beta.m[,-c(1,2,3,4,5,6,7,8)]),col="#66C3A5",lwd=2.5)

#DMPs#
myDMP<-champ.DMP(beta=beta.m,pheno=group_list$Sample_Group,arraytype = "EPIC",adjPVal = 0.05)
DMP.GUI(DMP=myDMP[[1]],beta = beta.m,pheno = group_list$Sample_Group)
logFoldChange=0.2
P=0.05
allDiff<-myDMP[[1]]
diffSig0.2 <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & adj.P.Val< P )), ]
diffUp0.2 <- allDiff[with(allDiff, (logFC>logFoldChange & adj.P.Val < P )), ]
diffDown0.2 <- allDiff[with(allDiff, (logFC<(-logFoldChange) & adj.P.Val < P )), ]
beta.diff0.2<-beta.m[rownames(diffSig0.2),]
beta.UP0.2<-beta.m[rownames(diffUp0.2),]
beta.DOWN0.2<-beta.m[rownames(diffDown0.2),]

#feature#
library(readxl)
table(diffUp0.2$feature)
table(diffDown0.2$feature)
Feature<- read_excel("Hyper-Hypo-Feature.xlsx")
ggplot(Feature,aes(x=group,y=num,fill = feature))+
  geom_bar(stat = 'identity', position = 'fill')+
  theme_classic()+ 
  scale_fill_brewer(palette = "Set2")

#cgi#
table(diffUp0.2$cgi)
table(diffDown0.2$cgi)
Cgi<- read_excel("Hyper-Hypo-Cgi.xlsx")
ggplot(Cgi,aes(x=group,y=num,fill = cgi))+
  geom_bar(stat = 'identity', position = 'fill')+
  theme_classic()+ 
  scale_fill_brewer(palette = "Set2")

#heatmap#
library(pheatmap)
cg<-beta.diff0.2
ac=data.frame(group=group_list$Sample_Group)
n=t(scale(t(cg)))
n[n>2]=2
n[n< -2]= -2
n[1:4,1:4]
rownames(ac)=colnames(n)
pheatmap(n,color = viridis(100),show_colnames =T,
         show_rownames = F,annotation_col = ac)




###################MethylCIBERSORT##########################
#install.packages(c("caret", "glmnet", "NMF"))
library("caret", "glmnet", "NMF")
library("MethylCIBERSORT")
library(tibble)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
Mat <- beta.m
data("StromalMatrix_V2")
dim(Stromal_v2)
Stromal_v2.pheno
Int <- intersect(rownames(Mat), rownames(Stromal_v2))
Mat <- Mat[match(Int, rownames(Mat)),]
Stromal_v2 <- Stromal_v2[match(Int, rownames(Stromal_v2)),]

RefData <- Stromal_v2
RefPheno <- Stromal_v2.pheno
Signature <- FeatureSelect.V4(CellLines.matrix = NULL,
                              Heatmap = F,
                              export = TRUE,
                              sigName = "MyReference",
                              Stroma.matrix = RefData,
                              deltaBeta = 0.2,
                              FDR = 0.01,
                              MaxDMRs = 100,
                              Phenotype.stroma = RefPheno)
Prep.CancerType(Beta = Mat, Probes = rownames(Signature$SignatureMatrix), fname = "MixtureMatrix")
a<-Signature$SignatureMatrix

source("./CIBERSORT.R")
results <-CIBERSORT(sig_matrix = "./MyReference_Signature.txt",
                    ## the signature from our study
                    mixture_file = "MixtureMatrix.txt",
                    ## i.e. make your own beta matrix
                    perm = 1000,
                    QN = F)

re <- results[,-(11:13)]


k <- apply(re,2,function(x) {sum(x == 0) < nrow(results)/2})
table(k)
re2 <- as.data.frame(t(re[,k]))
an = data.frame(group = pd$Sample_Group,
                row.names = pd$Sample_Name)
pheatmap(re2,scale = "row",
         show_colnames = T,
         annotation_col = an,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

mypalette <- colorRampPalette(brewer.pal(8,"Set2"))

dat <- re %>% as.data.frame() %>%
  rownames_to_column("Sample") %>%
  gather(key = Cell_type,value = Proportion,-Sample)

ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) +
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "bottom") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22))

ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) +
  geom_boxplot(outlier.shape = 21,color = "black") +
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = mypalette(22))

a = dat %>%
  group_by(Cell_type) %>%
  summarise(m = median(Proportion)) %>%
  arrange(desc(m)) %>%
  pull(Cell_type)

dat$Cell_type = factor(dat$Cell_type,levels = a)

ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) +
  geom_boxplot(outlier.shape = 21,color = "black") +
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = mypalette(22))

dat$Group = rep(c(an$group),10)

ggplot(dat,aes(Cell_type,Proportion,fill = Group)) +
  geom_boxplot(outlier.shape = 21,color = "black") +
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(4,1)])+ 
  stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")



#############	Single cell transcriptome ########################
library(scDblFinder)
library(Seurat)
library(dplyr)
library(cowplot)
library(sctransform)
library(SoupX)
library(glmGamPoi)
library(CellChat)
library(tidyverse)
library(ggalluvial)
library(SeuratObject)
library(showtext)
setwd("./GSE173205")

#data load#
M1 = Read10X("./GSE173205/GSE173205_RAW/P112/")
M1 <- CreateSeuratObject(counts = M1, min.cells = 10,min.features = 200,project="M1")
set.seed(1)
dbl<-scDblFinder(M1@assays$RNA@counts)
table(dbl$scDblFinder.class) 
M1$doublet_logic<-ifelse(dbl$scDblFinder.class=="doublet",T,F)
M1<-subset(M1,subset=doublet_logic==F)
M1[["percent.mt"]] <- PercentageFeatureSet(M1, pattern = "^MT-")
M1[["percent.ribo"]] <- PercentageFeatureSet(M1, pattern = "^RP[SL]")
VlnPlot(M1, features = c("nFeature_RNA","percent.mt","percent.ribo","nCount_RNA"), pt.size = 0.1) + NoLegend()
M1<- RenameCells(M1, add.cell.id = "M1")%>%subset(subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25 & nCount_RNA > 1000)
M1<-SCTransform(M1,vst.flavor = "v2", verbose = FALSE)

M2 = Read10X("./GSE173205/GSE173205_RAW/P115/")
M2 <- CreateSeuratObject(counts = M2, min.cells = 10,min.features = 200,project="M2")
set.seed(2)
dbl<-scDblFinder(M2@assays$RNA@counts)
table(dbl$scDblFinder.class) 
M2$doublet_logic<-ifelse(dbl$scDblFinder.class=="doublet",T,F)
M2<-subset(M2,subset=doublet_logic==F)
M2[["percent.mt"]] <- PercentageFeatureSet(M2, pattern = "^MT-")
M2[["percent.ribo"]] <- PercentageFeatureSet(M2, pattern = "^RP[SL]")
VlnPlot(M2, features = c("nFeature_RNA","percent.mt","percent.ribo","nCount_RNA"), pt.size = 0.1) + NoLegend()
M2<- RenameCells(M2, add.cell.id = "M2")%>%subset(subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25 & nCount_RNA > 1000)
M2<-SCTransform(M2,vst.flavor = "v2", verbose = FALSE)

M3 = Read10X("./GSE173205/GSE173205_RAW/P116/")
M3 <- CreateSeuratObject(counts = M3, min.cells = 10,min.features = 200,project="M3")
set.seed(3)
dbl<-scDblFinder(M3@assays$RNA@counts)
M3$doublet_logic<-ifelse(dbl$scDblFinder.class=="doublet",T,F)
M3<-subset(M3,subset=doublet_logic==F)
M3[["percent.mt"]] <- PercentageFeatureSet(M3, pattern = "^MT-")
M3[["percent.ribo"]] <- PercentageFeatureSet(M3, pattern = "^RP[SL]")
VlnPlot(M3, features = c("nFeature_RNA","percent.mt","percent.ribo","nCount_RNA"), pt.size = 0.1) + NoLegend()
M3<- RenameCells(M3, add.cell.id = "M3")%>%subset(subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25 & nCount_RNA > 1000)
M3<-SCTransform(M3,vst.flavor = "v2", verbose = FALSE)

M4 = Read10X("./GSE173205/GSE173205_RAW/P121/")
M4 <- CreateSeuratObject(counts = M4, min.cells = 10,min.features = 200,project="M4")
set.seed(4)
dbl<-scDblFinder(M4@assays$RNA@counts)
M4$doublet_logic<-ifelse(dbl$scDblFinder.class=="doublet",T,F)
M4<-subset(M4,subset=doublet_logic==F)
M4[["percent.mt"]] <- PercentageFeatureSet(M4, pattern = "^MT-")
M4[["percent.ribo"]] <- PercentageFeatureSet(M4, pattern = "^RP[SL]")
VlnPlot(M4, features = c("nFeature_RNA","percent.mt","percent.ribo","nCount_RNA"), pt.size = 0.1) + NoLegend()
M4<- RenameCells(M4, add.cell.id = "M4")%>%subset(subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25 & nCount_RNA > 1000)
M4<-SCTransform(M4,vst.flavor = "v2", verbose = FALSE)

M1$gender<-"F"
M1$year<-"51"
M2$gender<-"M"
M2$year<-"48"
M3$gender<-"F"
M3$year<-"57"
M4$gender<-"F"
M4$year<-"44"

#datasets integration#
features <- SelectIntegrationFeatures(object.list = list(M1,M2,M3,M4), nfeatures = 3000)
alldata.list <- PrepSCTIntegration(object.list = list(M1,M2,M3,M4), anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = alldata.list, normalization.method = "SCT",anchor.features = features)
alldata <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

library(ggplot2)
library(Seurat)
library(dplyr)
library(cowplot)
library(sctransform)
alldata <- CellCycleScoring(object = alldata,g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
VlnPlot(alldata, features = c("S.Score","G2M.Score"),pt.size = 0)
alldata <- RunPCA(alldata, npcs = 30, verbose = FALSE)
alldata <- RunUMAP(alldata, eduction = "pca", dims = 1:30, verbose = FALSE)
alldata <- FindNeighbors(alldata, reduction = "pca", dims = 1:30, verbose = FALSE)
alldata <- FindClusters(alldata, resolution = 1,verbose=FALSE)
alldata <- RunTSNE(alldata, reduction = "pca", dims = 1:30, verbose = FALSE)


p1 <- DimPlot(alldata, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(alldata, reduction = "umap", label = TRUE)
p1 | p2
DimPlot(alldata, reduction = "umap", split.by = "orig.ident",label = TRUE)
DimPlot(alldata, reduction = "umap", split.by = "gender",label = TRUE)
DimPlot(alldata, reduction = "umap", split.by = "year",label = TRUE)

p1 <- DimPlot(alldata, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(alldata, reduction = "tsne", label = TRUE)
p1 | p2
DimPlot(alldata, reduction = "tsne", split.by = "orig.ident",label = TRUE)
DimPlot(alldata, reduction = "tsne", split.by = "gender",label = TRUE)
DimPlot(alldata, reduction = "tsne", split.by = "year",label = TRUE)

DefaultAssay(alldata) <- "RNA"   
alldata <- NormalizeData(alldata, verbose=F)
alldata <- ScaleData(alldata, verbose=F)
save(alldata,file = "GSE173205_alldata.Rdata")
load(file = "GSE173205_alldata.Rdata")

alldata.markers <- FindAllMarkers(alldata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
markers_seurat <-alldata.markers%>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(alldata, features = markers$gene) + NoLegend()


unique_clusters <- unique(alldata$seurat_clusters)
for (cluster_label in unique_clusters) {
  markers <- subset(markers_seurat,cluster==cluster_label)
  vln_plot <- VlnPlot(alldata, features = markers$gene, pt.size = 0)
  pdf_file <- file.path("marker_plots", paste0("cluster_", cluster_label, "_markers.pdf"))
  pdf(pdf_file,height = 20,width = 20)
  print(vln_plot)
  dev.off()
}

# Cell Type#
new.cluster.ids <- c("Endothelial cells","Fibroblasts","Endothelial cells","Melanocyte","Endothelial cells"
                     ,"T_cells","Fibroblasts","Smooth muscle cells","T_cells","Fibroblasts","Keratinocytes"
                     ,"Endothelial cells","Fibroblasts","Keratinocytes","Myeloid cells","T_cells"
                     ,"Fibroblasts","Myeloid cells","Keratinocytes","Smooth muscle cells","Keratinocytes"
                     ,"Fibroblasts","Fibroblasts","Fibroblasts","Keratinocytes","Myeloid cells"
                     ,"Endothelial cells","Myeloid cells","Keratinocytes")
names(new.cluster.ids) <- levels(alldata)              
new.cluster.ids
alldata <- RenameIdents(alldata, new.cluster.ids)
DimPlot(alldata, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(alldata, reduction = "tsne", label = TRUE, pt.size = 0.5) 

#Cell select markers
{
  FeaturePlot(alldata,features = c("KRT5","KRT14","KRT15"))
  FeaturePlot(alldata,features = c("VWF","CCL14",'CAVIN2'))
  FeaturePlot(alldata,features = c("COL1A1","COL3A1","LUM"))
  FeaturePlot(alldata,features = c("CD3E","CD69","IL7R"))
  FeaturePlot(alldata,features = c("DCT","MLANA","TYRP1"))
  FeaturePlot(alldata,features = c("RERGL","MYH11","SORBS2"))
  FeaturePlot(alldata,features = c("TYROBP",'AIF1',"CD68"))
  
  features <- c("VWF","CCL14",'CAVIN2',
                "COL1A1","COL3A1","LUM",
                "DCT","MLANA","TYRP1",
                "CD3E","CD69","IL7R",
                "RERGL","MYH11","SORBS2",
                "KRT5","KRT14","KRT15",
                "TYROBP",'AIF1',"CD68"
  )           
  DotPlot(alldata, features = features) + RotatedAxis()+scale_color_gradient(high = "#E64B35FF", low = "#86CAB9")
  DoHeatmap(subset(alldata, downsample = 100), features = features, size = 3)
}

save(alldata,file="GSE173205_alldata_cellsort.Rdata")


#cell Composition#
{
  table(alldata$orig.ident)
  Cellratio <- prop.table(table(Idents(alldata)))
  table(Idents(alldata), alldata$orig.ident)
  Cellratio <- prop.table(table(Idents(alldata), alldata$orig.ident), margin = 2)
  Cellratio
  Cellratio <- as.data.frame(Cellratio)
  colourCount = length(unique(Cellratio$Var1))
  library(ggplot2)
  ggplot(Cellratio) + 
    geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
    theme_classic() +
    labs(x='Sample',y = 'Ratio')+
    coord_flip()+
    theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
  
  ggplot(Cellratio) + 
    geom_bar(aes(x =Var1, y= Freq,fill=Var1),stat = "identity",width = 0.7,size = 0.5)+ 
    theme_classic() +
    labs(x='Sample',y = 'Ratio')+
    coord_flip()+
    theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
  
}

#Two Group Cell_Composition_comparison
{
  load("/project/hlFan/Meningioma/Cell_Composition_comparison.Rdata")
  p = ggbarplot(a,x = "Var1",y="Freq",add = "mean_se",position = position_dodge(0.8),fill="group",color = "group",palette="Reds", add.params = list(color = "black"))+
    theme_classic()+
    stat_compare_means(aes(group=group), label = "p.signif",method = "t.test")
  
  p
  
}



######################EpiSCORE############################ 
library(EpiSCORE)
library(ggsci)
library(tidyr)
library(ggpubr)
library(RColorBrewer)
library(tibble)
library(pheatmap)
library(tidyverse)
library(EpiDISH)
library(SuperExactTest)
expref.o <- ConstExpRef(lungSS2mca1.m,celltypeSS2.idx,celltypeSS2.v,markspecTH=rep(3,4));
print(dim(expref.o$ref$med));
head(expref.o$ref$med)

celltypeSS2.idx<-alldata@meta.data$seurat_clusters
celltypeSS2.idx<-as.numeric(celltypeSS2.idx)
celltypeSS2.idx[which(celltypeSS2.idx %in% c(1,3,5,12,27))]<-"31" ##Endothelial cells
celltypeSS2.idx[which(celltypeSS2.idx %in% c(2,7,10,13,17,22,23,24))]<-"32" ##Fibroblasts
celltypeSS2.idx[which(celltypeSS2.idx =="4")]<-"33" ##Melanocyte
celltypeSS2.idx[which(celltypeSS2.idx %in% c(6,9,16))]<-"34" ##T_cells
celltypeSS2.idx[which(celltypeSS2.idx %in% c(8,20))]<-"35" ##Smooth muscle cells
celltypeSS2.idx[which(celltypeSS2.idx %in% c(11,14,19,21,25,29))]<-"36" ##Keratinocytes
celltypeSS2.idx[which(celltypeSS2.idx %in% c(15,18,26,28))]<-"37" ##Myeloid cells

celltypeSS2.idx[which(celltypeSS2.idx =="31")]<-"1" 
celltypeSS2.idx[which(celltypeSS2.idx =="32")]<-"2" 
celltypeSS2.idx[which(celltypeSS2.idx =="33")]<-"3" 
celltypeSS2.idx[which(celltypeSS2.idx =="34")]<-"4" 
celltypeSS2.idx[which(celltypeSS2.idx =="35")]<-"5" 
celltypeSS2.idx[which(celltypeSS2.idx =="36")]<-"6"
celltypeSS2.idx[which(celltypeSS2.idx =="37")]<-"7" 
celltypeSS2.idx <- factor(celltypeSS2.idx, levels=c('1','2',"3","4","5","6","7"))
table(celltypeSS2.idx)
levels(celltypeSS2.idx)
count<- as(alldata@assays$RNA@counts, "matrix")

save(celltypeSS2.idx,file="celltypeSS2.idx.Rdata")
load(file="celltypeSS2.idx.Rdata")
save(count,file="count.Rdata")
load(count,file="count.Rdata")


celltypeSS2.v<-c("Endothelial_cells","Fibroblasts","Melanocyte",
                 "T_cells","Smooth_muscle_cells","Keratinocytes",
                 "Myeloid_cells")
ncpct.v <- summary(factor(celltypeSS2.idx));
names(ncpct.v) <- celltypeSS2.v;
print(ncpct.v);


out.l <- ConstExpRef(count,celltypeSS2.idx,celltypeSS2.v)
class(celltypeSS2.idx)
class(celltypeSS2.v)
celltypeSS2.idx<-as.numeric(celltypeSS2.idx)
#save(out.l,file="EpiSCORE_out.l.Rdata")
refDNAm1.m <- ImputeDNAmRef(out.l$ref$med,db="SCM2",geneID="SYMBOL");
refDNAm2.m <- ImputeDNAmRef(out.l$ref$med,db="RMAP",geneID="SYMBOL");
refDNAm.m <- ConstMergedDNAmRef(refDNAm1.m,refDNAm2.m)

load("./myNorm.Rdata")
avDNAm.m <- constAvBetaTSS(myNorm,type="850k")
wRPC.o <- wRPC(avDNAm.m,refDNAm.m,useW=TRUE,wth=0.4,maxit=2000)
estF<-wRPC.o[["estF"]]
boxplot(estF)
heatmap(estF)

load("./pd_grouplist.Rdata")
an = data.frame(group = pd$Sample_Group,
                row.names = pd$Sample_Name)
pheatmap(t(estF),show_colnames =T,show_rownames = T,annotation_col = an,
         color =  rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)))


data(SkinRef)
estF.m <- epidish(exprefSkin.m,ref.m=out.l$ref$med,method="RPC",maxit=200)$estF
pheatmap(t(estF.m),show_colnames =T,show_rownames = T)


a <- estF
a<-as.data.frame(a)
rownames(pd)<-pd$Sample_Name
identical(rownames(a),rownames(pd))
b <- pd
class(b$Sample_Group)
a$group <- b$Sample_Group
a <- a %>% rownames_to_column("Sample_Name")

b <- gather(a,key=cell,value = Proportion,-c(group,Sample_Name))
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
ggboxplot(b, x = "cell", y = "Proportion",
          fill = "group",
          title = "Cell_Proportion")+
  scale_fill_manual(values = mypalette(22)[c(4,1)])+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=60, hjust=1),plot.title = element_text(hjust = 0.5,size=15)) 

c = b %>%
  group_by(cell) %>%
  summarise(m = median(Proportion)) %>%
  arrange(desc(m)) %>%
  pull(cell)
b$cell = factor(b$cell,levels = c)

ggplot(b,aes(cell,Proportion,fill = cell)) +
  geom_boxplot(outlier.shape = 21,color = "black") +
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")+
  scale_fill_manual(values = mypalette(22))

ggplot(b,aes(Sample_Name,Proportion,fill = cell)) +
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",y = "Estiamted Proportion") +
  #theme_bw() +
  theme(axis.text.x = element_text(angle=60, hjust=1),
        legend.position = "bottom") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22))


pheno.v <- ifelse(an$group == "cluster_1", 0, 1)
pheno.v <- as.numeric(pheno.v)
names(pheno.v) <- rownames(pd)
frac.m<-estF
save(beta,)
cdmc.o <- CellDMC(myNorm,pheno.v=pheno.v, frac.m=frac.m, 
                  adjPMethod = "fdr", adjPThresh = 0.05, 
                  cov.mod = NULL, sort = FALSE)
load(CellDMC.Rdata)

dmctLUSC.lv <- list();
for(ct in 1:6){
  dmctLUSC.lv[[ct]] <- rownames(cdmc.o$coe[[ct]][which(cdmc.o$dmct[,1+ct]!=0),]);
}
res.o <- supertest(dmctLUSC.lv,n=720903)

plot(res.o,"landscape",sort.by="size")
plot(res.o,"landscape")
plot(res.o)



