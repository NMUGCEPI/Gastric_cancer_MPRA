##########################################################################################################################
#--------------------------------------------------------------------------------------------------- crop seq pre-process
##########################################################################################################################

#1.transcript analysis
cd /Public/wtp/blj/singlecell/crop_seq/code
multi_rna\
 --mapfile /Public/wtp/blj/singlecell/crop_seq/code/mapfile.txt\
 --chemistry auto\
 --genomeDir /Public/wtp/blj/singlecell/CeleScope/ensembl_75\
 --thread 8\
 --mod shell 

cd /Public/wtp/blj/singlecell/crop_seq/work
celescope rna sample --outdir .//crop_hgc/00.sample --sample crop_hgc --thread 8 --chemistry auto --wells 384  --fq1 /Public/wtp/blj/singlecell/crop_seq/rawdata/LHC241024041/crop_hgc_L2_1.fq.gz 
celescope rna starsolo --outdir .//crop_hgc/01.starsolo --sample crop_hgc --thread 8 --chemistry auto --adapter_3p AAAAAAAAAAAA --genomeDir /Public/wtp/blj/singlecell/CeleScope/ensembl_75 --outFilterMatchNmin 50 --soloCellFilter "EmptyDrops_CR 30000 0.99 10 45000 90000 500 0.01 20000 0.001 10000" --starMem 32 --soloFeatures "Gene GeneFull_Ex50pAS"  --fq1 /Public/wtp/blj/singlecell/crop_seq/rawdata/LHC241024041/crop_hgc_L2_1.fq.gz --fq2 /Public/wtp/blj/singlecell/crop_seq/rawdata/LHC241024041/crop_hgc_L2_2.fq.gz 
celescope rna analysis --outdir .//crop_hgc/02.analysis --sample crop_hgc --thread 8 --genomeDir /Public/wtp/blj/singlecell/CeleScope/ensembl_75  --matrix_file .//crop_hgc/outs/filtered

#2.CROPseq splite sgRNA
multi_tag \
--mapfile /Public/wtp/blj/singlecell/crop_seq/code/mapfile_tag.txt \
--linker_fasta /Public/wtp/blj/singlecell/crop_seq/code/linker.fasta \
--barcode_fasta /Public/wtp/blj/singlecell/crop_seq/code/sgRNA_sequences.fasta \
--dim 1 \
--fq_pattern L26C20 \
--mod shell \
--split_matrix

cd /Public/wtp/blj/singlecell/crop_seq/tag
celescope tag sample --outdir .//crop_ges/00.sample --sample crop_ges --thread 10 --chemistry auto --wells 384  --fq1 /Public/wtp/blj/singlecell/crop_seq/rawdata/20241119/CUSLHC241024041/CUSLHC241024041_S2_L001_R1_001.fastq.gz 
celescope tag barcode --outdir .//crop_ges/01.barcode --sample crop_ges --thread 10 --chemistry auto --lowNum 2 --wells 384  --fq1 /Public/wtp/blj/singlecell/crop_seq/rawdata/20241119/CUSLHC241024041/CUSLHC241024041_S2_L001_R1_001.fastq.gz --fq2 /Public/wtp/blj/singlecell/crop_seq/rawdata/20241119/CUSLHC241024041/CUSLHC241024041_S2_L001_R2_001.fastq.gz 
celescope tag mapping_tag --outdir .//crop_ges/02.mapping_tag --sample crop_ges --thread 10 --fq_pattern L26C20 --barcode_fasta /Public/wtp/blj/singlecell/crop_seq/code/sgRNA_sequences.fasta --linker_fasta /Public/wtp/blj/singlecell/crop_seq/code/linker.fasta  --fq .//crop_ges/01.barcode/crop_ges_2.fq 
celescope tag count_tag --outdir .//crop_ges/03.count_tag --sample crop_ges --thread 10 --UMI_min auto --dim 1 --SNR_min auto --coefficient 0.1  --match_dir /Public/wtp/blj/singlecell/crop_seq/work/crop_ges --read_count_file .//crop_ges/02.mapping_tag/crop_ges_read_count.tsv 
celescope tag analysis_tag --outdir .//crop_ges/04.analysis_tag --sample crop_ges --thread 10  --match_dir /Public/wtp/blj/singlecell/crop_seq/work/crop_ges --tsne_tag_file .//crop_ges/03.count_tag/crop_ges_tsne_tag.tsv 
celescope tag split_tag --outdir .//crop_ges/05.split_tag --sample crop_ges --thread 10 --split_matrix  --match_dir /Public/wtp/blj/singlecell/crop_seq/work/crop_ges --umi_tag_file .//crop_ges/03.count_tag/crop_ges_umi_tag.tsv 

#3.process expression matrix
cd /Public/wtp/blj/singlecell/crop_seq/work
library(data.table)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
#QC
ges <- Read10X(data.dir = '/Public/wtp/blj/singlecell/crop_seq/work/crop_ges/outs/filtered')
ges <- CreateSeuratObject(counts = ges ,project = 'seurat', min.cells = 3, min.features = 200,names.delim = '_')
ges[["percent.mt"]]<-PercentageFeatureSet(object = ges, pattern = "^MT-")
ges@meta.data$orig.ident <- 'ges'
tag_ges <- fread('/Public/wtp/blj/singlecell/crop_seq/tag/crop_ges/03.count_tag/crop_ges_umi_tag.tsv')
a <- rownames(ges@meta.data)==tag_ges$V1
ges@meta.data$tag <- tag_ges$tag
save(ges,file='/Public/wtp/blj/singlecell/crop_seq/work/crop_ges/ges.Rdata')
hgc <- Read10X(data.dir = '/Public/wtp/blj/singlecell/crop_seq/work/crop_hgc/outs/filtered')
hgc <- CreateSeuratObject(counts = hgc ,project = 'seurat', min.cells = 3, min.features = 200,names.delim = '_')
hgc[["percent.mt"]]<-PercentageFeatureSet(object = hgc, pattern = "^MT-")
hgc@meta.data$orig.ident <- 'hgc'
tag_hgc <- fread('/Public/wtp/blj/singlecell/crop_seq/tag/crop_hgc/03.count_tag/crop_hgc_umi_tag.tsv')
a <- rownames(hgc@meta.data)==tag_hgc$V1
hgc@meta.data$tag <- tag_hgc$tag
save(hgc,file='/Public/wtp/blj/singlecell/crop_seq/work/crop_hgc/hgc.Rdata')
setwd('/Public/wtp/blj/singlecell/crop_seq/result_v2')
ALL <- merge(x=ges,y=hgc,merge.data=TRUE,project="SeuratProject")
mean(ALL$'percent.mt')
p<-VlnPlot(object=ALL,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3,pt.size=0,group.by='orig.ident')
ggsave('featureViolin.png', p, width = 7, height = 6)
ALL_SIN1<-subset(x=ALL,subset=nFeature_RNA>200&nFeature_RNA<5000)
ALL_SIN2<-subset(x=ALL_SIN1,subset=nCount_RNA>1000&nCount_RNA<20000) 
ALL_SIN<-subset(x=ALL_SIN2,subset=percent.mt<20)

##romove cells have multiple sgRNAs
Idents(ALL_SIN) <- "tag"
ALL_SIN_SG <- subset(ALL_SIN,idents=c('Undetermined','Multiplet'),invert = TRUE)
##cell PC
ALL_SIN_SG_nor <- NormalizeData(object = ALL_SIN_SG, normalization.method = 'LogNormalize',scale.factor = 10000)
ALL_SIN_SG_nor <- FindVariableFeatures(object = ALL_SIN_SG_nor,selection.method='vst',mean.function = ExpMean,dispersion.function = LogVMR,mean.cutoff=c(0.125,3),dispersion.cutoff=c(0.5,Inf))
top10<-head(x = VariableFeatures(object=ALL_SIN_SG_nor),10)
plot1<-VariableFeaturePlot(object = ALL_SIN_SG_nor)
plot2<-LabelPoints(plot = plot1,points = top10,repel =T)
p<-CombinePlots(plots=list(plot1,plot2))
ggsave('featureVar_MT50.png', p, width = 10, height = 6)
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
ALL_SIN_SG_nor <- CellCycleScoring(ALL_SIN_SG_nor, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ALL_SIN_SG_nor_PCA=ScaleData(ALL_SIN_SG_nor,vars.to.regress = c("percent.mt","nCount_RNA","orig.ident","S.Score", "G2M.Score"))
ALL_SIN_SG_nor_PCA=RunPCA(object=ALL_SIN_SG_nor_PCA,npcs=100,pc.genes=VariableFeatures(object=ALL_SIN_SG_nor_PCA))
p<-ElbowPlot(object=ALL_SIN_SG_nor_PCA,ndims=100)
ggsave('PCA_MT50_100.png', p, width = 8, height = 6)
p <- DimPlot(object = ALL_SIN_SG_nor_PCA, reduction = "pca", pt.size = 0.1, group.by = "orig.ident",raster=FALSE)
ggsave('epi_pca_batch.png', p, width = 6, height = 6)
pcSelect=20
ALL_SIN_SG_nor_PCA<-FindNeighbors(object=ALL_SIN_SG_nor_PCA,reduction = "pca",dims=1:pcSelect)
ALL_SIN_SG_nor_PCA<-FindClusters(object=ALL_SIN_SG_nor_PCA,resolution=c(0.1,0.2,0.3))
ALL_SIN_SG_nor_PCA <- RunUMAP(ALL_SIN_SG_nor_PCA, reduction = "pca", dims = 1:pcSelect,min.dist = 0.5, n.neighbors = 25L)
Idents(ALL_SIN_SG_nor_PCA) <- 'RNA_snn_res.0.1'
ALL_SIN_SG_nor_PCA$seurat_clusters <- ALL_SIN_SG_nor_PCA$RNA_snn_res.0.1
p<-DimPlot(ALL_SIN_SG_nor_PCA, reduction = "umap",label.size=4,pt.size=0.1,label=T,raster=FALSE)
ggsave('ALL_SIN_SG_nor_PCA_UMAP_cluster.png', p, width = 6.5, height = 5.5)
p<-DimPlot(ALL_SIN_SG_nor_PCA, reduction = "umap",pt.size=0.1,group.by='RNA_snn_res.0.2',raster=FALSE)
ggsave('ALL_SIN_SG_nor_PCA_UMAP_cluster2.png', p, width = 6.5, height = 5.5)
p<-DimPlot(ALL_SIN_SG_nor_PCA, reduction = "umap",pt.size=0.1,group.by='RNA_snn_res.0.3',raster=FALSE)
ggsave('ALL_SIN_SG_nor_PCA_UMAP_cluster3.png', p, width = 6.5, height = 5.5)
p<-DimPlot(ALL_SIN_SG_nor_PCA, reduction = "umap",pt.size=0.1,group.by='orig.ident',raster=FALSE)
ggsave('ALL_SIN_SG_nor_PCA_UMAP_cellline.png', p, width = 6.5, height = 5.5)
save(ALL_SIN_SG_nor_PCA,file='ALL_SIN_SG_nor_PCA.Rdata')
#remove outline cells
library(factoextra)
umap_type = ALL_SIN_SG_nor_PCA@reductions$umap@cell.embeddings %>%
      as.data.frame() %>% cbind(type = ALL_SIN_SG_nor_PCA@meta.data$orig.ident)
km <- kmeans(umap_type[,1:2], 2)
p<-fviz_cluster(km, data = umap_type[,1:2],
			 geom = c("point"),
			 ellipse=F,
             ggtheme = theme_bw())
ggsave('UMAP_KM2.png', p, width = 6.5, height = 6)
umap_type$km<-km$cluster
ges<-subset(umap_type,type=='ges'&km %in% c('2'))
hgc<-subset(umap_type,type=='hgc'&km %in% c('1'))
all<-rbind(ges,hgc)
ALL_SIN_SG_nor_filter<-subset(x=ALL_SIN_SG_nor_PCA,cell=rownames(all))
save(ALL_SIN_SG_nor_filter,file='ALL_SIN_SG_nor_filter.Rdata')
p<-DimPlot(ALL_SIN_SG_nor_filter, reduction = "umap",pt.size=0.1,group.by='orig.ident',raster=FALSE)
ggsave('ALL_SIN_SG_nor_PCA_UMAP_cellline_filter.png', p, width = 6.5, height = 5.5)

##Figure 5D
use_colors <- c(ges="#65c294",hgc="#f3715c")
p<-DimPlot(ALL_SIN_SG_nor_filter, reduction = "umap",pt.size=0.1,group.by='orig.ident',raster=FALSE)+
   scale_color_manual(values = use_colors)
ggsave('ALL_SIN_SG_nor_PCA_UMAP_cellline_filter.pdf', p, width = 6.5, height = 5.5)

#Figure_S13C
R1 <- fread('/data/blj/crop_seq/results_v2/outY1.txt')
R2 <- fread('/data/blj/crop_seq/results_v2/outY2.txt')
R3 <- fread('/data/blj/crop_seq/results_v2/outY3.txt')
R4 <- cbind(R1,R2,R3)
colnames(R4)[c(2,4,6)] <- c('f1','f2','f3')
R4$number <- R4$f1+R4$f2+R4$f3
R4$number <- log10(R4$number)
R4$group <- 'aa'
data <- data.frame(Value = R4$number)
data$group <- 1
p=ggboxplot(data,x='group',y='Value',fill = "grey")+
   labs(x = '', y = 'log10 read count per sgRNA')+guides(fill=F)
ggsave('/data/blj/crop_seq/results_v2/log10 read count per sgRNA.pdf', p, width = 3, height =5)

#Figure_S13G
meta_sg<-ALL_SIN_SG_nor_filter@meta.data
a<-as.data.frame(table(meta_sg$orig.ident,meta_sg$tag))
write.csv(a,"sample_number_sg.csv",quote=F)
a$group <- 'sg'
ges_sg <- subset(a,Var1=='ges')
hgc_sg <- subset(a,Var1=='hgc')
summary(ges$Freq)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   81.0   268.2   319.5   317.2   376.5   624.0
summary(hgc$Freq)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   73.0   212.8   254.0   258.7   303.0   480.0 
p<-ggplot(data=ges,aes(x=Freq)) + 
   geom_histogram(binwidth=30,fill="#94aad6",color="#e9ecef", alpha=0.9)+
   labs(x="Cells per sgRNA in GES-1",y="sgRNA")+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 14))
ggsave('cellnumber_sgRNA_ges.pdf', p, width = 5, height =4)

p<-ggplot(data=hgc,aes(x=Freq)) + 
   geom_histogram(binwidth=30,fill="#94aad6",color="#e9ecef", alpha=0.9)+
   labs(x="Cells per sgRNA in HGC-27",y="sgRNA")+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 14))
ggsave('cellnumber_sgRNA_hgc.pdf', p, width = 5, height =4)

##Figure_S13H
library(tidyr)
library(ggpubr)
dat <- read.csv('/Public/wtp/blj/singlecell/crop_seq/result_v2/sample_number_sg.csv',row.names=1)
data_wide <- pivot_wider(dat, 
                         names_from = Var1, 
                         values_from = Freq)
p<-ggplot(data_wide, aes(x = ges, y = hgc)) + 
  geom_point() +
  labs(x = 'Cells recovered by knockdown in GES-1', y = 'Cells recovered by knockdown in HGC-27') +
  theme(axis.title = element_text(size = 14))+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'))+
  stat_cor(data=data_wide, method = "pearson",size = 6)
ggsave('cellnumber_sgRNA_pearson.pdf', p, width = 5, height =5)
result <- cor.test(data_wide$ges, data_wide$hgc, method = "pearson")

###Figure 5E
a<-as.data.frame(table(meta_sg$orig.ident,meta_sg$group))
a$group <- 'sgregion'
ges_region <- subset(a,Var1=='ges')
hgc_region <- subset(a,Var1=='hgc')
summary(ges$Freq)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  179.0   815.5   965.0   935.0  1066.2  2377.0
summary(hgc$Freq)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  180.0   657.0   775.0   762.6   871.8  1982.0
ges <- rbind(ges_sg,ges_region)
hgc <- rbind(hgc_sg,hgc_region)
all <- rbind(ges,hgc)
p=ggboxplot(all,'group', 'Freq',fill = "group",palette = c("#C2C2C2","#707070"),add = "none")+
  facet_wrap(~Var1, nrow = 1)+
  labs(x="",y="Number of cells")+
ggsave('cellnumber_sg_and_target_boxplot.pdf', p, width = 6, height =4.5) 

##########################################################################################################################
#--------------------------------------------------------------------------------------------calculate sgRNA target genes
##########################################################################################################################

cd /Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts
library(data.table)
library(dplyr)
library(Seurat)
library(edgeR)
library(MatchIt)
load('/Public/wtp/blj/singlecell/crop_seq/result_v2/ALL_SIN_SG_nor_filter.Rdata')
meta <- ALL_SIN_SG_nor_filter@meta.data
meta_ges <- subset(meta,orig.ident=='ges')
meta_hgc <- subset(meta,orig.ident=='hgc')

count_data<-GetAssayData(ALL_SIN_SG_nor_filter,slot="counts")
count_data<-as.data.frame(count_data)
count_data<-as.data.frame(t(count_data))
count_data$id<-rownames(count_data)
meta <- ALL_SIN_SG_nor_filter@meta.data
meta$id<-rownames(meta)
meta <- meta[,c(1,14,15,17)]
count_data<-merge(count_data,meta,by='id')
count_data$ID <- paste(count_data$orig.ident,count_data$group,sep='_')
N=ncol(count_data)-4
exp_sum<-aggregate(count_data[,2:N], by=list(group=count_data$ID),sum)
rownames(exp_sum)<-exp_sum$group
exp_sum<-as.data.frame(t(exp_sum[,2:ncol(exp_sum)]))
fwrite(exp_sum, 'crop_all_sum_exp_v2.xls', row.names=T, col.names=T, sep='\t')

## neg control
#######ges
dat <- table(meta_ges$sgRNA)
dat <- as.data.table(dat)
dat[, "sgRNA" := tstrsplit(V1, "_")[1]]
dat1 <- dat %>%
  group_by(sgRNA) %>%
  summarise(N_sum = sum(N))
dat_list <- list()
dat_list$sgRNA_id <- dat1$sgRNA
dat_list$number <- round(dat1$N_sum)
for (j in 3:38) {
sgRNA <- dat_list$sgRNA_id[j]
print(sgRNA)
subsg <- c(sgRNA, "control")
meta_ges1 <- subset(meta_ges,group %in% subsg)
meta_ges1$info <- ifelse(meta_ges1$group == sgRNA,"case","control")
meta_ges1$info =recode(meta_ges1$info, case='1', control='0')
meta_ges1$info <- as.numeric(meta_ges1$info)
##1:2 match ref
extractSamples <- matchit(
info ~ nCount_RNA + nFeature_RNA + percent.mt +RNA_snn_res.0.2,
data = meta_ges1, 
method = "nearest",
distance = "glm",
link = "logit",
ratio = 2,
caliper = 0.1,
replace=F
)
summary(extractSamples)
cell_Match <- match.data(extractSamples)
cell_Match_ref <- subset(cell_Match,info==0)
cell_id <- rownames(cell_Match_ref)
sc_data_ref <- subset(ALL_SIN_SG_nor_filter,cell=cell_id)
##counts data
count_data<-GetAssayData(sc_data_ref,slot="counts")
count_data<-as.data.frame(count_data)

set.seed(123)
n_samples <- ncol(count_data) 
sample_size <- dat_list$number[j]
n_repeats <- 100 
#mean exp
sample_sums <- list()
#100 sample
     for (i in 1:n_repeats) {
         sample_indices <- sample(n_samples, size = sample_size)
         sampled_matrix <- count_data[, sample_indices]
         sum_expression <- rowSums(sampled_matrix)
         sample_sums[[i]] <- sum_expression
         }
sample_sums_df <- do.call(cbind, sample_sums)
colnames(sample_sums_df) <- paste("Sample", 1:ncol(sample_sums_df), sep="")
sample_sums_df <- as.data.frame(sample_sums_df)
output <- sprintf('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/ref_sgrna_v2/%s_ges_control_exp.xls',sgRNA)
fwrite(sample_sums_df, output, row.names=T, col.names=T, sep='\t')
}

########hgc
dat <- table(meta_hgc$sgRNA)
dat <- as.data.table(dat)
dat[, "sgRNA" := tstrsplit(V1, "_")[1]]
dat1 <- dat %>%
  group_by(sgRNA) %>%
  summarise(N_sum = sum(N))
dat_list <- list()
dat_list$sgRNA_id <- dat1$sgRNA
dat_list$number <- round(dat1$N_sum)
for (j in 3:38) {
sgRNA <- dat_list$sgRNA_id[j]
print(sgRNA)
subsg <- c(sgRNA, "control")
meta_hgc1 <- subset(meta_hgc,group %in% subsg)
meta_hgc1$info <- ifelse(meta_hgc1$group == sgRNA,"case","control")
meta_hgc1$info =recode(meta_hgc1$info, case='1', control='0')
meta_hgc1$info <- as.numeric(meta_hgc1$info)
##1:2 match ref
extractSamples <- matchit(
info ~ nCount_RNA + nFeature_RNA + percent.mt +RNA_snn_res.0.2,
data = meta_hgc1, 
method = "nearest",
distance = "glm",
link = "logit",
ratio = 2,
caliper = 0.1,
replace=F
)
summary(extractSamples)
cell_Match <- match.data(extractSamples)
cell_Match_ref <- subset(cell_Match,info==0)
cell_id <- rownames(cell_Match_ref)
sc_data_ref <- subset(ALL_SIN_SG_nor_filter,cell=cell_id)
##counts data
count_data<-GetAssayData(sc_data_ref,slot="counts")
count_data<-as.data.frame(count_data)
set.seed(123)
n_samples <- ncol(count_data) 
sample_size <- dat_list$number[j]
n_repeats <- 100 
#mean exp
sample_sums <- list()
#100 sample
     for (i in 1:n_repeats) {
         sample_indices <- sample(n_samples, size = sample_size)
         sampled_matrix <- count_data[, sample_indices]
         sum_expression <- rowSums(sampled_matrix)
         sample_sums[[i]] <- sum_expression
         }
sample_sums_df <- do.call(cbind, sample_sums)
colnames(sample_sums_df) <- paste("Sample", 1:ncol(sample_sums_df), sep="")
sample_sums_df <- as.data.frame(sample_sums_df)
output <- sprintf('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/ref_sgrna_v2/%s_hgc_control_exp.xls',sgRNA)
fwrite(sample_sums_df, output, row.names=T, col.names=T, sep='\t')
}

##calculate DEG
library(limma)
library(edgeR)
library(data.table)
exp_sum <- read.table('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/crop_all_sum_exp_v2.xls',sep='\t',row.names=1,h=T)
exp_sum_ges<-exp_sum[,c(1:38)]
exp_sum_hgc<-exp_sum[,c(39:76)]
#ges 
setwd('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/ges')
meta <- as.data.table(colnames(exp_sum_ges))
colnames(meta) <- 'sample'
meta[, "sgRNA" := tstrsplit(sample, "_")[2]]
sg <- unique(meta$sgRNA)
#DEG vs ref
for(i in 3:38){
sgid<-sg[i]
print(sgid)
meta1 <- subset(meta,sgRNA %in% c("control",sgid))
exp_sum_sub1 <- exp_sum_ges[,meta1$sample]
ref <- read.table(sprintf('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/ref_sgrna_v2/%s_ges_control_exp.xls',sgid),sep='\t',row.names=1,h=T)
meta_ref <- as.data.table(colnames(ref))
colnames(meta_ref) <- 'sample'
meta_ref$sgRNA <- 'ref_sample'
meta2 <- rbind(meta1,meta_ref)
meta2 <- subset(meta2,sgRNA != 'control')
exp_sum_sub2 <- cbind(exp_sum_sub1,ref)
exp_sum_sub2 <- exp_sum_sub2[,meta2$sample]
group <- meta2$sgRNA
y <- DGEList(counts=exp_sum_sub2,group = group)
keep <- filterByExpr(y,min.count = 5, min.total.count = 10)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y,method = 'TMM')
group=factor(y$samples$group)
design <- model.matrix(~group)
rownames(design) <- colnames(y)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmFit(y, design, robust=T)
qlf <- glmLRT(fit)
res = as.data.frame(topTags(qlf,n = nrow(y)))
res = cbind(gene_id=rownames(res),res)
write.table(res,sprintf("%s_vs_ref_sample_ges.xls",sgid),sep="\t",quote=FALSE,row.names=F)
}

#hgc
setwd('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/hgc')
meta <- as.data.table(colnames(exp_sum_hgc))
colnames(meta) <- 'sample'
meta[, "sgRNA" := tstrsplit(sample, "_")[2]]
sg <- unique(meta$sgRNA)
#DEG vs ref
for(i in 3:38){
sgid<-sg[i]
print(sgid)
meta1 <- subset(meta,sgRNA %in% c("control",sgid))
exp_sum_sub1 <- exp_sum_hgc[,meta1$sample]
ref <- read.table(sprintf('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/ref_sgrna_v2/%s_hgc_control_exp.xls',sgid),sep='\t',row.names=1,h=T)
meta_ref <- as.data.table(colnames(ref))
colnames(meta_ref) <- 'sample'
meta_ref$sgRNA <- 'ref_sample'
meta2 <- rbind(meta1,meta_ref)
meta2 <- subset(meta2,sgRNA != 'control')
exp_sum_sub2 <- cbind(exp_sum_sub1,ref)
exp_sum_sub2 <- exp_sum_sub2[,meta2$sample]
group <- meta2$sgRNA
y <- DGEList(counts=exp_sum_sub2,group = group)
keep <- filterByExpr(y,min.count = 5, min.total.count = 10)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y,method = 'TMM')
group=factor(y$samples$group)
design <- model.matrix(~group)
rownames(design) <- colnames(y)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmFit(y, design, robust=T)
qlf <- glmLRT(fit)
res = as.data.frame(topTags(qlf,n = nrow(y)))
res = cbind(gene_id=rownames(res),res)
write.table(res,sprintf("%s_vs_ref_sample_hgc.xls",sgid),sep="\t",quote=FALSE,row.names=F)
}

#cis-eGene
target <- fread('/Public/wtp/blj/singlecell/crop_seq/result_v2/snp_gene_1mb.txt')
library(data.table)
library(dplyr)
#ges
setwd('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/ges')
# get the filenames
filenames1 <- list.files(path=".", pattern="_vs_ref_sample_ges.xls")
# read all the files
df_ges <- lapply(filenames1, function (x) { tryCatch(fread(x), error=function(e) NULL)})
length(df_ges)
name_lst <- sub("_vs_ref_sample_ges.xls","",filenames1)
names(df_ges) <- name_lst
df_ges <- bind_rows(df_ges, .id="sgRNA")
colnames(df_ges)
df_ges$sg_gene <- paste(df_ges$sgRNA,df_ges$gene_id,sep='_')
target$sg_gene <- paste(target$RS,target$Symbol,sep='_')
df_ges1 <- subset(df_ges,sg_gene %in% c(target$sg_gene,'Essential_MVD'))
df_ges1$threshold = factor(ifelse(df_ges1$PValue < 0.05 & df_ges1$logFC < 0, ifelse(df_ges1$FDR < 0.05,'Down','Suggestive'),'Insignificance'),levels=c('Down','Suggestive','Insignificance'))
write.table(df_ges1,'/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/ges/deg_1mb_target_ges.xls',sep="\t",quote=FALSE,row.names=F)
#hgc
setwd('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/hgc')
# get the filenames
filenames1 <- list.files(path=".", pattern="_vs_ref_sample_hgc.xls")
# read all the files
df_hgc <- lapply(filenames1, function (x) { tryCatch(fread(x), error=function(e) NULL)})
length(df_hgc)
name_lst <- sub("_vs_ref_sample_hgc.xls","",filenames1)
names(df_hgc) <- name_lst
df_hgc <- bind_rows(df_hgc, .id="sgRNA")
colnames(df_hgc)
df_hgc$sg_gene <- paste(df_hgc$sgRNA,df_hgc$gene_id,sep='_')
target$sg_gene <- paste(target$RS,target$Symbol,sep='_')
df_hgc1 <- subset(df_hgc,sg_gene %in% c(target$sg_gene,'Essential_MVD'))
df_hgc1$threshold = factor(ifelse(df_hgc1$PValue < 0.05 & df_hgc1$logFC < 0, ifelse(df_hgc1$FDR < 0.05,'Down','Suggestive'),'Insignificance'),levels=c('Down','Suggestive','Insignificance'))
write.table(df_hgc1,'/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/hgc/deg_1mb_target_hgc.xls',sep="\t",quote=FALSE,row.names=F)

df_ges1 <- read.table('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/ges/deg_1mb_target_ges.xls',h=T)
df_hgc1 <- read.table('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/hgc/deg_1mb_target_hgc.xls',h=T)
df_ges1 <- subset(df_ges1,sgRNA!='rs57214277')
df_hgc1 <- subset(df_hgc1,sgRNA!='rs57214277')
write.table(df_ges1,'/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/ges/deg_1mb_target_ges.xls',sep="\t",quote=FALSE,row.names=F)
write.table(df_hgc1,'/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/hgc/deg_1mb_target_hgc.xls',sep="\t",quote=FALSE,row.names=F)

#Figure_5F
df_ges2 <- subset(df_ges1,threshold=='Down') 
df_hgc2 <- subset(df_hgc1,threshold=='Down')
both1 <- intersect(df_ges2$sg_gene,df_hgc2$sg_gene)
both <- unique(c(both1,df_ges2$sg_gene,df_hgc2$sg_gene))
both <- as.data.table(both)
colnames(both) <- 'sg_gene'
both[, "sg" := tstrsplit(sg_gene, "_")[1]]
both[, "gene" := tstrsplit(sg_gene, "_")[2]]
loci <- read.csv('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/loci_id.csv')
both <- left_join(both,loci,by=c('sg'='rsid'))
both <- as.data.frame(both)
both<-both[order(both[,'order']),]
ges_plot <- left_join(both,df_ges1,by = 'sg_gene')
hgc_plot <- left_join(both,df_hgc1,by = 'sg_gene')
ges_plot$sg_gene <- paste(ges_plot$loci,ges_plot$sg_gene,sep='_')
hgc_plot$sg_gene <- paste(hgc_plot$loci,hgc_plot$sg_gene,sep='_')
ges_plot1 <- subset(ges_plot,gene!='MVD')
hgc_plot1 <- subset(hgc_plot,gene!='MVD')
both2 <- both1[-1]
ges_plot1$sg_gene1 <- paste(ges_plot1$sg,ges_plot1$gene,sep='_')
ges_plot2 <- subset(ges_plot1,sg_gene1%in%both2)
ges_plot3 <- subset(ges_plot1,!sg_gene1%in%both2)
ges_plot4 <- rbind(ges_plot2,ges_plot3)
hgc_plot1$sg_gene1 <- paste(hgc_plot1$sg,hgc_plot1$gene,sep='_')
hgc_plot2 <- subset(hgc_plot1,sg_gene1%in%both2)
hgc_plot3 <- subset(hgc_plot1,!sg_gene1%in%both2)
hgc_plot4 <- rbind(hgc_plot2,hgc_plot3)
ges_plot4<-ges_plot4[order(ges_plot4[,'order']),]
hgc_plot4<-hgc_plot4[order(hgc_plot4[,'order']),]
library(ggplot2)
ges_plot4$sg_gene <-  factor(ges_plot4$sg_gene , levels=ges_plot4$sg_gene, ordered=TRUE)
p<-ggplot(ges_plot4,aes(x=sg_gene,y=logFC,color=threshold))+
  geom_point(size=1.6)+
  scale_color_manual(values=c("Down"="#AE0000","Suggestive" = "#FF9797"))+
  theme_bw()+
  theme(legend.position="right",
  axis.text = element_text(size = 11),axis.title = element_text(size = 14),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  geom_hline(yintercept = 0,lty=3,col="black",lwd=0.5)+
  ylab('Log2 fold change')+
  xlab('')
ggsave('ges_1mb_both_fdr1.pdf', p, width = 12, height = 6.5)
hgc_plot4$sg_gene <-  factor(hgc_plot4$sg_gene , levels=hgc_plot4$sg_gene, ordered=TRUE)
p<-ggplot(hgc_plot4,aes(x=sg_gene,y=logFC,color=threshold))+
  geom_point(size=1.6)+
  scale_color_manual(values=c("Down"="#AE0000","Suggestive" = "#FF9797"))+
  theme_bw()+
  theme(legend.position="right",
  axis.text = element_text(size = 11),axis.title = element_text(size = 14),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  geom_hline(yintercept = 0,lty=3,col="black",lwd=0.5)+
  ylab('Log2 fold change')+
  xlab('')
ggsave('hgc_1mb_both_fdr1.pdf', p, width = 12, height = 6.5)


#Figure_S13J
df_ges2 <- subset(df_ges1,threshold!='Insignificance')
df_hgc2 <- subset(df_hgc1,threshold!='Insignificance')
both1 <- intersect(df_ges2$sg_gene,df_hgc2$sg_gene)
both <- unique(c(both1,df_ges2$sg_gene,df_hgc2$sg_gene))
both <- as.data.table(both)
colnames(both) <- 'sg_gene'
both[, "sg" := tstrsplit(sg_gene, "_")[1]]
both[, "gene" := tstrsplit(sg_gene, "_")[2]]
loci <- read.csv('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/loci_id.csv')
both <- left_join(both,loci,by=c('sg'='rsid'))
both <- as.data.frame(both)
both<-both[order(both[,'order']),]
ges_plot <- left_join(both,df_ges1,by = 'sg_gene')
hgc_plot <- left_join(both,df_hgc1,by = 'sg_gene')
ges_plot$sg_gene <- paste(ges_plot$loci,ges_plot$sg_gene,sep='_')
hgc_plot$sg_gene <- paste(hgc_plot$loci,hgc_plot$sg_gene,sep='_')
ges_plot1 <- subset(ges_plot,gene!='MVD')
hgc_plot1 <- subset(hgc_plot,gene!='MVD')
both2 <- both1[-1]
ges_plot1$sg_gene1 <- paste(ges_plot1$sg,ges_plot1$gene,sep='_')
ges_plot2 <- subset(ges_plot1,sg_gene1%in%both2)
ges_plot3 <- subset(ges_plot1,!sg_gene1%in%both2)
ges_plot4 <- rbind(ges_plot2,ges_plot3)
hgc_plot1$sg_gene1 <- paste(hgc_plot1$sg,hgc_plot1$gene,sep='_')
hgc_plot2 <- subset(hgc_plot1,sg_gene1%in%both2)
hgc_plot3 <- subset(hgc_plot1,!sg_gene1%in%both2)
hgc_plot4 <- rbind(hgc_plot2,hgc_plot3)
ges_plot4<-ges_plot4[order(ges_plot4[,'order']),]
hgc_plot4<-hgc_plot4[order(hgc_plot4[,'order']),]
library(ggplot2)
ges_plot4$sg_gene <-  factor(ges_plot4$sg_gene , levels=ges_plot4$sg_gene, ordered=TRUE)
p<-ggplot(ges_plot4,aes(x=sg_gene,y=logFC,color=threshold))+
  geom_point(size=1.6)+
  scale_color_manual(values=c("Down"="#AE0000","Suggestive" = "#FF9797"))+
  theme_bw()+
  theme(legend.position="right",
  axis.text = element_text(size = 11),axis.title = element_text(size = 14),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  geom_hline(yintercept = 0,lty=3,col="black",lwd=0.5)+
  ylab('Log2 fold change')+
  xlab('')
ggsave('ges_1mb_both_fdr2.pdf', p, width = 12, height = 6.5)
hgc_plot4$sg_gene <-  factor(hgc_plot4$sg_gene , levels=hgc_plot4$sg_gene, ordered=TRUE)
p<-ggplot(hgc_plot4,aes(x=sg_gene,y=logFC,color=threshold))+
  geom_point(size=1.6)+
  scale_color_manual(values=c("Down"="#AE0000","Suggestive" = "#FF9797"))+
  theme_bw()+
  theme(legend.position="right",
  axis.text = element_text(size = 11),axis.title = element_text(size = 14),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  geom_hline(yintercept = 0,lty=3,col="black",lwd=0.5)+
  ylab('Log2 fold change')+
  xlab('')
ggsave('hgc_1mb_both_fdr2.pdf', p, width = 12, height = 6.5)

##Figure_S13I
ges_plot_mvd <- subset(ges_plot,gene_id=='MVD')
hgc_plot_mvd <- subset(hgc_plot,gene_id=='MVD')
ges_plot_mvd$group <- 'GES-1' 
hgc_plot_mvd$group <- 'HGC-27' 
all_mvd <- rbind(ges_plot_mvd,hgc_plot_mvd)
all_mvd$signp <- all_mvd$logFC * -log10(all_mvd$PValue)
p <- ggplot(all_mvd) +
   theme_bw()+
  geom_segment(aes(x=group, xend=group, y=0, yend=signp), color="grey",size=1) +
  geom_point(aes(x=group, y=signp, color=threshold), size=4) +
  scale_color_manual(values = c("#f391a9"),
                     labels = c("FDR < 0.05"),
                     name = "Significance")+
  geom_hline(yintercept = 0, lty=2,color = 'grey', lwd=0.8)+
  labs(title = "Downregulation of MVD",
       x = "",
       y = "Signed -log10P")+
  theme(axis.text.x = element_text(size=14))
ggsave('Downregulation of MVD.pdf', p, width = 4.5, height = 5) 

##Figure_5I
cd /Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2
library(data.table)
library(dplyr)
library(ggplot2)
library(ggprism)
all_hgc <- fread('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/hgc/deg_1mb_target_hgc.xls')
all_ges <- fread('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/ges/deg_1mb_target_ges.xls')
res_hgc <- subset(all_hgc,sgRNA=='rs7927406')
res_ges <- subset(all_ges,sgRNA=='rs7927406')
dat <- intersect(res_hgc$gene_id,res_ges$gene_id)
#dat <- dat[dat!='VWA5A']
res_hgc <- subset(res_hgc,gene_id%in%dat)
res_ges <- subset(res_ges,gene_id%in%dat)
res_hgc$signp <- res_hgc$logFC * -log10(res_hgc$PValue)
res_ges$signp <- res_ges$logFC * -log10(res_ges$PValue)
res_hgc$gene_id <- factor(res_hgc$gene_id, levels=res_hgc$gene_id, ordered=TRUE)
res_hgc$threshold <- factor(res_hgc$threshold, levels=c('Down','Suggestive','Insignificance'), ordered=TRUE)
p <- ggplot(res_hgc) +
   theme_bw()+
  geom_segment(aes(x=gene_id, xend=gene_id, y=0, yend=signp), color="grey",size=1) +
  geom_point(aes(x=gene_id, y=signp, color=threshold), size=4) +
  scale_color_manual(values = c("Down"="#AE0000","Suggestive" = "#FF9797","Insignificance"="gray50"),
                     name = "Significance")+
  geom_hline(yintercept = 0, lty=2,color = 'grey', lwd=0.8)+
  labs(title = "Downregulation of genes by perturbation of rs7927406 in HGC-27",
       x = "",
       y = "Signed -log10P")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
ggsave('Downregulation of genes by perturbation of rs7927406_hgc.pdf', p, width = 6, height = 5)  

##Figure_S13M
cd /Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2
library(data.table)
library(dplyr)
library(ggplot2)
library(ggprism)
all_hgc <- fread('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/hgc/deg_1mb_target_hgc.xls')
all_ges <- fread('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/ges/deg_1mb_target_ges.xls')
res_hgc <- subset(all_hgc,sgRNA=='rs56747346')
#res_ges <- subset(all_ges,sgRNA=='rs56747346')
res_ges <- read.csv('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/deg_1mb_target_ges_muc1.csv')
res_ges1 <- res_ges[,c(1,2)]
res_hgc <- left_join(res_ges1,res_hgc,by='gene_id')
res_hgc$signp <- res_hgc$logFC * -log10(res_hgc$PValue)
res_ges$signp <- res_ges$logFC * -log10(res_ges$PValue)
res_hgc$gene_id <- factor(res_hgc$gene_id, levels=res_hgc$gene_id, ordered=TRUE)
res_hgc$threshold <- factor(res_hgc$threshold, levels=c('Down','Suggestive','Insignificance'), ordered=TRUE)
p <- ggplot(res_hgc) +
   theme_bw()+
  geom_segment(aes(x=gene_id, xend=gene_id, y=0, yend=signp), color="grey",size=1) +
  geom_point(aes(x=gene_id, y=signp, color=threshold), size=3) +
  scale_color_manual(values = c("Down"="#AE0000","Suggestive" = "#FF9797","Insignificance"="gray50"),
                     name = "Significance")+
  geom_hline(yintercept = 0, lty=2,color = 'grey', lwd=0.8)+
  labs(title = "Downregulation of genes by perturbation of rs56747346 in HGC-27",
       x = "",
       y = "Signed -log10P")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('Downregulation of genes by perturbation of rs56747346_hgc.pdf', p, width = 8, height = 5)  

res_ges$gene_id <- factor(res_ges$gene_id, levels=res_ges$gene_id, ordered=TRUE)
res_ges$threshold <- factor(res_ges$threshold, levels=c('Down','Suggestive','Insignificance'), ordered=TRUE)
p <- ggplot(res_ges) +
   theme_bw()+
  geom_segment(aes(x=gene_id, xend=gene_id, y=0, yend=signp), color="grey",size=1) +
  geom_point(aes(x=gene_id, y=signp, color=threshold), size=3) +
  scale_color_manual(values = c("Down"="#AE0000","Suggestive" = "#FF9797","Insignificance"="gray50"),
                     name = "Significance")+
  geom_hline(yintercept = 0, lty=2,color = 'grey', lwd=0.8)+
  labs(title = "Downregulation of genes by perturbation of rs56747346 in GES-1",
       x = "",
       y = "Signed -log10P")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('Downregulation of genes by perturbation of rs56747346_ges.pdf', p, width = 8, height = 5)  

####Figure_S13N
cd /Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2
library(data.table)
library(dplyr)
library(ggplot2)
library(ggprism)
all_hgc <- fread('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/hgc/deg_1mb_target_hgc.xls')
all_ges <- fread('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/ges/deg_1mb_target_ges.xls')
res_hgc <- subset(all_hgc,sgRNA=='rs145938518')
#res_ges <- subset(all_ges,sgRNA=='rs145938518')
res_ges <- read.csv('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/deg_1mb_target_ges_KRTCAP2.csv')
res_ges1 <- res_ges[,c(1,2)]
res_hgc <- left_join(res_ges1,res_hgc,by='gene_id')
res_hgc$signp <- res_hgc$logFC * -log10(res_hgc$PValue)
res_ges$signp <- res_ges$logFC * -log10(res_ges$PValue)
res_hgc$gene_id <- factor(res_hgc$gene_id, levels=res_hgc$gene_id, ordered=TRUE)
res_hgc$threshold <- factor(res_hgc$threshold, levels=c('Down','Suggestive','Insignificance'), ordered=TRUE)
p <- ggplot(res_hgc) +
   theme_bw()+
  geom_segment(aes(x=gene_id, xend=gene_id, y=0, yend=signp), color="grey",size=1) +
  geom_point(aes(x=gene_id, y=signp, color=threshold), size=3) +
  scale_color_manual(values = c("Down"="#AE0000","Suggestive" = "#FF9797","Insignificance"="gray50"),
                     name = "Significance")+
  geom_hline(yintercept = 0, lty=2,color = 'grey', lwd=0.8)+
  labs(title = "Downregulation of genes by perturbation of rs145938518 in HGC-27",
       x = "",
       y = "Signed -log10P")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('Downregulation of genes by perturbation of rs145938518_hgc.pdf', p, width = 8, height = 5)  
res_ges$gene_id <- factor(res_ges$gene_id, levels=res_ges$gene_id, ordered=TRUE)
res_ges$threshold <- factor(res_ges$threshold, levels=c('Down','Suggestive','Insignificance'), ordered=TRUE)
p <- ggplot(res_ges) +
   theme_bw()+
  geom_segment(aes(x=gene_id, xend=gene_id, y=0, yend=signp), color="grey",size=1) +
  geom_point(aes(x=gene_id, y=signp, color=threshold), size=3) +
  scale_color_manual(values = c("Down"="#AE0000","Suggestive" = "#FF9797","Insignificance"="gray50"),
                     name = "Significance")+
  geom_hline(yintercept = 0, lty=2,color = 'grey', lwd=0.8)+
  labs(title = "Downregulation of genes by perturbation of rs145938518 in GES-1",
       x = "",
       y = "Signed -log10P")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('Downregulation of genes by perturbation of rs145938518_ges.pdf', p, width = 8, height = 5)  

####Figure_S13O
cd /Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2
library(data.table)
library(dplyr)
library(ggplot2)
library(ggprism)
all_ges <- fread('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/ges/deg_1mb_target_ges.xls')
res_ges <- subset(all_ges,sgRNA=='rs28734952')
res_ges <- subset(all_ges,sgRNA=='rs2976385')
res_ges <- subset(all_ges,sgRNA=='rs2978980')
res_ges <- subset(all_ges,sgRNA=='rs2572874 ')
res_ges$signp <- res_ges$logFC * -log10(res_ges$PValue)
res_ges$gene_id <- factor(res_ges$gene_id, levels=res_ges$gene_id, ordered=TRUE)
res_ges$threshold <- factor(res_ges$threshold, levels=c('Down','Suggestive','Insignificance'), ordered=TRUE)
p <- ggplot(res_ges) +
   theme_bw()+
  geom_segment(aes(x=gene_id, xend=gene_id, y=0, yend=signp), color="grey",size=1) +
  geom_point(aes(x=gene_id, y=signp, color=threshold), size=3) +
  scale_color_manual(values = c("Down"="#AE0000","Suggestive" = "#FF9797","Insignificance"="gray50"),
                     name = "Significance")+
  geom_hline(yintercept = 0, lty=2,color = 'grey', lwd=0.8)+
  labs(title = "Downregulation of genes by perturbation of rs145938518 in GES-1",
       x = "",
       y = "Signed -log10P")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('Downregulation of genes by perturbation of rs145938518_ges.pdf', p, width = 8, height = 5)  

##Figure_5G
df_ges <- fread('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/ges/deg_1mb_target_ges.xls')
df_hgc <- fread('/Public/wtp/blj/singlecell/crop_seq/result_v2/sum_counts/control_deg_v2/hgc/deg_1mb_target_hgc.xls')
df_ges1 <- subset(df_ges,sgRNA!='rs57214277')
df_hgc1 <- subset(df_hgc,sgRNA!='rs57214277')
df_ges1$group <- 'ges'
df_hgc1$group <- 'hgc'
df_all <- rbind(df_ges1,df_hgc1)
df_all$FC <- 2^df_all$logFC
df_all$percentage <- df_all$FC-1
df_all$group1[df_all$threshold %in% c("Down")]<-"Hit"
df_all$group1[df_all$threshold %in% c("Suggestive")]<-"Hit-Suggestive"
df_all$group1[df_all$threshold %in% c("Insignificance")]<-"Non-Hit"
df_all$percentage[df_all$percentage > 1]<-1
p=ggboxplot(df_all,'group1', 'percentage',fill = "group1",palette = c("#f47920","#1d953f","#4e72b8"),add = "none")+
   labs(x = '', y = 'percentage downregulation')+guides(fill=F)
ggsave('percentage downregulation.pdf', p, width = 4, height =5)

##Figure_S14E
setwd('/data/gc_sceqtl/sceqtl/seuratall/allfilter_orig/final')
load('MIXALL_final.Rdata')
library(Seurat)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(Nebulosa)
DefaultAssay(MIXALL_final) <- "RNA"
markers <- c('SPA17')
p <- plot_density(MIXALL_final, markers, reduction = 'umap',joint = FALSE, combine = FALSE, direction = -1)+
     scale_color_continuous_sequential(palette = "Sunset")
ggsave('density_markers_SPA17.png', p, width = 6, height = 6)