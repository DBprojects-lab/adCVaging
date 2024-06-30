library(Seurat)
library(dplyr)
library(patchwork)
library(SingleR)
library(patchwork)
library(DoubletFinder)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggpubr)

setwd("")
count <- readRDS("ROSMAP.VascularCells.counts_full.rds")
meta <- readRDS("ROSMAP.VascularCells.meta_full.rds")
dim(count)
dim(meta)
table(meta$brain_region)

meta <- meta[meta$brain_region == "Prefrontal_cortex", ]
#Calculating all brain regions is not necessary for this step.

count <- count[, rownames(meta)]
identical(colnames(count), rownames(meta))
table(meta$cellsubtype)

pbmc_ctx <- CreateSeuratObject(counts = count)
colnames(meta)
pbmc_ctx@meta.data$batch <- meta$batch
pbmc_ctx@meta.data$apoe <- meta$apoe
pbmc_ctx@meta.data$nFeature_RNA <- meta$nFeature_RNA
pbmc_ctx@meta.data$nCount_RNA <- meta$nCount_RNA
pbmc_ctx@meta.data$race <- meta$race
pbmc_ctx@meta.data$pmi <- meta$pmi
pbmc_ctx@meta.data$age_bl <- meta$age_bl
pbmc_ctx@meta.data$age_death <- meta$age_death
pbmc_ctx@meta.data$orig.ident <- meta$id
pbmc_ctx@meta.data$ADdiag2types <- meta$ADdiag2types
pbmc_ctx@meta.data$ADdiag3types <- meta$ADdiag3types
pbmc_ctx@meta.data$cogdx <- meta$cogdx
pbmc_ctx@meta.data$braaksc <- meta$braaksc
pbmc_ctx@meta.data$nFeature_RNA <- meta$nFeature_RNA
pbmc_ctx@meta.data$percent.mt <- meta$percent.mt
pbmc_ctx@meta.data$percent.rp <- meta$percent.rp

pbmc_ctx@meta.data$celltype <- meta$celltype
pbmc_ctx@meta.data$brain_region <- meta$brain_region
pbmc_ctx@meta.data$msex <- meta$msex
pbmc_ctx@meta.data$msex[pbmc_ctx@meta.data$msex == 0] <- "female"
pbmc_ctx@meta.data$msex[pbmc_ctx@meta.data$msex == 1] <- "male"


pbmc_ctx@meta.data$orig.ident2 <- 0
pbmc_ctx@meta.data$orig.ident2[(pbmc_ctx@meta.data$msex == "female" &
                                  pbmc_ctx@meta.data$ADdiag2types ==
                                  "AD")] <- "Ctx_AD_F"

pbmc_ctx@meta.data$orig.ident2[(pbmc_ctx@meta.data$msex == "male" &
                                  pbmc_ctx@meta.data$ADdiag2types ==
                                  "AD")] <- "Ctx_AD_M"
pbmc_ctx@meta.data$orig.ident2[(pbmc_ctx@meta.data$msex == "female" &
                                  pbmc_ctx@meta.data$ADdiag2types ==
                                  "nonAD")] <- "Ctx_C_F"
pbmc_ctx@meta.data$orig.ident2[(pbmc_ctx@meta.data$msex == "male" &
                                  pbmc_ctx@meta.data$ADdiag2types ==
                                  "nonAD")] <- "Ctx_C_M"

table(pbmc_ctx$orig.ident2)
table(pbmc_ctx$msex, pbmc_ctx$ADdiag2types)

pbmc_ctx <-
  subset(pbmc_ctx,
         subset = nFeature_RNA > 200 &
           nFeature_RNA < 5000 & percent.mt < 5 & percent.rp < 5)
pbmc_ctx <-
  NormalizeData(pbmc_ctx,
                normalization.method = "LogNormalize",
                scale.factor = 10000)
pbmc_ctx <-
  FindVariableFeatures(pbmc_ctx, selection.method = "vst", nfeatures = 2000)
pbmc_ctx <- ScaleData(pbmc_ctx, features = rownames(pbmc_ctx))
pbmc_ctx <-
  RunPCA(pbmc_ctx, features = VariableFeatures(object = pbmc_ctx))
k = 1:30
pbmc_ctx <- FindNeighbors(pbmc_ctx, dims = k)
pbmc_ctx <- FindClusters(pbmc_ctx, resolution = 0.5)
head(Idents(pbmc_ctx), 5)
pbmc_ctx <- RunUMAP(pbmc_ctx, dims = k)

library(harmony)
pbmc_ctx <- RunHarmony(pbmc_ctx, c("batch"), plot_convergence = T)
harmony_embeddings <- Embeddings(pbmc_ctx, 'harmony')


pbmc_ctx <- RunUMAP(pbmc_ctx, reduction = "harmony", dims = k)
pbmc_ctx <- FindNeighbors(pbmc_ctx, reduction = "harmony", dims = k)
pbmc_ctx <- FindClusters(pbmc_ctx, resolution = 0.5)
identity(pbmc_ctx)

sweep.res.list <- paramSweep_v3(pbmc_ctx, PCs = 1:30, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)


colnames(sweep.stats)
pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK <- as.numeric(as.character(pK[[1]]))
pK = 0.01
annotations <- pbmc_ctx@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(ncol(pbmc_ctx) * 8 * 1e-6 * nrow(pbmc_ctx@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
pbmc_ctx <- doubletFinder_v3(
  pbmc_ctx,
  PCs = 1:30,
  pN = 0.2,
  pK = pK,
  nExp = nExp_poi.adj,
  reuse.pANN = F,
  sct = F
)

names(pbmc_ctx@meta.data)
names(pbmc_ctx@meta.data)[23] <- "DF.classifications"
names(pbmc_ctx@meta.data)[22] <- "pANN"
pbmc_ctx <- subset(pbmc_ctx, subset = DF.classifications == "Singlet")

DimPlot(
  pbmc_ctx,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  raster = FALSE,
  group.by = "celltype"
)

#save(pbmc_ctx,file = "/boot3/LPH/111data/ROSMAP/pbmc_ctx_new.Rdata")

load("/boot3/LPH/111data/ROSMAP/pbmc_ctx_new.Rdata")

data <- pbmc_ctx[["RNA"]]@data

table(pbmc_ctx$brain_region)

DimPlot(
  pbmc_ctx,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  raster = F,
  group.by = "celltype"
)

cell1 <-
  colnames(pbmc_ctx)[intersect(which(pbmc_ctx$celltype == "Endo"),
                               which(pbmc_ctx$orig.ident2 == "Ctx_AD_M"))]
cell2 <-
  colnames(pbmc_ctx)[intersect(which(pbmc_ctx$celltype == "Endo"),
                               which(pbmc_ctx$orig.ident2 == "Ctx_C_M"))]
earlyAD_male_cell <-
  colnames(pbmc_ctx)[intersect(intersect(
    which(pbmc_ctx$celltype == "Endo"),
    which(pbmc_ctx$ADdiag3types == "earlyAD")
  ),
  which(pbmc_ctx$msex == "male"))]
lateAD_mele_cell <-
  colnames(pbmc_ctx)[intersect(intersect(
    which(pbmc_ctx$celltype == "Endo"),
    which(pbmc_ctx$ADdiag3types ==
            "lateAD")
  ),
  which(pbmc_ctx$msex == "male"))]
nonAD_male_cell <-
  colnames(pbmc_ctx)[intersect(which(pbmc_ctx$celltype == "Endo"),
                               which(pbmc_ctx$orig.ident2 == "Ctx_C_M"))]


#Calculate differentially DEGs between different cell groups
male_DEGs <- FindMarkers(
  pbmc_ctx,
  ident.1 = cell1,
  #replaceable
  ident.2 = cell2,
  #replaceable
  logfc.threshold = 0,
  test.use = "MAST",
  latent.vars = c(
    "race",
    "pmi",
    "nCount_RNA",
    "nFeature_RNA",
    "age_bl",
    "age_death",
    "batch"
  )
)
#save(male_DEGs,file = "")



cell3 <-
  colnames(pbmc_ctx)[intersect(which(pbmc_ctx$celltype == "Endo"),
                               which(pbmc_ctx$orig.ident2 == "Ctx_AD_F"))]
cell4 <-
  colnames(pbmc_ctx)[intersect(which(pbmc_ctx$celltype == "Endo"),
                               which(pbmc_ctx$orig.ident2 == "Ctx_C_F"))]
nonAD_female_cell <-
  colnames(pbmc_ctx)[intersect(which(pbmc_ctx$celltype == "Endo"),
                               which(pbmc_ctx$orig.ident2 == "Ctx_C_F"))]
earlyAD_female_cell <-
  colnames(pbmc_ctx)[intersect(intersect(
    which(pbmc_ctx$celltype == "Endo"),
    which(pbmc_ctx$ADdiag3types ==
            "earlyAD")
  ),
  which(pbmc_ctx$msex ==
          "female"))]
lateAD_femele_cell <-
  colnames(pbmc_ctx)[intersect(intersect(
    which(pbmc_ctx$celltype == "Endo"),
    which(pbmc_ctx$ADdiag3types ==
            "lateAD")
  ),
  which(pbmc_ctx$msex ==
          "female"))]

#Calculate differentially DEGs between different cell groups
female_DEGs <- FindMarkers(
  pbmc_ctx,
  ident.1 = cell3,
  #replaceable
  ident.2 = cell4,
  #replaceable
  logfc.threshold = 0,
  test.use = "MAST",
  latent.vars = c(
    "race",
    "pmi",
    "nCount_RNA",
    "nFeature_RNA",
    "age_bl",
    "age_death",
    "batch"
  )
)

#save(female_DEGs,file = "")
#Select DEGs that meet the threshold
#F_UP <- rownames(female_DEGs[female_DEGs$p_val<0.05&female_DEGs$avg_log2FC>0.1,])
#F_DOWN <- rownames(female_DEGs[female_DEGs$p_val<0.05&female_DEGs$avg_log2FC<(-0.1),])
#M_UP <- rownames(male_DEGs[male_DEGs$p_val<0.05&male_DEGs$avg_log2FC>0.1,])
#M_DOWN <- rownames(male_DEGs[male_DEGs$p_val<0.05&male_DEGs$avg_log2FC<(-0.1),])


#DEGs of other cell types
table(pbmc_ctx$celltype)
cell1 <-
  colnames(pbmc_ctx)[intersect(which(pbmc_ctx$celltype == "Per"),
                               #replaceable
                               which(pbmc_ctx$orig.ident2 == "Ctx_AD_M"))]
cell2 <-
  colnames(pbmc_ctx)[intersect(which(pbmc_ctx$celltype == "Per"),
                               #replaceable
                               which(pbmc_ctx$orig.ident2 == "Ctx_C_M"))]
table(pbmc_ctx$orig.ident2, pbmc_ctx$celltype)
length(cell1)#1001
length(cell2)#1473
table(pbmc_ctx$celltype, pbmc_ctx$orig.ident2)
M <- FindMarkers(
  pbmc_ctx,
  ident.1 = cell1,
  ident.2 = cell2,
  #features = features,
  logfc.threshold = 0,
  test.use = "MAST",
  latent.vars = c(
    "race",
    "pmi",
    "nCount_RNA",
    "nFeature_RNA",
    "age_bl",
    "age_death",
    "batch"
  )
)

cell3 <-
  colnames(pbmc_ctx)[intersect(which(pbmc_ctx$celltype == "Per"),
                               #replaceable
                               which(pbmc_ctx$orig.ident2 == "Ctx_AD_F"))]
cell4 <-
  colnames(pbmc_ctx)[intersect(which(pbmc_ctx$celltype == "Per"),
                               #replaceable
                               which(pbmc_ctx$orig.ident2 == "Ctx_C_F"))]

F <- FindMarkers(
  pbmc_ctx,
  ident.1 = cell3,
  ident.2 = cell4,
  logfc.threshold = 0,
  test.use = "MAST",
  latent.vars = c(
    "race",
    "pmi",
    "nCount_RNA",
    "nFeature_RNA",
    "age_bl",
    "age_death",
    "batch"
  )
)

#save(M,file = "")
#save(F,file = "")

