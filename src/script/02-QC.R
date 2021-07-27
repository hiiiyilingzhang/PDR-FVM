## Data examining ----
# cell number:
dim(integData)

# Counts per cell:
counts_per_cell <- Matrix::colSums(integData)
hist(log10(counts_per_cell+1),main='counts per cell',col='Thistle')

## Genes per cell:
genes_per_cell <- Matrix::colSums(integData@assays$integrated@data > 0)
hist(log10(genes_per_cell+1), main='genes per cell', col='Thistle')

## Remove Doublets ----
# DoubletDecon 
# https://github.com/EDePasquale/DoubletDecon
# Fail to install

# DoubletFinder
# https://github.com/chris-mcginnis-ucsf/DoubletFinder
library(DoubletFinder)
library(patchwork)
library(clustree)

inteData <- ScaleData(inteData, verbose = FALSE) %>% RunPCA(npcs = 30)

ElbowPlot(inteData,ndims = 30)
pc.num=1:20 # Set dim 1:20

inteData <- RunUMAP(inteData, dims=pc.num,verbose = F)
inteData <- FindNeighbors(inteData, reduction = "pca", dims = pc.num)

checkRes <- FindClusters(object = inteData, resolution = c(seq(.1,1.6,.2)), verbose = F)
clustree(checkRes@meta.data,prefix = "integrated_snn_res.")
inteData <- FindClusters(inteData, resolution = 0.7)

# Optimize the parameters
sweep.res.list <- paramSweep_v3(inteData, PCs = pc.num, sct = T)
# Use log transform
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
# Show the best parameter
bcmvn <- find.pK(sweep.stats)
# Extract the best pK
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
pK_bcmvn <- 0.06

# Run DoubletFinder
# Exclude doublets
# doublet rate: https://singleronbio.com/product/?type=detail&id=8
DoubletRate = 0.0268
# Estimate the percentage of homotypic doublets
homotypic.prop <- modelHomotypic(inteData$seurat_clusters)
nExp_poi <- round(DoubletRate*ncol(inteData)) 
# Adjust for homotypic doublets
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
# Run DoubletFinder with varying classification stringencies
inteData <- doubletFinder_v3(inteData, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                             nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
# Present the res, classification info is saved in meta.data
DF.name <- colnames(inteData@meta.data)[grepl("DF.classification", colnames(inteData@meta.data))]
# Visualization
pdf(file = "results/QC-DoubletFinder/res-plot.pdf", width = 8, height = 4)
cowplot::plot_grid(ncol = 2, DimPlot(inteData, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(inteData, group.by = DF.name) + NoAxes())
dev.off()

## Calculate QC ----

# ---- troubleshooting ----
# https://github.com/satijalab/seurat/issues/3596
# Error in h(simpleError(msg, call)) : 
#   error in evaluating the argument 'x' in selecting a method for function 'colSums': no 'dimnames[[.]]': cannot use character indexing
# 
# PercentageFeatureSet uses the raw counts matrix, which is not present in the integrated assay. 
# If you pass assay = "RNA"/"SCT", then PercentageFeatureSet will run correctly
# -------------------------

# mt
inteData[["percent.mt"]] <- PercentageFeatureSet(object = inteData, pattern = "^MT-", assay = "RNA")
inteData$percent.mt <- inteData@meta.data$percent.mt / 100
# ribosome
inteData[["percent.ribo"]] <- PercentageFeatureSet(object = inteData, pattern = "^RP[SL]", assay = "RNA")
# hb
inteData[["percent.hb"]] <- PercentageFeatureSet(object = inteData, pattern = "^HB[^(P)]", assay = "RNA")
# genes per UMI for each cell
inteData$log10GenesPerUMI <- log10(inteData$nFeature_RNA) / log10(inteData$nCount_RNA)

## Plot QC ----
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb")
pdf(file = "results/QC-Calculate/vinplot.pdf", height = 6, width = 10)
VlnPlot(inteData, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
dev.off()

# VlnPlot(inteData, group.by = "orig.ident", features = "percent.hb", pt.size = 0.1, y.max = 1)
# hb???? weird

## Filtering ----
# Detection-based filtering
par(mar = c(4, 8, 2, 1))
C <- inteData@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed,])), cex = 0.1, las = 1, xlab = "% total count per cell", col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

# Mito/Ribo/hb filtering
selected_mito <- WhichCells(inteData, expression = percent.mt < 20)
selected_ribo <- WhichCells(inteData, expression = percent.ribo > 5)
