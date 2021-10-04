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
library(ggraph)
library(clustree)
library(ggsci)

# Pre-process Seurat object
raw.Singleron <- SCTransform(raw.Singleron, verbose = F)
raw.Singleron <- RunPCA(raw.Singleron, npcs = 30, verbose = F)
ElbowPlot(raw.Singleron,ndims = 30)

# Set dim 1:20
pc.num=1:25

raw.Singleron <- RunUMAP(raw.Singleron, dims=pc.num,verbose = F)
raw.Singleron <- FindNeighbors(raw.Singleron, reduction = "pca", dims = pc.num)

checkRes <- FindClusters(object = raw.Singleron, resolution = c(seq(.1,1.4,.2)), verbose = F)
clustree(checkRes@meta.data,prefix = "SCT_snn_res.")

raw.Singleron <- FindClusters(raw.Singleron, resolution = 0.7)

# pK Identification

# Optimize the parameters
sweep.res.list <- paramSweep_v3(raw.Singleron, PCs = pc.num, sct = T)
# Use log transform
# Show the best parameter
bcmvn <- find.pK(sweep.stats)

# Extract the best pK
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
pK_bcmvn

# Homotypic Doublet Proportion Estimate
# doublet rate: https://www.singleronbio.com/resources/list-19.html. Because of the various distribution of cell number detected, we set 2.68% based on the Standard Chip, just in case we miss important data especially p2-ERM(589 cells). 

# Estimate the percentage of homotypic doublets
homotypic.prop <- modelHomotypic(raw.Singleron$seurat_clusters)
nExp_poi <- round(DoubletRate*ncol(raw.Singleron)) 
# Adjust for homotypic doublets
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder
raw.Singleron <- doubletFinder_v3(raw.Singleron, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = T)

# Present the res, classification info is saved in meta.data
DF.name <- colnames(raw.Singleron@meta.data)[grepl("DF.classification", colnames(raw.Singleron@meta.data))]
# Visualization
cowplot::plot_grid(ncol = 2, DimPlot(raw.Singleron, group.by = "orig.ident") + NoAxes(),DimPlot(raw.Singleron, group.by = DF.name) + NoAxes())


raw.Singleron$DF.classifications_0.25_0.22_1336 <- factor(raw.Singleron$DF.classifications_0.25_0.22_1336)

raw.Singleron@meta.data %>% 
  ggplot(aes(x=orig.ident, fill=DF.classifications_0.25_0.22_1336)) + 
  geom_bar(position = "stack") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")+
  scale_fill_npg()

raw.Singleron.singlet <- raw.Singleron[, raw.Singleron@meta.data[, DF.name] == "Singlet"]
raw.Singleron.singlet

## Calculate QC----
# Generating quality metrics

# mitochondrial ratio
raw.Singleron.singlet[["percent.mt"]] <- PercentageFeatureSet(object = raw.Singleron.singlet, pattern = "^MT-", assay = "RNA")
raw.Singleron.singlet$percent.mt <- raw.Singleron.singlet@meta.data$percent.mt / 100

# ribosome ratio
raw.Singleron.singlet[["percent.ribo"]] <- PercentageFeatureSet(object = raw.Singleron.singlet, pattern = "^RP[SL]", assay = "RNA")
raw.Singleron.singlet$percent.ribo <- raw.Singleron.singlet@meta.data$percent.ribo / 100

# hemoglobin ratio
raw.Singleron.singlet[["percent.hb"]] <- PercentageFeatureSet(object = raw.Singleron.singlet, pattern = "^HB[^(P)]", assay = "RNA")
raw.Singleron.singlet$percent.ribo <- raw.Singleron.singlet@meta.data$percent.ribo / 100

# number of genes detected per UMI
raw.Singleron.singlet$log10GenesPerUMI <- log10(raw.Singleron.singlet$nFeature_RNA) / log10(raw.Singleron.singlet$nCount_RNA)

# Assessing the quality metrics----

# UMI counts (transcripts) per cell
# Visualize the number UMIs/transcripts per cell    
raw.Singleron.singlet@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Genes detected per cell 
# Visualize the distribution of genes detected per cell via histogram
raw.Singleron.singlet@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 200)

# UMIs vs. genes detected
raw.Singleron.singlet@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 200) +
  facet_wrap(~orig.ident)

# Complexity
raw.Singleron.singlet@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# Mitochondrial counts ratio

raw.Singleron.singlet@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)


# Ribosome counts ratio

raw.Singleron.singlet@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent.ribo, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.01)

# Hemoglobin counts ratio

raw.Singleron.singlet@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent.hb, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.25)


## Filtering----
# Cell-level filtering

# * nCount > 500  
# * nFeature > 200 (p2-ERM 589 cells, preserve the data)  
# * log10GenesPerUMI > 0.8  
# * percent.mt < 0.2  
# * percent.ribo < 0.01  
# * percent.hb < 0.25  


filteredData <- subset(x = raw.Singleron.singlet, 
                       subset= (nCount_RNA >= 500) & 
                         (nFeature_RNA >= 200) & 
                         (log10GenesPerUMI > 0.80) & 
                         (percent.mt < 0.20) & 
                         (percent.ribo < 0.01) & 
                         (percent.hb < 0.25))
filteredData

# Gene-level filtering
## Extract counts
counts <- GetAssayData(object = filteredData, slot = "counts",assay = "RNA")

## Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

## Sums all TRUE values and returns TRUE if more than 3 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 3

## Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

## Reassign to filtered Seurat object
clean.Singleron.singlet <- CreateSeuratObject(filtered_counts, meta.data = filteredData@meta.data)

clean.Singleron.singlet
