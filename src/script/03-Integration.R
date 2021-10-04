
# Harmony 
# https://github.com/immunogenomics/harmony  
# Overview of Harmony algorithm: PCA embeds cells into a space with reduced dimensionality. Harmony accepts the cell coordinates in this reduced space and runs an iterative algorithm to adjust for dataset specific effects.  

# Cell cycle scoring check ----

cell_cycle_genes <- read.csv("../data/cell-cycle/Homo_sapiens.csv")
head(cell_cycle_genes)

library(AnnotationHub)
library(ensembldb)

# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)


# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")


## SCTransform ----

mergedData <- merge(x = clean10X, y = cleanSingleron)

# Perform the cell cycle scoring and sctransform on all samples.
mergedData <- NormalizeData(mergedData, verbose = F) %>% 
  CellCycleScoring(g2m.features=g2m_genes, s.features=s_genes)

VlnPlot(mergedData, features = "S.Score", group.by = "orig.ident", pt.size = 0.1)
VlnPlot(mergedData, features = "G2M.Score", group.by = "orig.ident", pt.size = 0.1)

mergedData <- SCTransform(mergedData, vars.to.regress = c("percent.mt"), variable.features.n = 3000, verbose = F, seed.use = 27)
mergedData <- RunPCA(mergedData, npcs = 30, verbose = FALSE)

## Integration ----

sampledata <- read.table("../data/tmp.data/sampleData.txt", header = T, sep = "\t")
head(sampledata)

metaNow <- mergedData@meta.data
metaNow$barcode <- rownames(metaNow)

usedCol <- metaNow[,c("barcode","orig.ident")]
usedCol <- merge(usedCol, sampledata, by = "orig.ident")
rownames(usedCol) <- usedCol$barcode

mergedData <- AddMetaData(object = mergedData, metadata = usedCol)

integrated <- RunHarmony(mergedData, "Platform", assay.use = "SCT", plot_convergence = TRUE)
head(integrated@meta.data)


## Check integration ----
# Run tSNE
integrated <- RunTSNE(object = integrated, dims.use = 25,reduction = "harmony")

# Plot tSNE
DimPlot(object = integrated, group.by = "Status", reduction = "tsne")
FeaturePlot(object = integrated, features = c("CCR5"), cols = c("grey", "blue"), reduction = "tsne", split.by = "Status")
FeaturePlot(object = integrated, features = c("ENO2"), cols = c("grey", "blue"), reduction = "tsne",split.by = "Status")
