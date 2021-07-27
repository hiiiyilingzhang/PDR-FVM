## Library needed packages ----
library(Seurat)
library(dplyr)

## Read in patient 1&2 data and Initialize Seurat Object ----
raw.p1ERM <- read.table("data/rawMatix/patient_1/PDR-ERM-0426_matrix.tsv",header = T) %>% CreateSeuratObject(project = "p1-ERM")
raw.p1Blood <- read.table("data/rawMatix/patient_1/PDR-Blood-0426_matrix.tsv", header = T) %>% CreateSeuratObject(project = "p1-Blood")

raw.p2ERM <- read.table("data/rawMatix/patient_2/PDR-ERM_matrix.tsv", header = T) %>% CreateSeuratObject(project = "p2-ERM")
raw.p2Blood <- read.table("data/rawMatix/patient_2/PDR-Blood_matrix.tsv", header = T) %>% CreateSeuratObject(project = "p2-Blood")

## SCTransform ----
raw.p1ERM <- SCTransform(raw.p1ERM, verbose = FALSE)
raw.p1Blood <- SCTransform(raw.p1Blood, verbose = FALSE)
raw.p2ERM <- SCTransform(raw.p2ERM, verbose = FALSE)
raw.p2Blood <- SCTransform(raw.p2Blood, verbose = FALSE)

### Integration ----
features <- SelectIntegrationFeatures(object.list = list(raw.p1ERM,raw.p1Blood,raw.p2ERM,raw.p2Blood))
all.anchors <- FindIntegrationAnchors(object.list = list(raw.p1ERM,raw.p1Blood,raw.p2ERM,raw.p2Blood), anchor.features = features)
integData <- IntegrateData(anchorset = all.anchors)
DefaultAssay(integData) <- "integrated"

# integData
# integData@assays$integrated[1:5,1:5]
# dim(integData[["RNA"]]@counts)
# dim(integData@assays$integrated@data)
# integData@meta.data[1:5,1:5]
# table(integData@meta.data["orig.ident"])


