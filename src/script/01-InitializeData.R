## Library needed packages ----
library(Seurat)
library(dplyr)

## Read in and Initialize Seurat Object ----
# Singleron
raw.p1ERM <- read.table("data/rawMatix/patient_1/PDR-ERM-0426_matrix.tsv",header = T) %>% CreateSeuratObject(project = "p1-ERM")
raw.p1Blood <- read.table("data/rawMatix/patient_1/PDR-Blood-0426_matrix.tsv", header = T) %>% CreateSeuratObject(project = "p1-Blood")

raw.p2ERM <- read.table("data/rawMatix/patient_2/PDR-ERM_matrix.tsv", header = T) %>% CreateSeuratObject(project = "p2-ERM")
raw.p2Blood <- read.table("data/rawMatix/patient_2/PDR-Blood_matrix.tsv", header = T) %>% CreateSeuratObject(project = "p2-Blood")

raw.p3ERM <- read.table("data/rawMatix/patient_3/PDR-ERM-0609_matrix.tsv.gz",header = T) %>% CreateSeuratObject(project = "p3-ERM")
raw.p3Blood <- read.table("data/rawMatix/patient_3/PDR-blood-0609_matrix.tsv.gz",header = T) %>% CreateSeuratObject(project = "p3-Blood")

raw.p4ERM <- read.table("data/rawMatix/patient_4/PDR-ERM-210630_matrix.tsv.gz",header = T) %>% CreateSeuratObject(project = "p4-ERM")
raw.p4Blood <- read.table("data/rawMatix/patient_4/PDR-B-210630_matrix.tsv.gz",header = T) %>% CreateSeuratObject(project = "p4-Blood")

raw.ctl1.Singleron <- read.table("data/rawMatix/control_PBMC/mid-aged/DZ0804_matrix.tsv.gz",header = T) %>% CreateSeuratObject(project = "p4-ERM")

# 10X
raw.ctl2.10X <- Read10X(data.dir = "data/rawMatix/control_PBMC/mid-aged/GSM5335490_C1/C1/") %>% CreateSeuratObject(project = "10X_ctr2")
raw.ctl3.10X <- Read10X(data.dir = "data/rawMatix/control_PBMC/mid-aged/GSM5335491_C2/C2/") %>% CreateSeuratObject(project = "10X_ctr3")
raw.ctl4.10X <- Read10X(data.dir = "data/rawMatix/control_PBMC/young/GSM4905216_YYZ/YYZ/") %>% CreateSeuratObject(project = "10X_ctr4")
raw.ctl5.10X <- Read10X(data.dir = "data/rawMatix/control_PBMC/young/GSM4905217_YTK/YTK/") %>% CreateSeuratObject(project = "10X_ctr5")
raw.ctl6.10X <- Read10X(data.dir = "data/rawMatix/control_PBMC/young/GSM4954813_C-1/") %>% CreateSeuratObject(project = "10X_ctr6")


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


