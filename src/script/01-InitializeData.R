## Library needed packages ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

## Read in and Initialize Seurat Object ----
# Singleron
raw.p1ERM <- read.table("data/rawMatix/patient_1/PDR-ERM-0426_matrix.tsv",header = T) %>% CreateSeuratObject(project = "p1-ERM")
raw.p1Blood <- read.table("data/rawMatix/patient_1/PDR-Blood-0426_matrix.tsv", header = T) %>% CreateSeuratObject(project = "p1-Blood")
raw.p2ERM <- read.table("data/rawMatix/patient_2/PDR-ERM_matrix.tsv", header = T) %>% CreateSeuratObject(project = "p2-ERM")
raw.p2Blood <- read.table("data/rawMatix/patient_2/PDR-Blood_matrix.tsv", header = T) %>% CreateSeuratObject(project = "p2-Blood")
raw.p3ERM <- read.table("data/rawMatix/patient_3/PDR-ERM-0609_matrix.tsv",header = T) %>% CreateSeuratObject(project = "p3-ERM")
raw.p3Blood <- read.table("data/rawMatix/patient_3/PDR-blood-0609_matrix.tsv",header = T) %>% CreateSeuratObject(project = "p3-Blood")
raw.p4ERM <- read.table("data/rawMatix/patient_4/PDR-ERM-210630_matrix.tsv",header = T) %>% CreateSeuratObject(project = "p4-ERM")
raw.p4Blood <- read.table("data/rawMatix/patient_4/PDR-B-210630_matrix.tsv.gz",header = T) %>% CreateSeuratObject(project = "p4-Blood")
raw.ctl1.Singleron <- read.table("data/rawMatix/control_PBMC/mid-aged/DZ0804_matrix.tsv.gz",header = T) %>% CreateSeuratObject()

# 10X
raw.ctl2.10X <- Read10X(data.dir = "data/rawMatix/control_PBMC/mid-aged/GSM5335490_C1/C1/") %>% CreateSeuratObject(project = "10X_ctr2")
raw.ctl3.10X <- Read10X(data.dir = "data/rawMatix/control_PBMC/mid-aged/GSM5335491_C2/C2/") %>% CreateSeuratObject(project = "10X_ctr3")
raw.ctl4.10X <- Read10X(data.dir = "data/rawMatix/control_PBMC/young/GSM4905216_YYZ/YYZ/") %>% CreateSeuratObject(project = "10X_ctr4")
raw.ctl5.10X <- Read10X(data.dir = "data/rawMatix/control_PBMC/young/GSM4905217_YTK/YTK/") %>% CreateSeuratObject(project = "10X_ctr5")
raw.ctl6.10X <- Read10X(data.dir = "data/rawMatix/control_PBMC/young/GSM4954813_C-1/") %>% CreateSeuratObject(project = "10X_ctr6")

## Merge ----
# Singleron
raw.Singleron <- merge(x = raw.p1ERM, y = c(raw.p1Blood, raw.p2ERM, raw.p2Blood,raw.p3ERM, raw.p3Blood, raw.p4ERM, raw.p4Blood, raw.ctl1.Singleron))

# 10X
raw.10X <- merge(x = raw.ctl2.10X, y = c(raw.ctl3.10X,raw.ctl4.10X, raw.ctl5.10X, raw.ctl6.10X))

## Check Data ----

raw.Singleron@assays$RNA[1:5,1:5]
raw.Singleron@meta.data[1:5,]

raw.Singleron@meta.data %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

raw.10X@meta.data %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

