library(SingleR)
library(celldex)
library(scater)

# https://github.com/dviraran/SingleR
# https://bioconductor.org/packages/3.13/data/experiment/vignettes/celldex/inst/doc/userguide.html


hpca.se <- HumanPrimaryCellAtlasData()
save(hpca.se,file="data/annotation/hpca-se.RData")

common_hpca <- intersect(rownames(integrated), rownames(hpca.se))
hpca.se <- hpca.se[common_hpca,]

dat_forhpca.se <- as.data.frame(integrated[["RNA"]]@data)

dat_forhpca.se <- SummarizedExperiment(assays=list(counts=dat_forhpca.se[common_hpca,])) %>% logNormCounts()

dat_forhpca.se@metadata <- integrated@meta.data

pred.main.hpca <- SingleR(test = dat_forhpca.se, ref = hpca.se,
                          labels = hpca.se$label.main,de.method="wilcox",
                          clusters = dat_forhpca.se$seurat_clusters)
table(pred.main.hpca$pruned.labels)

plotScoreHeatmap(pred.main.hpca)

result_main_hpca <- as.data.frame(pred.main.hpca$labels)
result_main_hpca$CB <- rownames(pred.main.hpca)
colnames(result_main_hpca) <- c('HPCA_Main', 'seurat_clusters')
write.table(result_main_hpca, file = "output/Annotation/HPCA_Main.txt", sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)


integrated@meta.data <- merge(integrated@meta.data,result_main_hpca,by="seurat_clusters")
rownames(set1.obj@meta.data) <- set1.obj@meta.data$CB