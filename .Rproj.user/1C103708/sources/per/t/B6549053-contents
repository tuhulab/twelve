# install.packages("BiocManager")
# BiocManager::install("GEOquery")

library(GEOquery)
gse_121212 <- getGEO("GSE121212", GSEMatrix =TRUE, AnnotGPL=FALSE)
gse_121212 <- gse_121212[[1]]
exprs(gse_121212)
