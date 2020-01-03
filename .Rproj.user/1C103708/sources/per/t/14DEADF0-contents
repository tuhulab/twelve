####packages####
library(GEOquery)


####Download raw data####
geofile <- getGEO(GEO = "GSE68592")




#####
################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)
library(tidyverse)

# load series and platform data from GEO

gset <- getGEO("GSE68592", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) 
  idx <- grep("GPL15084", attr(gset, "names")) 
else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("XXXXXXXXXX00000XXXXXXXXXX00000XXXXXXXXXXXXXXXXXXXX",
               "XXXXX11111XXXXXXXXXXXXXXXXXXXXXXXX11111XXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXX")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GB_ACC","SEQUENCE","SPOT_ID"))

featuredata <- gset[["GSE68592_series_matrix.txt.gz"]]@featureData@data %>% as.tibble()

tT %>% as.tibble() %>% filter(logFC > 1,
                              adj.P.Val < .05, GB_ACC != "") %>% 
  left_join(featuredata) %>% 
  select(ID, adj.P.Val,logFC,SEQUENCE, GENE_SYMBOL, GENE_NAME, DESCRIPTION) %>% 
  write_csv("90d_up.csv")

tT %>% as.tibble() %>% filter(logFC < -1,
                              adj.P.Val < .05, GB_ACC != "") %>% 
  left_join(featuredata) %>% 
  select(ID, adj.P.Val,logFC,SEQUENCE, GENE_SYMBOL, GENE_NAME, DESCRIPTION) %>% 
  write_csv("90d_down.csv")

