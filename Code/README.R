# Pakages
library(dplyr)
library(tidyverse)
library(GEOquery)
library(Biobase)
library(ggplot2)
library(oligo)
library(DESeq2)
library(pheatmap)
library(WGCNA)
library(gridExtra)
library(CorLevelPlot)

# Downloading raw data
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE203206", "file=GSE203206_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
raw_data <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

# Downloading Metadata
meta_data <- getGEO('GSE203206', GSEMatrix = TRUE)

# Downloading Annotation data
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID

# Phenotype Data 
phenoData <- pData(phenoData(meta_data[[1]]))
View(phenoData)

## Extract
phenoData <- phenoData[,c(1,2,40:42)]

# Preparation of data 
raw_data[1:10,1:10]

# Quality control
## Detection of Outliers 
gsg <- goodSamplesGenes(t(raw_data)) # goodSamplesGenes (gsg)
summary(gsg)

#
gsg$allOK
table(gsg$goodGenes) #
table(gsg$goodSamples) #
gsg$allOK # True: all genes passed the quality control; False: some genes failed the QC and should be removed
table(gsg$goodGenes) # True: the genes can be used for further analysis; False: genes need to be removed  
table(gsg$goodSamples) # True: samples are usable for the analysis; False: some samples should be removed

### Detection of outliers by clustering 
htree <- hclust(dist(t(raw_data)), method = "average")
plot(htree)

### PCA 
pca <- prcomp(t(raw_data))
pca_data <- pca$x
pca_var <- pca$sdev^2
pca_var_percentage <- round(pca_var/sum(pca_var)*100, digits = 2)
pca_data <- as.data.frame(pca_data)

pca_data |> 
  ggplot(aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca_data)) +
  labs(x = paste0("PC1: ", pca_var_percentage[1], "%"), 
       y = paste0("PC2: ", pca_var_percentage[2], "%"))

# Elimination of outliers
outliers <-  c('GSM6164667','GSM6164668','GSM6164676','GSM6164674','GSM6164669','GSM6164692','GSM6164681')
raw_data_subset <- raw_data[,!(colnames(raw_data) %in% outliers)]

# Data normalization 
## Elimination of outliers from the phenotype data 
colData <- phenoData %>%
  filter(!row.names(.) %in% outliers)
## Data manipulation 
names(colData)
names(colData) <- gsub(":ch1", "", names(colData)) # removal of special characters
names(colData) # verification
names(colData) <- gsub("\\s", "_", names(colData))

## Making the column names and row names identical 
all(rownames(colData) %in% colnames(raw_data_subset))
all(rownames(colData) == colnames(raw_data_subset))

## Create a DESeq Matrix 
dds <- DESeqDataSetFromMatrix(countData = raw_data_subset,
                              colData = colData,
                              design = ~1)

### filtering low-count genes in more than 75% of the samples 
dds_filter <- dds[rowSums(counts(dds) >= 15) >= 24,]  
nrow(dds_filter) # number of genes left: 18401

### Normalization
dds_norm <- vst(dds_filter)

### Normalised counts
### Normalized counts
norm_count <- assay(dds_norm) |> 
  t()

# Network construction  
## Tresholding powers 
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

## Calling the network topology function
pst <- pickSoftThreshold(norm_count, 
                         powerVector = power, 
                         networkType = "signed",
                         verbose = 5)

## Storing indices data
pst_data <- pst$fitIndices

## Choosing power for the analysis 
### Power vs Scale free topology model fit
a1 <- ggplot(pst_data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.85, color = "red") +
  labs(x = "Power", y = "Scale free topology model fit, signed R^2") +
  theme_classic()

### Power vs Mean Conectivity
a2 <- pst_data |> 
  ggplot(aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = "Power", y = "Mean Conectivity") +
  theme_classic()

## Power threshold
grid.arrange(a1, a2, nrow = 2)

## Convert metric to numeric 
norm_count[] <- sapply(norm_count, as.numeric)
soft_power <- 26
temp_cor <- cor
cor <- WGCNA::cor

## Memory estimate with respect to blocksite
bwnet <- blockwiseModules(norm_count,
                          maxBlockSize = 19000,
                          TOMType = "signed", 
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

cor <- temp_cor

## Module Eigengens
module_eigengens <- bwnet$MEs
head(module_eigengens)

##
table(bwnet$colors)

## Dendogram analysis 
plotDendroAndColors(bwnet$dendrograms[[1]],
                    cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)

# Module trait association 
## Create traits fields, binarise categorical variables
traits <- colData |> 
  mutate(disease_state_bin = ifelse(grepl("Alzheimer's disease", disease_state), 1, 0)) |> 
  dplyr::select(6)

## Binarise categorical variables
colData$group <-  factor(colData$group,
                         levels = c("control", "Early-onset sporadic AD", "Late-onset sporadic AD"))

group_out <- binarizeCategoricalColumns(colData$group,
                                        includePairwise = FALSE, 
                                        includeLevelVsAll = TRUE,
                                        minCount = 1)

traits <- cbind(traits, group_out) #*

## Define number of genes 
nSamples <- nrow(norm_count)
nGenes <-  ncol(norm_count)

## Module-trait correlation 
module_trait_corr <- cor(module_eigengens, traits, use = "p")

## Module-trait correlation p values
module_trait_corr_pvalues <- corPvalueStudent(module_trait_corr, nSamples)

# Visualization of module-trait association as a heatmap

heatmap_data <- merge(module_eigengens, traits, by = "row.names")
head(heatmap_data)

heatmap_data <- heatmap_data |> 
  column_to_rownames(var = "Row.names")

colnames(heatmap_data)[12:14] <- c("Disease vs Control", "EOAD vs ALL", "LAD vs ALL") # EOAD: Early-onset sporadic AD; LAD: Late-onset sporadic AD 

CorLevelPlot(heatmap_data,
             x = names(heatmap_data)[12:14],
             y = names(heatmap_data)[1:11],
             col = c("blue3", "skyblue", "white", "pink", "red"),
             main = "Module-Traits association"
             )

# Gene mapping 
module_genes_mapping <- merge(as.data.frame(bwnet$colors), annot, by = 0) # annotation data merged

## Gene extraction from each significance module trait based on the heatmap 
gene_MEblue <- module_genes_mapping |> 
  filter(`bwnet$colors` == "blue") |> ## "blue" can be change for any other variable depending on the results
  dplyr::select(Symbol)

# Module membership measurements (mod_mem_mea)
mod_mem_mea <- cor(module_eigengens, norm_count, use = "p")

## Module membership measurements and associated p-values (mmm_pvals)
mmm_pvals <- corPvalueStudent(mod_mem_mea, nSamples)

# Gene significance analysis (gsa) and associated p-values (pvals)
gsa <- cor(norm_count, traits$`data.Early-onset sporadic AD.vs.all`, use = "p") 

## Associated p-values
gsa_pvals <- corPvalueStudent(gsa, nSamples)

results_gsa_pvals <- merge(as.data.frame(gsa_pvals), annot, by = 0)

results_gsa_pvals |> 
  dplyr::select(V1, GeneID, Symbol)

