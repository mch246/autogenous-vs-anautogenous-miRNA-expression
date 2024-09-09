# This script runs DESeq2 (specifically starting with a count matrix file)

##Install DESeq2
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install('EnhancedVolcano')



#Load Libraries
library(DESeq2); library(ggplot2); library(EnhancedVolcano); library(tidyverse)


# Set working directory to source file location
setwd('/Users/mch284/Documents/Armbruster Lab/Auto_biting/auto_miRNA/')
# Read in the count data matrix
countData <- read.csv("/Users/mch284/Documents/Armbruster Lab/Auto_biting/auto_miRNA/AUTO_counts_miRNA_DESeq_set.csv", header = TRUE)
head(countData)

# IF removing a column--I removed outlier MAN_15H due to low read count from input csv 
countData <- countData %>%
  select(-columnID)

#Create the metaData table (this could alternatively be made externally and read in)
sampleNames <- colnames(countData)[2:8] #Here I am calling the column names of the Count matrix, minus the first column which is the name for the miRNA annotations
sampleConditions <- substr(sampleNames, 1, 2) #to get conditions I'm pulling the second letter from the sample names, which is either AU(Auto) or MA (MAN)

metaData <- data.frame(sampleName = sampleNames,
                        condition = sampleConditions)

# Make the DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = metaData, 
                              design= ~ condition, tidy = TRUE)

dds


#DESeq recommends a pre-filtering step to reduce memory size and increase speed. They suggest keeping only rows which have 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Relevel the condition for what to compare to. Likely not important in our case with two conditions, but you would want everything compared to the control. (Default is first condition alphabetically)
dds$condition <- relevel(dds$condition, ref = "AU")

# A bit of quality control
# Look at the distribution of the count data across the samples
librarySizes <- colSums(counts(dds))

barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, ylim = c(0,1.6e+07),
        main="Barplot of library sizes")

logcounts <- log2(counts(dds) + 1)
head(logcounts)

#Is there any difference between per gene counts for each of the sample groups?
statusCol <- as.numeric(factor(dds$condition)) + 1  # make a colour vector

boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=statusCol)

#Adding median log counts
abline(h=median(as.matrix(logcounts)), col="blue")

# Transform normalized counts using the rlog function
rld <- rlog(dds, blind=TRUE)
plotPCA(rld, intgroup="condition")+
  scale_color_manual(values=c("#00BFC4","#F8766D"))

plotPCA(rld)+
  scale_color_manual(values=c("#29788e", "#b93555"), 
                     name="",
                     labels=c("AUTO","NON-AUTO-LAB"))+


#try adding labels
plotPCA(rld, intgroup="condition") + geom_text(aes(label=sampleNames), vjust=2)

########################## Quality control then PCA Visualization ##########################
cds <- estimateSizeFactors(dds)
cds <- estimateDispersions(cds)
vsd = varianceStabilizingTransformation(cds, blind=TRUE)
theme_set(theme_bw())
meanSdPlot(assay(vsd))
plotDispEsts(cds)
nudge <- position_nudge(y = 10)
PCA_data <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE,ntop=500)
View(PCA_data)
percentVar <- round(100 * attr(PCA_data, "percentVar"))
nudge <- position_nudge(y = 4)
miRNA_pca_plot <- ggplot(PCA_data, aes(x = PC1, y = PC2, color = condition, position_nudge(y=10))) +
  geom_point(aes(shape = condition, color = condition),size=3)+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = c("#046C9A", "#A42820"),
                     labels=c("AUTO","NON-AUTO-LAB"))+ 
  scale_shape_manual(values = c(15,17),
                     labels=c("AUTO","NON-AUTO-LAB")) +
  coord_fixed()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.caption = element_text(face = "italic"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.box.just = "center",
        legend.box.background = element_blank(),
        plot.margin = margin(5.5, 30, 5.5, 5.5, "points"),
        legend.margin = margin(10, 10, 10, 10))+
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=12),
        legend.title= element_text(size=12))+
  labs(color = "Population")+
  guides(color = guide_legend(title = "Population", order = 1, ncol = 2,override.aes = list(size=4)),
         shape = guide_legend(title = "Population", order = 1, ncol = 2))
miRNA_pca_plot
ggsave("miRNA_pca_plot.png",plot=miRNA_pca_plot,dpi=600,units='in',width=6,height=5)


# PCA plot with other PCs besides 1 & 2
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))

view(sampleNames)

pca_df <- cbind(metaData, pca$x)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
ggplot(pca_df) + 
  geom_point(aes(x=PC6, y=PC7, color = condition), size = 3) +
  xlab(paste0("PC6: ",round(percentVar[6] * 100),"% variance")) +
  ylab(paste0("PC7: ",round(percentVar[7] * 100),"% variance"))

# look at loadings of PCA
head(pca$rotation)

# Differential Expression Analysis
#Run DESeq
dds <- DESeq(dds)
res <- results(dds)
res

# Apparently, it is common to shrink the log fold change estimate for better visualization and ranking of genes. Often, lowly expressed genes tend to have relatively high levels of variability so shrinking can reduce this.
resultsNames(dds)
res_LFC <- lfcShrink(dds, coef="condition_MA_vs_AU", type="apeglm")
res_LFC

# Volcano plot of the lfcShrink data. For various setting try: `browseVignettes("EnhancedVolcano")`
# Note for selecting fold change cutoff: log2foldchange 0.58 is equal to a 1.5 fold change
EV_AUTO_miRNA<-EnhancedVolcano(res_LFC,
                lab = rownames(res_LFC),
                x = 'log2FoldChange',
                y = 'padj',
                ylab = "-Log10(p-adjusted)",
                selectLab = NA,
                title = NULL,
                subtitle = NULL,
                xlim = c(-2.5, 2.5),
                ylim = c(0,5),
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 2.0,
                labSize = 5.0,
                caption = NULL)+
  title = "",
subtitle="")+
  theme(plot.caption = element_text(face = "italic"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.box.just = "center",
        legend.box.background = element_blank(),
        plot.margin = margin(5.5, 30, 5.5, 5.5, "points"),
        legend.margin = margin(10, 10, 10, 10))+
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=12),
        legend.title= element_text(size=12))
EV_AUTO_miRNA
ggsave(
 filename="EV_AUTO_miRNA.tiff",
  plot = EV_AUTO_miRNA,
  device = tiff,
  dpi=600,
 units='in',width=6,height=5)

##PCA for inset
miRNA_pca_plot_inset <- ggplot(PCA_data, aes(x = PC1, y = PC2, color = condition, position_nudge(y=10))) +
  geom_point(aes(shape = condition, color = condition),size=3)+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = c("#046C9A", "#A42820"),
                     labels=c("AUTO","NON-AUTO-LAB"))+ 
  scale_shape_manual(values = c(15,17),
                     labels=c("AUTO","NON-AUTO-LAB")) +
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(tag = "B") +
  theme(plot.tag.position=c(.01,.91),
        plot.tag=element_text(face="bold",size = 16))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(plot.title = element_text(hjust = 0.3))+
  theme(plot.caption = element_text(face = "italic"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.box.just = "center",
        legend.box.background = element_blank(),
        plot.margin = margin(5.5, 30, 5.5, 5.5, "points"),
        legend.margin = margin(10, 10, 10, 10))+
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=12),
        legend.title= element_text(size=12))+
  labs(color = "Population")+
  guides(color = guide_legend(title = "Population", order = 1, ncol = 2,override.aes = list(size=4)),
         shape = guide_legend(title = "Population", order = 1, ncol = 2))

miRNA_pca_plot_inset


###Enhanced volcano for inset
EV_AUTO_miRNA_inset<-EnhancedVolcano(res_LFC,
                               lab = rownames(res_LFC),
                               legendPosition = "none",
                               x = 'log2FoldChange',
                               y = 'padj',
                               ylab = "-Log10(p-adjusted)",
                               selectLab = NA,
                               title = NULL,
                               subtitle = NULL,
                               xlim = c(-5.5, 5.5),
                               ylim = c(0,20),
                               pCutoff = 0.05,
                               FCcutoff = 1,
                               pointSize = 2.0,
                               labSize = 5.0,
                               caption = NULL,
gridlines.major = FALSE,
gridlines.minor = FALSE)+
  labs(tag = "A") +
  theme(plot.tag.position=c(.14,.98),
        plot.tag=element_text(face="bold",size = 18))

EV_AUTO_miRNA_inset

miRNA_pca_plot_inset
AUTO_miRNA_inset_fig<-EV_AUTO_miRNA_inset+annotation_custom(ggplotGrob(miRNA_pca_plot_inset), xmin = -0.5, xmax = 3.6, 
                                            ymin = 3.1, ymax = 7.8)

AUTO_miRNA_inset_fig
ggsave(
  filename="AUTO_miRNA_inset_fig.tiff",
  plot = AUTO_miRNA_inset_fig,
  device = tiff,
  dpi=600,
  limitsize = TRUE,
  bg = NULL)
# Full result print out
full_res <- res_LFC %>%
  data.frame() %>%
  rownames_to_column(var="miRNA") %>% 
  as.data.frame()

normalized_counts <- counts(dds, normalized=TRUE)
view(normalized_counts)
write.csv((normalized_counts), file="/Users/mch284/Documents/Armbruster Lab/Auto_biting/auto_miRNA/AUTO_MAN_miRNA_DEseq2_normalized_counts.csv" , row.names=T)
# Write out a table of these results
#write.csv(dplyr::select(full_res, miRNA, log2FoldChange, padj), 
          file="/Users/mch284/Documents/Armbruster Lab/Auto_biting/auto_miRNA/MANvAUTO_miRNA_LFCshrink_padj_FULL_re.csv", row.names = F)
