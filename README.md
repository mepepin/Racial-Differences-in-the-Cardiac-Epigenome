---
title: "Supplemental Analysis Workflow - Epigenenomics of Race in Heart Failure"
author: "Mark E. Pepin"
date: "01/19/2021"
output:
  html_document:
    keep_md: yes
    code_folding: hide
    toc: yes
    toc_float: yes
  pdf_document:
    toc: yes
  word_document:
    toc: yes
geometry: margin=1 in
header-includes:
- \usepackage{booktabs}
- \usepackage{longtable}
- \usepackage{array}
- \usepackage{multirow}
- \usepackage[table]{xcolor}
- \usepackage{wrapfig}
- \usepackage{float}
- \usepackage{colortbl}
- \usepackage{pdflscape}
- \usepackage{tabu}
- \usepackage{threeparttable}
mainfont: Times
fontsize: 10pt
always_allow_html: yes
---



**Code Authors**: Mark E. Pepin
**Contact**: pepinme@uab.edu
**Institution**: University of Alabama at Birmingham  
**Location**: 542 Biomedical Research Building 2, Birmingham, AL 35294  

# Analysis Parameters

## Parameters

Define the parameters used, along with the conditions required for the current analysis. This segment must be modified for each analysis performed.


```r
##Set the experimental conditions [MUST DO THIS MANUALLY]
DIABETES=c("ND", "T2D") #Should define the reference group first
ISCHEMIA=c("NICM", "ICM")
RACE=c("Caucasian_American", "African_American") #, "Asian_American", "Eastern_African"
ISO=c("NONE", "DOB")
STATISTIC = 0.05 #P statistic threshold used in this combination.

VARIABLE = RACE

COMPARISON= paste0(VARIABLE[2], ".vs.", VARIABLE[1])

# Candidate Gene Selection (RNA-sequencing) EDIT THIS LIST BASED ON INTERESTS.
GENES=c("DNMT3A", "GADD45B", "HADH", "HK1", "HK2", "PFKM", "ALDOA", "GAPDH", "ENO3", "TP53", "MYC", "MDH1", "NDUFA4")
VAR1="Race"
VAR2="Diabetes"
# Single Bar Graph
GENES=c("DNMT3A", "MYH6", "SDHC", "GADD45B")
# Pathway Gene Screening
GENE.Pathway<-read.csv("../2_Input/Gene.Sets/Oxphos.KEGG.csv", header = F) #Can alter the gene pathway (just change path)
colnames(GENE.Pathway)<-"external_gene_name"

# Create Output Folder Structure
ifelse(!dir.exists(file.path(paste0("../3_Output/"))), dir.create(file.path(paste0("../3_Output/"))), FALSE)
```

```
## [1] FALSE
```

```r
ifelse(!dir.exists(file.path(paste0("../3_Output/1_RNA/"))), dir.create(file.path(paste0("../3_Output/1_RNA/"))), FALSE)
```

```
## [1] FALSE
```

```r
ifelse(!dir.exists(file.path(paste0("../3_Output/2_Methyl/"))), dir.create(file.path(paste0("../3_Output/2_Methyl/"))), FALSE)
```

```
## [1] FALSE
```

```r
ifelse(!dir.exists(file.path(paste0("../3_Output/3_Combined/"))), dir.create(file.path(paste0("../3_Output/3_Combined/"))), FALSE)
```

```
## [1] FALSE
```

## Packages


```r
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, Hmisc, openxlsx, corrplot, RColorBrewer, kableExtra, ggplot2, gridExtra, ggpubr, ggsignif, DESeq2, data.table, GenomicFeatures, biomaRt, Haplin, pheatmap, calibrate, ggrepel, tidyr, gtools)
```

# Patient Information
## Correlation and Covariate Network Analysis

The first task in this project was to verify that the heart failure patients from which the cardiac tissues were obtained allowed for a simple comparison (i.e. not requiring multiple regression) by both racial and etiologic factors, as a means of examining the experimental subjects for confounding variables.


```r
library(openxlsx)
library(dplyr)
library(Hmisc)
library(corrplot)
library(RColorBrewer)
options(kableExtra.latex.load_packages = FALSE)
library(kableExtra)
Patient_Data <- openxlsx::read.xlsx("../2_Input/1_Patient/colData_CHAAMPS.FINAL.xlsx", sheet = "Matrix", rowNames = T)
Patient_Data<-as.data.frame(Patient_Data)
# Patient_Data<-dplyr::select(Patient_Data, -Syst.PA, -PCWP,-RA)
Patient_Data<-Patient_Data[complete.cases(Patient_Data),]

# Print the table for use
Patient_Data %>% kable(format="latex", 
                       align="c", 
                       booktabs=T, 
                       caption="Patient Characteristics") %>% 
  kable_styling(latex_options=c("striped", 
                                "condensed", 
                                "scale_down"))
```

```r
# Format for Correlation
Patient_corr<-data.matrix(Patient_Data)
cor.r<-cor(Patient_corr)
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "brown4"))(paletteLength)
p.mat<-cor.mtest(cor.r)$p
rownames(p.mat)<-rownames(cor.r)
colnames(p.mat)<-colnames(cor.r)
corrplot(cor.r, 
         order="original",
         type="full",
         method = "square",
         outline=FALSE,
         col = myColor,
         tl.cex=0.7,
         tl.col="black",
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         addgrid.col = NA)
```

![\label{correlation_patients}Correlation Matrix of Patient Characteristics. All data were obtained from electronic medical records and de-identified in accordance with the UAB Institutional Review Board requirements.](ComprehensiveAnalysis_2021_files/figure-html/correlation_patients-1.png)

```r
pdf(file=paste0("../3_Output/_Patient/Patient.Correlation.pdf"))
corrplot(cor.r, 
         order="original",
         type="full",
         method = "square",
         outline=FALSE,
         col = myColor,
         tl.cex=0.7,
         tl.col="black",
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         addgrid.col = NA)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
library(igraph)
#convert to an adjacency matrix
adj.mat<-p.mat
rownames(adj.mat)<-rownames(cor.r)
colnames(adj.mat)<-colnames(cor.r)
adj.mat[adj.mat>0.05]=NA
adj.mat[adj.mat<0.05]=1
adj.mat[is.na(adj.mat)]=0

patient.igraph<-graph_from_adjacency_matrix(adj.mat, diag = F)
patient.igraph<-simplify(patient.igraph)
plot(patient.igraph, edge.arrow.size=0.0, cex=0.5)
```

![\label{correlation_patients}Correlation Matrix of Patient Characteristics. All data were obtained from electronic medical records and de-identified in accordance with the UAB Institutional Review Board requirements.](ComprehensiveAnalysis_2021_files/figure-html/correlation_patients-2.png)

```r
# library("threejs")
# library("htmlwidgets")
# library("igraph")
# patient.igraph.js<-patient.igraph
# graph_attr(patient.igraph.js, "layout") <- NULL
# gjs <- graphjs(patient.igraph.js, main="Patient Network", bg="white", attraction=0.9, repulsion=0.8, opacity=0.9, edge.color ="gray50", vertex.color = "gray50", edge.width = 4)
# print(gjs)
# saveWidget(gjs, file="Media-Network-gjs.html")
# browseURL("Media-Network-gjs.html")
```

## Patient Variables - Bar Graphs

To examine the changes of individual variables, we used bar graphs based on timepoint and patient response to LVAD-induced unloading.


```r
#####################################################
# Bar Graphs of the PHI
######################################################
#Location where graphs will be saved
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
VAR2="Race"
PARAMS=c("Age", "Creatinine", "Blood.Glucose", "LV.EF")

ifelse(!dir.exists(file.path(paste0("../3_Output/_Patient/Candidates/"))), dir.create(file.path(paste0("../3_Output/_Patient/Candidates/"))), FALSE)
```

```
## [1] FALSE
```

```r
# Patient_Data<-openxlsx::read.xlsx("../2_Input/_Patient/Patients_Health.Info.xlsx", sheet = "Matrix_Pre.v.CON", rowNames = T)
# Patient_Data$Sample_ID<-rownames(Patient_Data)
# Patient_Data<-dplyr::select(Patient_Data, -Sex)
# Patient_ALL<-reshape(Patient_Data, varying = 14:29, sep = "_", direction = "long")
# Patient_ALL<-dplyr::rename(Patient_ALL, Timing=time)

colData_all<-openxlsx::read.xlsx("../2_Input/1_Patient/colData_CHAAMPS.FINAL.xlsx", sheet = "Summary", rowNames = T, startRow = 2)
graph_info<-subset(colData_all, Race %in% RACE)

## create a function that converts factors to numeric
asNumeric=function(x){as.numeric(as.character(x))}
factorsNumeric=function(d){modifyList(d, lapply(d[, sapply(d, is.factor)], asNumeric))}

##
graph_info<-factorsNumeric(graph_info)
# graph_info$Response<-factor(graph_info$Response, levels = c("CON", "NR", "R"))
graph_info$Ischemia<-factor(graph_info$Ischemia, levels = c("NICM", "ICM"))
graph_info$Diabetes<-factor(graph_info$Diabetes, levels = c("ND", "T2D"))
graph_info$HTN<-factor(graph_info$HTN, levels = c("N", "Y"))

groupsize<-graph_info %>% group_by(Race) %>% tally()

#Define function for mean and standard deviation for each group
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}
########################################################
plotlist = list()
for (i in seq_along(PARAMS)){
PARAM<-as.name(PARAMS[i]) #define a "name" for the variable
ds <- data_summary(graph_info, varname=PARAM, groupnames=c(VAR2))
ds<-dplyr::left_join(ds, groupsize)
ds<-dplyr::mutate(ds, upper=ds[,2]+(sd/sqrt(n-1)))
patient<-ggplot(ds, aes_string(x=VAR2, y=PARAMS[i], fill=VAR2)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes_string(ymin=PARAMS[i], ymax="upper"), width=.2,
                 position=position_dodge(.9))
patient_plot<-patient+labs(title=paste0("Patient Data - ", PARAMS[i]), x=VAR2, fill=VAR2)+
    theme_classic2(base_size = 10) +
    scale_fill_manual(values=c('black', "white", "#999999")) + 
    stat_compare_means(mapping = aes_string(x=VAR2, y=PARAMS[i]), data = graph_info,label = "p.signif")
pdf(file=paste0("../3_Output/_Patient/Candidates/", PARAMS[i], "_", VAR2, "_barplot.pdf"), width = 3.5, height = 2.5)
print(patient_plot)
dev.off()
plotlist[[i]] = patient_plot
}
t<-marrangeGrob(grobs = plotlist, legend, nrow=3, ncol=2)
ggsave(paste0("../3_Output/_Patient/Candidates/_barplots_", VAR2, ".pdf"), t, width = 6, height = 7)

##########################################################
#Boxplot
##########################################################

plotlist = list()
for (i in seq_along(PARAMS)){
PARAM<-as.name(PARAMS[i]) #define a "name" for the variable
ds <- data_summary(graph_info, varname=PARAM, groupnames=c(VAR2))
ds<-dplyr::left_join(ds, groupsize)
ds<-dplyr::mutate(ds, upper=ds[,2]+(sd/sqrt(n-1)))
g<-ggplot(graph_info, aes_string(x=VAR2, y=PARAMS[i], fill=VAR2)) + 
  geom_boxplot() + stat_compare_means(label = "p.signif")
g_plot<-g+labs(title=paste0("Patient Data - ", PARAMS[i]), x=VAR2)+
   theme_classic2(base_size = 10) +
   scale_fill_manual(values=c('black', "white", "#999999")) + ylim(NA, max(ds$upper+.25*(ds$upper)))
pdf(file=paste0("../3_Output/_Patient/Candidates/", PARAMS[i], "_", VAR2, "_boxplot.pdf"), width = 3.5, height = 2.5)
print(g_plot)
dev.off()
plotlist[[i]] = g_plot
}
t<-marrangeGrob(grobs = plotlist, legend, nrow=3, ncol=2)
ggsave(paste0("../3_Output/_Patient/Candidates/_boxplots_", VAR2, ".pdf"), t, width = 6, height = 7)
```

\pagebreak

# RNA-Sequencing Analysis

## RNA-Sequencing Alignment using STAR

Once the proper patient samples were selected for the analysis, RNA was isolated from the left ventricle endocardial tissue using the RNeasy Lipid Mini-Kit according to the manufacturer's instructions (Qiagen, Valencia, CA). High-throughput RNA sequencing was performed at the Heflin Genomics Core at the University of Alabama at Birmingham. Once sample read quality was checked (multiQC analysis), the paired-end fastq files were then aligned to the reference genome, which was created using Gencode human sequence (GRCh38.p10.genome.fa) and annotation (gencode.v27.chr_patch_hapl_scaff.annotation.gtf). STAR aligner is the current gold-standard for this, which we used for the current analysis. Before aligning each fastq file to the genome, an annotated reference genome must first be assembled. This was performed as follows (this was performed in Cheaha as `bash GenomeReference.sh':

STAR=../../Tools/STAR-2.5.3a/bin/Linux_x86_64/STAR

$STAR \
--runThreadN 12 \
--runMode genomeGenerate \
--genomeDir ./ \
--genomeFastaFiles /data/scratch/pepinme/huHrt/Input/Genome/GRCh38.p10.genome.fa \

Alignment of short reads to this annotated genome could then proceed, using the following SLURM batch script which was submitted to the UAB *Cheaha* compute cluster (See **Appendix**). This shell script contains the following STAR alignment run settings:

$STAR_RUN \
--genomeDir $GENOME_DIR \
--readFilesCommand zcat \
--readFilesIn $INPUT_DIR/fastq/${VAR}_R1_001.fastq.gz \
    $INPUT_DIR/fastq/${VAR}_R2_001.fastq.gz \
--sjdbGTFfile $GENOME_DIR/gencode.v27.chr_patch_hapl_scaff.annotation.gtf \
--sjdbOverhang 99 \
--quantMode GeneCounts \
--runThreadN 12 \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${RESULTS_DIR}/Alignment/${VAR}_

## Read Count Compiling

Before the DESeq2-based differential expression can be computed, the counts generated by STAR need to be compiled, since the .tab file contains count information based on forward, reverse, and combined reads. Therefore, we will take the fourth column in each table and merge them.


```r
Count.files <- list.files(path = "../1_Cheaha/counts/", pattern = "*_ReadsPerGene.out.tab", full.names = TRUE, all.files = TRUE)
Counts <- lapply(Count.files, read.table, skip = 4) #skip the first 4 rows, since these are summary data.
#Create a data.frame containing the raw counts
countData.raw <- as.data.frame(sapply(Counts, function(x) x[,4])) #selects only the 4th column as the raw counts.
#Generate Column names and Row names for the counts (remove the extra nonsense from the path names)
colnames <- gsub( "_ReadsPerGene[.]out[.]tab", "", Count.files)
colnames <- gsub( "[.][.]/1_Cheaha/counts//", "", colnames)
colnames(countData.raw) <- colnames
row.names(countData.raw) <- Counts[[1]][,1]
```

## Data Pre-Processing

After alignment of the fastq files to the annotated genome assembly (hg38), the first step in the analysis is to consolidate the raw data from the provided files into data matrix that can be used to generate a normalized count matrix and differential expression dataset.


```
## [1] FALSE
```

### Count Normalization

DESeq2 (version 1.18.1) was used to perform the raw count normalization within R (version 3.4.2)


```r
######### RUN DESeq2
dds<-DESeqDataSetFromMatrix(countData=countData, colData = colData, design= ~Race)
dds
```

```
## class: DESeqDataSet 
## dim: 58381 29 
## metadata(1): version
## assays(1): counts
## rownames(58381): ENSG00000223972.5 ENSG00000227232.5 ...
##   ENSG00000210195.2 ENSG00000210196.2
## rowData names(0):
## colnames(29): LVAD001 LVAD006 ... LVAD153 LVAD157
## colData names(35): LVAD_ID Zip9 ... RA Color
```

```r
# dds$ICM<-relevel(dds$ICM, ref = "NICM") # setting the reference to wild-type genotype. only relevant for factors.
#Determine the Dispersion Relationship (determines which distribution to use for the differential analysis) - should take about 2 minutes
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
plotDispEsts(dds)
```

![](ComprehensiveAnalysis_2021_files/figure-html/DESeq2-1.png)<!-- -->

```r
png(file=paste0("../3_Output/1_RNA/", COMPARISON, "/", COMPARISON, "_Dispersion.png"))
plotDispEsts(dds)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

There appears to be a linear negative correlation between the mean and dispersion estimates, so the parametric "Wald" model should be an appropriate fit for differential expression analysis. Furthermore, we could get away with the parametric fit-type, but the run-time is not significantly impaired, allowing us to use the 'local' fit-type. NOTE: If it were nonlinear throughout, we would require a 'local' nonparametric fit-type.

### Differential Expression Analysis


```r
##Pre-Filter to reduce the size of this dataset (according to the DESeq2 document reccomendations)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds
```

```
## class: DESeqDataSet 
## dim: 32903 29 
## metadata(1): version
## assays(2): counts mu
## rownames(32903): ENSG00000223972.5 ENSG00000227232.5 ...
##   ENSG00000210195.2 ENSG00000210196.2
## rowData names(10): baseMean baseVar ... dispOutlier dispMAP
## colnames(29): LVAD001 LVAD006 ... LVAD153 LVAD157
## colData names(36): LVAD_ID Zip9 ... Color sizeFactor
```

```r
################Run DESeq2 differential quantification (Likelihood ratio test (LRT) or Wald-test)
dds<-DESeq(dds, test="Wald", fitType="parametric")
#compile the results tables
resultsNames(dds)
```

```
## [1] "Intercept"                                  
## [2] "Race_African_American_vs_Caucasian_American"
```

```r
resdf<-as.data.frame(results(dds, format = "DataFrame"))
resdf$ensembl_gene_id<-as.character(row.names(resdf))
```

Once the differential Expression analysis was performed, the following were compiled into a results data matrix: Log2FoldChange, P-value, Bonferroni-Adjusted P-Value (Q-value), and normalized counts for each sample.


```r
####Add Annotation to the results file (this will take some time, about 5 minutes...)
##Add Gene Information
library(biomaRt)
hsapiens <- biomaRt::useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                                     host = "www.ensembl.org",
                                     dataset = "hsapiens_gene_ensembl",
                                     version = 94)
bm <- getBM(attributes=c("ensembl_gene_id_version", "external_gene_name", "chromosome_name", "start_position", "end_position"),  mart=hsapiens)
write.csv(bm, "../2_Input/BiomaRt_Annotation.csv")
bm<-read.csv("../2_Input/BiomaRt_Annotation.csv", row.names = 1)
results<-merge(resdf, bm, by.x="ensembl_gene_id", by.y="ensembl_gene_id_version")

####Add normalized count data (for heatmap and sota)
normcount<-as.data.frame(counts(dds, normalized=TRUE))
normcount$ensembl_gene_id<-rownames(normcount)
results<-dplyr::left_join(results, normcount, by="ensembl_gene_id")
results<-results[order(results$pvalue),] # order table by pvalue
##Create a Counts table with annotated Gene name
# Counts_table<-dplyr::select(results, external_gene_name, contains("LVAD"))
# write.xlsx(Counts_table, "../2_Input/CountData.xlsx", row.names = FALSE)
#Create filters as tabs
results_p05<-dplyr::filter(results, pvalue<0.05) 
results_q05<-dplyr::filter(results, padj<0.05)
library(openxlsx)
wb_DESeq<-createWorkbook()
#Unfiltered
  addWorksheet(wb_DESeq, "Unfiltered")
  writeData(wb_DESeq, "Unfiltered", results, startCol = 1)
#P-value Significant (0.05)
  addWorksheet(wb_DESeq, "P_0.05")
  writeData(wb_DESeq, "P_0.05", results_p05, startCol = 1)
#Q-value Significant (0.05)
  addWorksheet(wb_DESeq, "Q_0.05")
  writeData(wb_DESeq, "Q_0.05", results_q05, startCol = 1)
saveWorkbook(wb_DESeq, file = paste0("../3_Output/1_RNA/", COMPARISON, "/", COMPARISON,"_DESeq2.xlsx"), overwrite = TRUE)
```

### MDS Plot - Transcriptomics


```r
library(limma)
#Rewrite the plotMDS function to output gene index.
plotMDS.default<-
function (x, top = 500, labels = colnames(x), col = NULL, cex = 1, 
    dim.plot = c(1, 2), ndim = max(dim.plot), gene.selection = "pairwise", 
    xlab = paste("Dimension", dim.plot[1]), ylab = paste("Dimension", 
        dim.plot[2]), ...) 
{
    x <- as.matrix(x)
    ok <- is.finite(x)
    if (!all(ok)) 
        x <- x[apply(ok, 1, all), ]
    if (is.null(labels)) 
        labels <- 1:dim(x)[2]
    nprobes <- nrow(x)
    nsamples <- ncol(x)
    if (ndim < 2) 
        stop("Need at least two dim.plot")
    if (nsamples < ndim) 
        stop("Two few samples")
    gene.selection <- match.arg(gene.selection, c("pairwise", 
        "common"))
    cn <- colnames(x)
    dd <- matrix(0, nrow = nsamples, ncol = nsamples, dimnames = list(cn, 
        cn))
    topindex <- nprobes - top + 1
    if (gene.selection == "pairwise") {
        for (i in 2:(nsamples)) for (j in 1:(i - 1)) dd[i, j] = sqrt(mean(sort.int((x[, 
            i] - x[, j])^2, partial = topindex)[topindex:nprobes]))
    }
    else {
    
    #		Same genes used for all comparisons ,"common"
   
        s <- rowMeans((x - rowMeans(x))^2)
        q <- quantile(s, p = (topindex - 1.5)/(nprobes - 1))

        x <- x[s >= q, ]
	
	#	 an extra line
        ind.top.genes<-which(s >= q)
        
        for (i in 2:(nsamples)) dd[i, 1:(i - 1)] = sqrt(colMeans((x[, 
            i] - x[, 1:(i - 1), drop = FALSE])^2))
    }
    a1 <- cmdscale(as.dist(dd), k = ndim)
    mds <- new("MDS", list(dim.plot = dim.plot, distance.matrix = dd, 
        cmdscale.out = a1, top = top, gene.selection = gene.selection))
    mds$x <- a1[, dim.plot[1]]
    mds$y <- a1[, dim.plot[2]]
     mdsPlot<-plotMDS(mds, labels = labels, col = col, cex = cex, xlab = xlab, 
        ylab = ylab, ...)
list       (mds=mds,  ind.top.genes=ind.top.genes)  
}

# Perform MDS on RNA-sequencing data
MDS_data<-read.xlsx("../2_Input/CountData.xlsx")
rownames(MDS_data)<-make.unique(MDS_data$external_gene_name, sep = ".")
MDS_data<-dplyr::select(MDS_data, -external_gene_name)
## Create color based on race
ann_colors = list(Race = c(African_American="#1b9e77", Asian_American = "#d95f02", Caucasian_American="#7570b3", Eastern_African = "#e7298a"))
## MDS PLot
NGENES = 1000
mds<-plotMDS.default(MDS_data, gene.selection = "common", top = NGENES, col = colData_all$Color)
```

![](ComprehensiveAnalysis_2021_files/figure-html/MDS-1.png)<!-- -->

```r
## Extract genes from this
MDS_GeneList<-rownames(as.data.frame(mds$ind.top.genes))
TopGenes_MDS<-subset(results, ensembl_gene_id %in% MDS_GeneList)
```

### QQ Plot

Before we examined the gene networks and pathways differentially regulated by NRF2 knockout, the first task was to determine whether transgene induction resulted in global changes. An effective way of determining this is the QQ plot, which compares the P-value distribution produced by the pairwise comparison (transgenic vs. WT mouse) to that of a random normal distribution. Below, it is evident that the two experimental groups produce robustly divergent expression patterns consistent with a true population difference worthy of differential expression analysis.


```r
#Create Q-Q plot
test<-results
test<-test[complete.cases(test),]
pQQ(test$pvalue, lim=c(0,10))
```

![](ComprehensiveAnalysis_2021_files/figure-html/QQ-Plot-1.png)<!-- -->

```r
png(file=paste0("../3_Output/1_RNA/", COMPARISON,  "/", COMPARISON,"_QQ.Plot.png"))
pQQ(test$pvalue, lim=c(0,10))
dev.off()
```

```
## quartz_off_screen 
##                 2
```

## Principal Components Analysis

Once we established that the populations under consideration truly display divergene expression patterns, we sought to determine whether unbiased global gene expression patterns recapitulate the described phenotypes within each heart failure group. To accomplish this, an unsupervised Principal Components Analysis (PCA) was initially used with normalized counts.

### PCA Features

Before running the principal components analysis, it was necessary to first determine the number of PC's required to account for 80% of the variance, a machine-learning algorithmm benchmark that provides sufficient confidence in the analysis.


```r
#Plot Features of the PCA
library(dplyr)
library(plotly)
##Import the data to be used for PCA
normCount_all<-data.matrix(read.xlsx("../2_Input/RNA_input.xlsx", sheet = "CountData", rowNames = T))
#transpose the dataset (required for PCA)
data.pca<-t(normCount_all)
data.pca<-as.data.frame(data.pca)
##Import the data to be used for annotation
rownames(colData_all)<-colData_all$LVAD_ID
Index<-colData_all
Index<-as.data.frame(Index)
##merge the file
data.pca_Final<-merge(Index, data.pca, by=0)
rownames(data.pca_Final)<-data.pca_Final$Row.names
pca.comp<-prcomp(data.pca_Final[,(ncol(Index)+2):ncol(data.pca_Final)])

pcaCharts=function(x) {
    x.var <- x$sdev ^ 2
    x.pvar <- x.var/sum(x.var)
    par(mfrow=c(2,2))
    plot(x.pvar,xlab="Principal component", 
         ylab="Proportion of variance", ylim=c(0,1), type='b')
    plot(cumsum(x.pvar),xlab="Principal component", 
         ylab="Cumulative Proportion of variance", 
         ylim=c(0,1), 
         type='b')
    screeplot(x)
    screeplot(x,type="l")
    par(mfrow=c(1,1))
}
pcaCharts(pca.comp)
```

![](ComprehensiveAnalysis_2021_files/figure-html/PCA_Features-1.png)<!-- -->

```r
png(file=paste0("../3_Output/1_RNA/", COMPARISON,  "/", COMPARISON, "_PCA.Charts.png"))
pcaCharts(pca.comp)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

### 3-Dimensional PCA

From the previous calculations, it is seens that only 2 principal components are necessary (accounting for >80% cumulative variance). Nonetheless, below is a 3-D PCA to ensure that all groups are characterize to higher-degree of stringency.


```r
##Create a 3D-PCA for Inspection
library(plotly)
##Index
Index_PCA<-openxlsx::read.xlsx("../2_Input/1_Patient/colData_CHAAMPS.FINAL.xlsx", sheet="Summary", startRow=2)
# Index_PCA<-dplyr::filter(Index_PCA, Race %in% c("Caucasian_American", "African_American"))
Index_PCA$Race<-factor(Index_PCA$Race, levels = c("Caucasian_American","African_American", "Asian_American", "Eastern_African"))
rownames(Index_PCA)<-Index_PCA$LVAD_ID
PCs<-merge(pca.comp$x, Index_PCA, by=0)
rownames(PCs)<-PCs$Row.names
# PCs$Race[which(PCs$Race == 1)] <- 'Caucasian American'
# PCs$Race[which(PCs$Race == 2)] <- 'African American'
# PCs$Race[which(PCs$Race == 3)] <- 'Asian American'
# PCs$Race[which(PCs$Race == 4)] <- 'Eastern African'
PCs$Race <- as.factor(PCs$Race)
fig <- plot_ly(PCs, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Race, colors = c("#7570b3", "#1b9e77", "#d95f02", "#e7298a"))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                     yaxis = list(title = 'PC2'),
                     zaxis = list(title = 'PC3')))
fig
```

```{=html}
<div id="htmlwidget-e54acf2131c7178d98e5" style="width:672px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-e54acf2131c7178d98e5">{"x":{"visdat":{"10942eb245e4":["function () ","plotlyVisDat"]},"cur_data":"10942eb245e4","attrs":{"10942eb245e4":{"x":{},"y":{},"z":{},"color":{},"colors":["#7570b3","#1b9e77","#d95f02","#e7298a"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter3d","mode":"markers","inherit":true}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"scene":{"xaxis":{"title":"PC1"},"yaxis":{"title":"PC2"},"zaxis":{"title":"PC3"}},"hovermode":"closest","showlegend":true},"source":"A","config":{"showSendToCloud":false},"data":[{"x":[-1529671.24100752,-334757.29228236,-612190.108321944,-954148.72411323,-902822.044248369,609478.534018162,838302.877830159,191569.728652581,963619.325902555,602660.703674678,471261.11287434,527877.08254015,-1211596.0775955,64433.4835743343,1688571.30867338],"y":[-187407.202959543,549287.304917815,85946.6639209321,89809.8055439567,304197.314381109,-250630.799549259,-370973.902133611,277785.29998333,126177.293778808,-484546.370347696,129179.468069339,-110589.928830888,-120541.782049278,-558287.285069333,339103.68268931],"z":[453288.274958443,-164874.661364243,105642.044912724,140990.272939514,228914.923366101,-11370.7903122676,-230390.086561291,-320132.303509603,278143.428539124,-109886.089499411,105435.023751819,65445.2533584074,-42602.8770476994,15142.060083232,259427.094967338],"type":"scatter3d","mode":"markers","name":"Caucasian_American","marker":{"color":"rgba(117,112,179,1)","line":{"color":"rgba(117,112,179,1)"}},"textfont":{"color":"rgba(117,112,179,1)"},"error_y":{"color":"rgba(117,112,179,1)"},"error_x":{"color":"rgba(117,112,179,1)"},"line":{"color":"rgba(117,112,179,1)"},"frame":null},{"x":[-954318.863345667,-172032.536883695,184025.483503864,-140766.536258928,-597065.572068181,-222599.692905327,1494948.31527561,-59856.2869103183,547321.253483613,437165.702451831,217896.226535602,-450259.488661536,-292387.13896002,-1216230.04964235],"y":[138465.850399628,502639.269282457,786839.199015143,-87282.0808489276,-179298.576552183,-96413.730363207,-247830.756351568,-318168.678997766,-457675.880711732,-356748.841614771,488565.790725704,296161.608119268,-221533.183279004,-433821.539547759],"z":[135648.919356215,155339.153463123,-157590.214381606,259209.307112512,-342666.93358393,-631791.414125601,-384415.853753897,-533452.889839404,718518.110445435,407934.387426665,-13489.339280251,-98562.5026960845,-9955.60165754406,-299699.164963146],"type":"scatter3d","mode":"markers","name":"African_American","marker":{"color":"rgba(27,158,119,1)","line":{"color":"rgba(27,158,119,1)"}},"textfont":{"color":"rgba(27,158,119,1)"},"error_y":{"color":"rgba(27,158,119,1)"},"error_x":{"color":"rgba(27,158,119,1)"},"line":{"color":"rgba(27,158,119,1)"},"frame":null},{"x":[659974.244667412],"y":[542144.516070678],"z":[-70433.1788337212],"type":"scatter3d","mode":"markers","name":"Asian_American","marker":{"color":"rgba(217,95,2,1)","line":{"color":"rgba(217,95,2,1)"}},"textfont":{"color":"rgba(217,95,2,1)"},"error_y":{"color":"rgba(217,95,2,1)"},"error_x":{"color":"rgba(217,95,2,1)"},"line":{"color":"rgba(217,95,2,1)"},"frame":null},{"x":[151596.269546688],"y":[-174552.527690999],"z":[92235.6467290363],"type":"scatter3d","mode":"markers","name":"Eastern_African","marker":{"color":"rgba(231,41,138,1)","line":{"color":"rgba(231,41,138,1)"}},"textfont":{"color":"rgba(231,41,138,1)"},"error_y":{"color":"rgba(231,41,138,1)"},"error_x":{"color":"rgba(231,41,138,1)"},"line":{"color":"rgba(231,41,138,1)"},"frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

## Heatmap and Clustering of DEGs (P < 0.01)

In order to visualize the distribution of differentially expressed genes, as well as determine the effect various heart failure etiologies on transcription, hierarchical clustering and heatmap visualization were performed at the Q < 0.05 statistical level. This analysis reveals that P < 0.05 is sufficient to separate all samples by genotype.


```r
library(pheatmap)
results_p05<-filter(results, pvalue<0.01)
normCount_all<-read.xlsx("../2_Input/CountData.xlsx", rowNames = F)
rownames(normCount_all)<-make.unique(normCount_all$external_gene_name, sep = ".")
normCount_all<-data.matrix(normCount_all[,2:ncol(normCount_all)])
hm_data<-subset(normCount_all, rownames(normCount_all) %in% results_p05$external_gene_name)
hm_data<-data.matrix(hm_data)
##
##Index file for annotating samples
rownames(colData_all)<-colData_all$LVAD_ID
Index<-dplyr::select(colData_all, Race, Age, Cardiac.Index, BMI)
Index<-as.data.frame(Index)
ann_colors = list(Race = c(African_American="#1b9e77", Asian_American = "#d95f02", Caucasian_American="#7570b3", Eastern_African = "#e7298a"))
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "gold2"))(paletteLength)
pheatmap(hm_data, scale="row", 
         cluster_cols = TRUE, 
         cluster_rows = TRUE,
         #breaks = myBreaks,
         cutree_cols = 2,
         cutree_rows = 2,
         angle_col = 45,
         fontsize_col = 8,
         color = myColor, 
         show_rownames = FALSE, 
         border_color = NA, 
         annotation_colors = ann_colors,
         annotation_col = Index,
         filename=paste0("../3_Output/1_RNA/", COMPARISON,  "/", COMPARISON,"_Heatmap_Normcount.P01.pdf"))
vst<-varianceStabilizingTransformation(data.matrix(read.xlsx("../2_Input/RNA_input.xlsx", sheet = "CountData", rowNames = T)))
dists <- dist(t(vst))
hc<-hclust(dists, method = "ward.D2")
plot(hc)
pdf(paste0("../3_Output/1_RNA/", COMPARISON,  "/", COMPARISON,"_Unsupervised Dendrogram.pdf"))
plot(hc)
dev.off()
```

```
## pdf 
##   3
```

```r
normhm<-vst[row.names(resdf[which(resdf$pvalue<0.05),]),]
normhm<-scale(t(normhm))
normhm<-t(normhm)
pheatmap(normhm, #Variance stabilized transformation
         cluster_cols=T, 
         clustering_method = "ward.D2",
         border_color=NA, 
         cluster_rows=T, 
         scale = 'row', 
         show_colnames = T, 
         show_rownames = F, 
         annotation_colors = ann_colors,
         color = myColor,
         annotation_col = Index,
         filename=paste0("../3_Output/1_RNA/", COMPARISON,  "/", COMPARISON,"_VST.Heatmap.P01.pdf"))
```
## DEGs involved in Epigenetic Regulation

We used the following script to identify all epigenetic regulators that are differentially-expressed in cardiac samples of African Americans relative to their age-adjusted Caucasian American counterparts.


```r
##########Find all Epigeneric Regulators 
EPIs<-openxlsx::read.xlsx("../2_Input/Miscelleneous/EpiGenes_main.xlsx")
EPIs_Genes<-EPIs$HGNC_symbol
DEG_List<-read.xlsx(paste0("../3_Output/1_RNA/", COMPARISON, "/", COMPARISON, "_DESeq2.xlsx"), sheet = "Unfiltered")
Epi_DEGs<-merge(DEG_List, EPIs, by.x="external_gene_name", by.y="HGNC_symbol")
openxlsx:::write.xlsx(Epi_DEGs,paste0("../3_Output/1_RNA/", COMPARISON,"/",COMPARISON, "_Epigenetic.DEGs.xlsx"))
#############################
```

## Differentially-expressed Genes in heart failure.
### Gene Expression Bar Graphs


```r
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(dplyr)
library(gtools)
library(openxlsx)
ifelse(!dir.exists(file.path(paste0("../3_Output/1_RNA/Candidates"))), dir.create(file.path(paste0("../3_Output/1_RNA/Candidates"))), FALSE)
```

```
## [1] FALSE
```

```r
Counts<-read.xlsx("../2_Input/CountData.xlsx", rowNames = F)
colData<-openxlsx::read.xlsx("../2_Input/1_Patient/colData_CHAAMPS.FINAL.xlsx", sheet = "Summary", startRow = 2)
colData<-dplyr::filter(colData, Race %in% RACE, Ischemia %in% ISCHEMIA, Diabetes %in% DIABETES)
#Filter results by the gene vector
DEGs<-subset(Counts, external_gene_name %in% GENES)
rownames(DEGs)<-make.unique(as.character(DEGs$external_gene_name, sep = "."))
tDEGs<-as.data.frame(t(DEGs))
tDEGs = tDEGs[-1, ]
## convert all genes to numeric (from factors)
asNumeric=function(x){as.numeric(as.character(x))}
factorsNumeric=function(d){modifyList(d, lapply(d[, sapply(d, is.character)], asNumeric))}
##
tDEGs<-factorsNumeric(tDEGs)
tDEGs$LVAD_ID<-rownames(tDEGs)
colData.ex<-dplyr::inner_join(tDEGs, colData)
colData.ex$Race<-factor(colData.ex$Race, levels = c("Caucasian_American", "African_American"))
colData.ex$Ischemia<-factor(colData.ex$Ischemia, levels = c("NICM", "ICM"))
colData.ex<-dplyr::group_by_(colData.ex, VAR1) #define groups for statistics
groupsize<-dplyr::tally(colData.ex) #define sample size (for SEM)
#Define function for mean and standard deviation for each group
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}
#Create a p-value for every summarized data (relative to CON)
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
    if( equal.variance==FALSE ) 
    {
        se <- sqrt( (s1^2/n1) + (s2^2/n2) )
        # welch-satterthwaite df
        df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    } else
    {
        # pooled standard deviation, scaled by the sample sizes
        se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
        df <- n1+n2-2
    }      
    t <- (m1-m2-m0)/se 
    dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
    names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
    return(dat) 
}
########################################################
plotlist = list()
p<-1
if(any(grepl("package:plyr", search()))) detach("package:plyr") else message("plyr not loaded")
for (i in GENES){
  GENE<-as.name(i) #define a "name" for the variable
  ds <- data_summary(colData.ex, varname=i, groupnames=c(VAR1))
  ds<-dplyr::left_join(ds, groupsize)
  ds<-dplyr::mutate(ds, upper=ds[,2]+(sd))
  for(j in seq_along(ds[,4])){
    p[j]=round(t.test2(m1=ds[j,2], m2=ds[1,i], s1=ds[j,"sd"], s2=ds[1,"sd"], n1=ds[j,"n"], n2=ds[1,"n"], equal.variance = FALSE)[4], 2)
  } ### THIS MUST BE CORRECTED BASED ON THE ANALYSIS (change [,THIS] as needed to create: m = mean, s = standard deviation, n = sample size.
  ds$p_con<-stars.pval(p)
g<-ggplot(ds, aes_string(x=VAR1, y=i, fill = VAR1)) +
    geom_bar(stat="identity", color="black",
           position=position_dodge()) +
    # geom_point(shape="O",stat="identity", mapping = aes_string(x=VAR2, y=i),data = colData.ex, position=position_dodge(0.9), stroke = 2) +
    geom_errorbar(aes_string(ymin=i, ymax="upper"), width=.2,
                 position=position_dodge(.9)) +
    geom_text(aes(label=p_con), position=position_dodge(width=0.9), vjust=-1, size=5) +
    scale_fill_manual(values=c("white", "black", "darkcyan")) +
    ylim(NA, 1.1*max(ds$upper))
g_plot<-g+labs(title=paste0("Expression - ", i), x=VAR1, y = "Normalized Counts")+
   theme_classic2(base_size = 10)
pdf(file=paste0("../3_Output/1_RNA/Candidates/", i, "_", VAR1, ".and.", VAR1, "_Expression.pdf"), width = 3, height = 1.5)
print(g_plot)
dev.off()
plotlist[[i]] = g_plot
}
t<-marrangeGrob(grobs = plotlist, legend, nrow=3, ncol=2)
ggsave(paste0("../3_Output/1_RNA/Candidates/_DEGs_", VAR1, ".v.", VAR2, ".pdf"), t, width = 6, height = 7)
########################################################
```

## Differential Expression bar graphs


```r
#Plot as single bar graph
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(dplyr)
library(tidyr)
library(openxlsx)
#Import dataset
Counts<-read.xlsx("../2_Input/CountData.xlsx")
Counts[,2:ncol(Counts)]<-sapply(Counts[,2:ncol(Counts)], as.numeric) #convert the count data to numeric (somehow becomes character in excel)
colData<-openxlsx::read.xlsx("../2_Input/1_Patient/colData_CHAAMPS.FINAL.xlsx", sheet = "Summary", startRow = 2)
colData<-dplyr::filter(colData, Race %in% RACE, Ischemia %in% ISCHEMIA, Diabetes %in% DIABETES)
#Filter results by the gene vector
DEGs<-dplyr::filter(Counts, external_gene_name %in% GENES)
rownames(DEGs)<-DEGs$external_gene_name
tDEGs<-as.data.frame(t(DEGs))
tDEGs = tDEGs[-1,]
## convert all genes to numeric (from factors)
asNumeric=function(x){as.numeric(as.character(x))}
factorsNumeric=function(d){modifyList(d, lapply(d[, sapply(d, is.character)], asNumeric))}
##
tDEGs<-factorsNumeric(tDEGs)
tDEGs$LVAD_ID<-rownames(tDEGs)
colData.ex<-dplyr::inner_join(tDEGs, colData)
colData.ex$Race<-factor(colData.ex$Race, levels = c("Caucasian_American", "African_American"))
colData.ex$Diabetes<-factor(colData.ex$Diabetes, levels = c("ND", "T2D"))
colData.ex$Ischemia<-factor(colData.ex$Ischemia, levels = c("NICM", "ICM"))
colData.ex<-colData.ex %>% group_by_(VAR1) #define groups for statistics
groupsize<-colData.ex %>% tally() #define sample size (for SEM)
#Define function for mean and standard deviation for each group
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}
colData.exg<-tidyr::gather(colData.ex, "Gene.Symbol", "Gene.Expression", 1:length(DEGs$external_gene_name))
colData.exg$Gene.Symbol<-factor(colData.exg$Gene.Symbol, levels = GENES)
detach("package:dplyr", unload=TRUE)
ds_timing<-data_summary(colData.exg, varname="Gene.Expression", groupnames=c("Gene.Symbol", VAR1))
groupsize<-colData.ex %>% group_by(Race) %>% dplyr::tally() #Calculate sample sizes
ds_timing<-dplyr::left_join(ds_timing, groupsize) #add sample size to summary table
ds_timing<-dplyr::mutate(ds_timing, #make upper and lower bounds based on SEM
                         lower=Gene.Expression-(sd/sqrt(n-1)), 
                         upper=Gene.Expression+(sd/sqrt(n-1))) 
g<-ggplot(ds_timing, aes_string(x="Gene.Symbol", y="Gene.Expression", fill = VAR1)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=Gene.Expression, ymax=upper), width=.2, position=position_dodge(.9)) + 
    labs(y = "Normalized Counts", x = "Gene Name")+
    theme_classic() +
    scale_fill_manual(values=c('white', "black", "darkcyan")) +
    stat_compare_means(mapping = aes_string(x="Gene.Symbol", y="Gene.Expression"), 
                       data = colData.exg, label = "p.signif")
pdf(file="../3_Output/1_RNA/Candidates/_DEGs_singlebar.pdf", width = 5, height = 2)
print(g)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
########################################################
```

## Pathway-based Gene Expression


```r
#GENE SET PLOTTING
####This code is useful in determining DEGs that associate with a given pathway
########################################################
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(dplyr)
library(gtools)
######################################################
#Import dataset
Counts<-read.xlsx("../2_Input/CountData.xlsx")
#Import the Pathway Information
GENE.Pathway<-read.csv("../2_Input/Gene.Sets/Oxphos.KEGG.csv", header = F)
colnames(GENE.Pathway)<-"external_gene_name"
#Filter results by the gene vector
DEGSET<-Counts %>% inner_join(., GENE.Pathway)
rownames(DEGSET)<-DEGSET$external_gene_name
GENESET_FIL=as.character(DEGSET$external_gene_name)
colData<-openxlsx::read.xlsx("../2_Input/1_Patient/colData_CHAAMPS.FINAL.xlsx", sheet = "Summary", startRow = 2)
colData<-dplyr::filter(colData, Race %in% RACE, Ischemia %in% ISCHEMIA, Diabetes %in% DIABETES)
#Filter results by the gene vector
tDEGs<-as.data.frame(t(DEGSET))
tDEGs = tDEGs[-1, ]
## convert all genes to numeric (from factors)
asNumeric=function(x){as.numeric(as.character(x))}
factorsNumeric=function(d){modifyList(d, lapply(d[, sapply(d, is.character)], asNumeric))}
##
tDEGs<-factorsNumeric(tDEGs)
tDEGs$LVAD_ID<-rownames(tDEGs)
colData.ex<-dplyr::inner_join(tDEGs, colData)
colData.ex$Race<-factor(colData.ex$Race, levels = c("Caucasian_American", "African_American"))
colData.ex$Ischemia<-factor(colData.ex$Ischemia, levels = c("NICM", "ICM"))
colData.ex<-dplyr::group_by_(colData.ex, VAR1) #define groups for statistics
groupsize<-dplyr::tally(colData.ex) #define sample size (for SEM)
#Define function for mean and standard deviation for each group
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}
#Create a p-value for every summarized data (relative to CON)
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
    if( equal.variance==FALSE ) 
    {
        se <- sqrt( (s1^2/n1) + (s2^2/n2) )
        # welch-satterthwaite df
        df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    } else
    {
        # pooled standard deviation, scaled by the sample sizes
        se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
        df <- n1+n2-2
    }      
    t <- (m1-m2-m0)/se 
    dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
    names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
    return(dat) 
}
########################################################
plotlist = list()
p<-1
if(any(grepl("package:plyr", search()))) detach("package:plyr") else message("plyr not loaded")
for (i in GENESET_FIL){
  GENE<-as.name(i) #define a "name" for the variable
  ds <- data_summary(colData.ex, varname=i, groupnames=c(VAR1))
  ds<-dplyr::left_join(ds, groupsize)
  ds<-dplyr::mutate(ds, upper=ds[,2]+(sd))
  for(j in seq_along(ds[,4])){
    p[j]=round(t.test2(m1=ds[j,2], m2=ds[1,i], s1=ds[j,"sd"], s2=ds[1,"sd"], n1=ds[j,"n"], n2=ds[1,"n"], equal.variance = FALSE)[4], 2)
  } ### THIS MUST BE CORRECTED BASED ON THE ANALYSIS (change [,THIS] as needed to create: m = mean, s = standard deviation, n = sample size.
  ds$p_con<-stars.pval(p)
g<-ggplot(ds, aes_string(x=VAR1, y=i, fill = VAR1)) +
    geom_bar(stat="identity", color="black",
           position=position_dodge()) +
    # geom_point(shape="O",stat="identity", mapping = aes_string(x=VAR2, y=i),data = colData.ex, position=position_dodge(0.9), stroke = 2) +
    geom_errorbar(aes_string(ymin=i, ymax="upper"), width=.2,
                 position=position_dodge(.9)) +
    geom_text(aes(label=p_con), position=position_dodge(width=0.9), vjust=-1, size=5) +
    scale_fill_manual(values=c("white", "black", "darkcyan")) +
    ylim(NA, 1.1*max(ds$upper))
g_plot<-g+labs(title=paste0("Expression - ", i), x=VAR1, y = "Normalized Counts")+
   theme_classic2(base_size = 10)
pdf(file=paste0("../3_Output/1_RNA/Candidates/", i, "_", VAR1, ".and.", VAR1, "_Expression.pdf"), width = 3, height = 1.5)
print(g_plot)
dev.off()
plotlist[[i]] = g_plot
}
t<-marrangeGrob(grobs = plotlist, legend, nrow=3, ncol=2)
ggsave(paste0("../3_Output/1_RNA/Candidates/_DEGs_Pathway", VAR1, ".v.", VAR2, ".pdf"), t, width = 6, height = 7)
########################################################
```


<!-- ## Covariate Analysis of Gene Expression -->

<!-- ```{r Gene.Confound} -->
<!-- library(ggplot2) -->
<!-- library(gridExtra) -->
<!-- library(ggpubr) -->
<!-- library(dplyr) -->

<!-- VAR_GENE = "POSTN" -->

<!-- GENE<-filter(results, external_gene_name==VAR_GENE) -->
<!-- GENE.ex<-t(select(GENE, contains("LVAD"))) -->
<!-- colData.ex<-cbind(GENE.ex, colData) -->
<!-- # colData.ex$NPPB<-colData.ex$GENE.ex -->

<!-- ggplotRegression<-function(fit){ -->
<!-- require(ggplot2) -->
<!-- ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +  -->
<!--          geom_point() + -->
<!--          stat_smooth(method = "lm", col = "red") + -->
<!--          labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), -->
<!--         "Intercept =",signif(fit$coef[[1]],5 ), -->
<!--         " Slope =",signif(fit$coef[[2]], 5), -->
<!--         " P =",signif(summary(fit)$coef[2,4], 5))) -->
<!-- } -->
<!-- #color by groups -->
<!-- scatterPlot_EF <- ggplotRegression(lm(LV.EF~GENE.ex, data = colData.ex)) -->
<!-- # Marginal density plot of x (top panel) -->
<!-- scatterPlot_Age <- ggplotRegression(lm(Age~GENE.ex, data = colData.ex)) -->
<!-- # Marginal density plot of y (right panel) -->
<!-- scatterPlot_A1C <- ggplotRegression(lm(HbA1C~GENE.ex, data = colData.ex)) -->
<!-- ggarrange(scatterPlot_EF+rremove("legend"), -->
<!--           scatterPlot_Age+rremove("legend"), -->
<!--           scatterPlot_A1C, -->
<!--          labels=c("A", "B", "C"), -->
<!--          ncol=2, nrow=2, -->
<!--          font.label = list(size=12, face="bold", family="Times")) -->
<!-- ``` -->

# Pilot Analysis - Illumina(R) HumanMethylation 450k Array


```r
#   Differential expression analysis with limma
library(openxlsx)
library(magrittr)
library(dplyr)
library(ggpubr)
library(matrixStats)
library("ggrepel")
library(ggplot2)
library(cowplot) 
m450k<-read.xlsx("../2_Input/_Pilot/Pepin.Wende_2018.xlsx", sheet = "betaData", rowNames = T)
colData_m450k<-read.xlsx("../2_Input/_Pilot/Pepin.Wende_2018.xlsx", sheet = "colData")
# MDS in ggplot2
Ntop = 1000
RowVar<-rowVars(data.matrix(m450k)) #calculate variances for each row (vector)
MDS.set<-as.data.frame(cbind(m450k, RowVar))
MDS_matrix<-MDS.set %>% arrange(desc(RowVar)) %>% top_n(Ntop,RowVar) #Select top N rows by variance
# Compute MDS
mds <- MDS_matrix %>% select(-RowVar) %>% t(.) %>%
  dist() %>%          
  cmdscale() %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")
rownames(mds)<-colnames(m450k)
mds$Sample_Name<-rownames(mds)
mds<-dplyr::inner_join(mds, colData_m450k)

#K-means clustering
clust <- kmeans(mds[,1:2], 2)$cluster %>%
  as.factor()
mds <- mds %>%
  mutate(kmeans.2 = clust)
# Main plot
COLOR=c("#1b9e77", "#7570b3")
pmain <- ggplot(mds, aes(x = Dim.1, y = Dim.2, color = Race))+
  scale_color_manual(values = COLOR) +
  theme(panel.background = element_rect("white", colour = "black", size=2), 
      panel.grid.major = element_line(colour = "gray50", size=.75), 
      panel.grid.minor = element_line(colour = "gray50", size=0.4),
      legend.position="bottom",
      legend.key=element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face="bold")) +
  geom_hline(yintercept = 0, size = 1) + 
  geom_vline(xintercept=0, size=1) +
  geom_point()+
  stat_ellipse()+
  geom_text_repel(data=mds, aes(label=Sample_Name), show_guide  = F) +
  labs(x="Principal Component 1", 
       y="Principal Component 2")
# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = mds, aes(x = Dim.1, fill = Race),
              alpha = 0.7, size = 0.2)+
  scale_fill_manual(values = COLOR)
# Marginal densities along y axis
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+ #must set coord_flip = true if using coord_flip() below
  geom_density(data = mds, aes(x = Dim.2, fill = Race),
                alpha = 0.7, size = 0.2)+
  coord_flip()+
  scale_fill_manual(values = COLOR)
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
pdf(file="../3_Output/2_Methyl/Pilot.Scatterhist.pdf", height = 7, width = 7, onefile = F)
ggdraw(p2)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

# Genome-wide DNA methylation - Illumina(R) HumanMethylation EPIC Array


```
## [1] FALSE
```

```
## [1] "../2_Input/3_Methyl/Wende_CHAAMP_MethylationEPIC_Sample_Sheet.csv"
```

```
## quartz_off_screen 
##                 2
```

```
## class: RGChannelSet 
## dim: 1051815 32 
## metadata(0):
## assays(2): Green Red
## rownames(1051815): 1600101 1600111 ... 99810990 99810992
## rowData names(0):
## colnames(32): LVAD001 LVAD084 ... LVAD083 LVAD092
## colData names(17): Sample_Name Sample_Well ... Basename filenames
## Annotation
##   array: IlluminaHumanMethylationEPIC
##   annotation: ilm10b4.hg19
```

```
## quartz_off_screen 
##                 2
```

![](ComprehensiveAnalysis_2021_files/figure-html/EPIC_Import-1.png)<!-- -->![](ComprehensiveAnalysis_2021_files/figure-html/EPIC_Import-2.png)<!-- -->

```
## quartz_off_screen 
##                 2
```

![](ComprehensiveAnalysis_2021_files/figure-html/EPIC_Import-3.png)<!-- -->

```
## quartz_off_screen 
##                 2
```

![](ComprehensiveAnalysis_2021_files/figure-html/EPIC_Import-4.png)<!-- -->

```
## quartz_off_screen 
##                 2
```

```
## quartz_off_screen 
##                 2
```

## SNPs in Methylation Data

SNPs have been shown to confound the interpretation of methylation arrays when mutations exist within the CpG sites. To address this concern, a package called “MethylToSNP” exists to identify novel SNPs based on methylation array data, which has been proposed to reduce the number of CpGs that are filtered. Using this approach, 1,294 CpGs were identified as putative SNPs which are likely influenced by minor allele fractions (MAFs); among these, 1,076 (83%) have been idenbtified as SNPs via published genomic sequencing analyses.

To determine whether the SNP frequency in our dataset contributes to racial differences, the beta values within these SNPs were clustered via hierarchical clustering, resulting in a complete separation according to patient race (see figure).


```r
#Remove SNPs
##Option 1: Built-in function that removes all known SNPs
gRatioSet_noSNPs<-dropLociWithSnps(gRatioSet.quantile, snps = c("SBE", "CpG"), maf = 0) #Doing this removes ~40% of all SNPs in the EPIC array
##Option 2: Identify putative SNPs using methylation barcoding (i.e. "gap hunting")
library("MethylToSNP")
library("RColorBrewer")
library("pheatmap")
Mvalues<-as.data.frame(getM(MSet))
x <- MethylToSNP(MSet, SNP=SNPs.147CommonSingle, verbose=TRUE)
x$CpG_ID<-rownames(x)
###
pdf(file="../3_Output/2_Methyl/Putative.SNPs.pdf", width = 10, height = 5) #Print Putative SNP Methylation
plotPotentialSNPs(x, MSet)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
plotPotentialSNPs(x, MSet)
```

![](ComprehensiveAnalysis_2021_files/figure-html/SNPs-1.png)<!-- -->

```r
SNPs<-merge(x, as.data.frame(getM(MSet)), by = "row.names")
write.csv(SNPs, "SNPs.csv")
SNPs_matrix<-SNPs %>% set_rownames(.$Row.names) %>% select(contains("LVAD"), -LVAD054)  %>% data.matrix()
SNPs_matrix<-SNPs_matrix[!is.infinite(rowSums(SNPs_matrix)),] #One infinite value exists!! (took ~2-3 days to troubleshoot)

# Index
Index<-openxlsx::read.xlsx("../2_Input/1_Patient/colData_CHAAMPS.FINAL.xlsx", sheet = "Summary", startRow = 2, rowNames = TRUE)
Index_SNPs<-Index[colnames(SNPs_matrix),] %>% select(Race)
ann_colors = list(Race = c(African_American="#1b9e77", Asian_American = "#d95f02", Caucasian_American="#7570b3", Eastern_African = "#e7298a"))
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "gold2"))(paletteLength)
ann_colors = list(Race = c(African_American="#1b9e77", 
                           Asian_American = "#d95f02", 
                           Caucasian_American="#7570b3", 
                           Eastern_African = "#e7298a"))

heatmap_SNP<-pheatmap(SNPs_matrix, scale="row", #Heatmap of SNPs
                      cluster_cols = TRUE, 
                      cluster_rows = TRUE,
                     cutree_cols = 3,
                     cutree_rows = 4,
                     angle_col = 45,
                     fontsize_col = 8,
                     color = myColor, 
                     show_rownames = FALSE, 
                     border_color = NA, 
                     annotation_colors = ann_colors,
                     annotation_col = Index_SNPs,
                    filename = paste0("../3_Output/2_Methyl/SNPS.heatmap.pdf"))
## Cluster Analysis
hc <-heatmap_SNP$tree_row
lbl <- cutree(hc, 4) # split gene dendrogram in 5 groups
cluster1<-which(lbl==1)
cluster2<-which(lbl==2)
cluster3<-which(lbl==3)
cluster4<-which(lbl==4) 
#
Cluster1_data<-SNPs_matrix[cluster1,]
Cluster2_data<-SNPs_matrix[cluster2,]
Cluster3_data<-SNPs_matrix[cluster3,]
Cluster4_data<-SNPs_matrix[cluster4,]
#
heatmap_c1<-pheatmap(Cluster1_data, scale="row", 
                      cluster_cols = TRUE, 
                      cluster_rows = TRUE,
                      #breaks = myBreaks,
                      cutree_cols = 2,
                      angle_col = 45,
                      fontsize_col = 8,
                      color = myColor, 
                      show_rownames = FALSE, 
                      border_color = NA, 
                      annotation_colors = ann_colors,
                     annotation_col = Index_SNPs,
                      filename = paste0("../3_Output/2_Methyl/SNPS_Cluster1.heatmap.pdf"))
#
heatmap_c2<-pheatmap(Cluster2_data, scale="row", 
                      cluster_cols = TRUE, 
                      cluster_rows = TRUE,
                      #breaks = myBreaks,
                      cutree_cols = 2,
                      angle_col = 45,
                      fontsize_col = 8,
                      color = myColor, 
                      show_rownames = FALSE, 
                      border_color = NA, 
                     annotation_colors = ann_colors,
                     annotation_col = Index_SNPs,
                      filename = paste0("../3_Output/2_Methyl/SNPS_Cluster2.heatmap.pdf"))
#
heatmap_c3<-pheatmap(Cluster3_data, scale="row", 
                      cluster_cols = TRUE, 
                      cluster_rows = TRUE,
                      #breaks = myBreaks,
                      cutree_cols = 2,
                      angle_col = 45,
                      fontsize_col = 8,
                      color = myColor, 
                      show_rownames = FALSE, 
                      border_color = NA, 
                     annotation_colors = ann_colors,
                     annotation_col = Index_SNPs,
                      filename = paste0("../3_Output/2_Methyl/SNPS_Cluster3.heatmap.pdf"))
#
heatmap_c4<-pheatmap(Cluster4_data, scale="row", 
                      cluster_cols = TRUE, 
                      cluster_rows = TRUE,
                      #breaks = myBreaks,
                      cutree_cols = 2,
                      angle_col = 45,
                      fontsize_col = 8,
                      color = myColor, 
                      show_rownames = FALSE, 
                      border_color = NA, 
                     annotation_colors = ann_colors,
                     annotation_col = Index_SNPs,
                      filename = paste0("../3_Output/2_Methyl/SNPS_Cluster4.heatmap.pdf"))
```

##MDS Plot - Unsupervised Clustering by Race

Because such a robust racial signature in cardiac DNA methylation was seen in the pilot analysis, we reproduced the unsupervised method in the larger cohort. This time, we continue to see a distinct racial difference. Furthermore, we found that this racially-determined clustering persisted to among the 500,000 most-variable CpG probes in the EPIC array.




# Differential Methylation of African American and Caucasian_American
##DMPs

Owing to the remarkable separation by patient race, we chose to identify the CpG sites responsible for a race-based epigenomic difference.


```r
library(minfi)
library(limma)
library(shinyMethyl)
library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
#########Part 1: Importing the Data
#Parameters
##Set the experimental conditions
DIABETES=c("ND", "T2D") #Should define the reference group first
ISCHEMIA=c("NICM", "ICM")
RACE=c("Caucasian_American", "African_American")
ISO=c("NONE", "DOB")
STATISTIC = 0.05 #P statistic threshold used in this combination.
VARIABLE = RACE
COMPARISON= paste0(VARIABLE[2], ".vs.", VARIABLE[1])
ifelse(!dir.exists(file.path("../3_Output/2_Methyl/", COMPARISON)), dir.create(file.path("../3_Output/2_Methyl/", COMPARISON)), FALSE)
```

```
## [1] FALSE
```

```r
##Get the array annotation
annoEPIC<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
annoEPIC<-dplyr::select(as.data.frame(annoEPIC), Name, chr, pos, Relation_to_Island, UCSC_RefGene_Name, UCSC_RefGene_Accession, UCSC_RefGene_Group, Regulatory_Feature_Group)
#Import the sample sheet
targets<-read.metharray.sheet(base="../2_Input/3_Methyl", pattern="csv$")
```

```
## [1] "../2_Input/3_Methyl/Wende_CHAAMP_MethylationEPIC_Sample_Sheet.csv"
```

```r
targets<-dplyr::filter(targets, Sample_Name!="LVAD054") #OUTLIER
# Filter targets
targets_filtered<-subset(targets, Race %in% RACE & Diabetes %in% DIABETES & Ischemia %in% ISCHEMIA)
#Import the array data from input directory (red and green .idat files)
RGSet<-read.metharray.exp(base = "../2_Input/3_Methyl/", targets = targets_filtered, verbose = TRUE)
sampleNames(RGSet)<-targets_filtered$Sample_Name
#Check the phenotype metadata to verify correct parsing
phenoData<-pData(RGSet)
#Get the manifest annotation (EPIC, methyl450k, etc...)
manifest<-getManifest(IlluminaHumanMethylationEPICmanifest)
typeIProbes <- getProbeInfo(manifest, type = "I")$Name
##Quality Control
#First step is to identify CpGs that failed to identify methylated positions (defined by expression intensity that reflects background levels)
detP<-detectionP(RGSet)
failed<-detP>0.01
#Determine the fraction of "failed" CpG probes (those which failed to identify a methylated CpG)
colMeans(failed)
```

```
##      LVAD001      LVAD084      LVAD058      LVAD014      LVAD122      LVAD006 
## 0.0002089850 0.0001916658 0.0001200798 0.0002678702 0.0001974388 0.0003175186 
##      LVAD013      LVAD133      LVAD111      LVAD075      LVAD066      LVAD067 
## 0.0001697281 0.0001766558 0.0003129001 0.0003486932 0.0003048179 0.0001731920 
##      LVAD125      LVAD142      LVAD060      LVAD127      LVAD079      LVAD087 
## 0.0002286134 0.0002678702 0.0003452293 0.0004018054 0.0003683216 0.0002771071 
##      LVAD141      LVAD132      LVAD121      LVAD032      LVAD130      LVAD106 
## 0.0001835835 0.0003025086 0.0002320772 0.0002112942 0.0002136034 0.0001177705 
##      LVAD157      LVAD153      LVAD140      LVAD083      LVAD092 
## 0.0002205311 0.0001951296 0.0002274588 0.0002147580 0.0001604912
```

```r
#Convert R/G to Methylated/Unmethylated in an object of class MethylSet
MSet<-preprocessRaw(RGSet)
#QC data
qc<-getQC(MSet)
plotQC(qc)
```

![](ComprehensiveAnalysis_2021_files/figure-html/African American.Caucasian_American-1.png)<!-- -->

```r
##Density plot
gRatioSet.quantile <- preprocessQuantile(RGSet, fixOutliers = TRUE, removeBadSamples = TRUE, badSampleCutoff = 10.5, quantileNormalize = TRUE, stratified = TRUE, mergeManifest = TRUE, sex = NULL)
#Quantification and Differential Expression Analysis
mVals<-getM(gRatioSet.quantile)
mVals<-mVals[,targets_filtered$Sample_Name]
phenoData<-phenoData[targets_filtered$Sample_Name,] #ensure that mVals has the same data as phenoData, otherwise dmpFinder won't work.
#Use Limma to perform differential methylation analysis using M-values (logit-transformed methylation data)
targets_filtered$Race<-factor(targets_filtered$Race)
targets_filtered$Race<-relevel(targets_filtered$Race, ref = "Caucasian_American")
targets_filtered$Diabetes<-factor(targets_filtered$Diabetes)
targets_filtered$ICM<-factor(targets_filtered$Ischemia)
design<-model.matrix(~0+Race, data = targets_filtered)
fit <- lmFit(mVals, design)
contMatrix <- makeContrasts(Race=RaceAfrican_American-RaceCaucasian_American,
                           levels=design)
contMatrix
```

```
##                         Contrasts
## Levels                   Race
##   RaceCaucasian_American   -1
##   RaceAfrican_American      1
```

```r
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))
```

```
##          Race
## Down      550
## NotSig 864898
## Up        411
```

```r
DMPs <- topTable(fit2, number = "all", adjust.method = "BH")
DMPs<-merge(DMPs, annoEPIC, by= 0)
rownames(DMPs)<-DMPs$Row.names
DMPs<-DMPs %>% select(-Row.names)
beta <- getBeta(gRatioSet.quantile)
#create an annotation table
beta.table<-beta
colnames(beta.table)<-phenoData$Sample_Name
write.csv(beta.table, paste0("../3_Output/2_Methyl/", COMPARISON, "/", COMPARISON, "_beta.table.csv"))
# annotated<-merge(beta.table, as.data.frame(Islands.UCSC), by = 0)
# annotated<-merge(annotated, as.data.frame(Locations), by.x = "Row.names", by.y = 0)
# annotated<-merge(annotated, as.data.frame(Other), by.x = "Row.names", by.y = 0)
# rownames(annotated)<-annotated$Row.names

#Merge annotation with differential methylation table
# Results_dmp<-merge(DMPs, annotated, by = 0) #
Results_dmp<-merge(DMPs, beta.table, by = 0)
#Calculate Average CpG Methylation by race - African_American
library(dplyr)
library(matrixStats)
Results<-Results_dmp %>% replace(is.na(.), 0) %>% mutate(
  African_American_SD = rowSds(as.matrix(Results_dmp[,targets_filtered$Sample_Name[targets_filtered$Race=="African_American"]])),
  African_American_Mean = rowMeans(as.matrix(Results_dmp[,targets_filtered$Sample_Name[targets_filtered$Race=="African_American"]])),
  Caucasian_American_SD = rowSds(as.matrix(Results_dmp[,targets_filtered$Sample_Name[targets_filtered$Race=="Caucasian_American"]])),
  Caucasian_American_Mean = rowMeans(as.matrix(Results_dmp[,targets_filtered$Sample_Name[targets_filtered$Race=="Caucasian_American"]])),
  Methylation.Diff=(African_American_Mean-Caucasian_American_Mean)*100
)
rownames(Results)<-Results$Row.names
Results_dmp_p05<-filter(Results, P.Value<0.05)
Results_dmp_q05<-filter(Results, adj.P.Val<0.05)
#########################################
#Identify Promoter-associated CpG Islands
library(tidyr)
PromCGI<-dplyr::filter(Results_dmp_p05, grepl("Island", Relation_to_Island), grepl("TSS", UCSC_RefGene_Group))
#Separate Gene Names into unique rows
PromCGI_sep<-PromCGI %>% mutate(UCSC_RefGene_Name = strsplit(as.character(UCSC_RefGene_Name), ";")) %>% unnest(UCSC_RefGene_Name) %>% distinct()

#Save a copy of the countData
library(openxlsx)
wb_countData<-createWorkbook()
addWorksheet(wb_countData, "Unfiltered")
  writeData(wb_countData, "Unfiltered", Results, startCol = 1)
addWorksheet(wb_countData, "P_0.05")
  writeData(wb_countData, "P_0.05", Results_dmp_p05, startCol = 1)
addWorksheet(wb_countData, "Q_0.05")
  writeData(wb_countData, "Q_0.05", Results_dmp_q05, startCol = 1)
addWorksheet(wb_countData, "Promoter.CGI")
  writeData(wb_countData, "Promoter.CGI", PromCGI_sep, startCol = 1)
saveWorkbook(wb_countData, file = paste0("../3_Output/2_Methyl/", COMPARISON, "/", COMPARISON, "_DMPs.xlsx"), overwrite = TRUE)
```

##DMRs


```r
library(DMRcate)
library(missMethyl)
library(biomaRt)
myAnnotation<-cpg.annotate(object = mVals, datatype = "array", what = "M", 
                           analysis.type = "differential", design = design, 
                           contrasts = TRUE, cont.matrix = contMatrix, 
                           coef = "Race", arraytype = "EPIC")
design_test<-model.matrix(~targets_filtered$Race)
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2) #calculates DMRs
results.ranges <- extractRanges(DMRs)
beta<-beta[,targets_filtered$Sample_Name] #ensure that columns are ordered exactly the same as the targets_filtered table.
groups <- c(African_American="#1b9e77", Caucasian_American="#7570b3")
type<-factor(targets_filtered$Race)
cols <- groups[as.character(type)] #creates a string of colors for the DMR.plot function
pdf(file = "../3_Output/2_Methyl/DMR.top.pdf")
DMR.plot(ranges = results.ranges, dmr = 3, CpGs = beta, what = "Beta", arraytype = "EPIC", phen.col = cols, genome = "hg38")
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
results.ranges.sig<-as.data.frame(results.ranges) %>% filter(Fisher < 0.05)
enrichment_GO <- goregion(results.ranges[(elementMetadata(results.ranges)[, "overlapping.genes"] %in% results.ranges.sig$overlapping.genes)], all.cpg = rownames(rownames(beta)), collection = "GO", array.type = "EPIC")

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('external_gene_name', 'ensembl_transcript_id', 'go_id'), filters = 'go', values = 'GO:0006103', mart = ensembl)
```
## Distribution of Methylation by Genomic and CpG Annotation

The following figure illustrates the enrichment of differential methylation within CpG "Open Seas" (non-Island associatd CpG sites) found within the associated gene body.  

**change the viewing angle of 3D-Histogram by dragging mouse across it.**


```r
library(plotly)
library(dplyr)
library(stringr)
library(reshape2)
library(readxl)
library(kableExtra)

##Create a regional annotation matrix of the EPIC array
Region_epic<-Results_dmp %>% 
  select(Row.names, UCSC_RefGene_Group, Relation_to_Island)
Stage1<-Region_epic %>% 
  mutate(UCSC_RefGene_Group = strsplit(as.character(UCSC_RefGene_Group), ";")) %>% 
  unnest(UCSC_RefGene_Group) 
Stage2<-Stage1 %>% 
  mutate(Relation_to_Island = strsplit(as.character(Relation_to_Island), ";")) %>%  
  unnest(Relation_to_Island)
Stage3<-distinct(Stage2)
Regional.Groups<-dplyr::group_by_(Stage3, "UCSC_RefGene_Group", "Relation_to_Island") %>% 
  tally()
Region.matrix<-Regional.Groups %>% 
  spread(Relation_to_Island, n)
Region.matrix<-Region.matrix %>% 
  select(UCSC_RefGene_Group, OpenSea, N_Shelf, N_Shore, Island, S_Shore, S_Shelf)
rownames(Region.matrix)<-Region.matrix$UCSC_RefGene_Group
Region.matrix<-Region.matrix[c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "ExonBnd", "3'UTR"),]
rownames(Region.matrix)<-Region.matrix$UCSC_RefGene_Group
Contour_3D_epic<-Region.matrix[,-1]
rownames(Contour_3D_epic)<-Region.matrix$UCSC_RefGene_Group
write.xlsx(Contour_3D_epic, "../3_Output/2_Methyl/Contour.3D_EPIC.xlsx")

##Create a regional annotation matrix of differentially-methylated positions
Region_p05<-Results_dmp_p05 %>% select(Row.names, UCSC_RefGene_Group, Relation_to_Island)
Stage1<-Region_p05 %>% mutate(UCSC_RefGene_Group = strsplit(as.character(UCSC_RefGene_Group), ";")) %>% unnest(UCSC_RefGene_Group) 
Stage2<-Stage1 %>% mutate(Relation_to_Island = strsplit(as.character(Relation_to_Island), ";")) %>%  unnest(Relation_to_Island)
Stage3<-distinct(Stage2)
Regional.Groups<-dplyr::group_by_(Stage3, "UCSC_RefGene_Group", "Relation_to_Island") %>% tally()
Region.matrix<-Regional.Groups %>% spread(Relation_to_Island, n)
Region.matrix<-Region.matrix %>% select(UCSC_RefGene_Group, OpenSea, N_Shelf, N_Shore, Island, S_Shore, S_Shelf)
rownames(Region.matrix)<-Region.matrix$UCSC_RefGene_Group
Region.matrix<-Region.matrix[c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "ExonBnd", "3'UTR"),]
rownames(Region.matrix)<-Region.matrix$UCSC_RefGene_Group
Contour_3D<-Region.matrix[,-1]
rownames(Contour_3D)<-Region.matrix$UCSC_RefGene_Group
write.xlsx(Contour_3D, "../3_Output/2_Methyl/Contour.3D_DMPs.p05.xlsx")
### Identify DMP Enrichment (IMPORTANT)
Enrichment_Region<-Contour_3D/Contour_3D_epic
Promoter_enr<-Enrichment_Region[rownames(Enrichment_Region)=="TSS200",] + Enrichment_Region[rownames(Enrichment_Region)=="TSS1500",]
rownames(Promoter_enr)<-"Promoter"
Enrichment_Region<-Enrichment_Region %>% 
  filter(!grepl("TSS", rownames(Enrichment_Region))) %>%
  rbind(Promoter_enr, .) %>%
  data.matrix()
paletteLength<-100
myColor <- colorRampPalette(c("dodgerblue4", "white", "brown4"))(paletteLength)
pheatmap(Enrichment_Region, color = myColor, cluster_rows = FALSE, cluster_cols = FALSE)
```

![](ComprehensiveAnalysis_2021_files/figure-html/Methylation Distribution-1.png)<!-- -->

```r
##Make a Table of the CpG Methylation Distribution
Enrichment_Region %>% kable( align="c", booktabs=T, 
                     caption="Methylation Distribution") %>% 
  kable_styling(latex_options=c("striped", "condensed", "repeat_header"))
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Methylation Distribution</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:center;"> OpenSea </th>
   <th style="text-align:center;"> N_Shelf </th>
   <th style="text-align:center;"> N_Shore </th>
   <th style="text-align:center;"> Island </th>
   <th style="text-align:center;"> S_Shore </th>
   <th style="text-align:center;"> S_Shelf </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Promoter </td>
   <td style="text-align:center;"> 0.1285693 </td>
   <td style="text-align:center;"> 0.1358225 </td>
   <td style="text-align:center;"> 0.1551988 </td>
   <td style="text-align:center;"> 0.1823682 </td>
   <td style="text-align:center;"> 0.1568289 </td>
   <td style="text-align:center;"> 0.1320130 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5'UTR </td>
   <td style="text-align:center;"> 0.0574760 </td>
   <td style="text-align:center;"> 0.0620397 </td>
   <td style="text-align:center;"> 0.0677872 </td>
   <td style="text-align:center;"> 0.0905560 </td>
   <td style="text-align:center;"> 0.0667659 </td>
   <td style="text-align:center;"> 0.0582758 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1stExon </td>
   <td style="text-align:center;"> 0.0564209 </td>
   <td style="text-align:center;"> 0.0434783 </td>
   <td style="text-align:center;"> 0.0705167 </td>
   <td style="text-align:center;"> 0.0847120 </td>
   <td style="text-align:center;"> 0.0715037 </td>
   <td style="text-align:center;"> 0.0453401 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Body </td>
   <td style="text-align:center;"> 0.0600375 </td>
   <td style="text-align:center;"> 0.0646733 </td>
   <td style="text-align:center;"> 0.0697822 </td>
   <td style="text-align:center;"> 0.0974912 </td>
   <td style="text-align:center;"> 0.0749590 </td>
   <td style="text-align:center;"> 0.0661616 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ExonBnd </td>
   <td style="text-align:center;"> 0.0463516 </td>
   <td style="text-align:center;"> 0.0348432 </td>
   <td style="text-align:center;"> 0.0659091 </td>
   <td style="text-align:center;"> 0.0752212 </td>
   <td style="text-align:center;"> 0.0555556 </td>
   <td style="text-align:center;"> 0.0375940 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3'UTR </td>
   <td style="text-align:center;"> 0.0630122 </td>
   <td style="text-align:center;"> 0.0531746 </td>
   <td style="text-align:center;"> 0.0837796 </td>
   <td style="text-align:center;"> 0.1099838 </td>
   <td style="text-align:center;"> 0.0596878 </td>
   <td style="text-align:center;"> 0.0493359 </td>
  </tr>
</tbody>
</table>

```r
write.xlsx(Enrichment_Region, paste0("../3_Output/2_Methyl/", COMPARISON, "/", COMPARISON, "_DMP.Enrichment_3D.xlsx"), rowNames = TRUE)
color <- colorRampPalette(c("grey", "orange", "red"))
t <- list(
  family = "times",
  size = 16,
  color = "black")
x_axis<-list(title = 'CpG Region', 
                     type="category", 
                     zeroline=TRUE, 
                     showline=TRUE, 
                     zerolinewidth = 4, 
            zerolinecolor="darkgrey", 
            linecolor="darkgrey", 
            linewidth=4, 
            titlefont=t, 
            tickfont=t)
y_axis<-list(title = 'Gene Region', 
                     type="category", 
                     zeroline=TRUE, 
                     showline=TRUE, 
                     zerolinewidth = 4, 
            zerolinecolor="darkgrey", 
            linecolor="darkgrey", 
            linewidth=4, 
            titlefont=t, 
            tickfont=t)
z_axis<-list(title = 'Number of DMPs', 
                     zerolinewidth = 4, 
                    zerolinecolor="darkgrey", 
                    linecolor="darkgrey", 
                    linewidth=4, 
                    titlefont=t, 
                    tickfont=t)
q<-plot_ly(z=~Enrichment_Region, colors=color(10), 
    text=as.character(rownames(Enrichment_Region))) %>% add_surface() %>% 
    layout(scene = list(xaxis = x_axis, yaxis = y_axis, zaxis = z_axis))

q #must comment out for PDF generation via knitr (Pandoc).
```

```{=html}
<div id="htmlwidget-4bb55bf61ae4e91cef1e" style="width:672px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-4bb55bf61ae4e91cef1e">{"x":{"visdat":{"10943b6d4374":["function () ","plotlyVisDat"]},"cur_data":"10943b6d4374","attrs":{"10943b6d4374":{"z":{},"text":["Promoter","5'UTR","1stExon","Body","ExonBnd","3'UTR"],"colors":["#BEBEBE","#CCB893","#DAB269","#E9AD3F","#F7A715","#FF9200","#FF6E00","#FF4900","#FF2400","#FF0000"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"surface","inherit":true}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"scene":{"xaxis":{"title":"CpG Region","type":"category","zeroline":true,"showline":true,"zerolinewidth":4,"zerolinecolor":"darkgrey","linecolor":"darkgrey","linewidth":4,"titlefont":{"family":"times","size":16,"color":"black"},"tickfont":{"family":"times","size":16,"color":"black"}},"yaxis":{"title":"Gene Region","type":"category","zeroline":true,"showline":true,"zerolinewidth":4,"zerolinecolor":"darkgrey","linecolor":"darkgrey","linewidth":4,"titlefont":{"family":"times","size":16,"color":"black"},"tickfont":{"family":"times","size":16,"color":"black"}},"zaxis":{"title":"Number of DMPs","zerolinewidth":4,"zerolinecolor":"darkgrey","linecolor":"darkgrey","linewidth":4,"titlefont":{"family":"times","size":16,"color":"black"},"tickfont":{"family":"times","size":16,"color":"black"}}},"hovermode":"closest","showlegend":false,"legend":{"yanchor":"top","y":0.5}},"source":"A","config":{"showSendToCloud":false},"data":[{"colorbar":{"title":"Enrichment_Region","ticklen":2,"len":0.5,"lenmode":"fraction","y":1,"yanchor":"top"},"colorscale":[["0","rgba(190,190,190,1)"],["0.0416666666666667","rgba(196,188,174,1)"],["0.0833333333333333","rgba(201,185,158,1)"],["0.125","rgba(206,183,142,1)"],["0.166666666666667","rgba(212,181,126,1)"],["0.208333333333333","rgba(217,179,110,1)"],["0.25","rgba(222,177,95,1)"],["0.291666666666667","rgba(228,175,80,1)"],["0.333333333333333","rgba(233,173,63,1)"],["0.375","rgba(238,171,50,1)"],["0.416666666666667","rgba(244,169,35,1)"],["0.458333333333333","rgba(248,164,19,1)"],["0.5","rgba(251,157,11,1)"],["0.541666666666667","rgba(254,149,3,1)"],["0.583333333333333","rgba(255,137,0,1)"],["0.625","rgba(255,124,0,1)"],["0.666666666666667","rgba(255,110,0,1)"],["0.708333333333333","rgba(255,97,0,1)"],["0.75","rgba(255,83,0,1)"],["0.791666666666667","rgba(255,69,0,1)"],["0.833333333333333","rgba(255,57,0,1)"],["0.875","rgba(255,42,0,1)"],["0.916666666666667","rgba(255,30,0,1)"],["0.958333333333333","rgba(255,19,0,1)"],["1","rgba(255,0,0,1)"]],"showscale":true,"z":[[0.128569289569748,0.135822515232267,0.155198834270766,0.182368236572356,0.156828904460773,0.13201295259128],[0.05747596205554,0.0620396600566572,0.0677872183228666,0.0905559931553359,0.066765873015873,0.0582757584860318],[0.0564208691863631,0.0434782608695652,0.070516717325228,0.0847120048829066,0.0715036803364879,0.0453400503778338],[0.0600375234521576,0.0646732964025896,0.069782159749517,0.0974912397618075,0.07495900679316,0.0661616161616162],[0.0463516204961072,0.0348432055749129,0.0659090909090909,0.0752212389380531,0.0555555555555556,0.037593984962406],[0.063012226252855,0.0531746031746032,0.0837796480489671,0.109983766233766,0.0596877869605142,0.0493358633776091]],"text":["Promoter","5'UTR","1stExon","Body","ExonBnd","3'UTR"],"type":"surface","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

## Heatmap and Hierarchical Clustering of Differential Methylation (P<0.05)


```r
library(pheatmap)
library(dplyr)
##Import Data Matrix
# betaHM<-read.csv("../2_Input/EPIC.betaValues.csv", row.names = 1)
## Filters to Apply to DMR
pvalue_threshold=0.05
DMP_location="Island"
Gene_region="TSS"

##Filter Differential Methylation Data
DMR.p05<-Results %>% filter(P.Value<pvalue_threshold)
DMR.p05<-DMR.p05 %>% select(Row.names, 
                            Methylation.Diff, 
                            P.Value, 
                            adj.P.Val, 
                            Relation_to_Island, 
                            UCSC_RefGene_Group, 
                            chr, 
                            pos, 
                            contains("LVAD"))
DMR.p05.Region<-DMR.p05 %>% 
  # filter(grepl(DMP_location, Relation_to_Island)) %>%
  filter(grepl(Gene_region, UCSC_RefGene_Group)) %>% 
  select(contains("LVAD")) %>%
  data.matrix()
#Import the Index File
LVAD_Counts_Data <- targets
rownames(LVAD_Counts_Data)<-targets$Sample_Name
Index<-LVAD_Counts_Data %>% select(Race, Age, BMI)
Index<-as.data.frame(Index)
paletteLength <- 100
ann_colors = list(Race = c(African_American="#1b9e77", Asian_American = "#d95f02", Caucasian_American="#7570b3", Eastern_African = "#e7298a"))
heatmap_DMC<-pheatmap(DMR.p05.Region, scale="row", 
                      cluster_cols = TRUE, 
                      cluster_rows = TRUE,
                      #breaks = myBreaks,
                      cutree_cols = 2,
                      cutree_rows = 2,
                      angle_col = 45,
                      fontsize_col = 8,
                      color = myColor, 
                      show_rownames = FALSE, 
                      border_color = NA, 
                      annotation_col = Index,
                      annotation_colors=ann_colors,
                      filename = paste0("../3_Output/2_Methyl/", COMPARISON, "/", COMPARISON, "_", DMP_location, ".", Gene_region,"_", "ALL_SAMPLES.heatmap.pdf"))
```

<!-- ## GEO Coding -->

<!-- ```{r GEO.Coding} -->
<!-- library(dplyr) -->
<!-- library("ggplot2") -->
<!-- theme_set(theme_bw()) -->
<!-- library("sf") -->
<!-- library("rnaturalearth") -->
<!-- library("rnaturalearthdata") -->
<!-- library("tools") -->
<!-- library("maps") -->
<!-- states <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) -->
<!-- states <- cbind(states, st_coordinates(st_centroid(states))) -->
<!-- states$ID <- toTitleCase(states$ID) -->
<!-- Patient_coord <- st_as_sf(colData_all, coords = c("Longitude", "Lattitude"), remove = FALSE,  -->
<!--     crs = 4326, agr = "constant") -->
<!-- pdf(file = "../3_Output/2_Methyl/GEO_Code.pdf", width = 10, height = 10) -->
<!-- ggplot(data = states) + -->
<!--     geom_sf() + -->
<!--     geom_sf(data = Patient_coord) + -->
<!--     geom_text_repel(data = Patient_coord, aes(x = Longitude, y = Lattitude, label = LVAD_ID), size = 3) + -->
<!--       geom_point(data = colData_all, aes(x = Longitude, y = Lattitude, colour=Race, label = LVAD_ID)) + -->
<!--     coord_sf(xlim = c(-90, -82), ylim = c(30, 36), expand = FALSE) -->
<!-- dev.off() -->
<!-- ``` -->

## Volcano Plot - DMPs


```r
# Load packages
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
#
# Results<-read.xlsx(paste0("../3_Output/2_Methyl/", COMPARISON, "/", COMPARISON, "_DMPs.xlsx"), sheet = "P_0.05")
#Read data from the web
Volcano_data = mutate(Results, sig=ifelse(Results$P.Value<0.05 & abs(Results$Methylation.Diff)>10, "P<0.05 and |Methylation| > 10%", "Not Sig"), minuslogpvalue = -log(P.Value), Methylation=Methylation.Diff)
max(Volcano_data$minuslogpvalue, na.rm = TRUE)
```

```
## [1] 36.50796
```

```r
# Results<-Results %>% filter(grepl("TSS", UCSC_RefGene_Group))
#plot the ggplot
p = ggplot(Volcano_data, aes(Methylation, minuslogpvalue, color = sig)) + 
  scale_color_manual(values=c('grey','#E69F00'))+
  theme(panel.background = element_rect("white", colour = "black", size=2), 
      panel.grid.major = element_line(colour = "gray50", size=.75), 
      panel.grid.minor = element_line(colour = "gray50", size=0.4),
      legend.position="bottom",
      legend.key=element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face="bold")) +
  geom_point(size = 0.08) + 
  labs(x="Percent Methylation", 
       y=expression(-Log[10](P-value))
       ) + xlim(-75,75)+ 
  ylim(0, 28) + 
  geom_hline(yintercept = 0, size = 1) + 
  geom_vline(xintercept=0, size=1) +
  geom_text_repel(data=filter(Volcano_data, minuslogpvalue>15, abs(Methylation) > 25), aes(label=Name), show_guide = F)
#add a repelling effect to the text labels.
p
```

![](ComprehensiveAnalysis_2021_files/figure-html/Volcano-1.png)<!-- -->

```r
pdf(file = paste0("../3_Output/2_Methyl/", COMPARISON, "/", COMPARISON, "Volcano.Plot.pdf"), height = 6, width = 5)
p
dev.off()
```

```
## quartz_off_screen 
##                 2
```

# Combining DNA methylation with RNA-sequencing Analysis


```r
library(dplyr)
library(openxlsx)
library(tidyr)
STATISTIC = 0.05 #P statistic threshold used in this combination.
ifelse(!dir.exists(file.path(paste0("../3_Output/3_Combined/", COMPARISON))), dir.create(file.path(paste0("../3_Output/3_Combined/", COMPARISON))), FALSE)
```

```
## [1] FALSE
```

```r
#Import the differential Methylation data
DEGs<-read.xlsx(paste0("../3_Output/1_RNA/", COMPARISON, "/", COMPARISON, "_DESeq2.xlsx"), sheet = "P_0.05")
colnames(DEGs)<-paste0("RNA_",colnames(DEGs))
DEGs$Gene.Symbol<-DEGs$RNA_external_gene_name
DMPs<-read.xlsx(paste0("../3_Output/2_Methyl/", COMPARISON, "/", COMPARISON, "_DMPs.xlsx"), sheet = "P_0.05")
colnames(DMPs)<-paste0("Methyl_", colnames(DMPs))
##
DMPs_separated<-DMPs %>% mutate(DMPs, chromStart=Methyl_pos, chromEnd=Methyl_pos+1, Gene.Symbol = strsplit(as.character(Methyl_UCSC_RefGene_Name), ";")) %>% unnest(Gene.Symbol) %>% distinct()
DMPs_separated$Gene.Symbol<-as.character(DMPs_separated$Gene.Symbol)
write.csv(DMPs_separated, "DMPs_separated.csv")
#Merge the datasets
DMPs_DEGs<-inner_join(DEGs, DMPs_separated, by="Gene.Symbol")
DMPs_promoter.DEGs<-dplyr::filter(DMPs_DEGs, grepl("TSS", Methyl_UCSC_RefGene_Group), RNA_pvalue<0.05) %>% filter((Methyl_Methylation.Diff>0 & RNA_log2FoldChange<0) | (Methyl_Methylation.Diff<0 & RNA_log2FoldChange>0))
Up.DEGs<-dplyr::filter(DMPs_promoter.DEGs, RNA_log2FoldChange<0)
Down.DEGs<-dplyr::filter(DMPs_promoter.DEGs,  RNA_log2FoldChange>0)

#Body Promoter DEGs
DMPs_body.DEGs<-DMPs_DEGs %>% filter(grepl("Body", Methyl_UCSC_RefGene_Group), RNA_pvalue<0.05) %>% filter((Methyl_Methylation.Diff>0 & RNA_log2FoldChange>0) | (Methyl_Methylation.Diff<0 & RNA_log2FoldChange<0))

library(openxlsx)
wb_combine<-createWorkbook()
#Unfiltered
  addWorksheet(wb_combine, "Inverse Promoter")
  writeData(wb_combine, "Inverse Promoter", DMPs_promoter.DEGs, startCol = 1)
  addWorksheet(wb_combine, "Direct Body")
  writeData(wb_combine, "Direct Body", DMPs_body.DEGs, startCol = 1)
saveWorkbook(wb_combine, file = paste0("../3_Output/3_Combined/", COMPARISON,"/",COMPARISON, "_Inverse_DMPs.Promoter.xlsx"), overwrite = TRUE)
```

##Combined Heatmap


```r
library(ComplexHeatmap)
library(openxlsx)
library(dplyr)
# detach("package:pheatmap", unload=TRUE)
#Index
Index_hm<-read.xlsx("../2_Input/1_Patient/colData_CHAAMPS.FINAL.xlsx", rowNames = T, startRow = 2)
Index_hm<-Index_hm %>% 
  subset(Race %in% RACE) %>% 
  filter(rownames(.)!="LVAD054") %>%
  select(Race, Age, BMI, Cardiac.Index)

#Gene Expression
DEG_hm.prom<-DMPs_promoter.DEGs %>% select(contains("RNA_LVAD")) %>% data.matrix()
rownames(DEG_hm.prom)<-make.unique(DMPs_promoter.DEGs$RNA_external_gene_name, sep = ".")
colnames(DEG_hm.prom)<-gsub('RNA_','',colnames(DEG_hm.prom)) #Remove RNA_ prefix from column names
DEG_hm.prom<-t(scale(t(DEG_hm.prom)))

DMP_hm.prom<-DMPs_promoter.DEGs %>% select(contains("Methyl_LVAD")) %>% data.matrix()
rownames(DMP_hm.prom)<-make.unique(DMPs_promoter.DEGs$RNA_external_gene_name, sep = ".")
colnames(DMP_hm.prom)<-gsub('Methyl_','',colnames(DMP_hm.prom)) #Remove Methyl_ prefix from column names
DMP_hm.prom<-t(scale(t(DMP_hm.prom)))

#Make heatmap for DMPs
paletteLength <- 100
myColor_DMC <- colorRampPalette(c("cyan4", "white", "firebrick2"))(paletteLength)
ann_colors = list(Race = c(African_American="#1b9e77", Asian_American = "#d95f02", Caucasian_American="#7570b3", Eastern_African = "#e7298a"),
                  BMI = brewer.pal(6, "Blues"),
                  Age = brewer.pal(6, "Purples"),
                  Cardiac.Index = brewer.pal(6, "Reds")
                    )
heatmap_DMC<-pheatmap(DMP_hm.prom, scale="none", 
                      cluster_cols = TRUE, 
                      cluster_rows = TRUE,
                      cutree_rows = 2,
                      fontsize_col = 8,
                      color = myColor_DMC, 
                      clustering_distance_cols = "correlation",
                      show_colnames = F,
                      show_rownames = FALSE, 
                      border_color = NA, 
                      annotation_col = Index_hm,
                      annotation_colors=ann_colors,
                      annotation_legend = F,
                      filename = paste0("../3_Output/3_Combined/", COMPARISON, "/", COMPARISON, "_", DMP_location, ".", Gene_region,"_", "DMP.promoters.heatmap.pdf"))

#DEG Heatmap
colnames_ordered<-colnames(DMP_hm.prom)[column_order(heatmap_DMC)] #Get column order from the DMC Heatmap
Testing<-DEG_hm.prom[,match(colnames_ordered, colnames(DEG_hm.prom))]
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "gold2"))(paletteLength) #DEG Colors
heatmap_DEG<-pheatmap(Testing, scale="row", 
                      cluster_cols = T, 
                      cluster_rows = TRUE,
                      cutree_rows = 2,
                      fontsize_col = 8,
                      color = myColor, 
                      clustering_distance_cols = "correlation",
                      show_colnames = F,
                      show_rownames = FALSE, 
                      border_color = NA, 
                      annotation_col = Index_hm,
                      annotation_colors=ann_colors,
                      filename = paste0("../3_Output/3_Combined/", COMPARISON, "/", COMPARISON, "_", DMP_location, ".", Gene_region,"_", "DEG.promoters.heatmap.pdf"))

pdf(file = paste0("../3_Output/3_Combined/", 
                  COMPARISON, "/", 
                  COMPARISON, "_", 
                  DMP_location, ".", 
                  Gene_region,"_", 
                  "DMC.DEG_Overlapping.Heatmap.pdf"), 
    height = 8, 
    width = 10
    )
heatmap_DMC + heatmap_DEG
dev.off()
```

```
## quartz_off_screen 
##                 2
```

##Dendrogram Link


```r
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms

##Methylation
res.dist_DMP <- dist(t(DMP_hm.prom), method = "euclidean")
hc1 <- hclust(res.dist_DMP, method = "ward.D2")
dend1 <- as.dendrogram (hc1)
colors_to_use <- as.numeric(targets_filtered$Race)
colors_to_use <- colors_to_use[order.dendrogram(dend1)]
labels_colors(dend1)<-colors_to_use

##DEG
res.dist_DEG <- dist(t(DEG_hm.prom), method = "euclidean")
hc2 <- hclust(res.dist_DEG, method = "ward.D2")
dend2 <- as.dendrogram (hc2)
colors_to_use2 <- as.numeric(targets_filtered$Race)
colors_to_use2 <- colors_to_use2[order.dendrogram(dend2)]
labels_colors(dend2)<-colors_to_use2
#
dend_list <- dendlist(dend1, dend2)
#
dendTangle<-dend_list %>% untangle(method = "step1side")
##

Test2<-tanglegram(dendTangle, sort = FALSE, common_subtrees_color_lines = TRUE, highlight_distinct_edges  = TRUE, highlight_branches_lwd = TRUE, lwd = 2, edge.lwd = 2, lab.cex=.7)
```

![](ComprehensiveAnalysis_2021_files/figure-html/dendrogram.link-1.png)<!-- -->

```r
pdf(file = "../3_Output/3_Combined/TANGLE.pdf", height = 8, width = 10)
tanglegram(dendTangle, sort = FALSE, common_subtrees_color_lines = TRUE, highlight_distinct_edges  = TRUE, highlight_branches_lwd = TRUE, lwd = 2, edge.lwd = 2, lab.cex=.7)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
entanglement(Test2)
```

```
## [1] 0.1390582
```

```r
cors <- cor.dendlist(Test2)

ha_index<-select(targets_filtered, Race, HbA1C, BMI)
rownames(ha_index)<-targets_filtered$Sample_Name
ha = HeatmapAnnotation(df = ha_index)

# Print correlation matrix
library(ComplexHeatmap)
p1=ComplexHeatmap::pheatmap(DMP_hm.prom, scale = "row")
p2=Heatmap(scale(DEG_hm.prom), name = "RNA-Seq", km = 2, top_annotation = ha)

p1+p2
```

##Circos Plot


```r
library(dplyr)
#create gene labels
Gene_labels_Body<-DMPs_body.DEGs%>% dplyr::filter(abs(RNA_log2FoldChange)>.585, grepl("Body", Methyl_UCSC_RefGene_Group)) %>% dplyr::filter(abs(Methyl_Methylation.Diff)>5) %>% dplyr::select(chrom=Methyl_chr, chromStart=RNA_start_position, chromEnd=RNA_end_position, GeneSymbol=RNA_external_gene_name, RNA_padj) %>% distinct() %>% dplyr::top_n(50, RNA_padj) %>% select(-RNA_padj)
Gene_labels_Body<-arrange(Gene_labels_Body, chromStart)
Gene_labels_Body$chrom<-factor(Gene_labels_Body$chrom, levels=c("chr1", "chr2", "chr3", "chr4", 
                                                      "chr5", "chr6", "chr7", "chr8", 
                                                      "chr9", "chr10", "chr11", "chr12", 
                                                      "chr13", "chr14", "chr15", "chr16", 
                                                      "chr17", "chr18", "chr19", "chr20", 
                                                      "chr21", "chr22", "chr23", "chrX", 
                                                      "chrY"))
Gene_labels_Body<-Gene_labels_Body[order(Gene_labels_Body$chrom),]
Gene_labels_Body<-Gene_labels_Body[!duplicated(Gene_labels_Body[,4]),]
Gene_labels_Body$Color<-"gray80" #Color the text of the corresponding circos
# Gene Labels based on Genomic Region
Gene_labels_Promoter<-DMPs_promoter.DEGs %>% dplyr::filter(abs(RNA_log2FoldChange)>.585, grepl("TSS", Methyl_UCSC_RefGene_Group)) %>% dplyr::filter(abs(Methyl_Methylation.Diff)>5) %>% dplyr::select(chrom=Methyl_chr, chromStart=RNA_start_position, chromEnd=RNA_end_position, GeneSymbol=RNA_external_gene_name, RNA_padj) %>% distinct() %>% dplyr::top_n(50, RNA_padj) %>% select(-RNA_padj)
Gene_labels_Promoter<-arrange(Gene_labels_Promoter, chromStart)
Gene_labels_Promoter$chrom<-factor(Gene_labels_Promoter$chrom, levels=c("chr1", "chr2", "chr3", "chr4", 
                                                      "chr5", "chr6", "chr7", "chr8", 
                                                      "chr9", "chr10", "chr11", "chr12", 
                                                      "chr13", "chr14", "chr15", "chr16", 
                                                      "chr17", "chr18", "chr19", "chr20", 
                                                      "chr21", "chr22", "chr23", "chrX", 
                                                      "chrY"))
Gene_labels_Promoter<-Gene_labels_Promoter[order(Gene_labels_Promoter$chrom),]
Gene_labels_Promoter<-Gene_labels_Promoter[!duplicated(Gene_labels_Promoter[,4]),]
Gene_labels_Promoter$Color<-"gray40" #Color the text of the corresponding circos
# Create a composite Labels frame with coloring based on CpG Region
Gene.labels<-rbind(Gene_labels_Body, Gene_labels_Promoter)
Gene.labels<-Gene.labels[order(Gene.labels$chrom),]


##Gene Expression
##Fold Change UP
Gene_FoldChange.UP<-dplyr::select(DMPs_DEGs, chrom=Methyl_chr, 
                                                      chromStart, chromEnd, FoldChange_DEG=RNA_log2FoldChange)
Gene_FoldChange.UP<-dplyr::filter(Gene_FoldChange.UP, FoldChange_DEG>0)
# Gene_FoldChange.UP<-dplyr::mutate(Gene_FoldChange.UP, chromEnd=chromStart+1)
Gene_FoldChange.UP<-dplyr::select(Gene_FoldChange.UP, chrom, chromStart, chromEnd, FoldChange_DEG)
Gene_FoldChange.UP<-arrange(Gene_FoldChange.UP, chromStart)
Gene_FoldChange.UP$chrom<-factor(Gene_FoldChange.UP$chrom, levels=c("chr1", "chr2", "chr3", "chr4", 
                                                      "chr5", "chr6", "chr7", "chr8", 
                                                      "chr9", "chr10", "chr11", "chr12", 
                                                      "chr13", "chr14", "chr15", "chr16", 
                                                      "chr17", "chr18", "chr19", "chr20", 
                                                      "chr21", "chr22", "chr23", "chrX", 
                                                      "chrY"))
Gene_FoldChange.UP<-Gene_FoldChange.UP[order(Gene_FoldChange.UP$chrom),]
##Fold Change DOWN
Gene_FoldChange.DOWN<-dplyr::select(DMPs_DEGs, chrom=Methyl_chr, 
                                                      chromStart, chromEnd, FoldChange_DEG=RNA_log2FoldChange)
Gene_FoldChange.DOWN<-dplyr::filter(Gene_FoldChange.DOWN, FoldChange_DEG<0)
# Gene_FoldChange.DOWN<-dplyr::mutate(Gene_FoldChange.DOWN, chromEnd=chromStart+1)
Gene_FoldChange.DOWN<-dplyr::select(Gene_FoldChange.DOWN, chrom, chromStart, chromEnd, FoldChange_DEG)
Gene_FoldChange.DOWN<-arrange(Gene_FoldChange.DOWN, chromStart)
Gene_FoldChange.DOWN$chrom<-factor(Gene_FoldChange.DOWN$chrom, levels=c("chr1", "chr2", "chr3", "chr4", 
                                                      "chr5", "chr6", "chr7", "chr8", 
                                                      "chr9", "chr10", "chr11", "chr12", 
                                                      "chr13", "chr14", "chr15", "chr16", 
                                                      "chr17", "chr18", "chr19", "chr20", 
                                                      "chr21", "chr22", "chr23", "chrX", 
                                                      "chrY"))
Gene_FoldChange.DOWN<-Gene_FoldChange.DOWN[order(Gene_FoldChange.DOWN$chrom),]
##Fold Change List
Gene_FoldChange_List<-list(Gene_FoldChange.UP, Gene_FoldChange.DOWN)

# Methylation Density - Islands
DMR_Islands<-DMPs %>% dplyr::filter(Methyl_Relation_to_Island=="Island") %>% select(chrom=Methyl_chr, Methyl_pos, perc.change=Methyl_Methylation.Diff)
DMR_Islands<-dplyr::filter(DMR_Islands, chrom!="chrM")
DMR_Islands<-dplyr::mutate(DMR_Islands, chromEnd=Methyl_pos+1) %>% select(chrom, chromStart=Methyl_pos, chromEnd, perc.change)
DMR_Islands$chrom<-factor(DMR_Islands$chrom, levels=c("chr1", "chr2", "chr3", "chr4", 
                                                      "chr5", "chr6", "chr7", "chr8", 
                                                      "chr9", "chr10", "chr11", "chr12", 
                                                      "chr13", "chr14", "chr15", "chr16", 
                                                      "chr17", "chr18", "chr19", "chr20", 
                                                      "chr21", "chr22", "chr23", "chrX", 
                                                      "chrY"))
DMR_Islands<-DMR_Islands[order(DMR_Islands$chrom),]
Methyl.UP.Island<-filter(DMR_Islands, perc.change>0)
Methyl.DOWN.Island<-filter(DMR_Islands, perc.change<0)
Methyl.List.island<-list(Methyl.DOWN.Island, Methyl.UP.Island)

# Methylation Density - Seas
DMR_Sea<-DMPs %>% dplyr::filter(Methyl_Relation_to_Island=="OpenSea") %>% select(chrom=Methyl_chr, Methyl_pos, perc.change=Methyl_Methylation.Diff)
DMR_Sea<-dplyr::filter(DMR_Sea, chrom!="chrM")
DMR_Sea<-dplyr::mutate(DMR_Sea, chromEnd=Methyl_pos+1) %>% select(chrom, chromStart=Methyl_pos, chromEnd, perc.change)
DMR_Sea$chrom<-factor(DMR_Sea$chrom, levels=c("chr1", "chr2", "chr3", "chr4", 
                                                      "chr5", "chr6", "chr7", "chr8", 
                                                      "chr9", "chr10", "chr11", "chr12", 
                                                      "chr13", "chr14", "chr15", "chr16", 
                                                      "chr17", "chr18", "chr19", "chr20", 
                                                      "chr21", "chr22", "chr23", "chrX", 
                                                      "chrY"))
DMR_Sea<-DMR_Sea[order(DMR_Sea$chrom),]
Methyl.UP.Sea<-filter(DMR_Sea, perc.change>0)
Methyl.DOWN.Sea<-filter(DMR_Sea, perc.change<0)
Methyl.List.Sea<-list(Methyl.DOWN.Sea, Methyl.UP.Sea)
# All DMPs
CpG_EPIC<-DMPs %>% dplyr::select(chrom=Methyl_chr, Methyl_pos, perc.change=Methyl_Methylation.Diff)
CpG_EPIC<-dplyr::filter(CpG_EPIC, chrom!="chrM")
CpG_EPIC<-dplyr::mutate(CpG_EPIC, chromEnd=Methyl_pos+1) %>% select(chrom, chromStart=Methyl_pos, chromEnd, perc.change)
CpG_EPIC$chrom<-factor(CpG_EPIC$chrom, levels=c("chr1", "chr2", "chr3", "chr4", 
                                                      "chr5", "chr6", "chr7", "chr8", 
                                                      "chr9", "chr10", "chr11", "chr12", 
                                                      "chr13", "chr14", "chr15", "chr16", 
                                                      "chr17", "chr18", "chr19", "chr20", 
                                                      "chr21", "chr22", "chr23", "chrX", 
                                                      "chrY"))
CpG_EPIC<-CpG_EPIC[order(CpG_EPIC$chrom),]
Methyl.UP.EPIC<-filter(CpG_EPIC, perc.change>0)
Methyl.DOWN.EPIC<-filter(CpG_EPIC, perc.change<0)
Methyl.List.EPIC<-list(Methyl.DOWN.EPIC, Methyl.UP.EPIC)

Methyl.List<-list(DMR_Sea,DMR_Islands)
#Plot the Circos
library(circlize)
library(gtools)
library(dplyr)

circos.genomicDensity1 = function (data, ylim.force = FALSE, window.size = NULL, overlap = TRUE, col = ifelse(area, "grey", "black"), lwd = par("lwd"), lty = par("lty"), type = "l", area = TRUE, area.baseline = NULL, baseline = 0, border = NA, ...) { if (!is.null(area.baseline)) 
data = normalizeToDataFrame(data)
if (!is.dataFrameList(data)) {
data = list(data)
}
if (length(col) == 1) {
col = rep(col, length(data))
}
if (length(lwd) == 1) {
lwd = rep(lwd, length(data))
}
if (length(lty) == 1) {
lty = rep(lty, length(data))
}
if (length(type) == 1) {
type = rep(type, length(data))
}
if (length(area) == 1) {
area = rep(area, length(data))
}

if (length(baseline) == 1) {
    baseline = rep(baseline, length(data))
}
if (length(border) == 1) {
    border = rep(border, length(data))
}
s = sapply(get.all.sector.index(), function(si) get.cell.meta.data("xrange", 
    sector.index = si))

if (is.null(window.size)) {
    window.size = 10^nchar(sum(s))/1000
}
df = vector("list", length = length(data))
for (i in seq_along(data)) {
    all.chr = unique(data[[i]][[1]])
    for (chr in all.chr) {
        region = data[[i]][data[[i]][[1]] == chr, 2:3, drop = FALSE]
        dn = genomicDensity(region, window.size = window.size, 
            overlap = overlap)
        dn = cbind(rep(chr, nrow(dn)), dn)
        df[[i]] = rbind(df[[i]], dn)
    }
}
if (ylim.force) {
    ymax = 1
}
else {
    ymax = max(sapply(df, function(gr) max(gr[[4]])))
}

circos.genomicTrackPlotRegion(df, ylim = c(-ymax,0), panel.fun = function(region, 
    value, ...) {
    i = getI(...)

    circos.genomicLines(region, -value, col = col[i], lwd = lwd[i], 
        lty = lty[i], type = type[i], border = border[i], 
        area = area[i], baseline = baseline[i])
}, ...)
}

environment(circos.genomicDensity1) <- asNamespace('circlize')

#to get error line number:

f <- function (data, ylim.force = FALSE, window.size = NULL, overlap = TRUE,
col = ifelse(area, "grey", "black"), lwd = par("lwd"), lty = par("lty"),
type = "l", area = TRUE, area.baseline = NULL, baseline = 0,
border = NA, ...)
{
circos.genomicDensity1(data, ylim.force = FALSE, window.size = NULL, overlap = TRUE,
col = ifelse(area, "grey", "black"), lwd = par("lwd"), lty = par("lty"),
type = "l", area = TRUE, area.baseline = NULL, baseline = 0,
border = NA, ...)
}
###########################################################
# pdf(file="trial_2.pdf",width=10,height=10)
# par(mar = c(1, 1, 1, 1))
# circos.par(gap.degree=1, track.margin = c(0, 0))
# circos.initializeWithIdeogram(track.height = 0.05)
# circos.genomicDensity(Methyl.UP, col = c("goldenrod3"), track.height = 0.1, baseline="bottom", bg.border ="white", track.margin = c(0, 0.0))
# circos.genomicDensity1(Methyl.DOWN, col = c("dodgerblue3"), track.height = 0.1, baseline="top", bg.border ="white", track.margin = c(0, 0.0))
# circos.clear()
# dev.off()
######################################################

om = circos.par("track.margin")
oc = circos.par("cell.padding")
circos.par(track.margin = c(0, 0), cell.padding = c(0, 0, 0, 0))
circos.par(start.degree = -250)
pdf(file=paste0("../3_Output/3_Combined/", COMPARISON, "/", COMPARISON, "_Circos.pdf"))
circos.initializeWithIdeogram(track.height = 0.05)

circos.genomicDensity(Methyl.List, col=add_transparency(c("gray80", "gray40"), transparency = 0.2) ,
                      track.height=0.2, bg.border=NA)
# ### Labels for DMRs in Islands
# circos.genomicDensity(Methyl.UP.Island, col = c("coral2"), track.height = 0.1, baseline="bottom", bg.border ="white", track.margin = c(0, 0.0))
# circos.genomicDensity1(Methyl.DOWN.Island, col = c("darkcyan"), track.height = 0.1, baseline="top", bg.border ="white", track.margin = c(0, 0.0))
# ### Labels for DMRs in Body
# circos.genomicDensity(Methyl.UP.Body, col = c("coral2"), track.height = 0.1, baseline="bottom", bg.border ="white", track.margin = c(0, 0.0))
# circos.genomicDensity1(Methyl.DOWN.Body, col = c("darkcyan"), track.height = 0.1, baseline="top", bg.border ="white", track.margin = c(0, 0.0))
#DEG with inverse GPI Islands Promoters
circos.genomicTrackPlotRegion(Gene_FoldChange_List,
                              ylim = c(-6, 6), bg.border=NA,
                              panel.fun = function(region, value, ...) {
 col = ifelse(value[[1]] > 0, "darkgoldenrod1", "dodgerblue2")
 circos.genomicPoints(region, value, col = add_transparency(col, 0.2), cex = 0.3, pch = 16)
 cell.xlim = get.cell.meta.data("cell.xlim")
 for(h in c(-3, -1.5, 0, 1.5, 3)) {
   circos.lines(cell.xlim, c(h, h), col ="#00000040")
 }
}, track.height = 0.2)

circos.genomicLabels(Gene.labels, labels.column=4, side='inside', col = Gene.labels$Color ,cex=0.6)
# circos.par(track.margin=om, cell.padding=oc)
# ## Add link for all DEGs with DMRs in promoter CGIs
# Link_Anchor <- read.csv("../1_Input/4_Combined/Circos/Link_Anchor.csv")
# Link<-read.csv("../1_Input/4_Combined/Circos/Link_DEG.DMR_Promoter.CGI_P<0.05.csv")
# Link$chrom<-factor(Link$chrom, levels=c("chr1", "chr2", "chr3", "chr4", 
#                                                       "chr5", "chr6", "chr7", "chr8", 
#                                                       "chr9", "chr10", "chr11", "chr12", 
#                                                       "chr13", "chr14", "chr15", "chr16", 
#                                                       "chr17", "chr18", "chr19", "chr20", 
#                                                       "chr21", "chr22", "chr23", "chrX", 
#                                                       "chrY"))
# Link<-Link[order(Link$chrom),]
# Link_Anchor<-Link_Anchor[1:nrow(Link),]
# circos.genomicLink(Link, Link_Anchor, col="black", lwd=0.5)
circos.clear()
dev.off()
```

```
## quartz_off_screen 
##                 2
```

##Ridgeline Plot


```r
# library
library(ggridges)
library(ggplot2)
library(dplyr)
DMPs<-openxlsx::read.xlsx(paste0("../3_Output/2_Methyl/", COMPARISON, "/", COMPARISON, "_DMPs.xlsx"), sheet = "Unfiltered")
colData<-openxlsx::read.xlsx("../2_Input/1_Patient/colData_CHAAMPS.FINAL.xlsx", sheet = "Summary", startRow=2)
GENE="POSTN"
DMPs_separated<-DMPs %>% mutate(DMPs, chromStart=pos, chromEnd=pos+1, Gene.Symbol = strsplit(as.character(UCSC_RefGene_Name), ";")) %>% unnest(Gene.Symbol) %>% distinct()
DMPs_separated$Gene.Symbol<-as.character(DMPs_separated$Gene.Symbol)
Gene_DMP<-subset(DMPs_separated, Gene.Symbol %in% GENE) 
Gene_DMP<-dplyr::select(Gene_DMP, chr,chromStart, chromEnd, Gene.Symbol, contains("LVAD"))
colnames(Gene_DMP)<-gsub("Methyl_", "", colnames(Gene_DMP))
gathered<-tidyr::gather(Gene_DMP, "Sample", "Methylation", 5:length(colnames(Gene_DMP)))
gathered_annot<-merge(colData, gathered, by.x = "LVAD_ID", by.y = "Sample")
gathered_annot$Methylation<-as.numeric(as.character(gathered_annot$Methylation))
gathered_annot$Race<-factor(gathered_annot$Race, levels = c("Causasian_American", "African_American", "Asian_American", "Eastern_African"))
gathered_annot$LVAD_ID<-factor(gathered_annot$LVAD_ID)
# basic example
pdf(file=paste0("../3_Output/3_Combined/", COMPARISON, "/", COMPARISON, "_", GENE, "_Methylation_gene.distribution.pdf"))
ggplot(gathered_annot, aes(x=chromStart, y=Methylation, group = LVAD_ID, color=Race))+
  theme_bw()+
  geom_line(size=1)+
  geom_point(size = 2, alpha = .6)+
  ggtitle(paste0("Methylation Distribution - ", GENE))
dev.off()
```

```
## quartz_off_screen 
##                 2
```

<!-- ## Supervised Machine Learning: Identifying a racial signature of cardiac DNA methylation -->

<!-- ```{r bglm} -->
<!-- library(BhGLM) -->
<!-- library(openxlsx) -->
<!-- library(tidyr) -->
<!-- library(dplyr) -->
<!-- CpGs<-read.csv("../2_Input/beta.table.csv", row.names=1) -->
<!-- tCpGs<-t(CpGs) -->
<!-- Annotation<-read.xlsx("../2_Input/1_Patient/Patient_Data.xlsx", startRow = 2) -->
<!-- rownames(Annotation)<-Annotation$LVAD_ID -->
<!-- Annotation<-subset(Annotation, Study=="CHAAMPS") -->
<!-- #Convert binary to boolean -->
<!-- Annotation$Race<-factor(Annotation$Race, levels = c("Caucasian_American", "African_American")) -->
<!-- Annotation<-dplyr::select(Annotation, LVAD_ID, Race, Age, Diabetes:RA) -->
<!-- Annotation$Diabetes<-factor(Annotation$Diabetes, levels = c("ND", "T2D")) -->
<!-- Annotation$Ischemia<-factor(Annotation$Ischemia, levels = c("NICM", "ICM")) -->
<!-- Annotation$HTN<-factor(Annotation$HTN, levels = c("Y", "N")) -->
<!-- Annotation$Obese<-factor(Annotation$Obese, levels = c("N", "Y")) -->
<!-- Annotation$Lipid.Lowering<-factor(Annotation$Lipid.Lowering, levels = c("0", "1")) -->
<!-- Annotation$ACE<-factor(Annotation$ACE, levels = c("0", "1")) -->
<!-- Annotation$B.Blocker<-factor(Annotation$B.Blocker, levels = c("0", "1")) -->
<!-- Annotation$Diuretic<-factor(Annotation$Diuretic, levels = c("0", "1")) -->
<!-- Annotation$`PDE3-Ant`<-factor(Annotation$`PDE3-Ant`, levels = c("0", "1")) -->
<!-- Annotation$`Adrenergic-Ag`<-factor(Annotation$`Adrenergic-Ag`, levels = c("0", "1")) -->
<!-- Annotation$Family.Hx<-factor(Annotation$Family.Hx, levels = c("0", "1")) -->
<!-- # Standardize only the continuous variables -->
<!-- Annotation_rescale <- Annotation %>% -->
<!-- 	mutate_if(is.numeric, funs(as.numeric(scale(.)))) -->
<!-- # Select the Training set using random selection of 50% sample population -->
<!-- Train_set<-dplyr::sample_frac(Annotation_rescale, 0.5) #randomly sample 50% of the samples for the training dataset. -->
<!-- rownames(Train_set)<-Train_set$LVAD_ID -->
<!-- Train_set<-select(Train_set, -LVAD_ID, -Diuretic, -Syst.PA, -PCWP, -RA) -->
<!-- # Perform GLM on the training set -->
<!-- Train_glm<-glm(Diabetes ~ ., data = Train_set, family = "binomial", na.action="na.exclude") -->
<!-- summary(Train_glm) -->
<!-- rankP <- predict(Train_glm, newdata = Annotation, type = "response") -->

<!-- Annotation<-dplyr::select(Annotation, Race) -->
<!-- #Join the EPIC data with the patient data -->
<!-- Merged<-merge(Annotation, tCpGs, by=0) -->
<!-- rownames(Merged)<-Merged$Row.names -->
<!-- Merged<-Merged[,!(names(Merged) %in% "Row.names")] -->
<!-- #Perform GLM -->
<!-- Test<-bglm(Race ~ ., data = Merged, family = "binomial", prior = mde(0, 0.04, 0.5), method.coef = 50) -->
<!-- ``` -->


<!-- #Combined Analysis: DEGs and DMPs -->
<!-- ##Correlation between DMRs and DEGs -->

<!-- ```{r Correlation.DMRs.DEGs} -->
<!-- library(openxlsx) -->
<!-- library(dplyr) -->
<!-- library(tidyr) -->
<!-- library(Hmisc) -->
<!-- COMPARISON = "African_American.vs.Caucasian_American_ALL" -->
<!-- #Import the differential data -->
<!-- RNA.Counts<-read.xlsx(paste0("../3_Output/1_RNA/", COMPARISON, "/", COMPARISON, "_DESeq2.xlsx"), sheet = "Unfiltered") -->
<!-- colnames(RNA.Counts)<-paste0("RNA_",colnames(RNA.Counts)) -->
<!-- CpGs<-read.xlsx(paste0("../3_Output/2_Methyl/", COMPARISON, "/", COMPARISON, "_DMPs.xlsx"), sheet = "Unfiltered") -->
<!-- colnames(CpGs)<-paste0("Methyl_", colnames(CpGs)) -->

<!-- # align methyl and RNAs, and separate into different tables -->
<!-- idx2 <- sapply(RNA.Counts$RNA_external_gene_name, grep, CpGs$Methyl_UCSC_RefGene_Name) -->
<!-- idx1 <- sapply(seq_along(idx2), function(i) rep(i, length(idx2[[i]]))) -->
<!-- test<-cbind(RNA.Counts[unlist(idx1),,drop=F], CpGs[unlist(idx2),,drop=F]) -->
<!-- rownames(test)<-make.unique(test$Methyl_Row.names.y) -->
<!-- RNAs<-test %>% select(contains("RNA_LVAD")) %>% select(-contains(".y"), -contains(".x")) -->
<!-- RNA_cor<-RNAs -->
<!-- Methyls<-test %>% select(contains("Methyl_LVAD")) %>% select(-contains(".y"), -contains(".x")) -->
<!-- Methyl_names<-gsub(pattern = "Methyl_", replacement = "", x  = colnames(Methyls)) -->
<!-- Methyl_cor<-Methyls -->
<!-- colnames(Methyl_cor)<-Methyl_names -->
<!-- RNA_names<-gsub(pattern = "RNA_", replacement = "", x  = colnames(RNAs)) -->
<!-- colnames(RNA_cor)<-RNA_names -->
<!-- #Ensure that the order and columns match exactly -->
<!-- col2keep = intersect(colnames(RNA_cor), colnames(Methyl_cor)) -->
<!-- RNAs = RNA_cor[, (match(col2keep, colnames(Methyl_cor)))] -->
<!-- # -->
<!-- tRNAs_norm<-scale(t(RNA_cor)) -->
<!-- tMethyls_norm<-scale(t(Methyl_cor)) -->
<!-- corr.spearman_p<-0 -->
<!-- corr.spearman_r<-0 -->
<!-- for (i in 1:ncol(tMethyls)) { -->
<!--   corr.test = rcorr(tRNAs_norm[,i], tMethyls_norm[,i], type = "pearson") -->
<!--   corr.spearman_p[i] = corr.test$P[1,2] -->
<!--   corr.spearman_r[i] = corr.test$r[1,2] -->
<!-- } -->
<!-- Final<-cbind(test, corr.spearman_p, corr.spearman_r) %>% select(-contains(".x"), -contains(".y")) %>% distinct() %>% arrange(corr.spearman_p) -->
<!-- Final_corr.sig<-Final %>% filter(corr.spearman_p<0.05) -->
<!-- Final_all.sig<-Final %>% filter(corr.spearman_p<0.1, RNA_pvalue<0.05, Methyl_pval<0.05) -->
<!-- # -->
<!-- wb_correlation<-createWorkbook() -->
<!-- #export correlation data -->
<!--   addWorksheet(wb_correlation, "Unfiltered") -->
<!--   writeData(wb_correlation, "Unfiltered", Final, startCol = 1) -->
<!--   addWorksheet(wb_correlation, "All Significant") -->
<!--   writeData(wb_correlation, "All Significant", Final_all.sig, startCol = 1) -->
<!--   addWorksheet(wb_correlation, "Correlation Significance") -->
<!--   writeData(wb_correlation, "Correlation Significance", Final_corr.sig, startCol = 1) -->
<!-- saveWorkbook(wb_correlation, file = paste0("../3_Output/3_Combined/", COMPARISON,"/",COMPARISON, "_correlation.xlsx"), overwrite = TRUE) -->

<!-- # Correlation Heatmap -->
<!-- library(corrplot) -->
<!-- M<-cor(tRNAs_norm, method = "spearman") -->
<!-- corrplot(M, ) -->
<!-- ``` -->

#Supplemental Table: R Session Information

All packages and setting are acquired using the following command: 


```r
options(kableExtra.latex.load_packages = FALSE)
library(kableExtra)
sinfo<-devtools::session_info()
sinfo$platform
```

```
##  setting  value                       
##  version  R version 4.0.3 (2020-10-10)
##  os       macOS Big Sur 10.16         
##  system   x86_64, darwin17.0          
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  ctype    en_US.UTF-8                 
##  tz       Europe/Berlin               
##  date     2021-01-19
```

```r
sinfo$packages %>% kable( 
                         align="c", 
                         longtable=T, 
                         booktabs=T,
                         caption="Packages and Required Dependencies") %>% 
    kable_styling(latex_options=c("striped", "repeat_header", "condensed"))
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Packages and Required Dependencies</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:center;"> package </th>
   <th style="text-align:center;"> ondiskversion </th>
   <th style="text-align:center;"> loadedversion </th>
   <th style="text-align:center;"> path </th>
   <th style="text-align:center;"> loadedpath </th>
   <th style="text-align:center;"> attached </th>
   <th style="text-align:center;"> is_base </th>
   <th style="text-align:center;"> date </th>
   <th style="text-align:center;"> source </th>
   <th style="text-align:center;"> md5ok </th>
   <th style="text-align:center;"> library </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> abind </td>
   <td style="text-align:center;"> abind </td>
   <td style="text-align:center;"> 1.4.5 </td>
   <td style="text-align:center;"> 1.4-5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/abind </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/abind </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2016-07-21 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> annotate </td>
   <td style="text-align:center;"> annotate </td>
   <td style="text-align:center;"> 1.68.0 </td>
   <td style="text-align:center;"> 1.68.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/annotate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/annotate </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AnnotationDbi </td>
   <td style="text-align:center;"> AnnotationDbi </td>
   <td style="text-align:center;"> 1.52.0 </td>
   <td style="text-align:center;"> 1.52.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/AnnotationDbi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/AnnotationDbi </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AnnotationFilter </td>
   <td style="text-align:center;"> AnnotationFilter </td>
   <td style="text-align:center;"> 1.14.0 </td>
   <td style="text-align:center;"> 1.14.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/AnnotationFilter </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/AnnotationFilter </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AnnotationHub </td>
   <td style="text-align:center;"> AnnotationHub </td>
   <td style="text-align:center;"> 2.22.0 </td>
   <td style="text-align:center;"> 2.22.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/AnnotationHub </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/AnnotationHub </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> askpass </td>
   <td style="text-align:center;"> askpass </td>
   <td style="text-align:center;"> 1.1 </td>
   <td style="text-align:center;"> 1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/askpass </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/askpass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-01-13 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> assertthat </td>
   <td style="text-align:center;"> assertthat </td>
   <td style="text-align:center;"> 0.2.1 </td>
   <td style="text-align:center;"> 0.2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/assertthat </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/assertthat </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-03-21 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> backports </td>
   <td style="text-align:center;"> backports </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/backports </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/backports </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-09 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> base64 </td>
   <td style="text-align:center;"> base64 </td>
   <td style="text-align:center;"> 2.0 </td>
   <td style="text-align:center;"> 2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/base64 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/base64 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2016-05-10 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> base64enc </td>
   <td style="text-align:center;"> base64enc </td>
   <td style="text-align:center;"> 0.1.3 </td>
   <td style="text-align:center;"> 0.1-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/base64enc </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/base64enc </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2015-07-28 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> beanplot </td>
   <td style="text-align:center;"> beanplot </td>
   <td style="text-align:center;"> 1.2 </td>
   <td style="text-align:center;"> 1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/beanplot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/beanplot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2014-09-19 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BiasedUrn </td>
   <td style="text-align:center;"> BiasedUrn </td>
   <td style="text-align:center;"> 1.7 </td>
   <td style="text-align:center;"> 1.07 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/BiasedUrn </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/BiasedUrn </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2015-12-28 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Biobase </td>
   <td style="text-align:center;"> Biobase </td>
   <td style="text-align:center;"> 2.50.0 </td>
   <td style="text-align:center;"> 2.50.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Biobase </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Biobase </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BiocFileCache </td>
   <td style="text-align:center;"> BiocFileCache </td>
   <td style="text-align:center;"> 1.14.0 </td>
   <td style="text-align:center;"> 1.14.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/BiocFileCache </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/BiocFileCache </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BiocGenerics </td>
   <td style="text-align:center;"> BiocGenerics </td>
   <td style="text-align:center;"> 0.36.0 </td>
   <td style="text-align:center;"> 0.36.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/BiocGenerics </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/BiocGenerics </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BiocManager </td>
   <td style="text-align:center;"> BiocManager </td>
   <td style="text-align:center;"> 1.30.10 </td>
   <td style="text-align:center;"> 1.30.10 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/BiocManager </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/BiocManager </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-11-16 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BiocParallel </td>
   <td style="text-align:center;"> BiocParallel </td>
   <td style="text-align:center;"> 1.24.1 </td>
   <td style="text-align:center;"> 1.24.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/BiocParallel </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/BiocParallel </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-06 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BiocVersion </td>
   <td style="text-align:center;"> BiocVersion </td>
   <td style="text-align:center;"> 3.12.0 </td>
   <td style="text-align:center;"> 3.12.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/BiocVersion </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/BiocVersion </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-05-14 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> biomaRt </td>
   <td style="text-align:center;"> biomaRt </td>
   <td style="text-align:center;"> 2.46.0 </td>
   <td style="text-align:center;"> 2.46.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/biomaRt </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/biomaRt </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Biostrings </td>
   <td style="text-align:center;"> Biostrings </td>
   <td style="text-align:center;"> 2.58.0 </td>
   <td style="text-align:center;"> 2.58.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Biostrings </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Biostrings </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> biovizBase </td>
   <td style="text-align:center;"> biovizBase </td>
   <td style="text-align:center;"> 1.38.0 </td>
   <td style="text-align:center;"> 1.38.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/biovizBase </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/biovizBase </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> bit </td>
   <td style="text-align:center;"> bit </td>
   <td style="text-align:center;"> 4.0.4 </td>
   <td style="text-align:center;"> 4.0.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/bit </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/bit </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-08-04 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> bit64 </td>
   <td style="text-align:center;"> bit64 </td>
   <td style="text-align:center;"> 4.0.5 </td>
   <td style="text-align:center;"> 4.0.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/bit64 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/bit64 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-08-30 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> bitops </td>
   <td style="text-align:center;"> bitops </td>
   <td style="text-align:center;"> 1.0.6 </td>
   <td style="text-align:center;"> 1.0-6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/bitops </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/bitops </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2013-08-17 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> blob </td>
   <td style="text-align:center;"> blob </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/blob </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/blob </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-01-20 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> broom </td>
   <td style="text-align:center;"> broom </td>
   <td style="text-align:center;"> 0.7.3 </td>
   <td style="text-align:center;"> 0.7.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/broom </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/broom </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-16 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BSgenome </td>
   <td style="text-align:center;"> BSgenome </td>
   <td style="text-align:center;"> 1.58.0 </td>
   <td style="text-align:center;"> 1.58.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/BSgenome </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/BSgenome </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> bsseq </td>
   <td style="text-align:center;"> bsseq </td>
   <td style="text-align:center;"> 1.26.0 </td>
   <td style="text-align:center;"> 1.26.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/bsseq </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/bsseq </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> bumphunter </td>
   <td style="text-align:center;"> bumphunter </td>
   <td style="text-align:center;"> 1.32.0 </td>
   <td style="text-align:center;"> 1.32.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/bumphunter </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/bumphunter </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cairo </td>
   <td style="text-align:center;"> Cairo </td>
   <td style="text-align:center;"> 1.5.12.2 </td>
   <td style="text-align:center;"> 1.5-12.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Cairo </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Cairo </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-07-07 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> calibrate </td>
   <td style="text-align:center;"> calibrate </td>
   <td style="text-align:center;"> 1.7.7 </td>
   <td style="text-align:center;"> 1.7.7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/calibrate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/calibrate </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-06-19 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> callr </td>
   <td style="text-align:center;"> callr </td>
   <td style="text-align:center;"> 3.5.1 </td>
   <td style="text-align:center;"> 3.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/callr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/callr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-13 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> car </td>
   <td style="text-align:center;"> car </td>
   <td style="text-align:center;"> 3.0.10 </td>
   <td style="text-align:center;"> 3.0-10 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/car </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/car </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-09-29 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> carData </td>
   <td style="text-align:center;"> carData </td>
   <td style="text-align:center;"> 3.0.4 </td>
   <td style="text-align:center;"> 3.0-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/carData </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/carData </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-05-22 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cellranger </td>
   <td style="text-align:center;"> cellranger </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/cellranger </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/cellranger </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2016-07-27 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> checkmate </td>
   <td style="text-align:center;"> checkmate </td>
   <td style="text-align:center;"> 2.0.0 </td>
   <td style="text-align:center;"> 2.0.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/checkmate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/checkmate </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-02-06 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> circlize </td>
   <td style="text-align:center;"> circlize </td>
   <td style="text-align:center;"> 0.4.12 </td>
   <td style="text-align:center;"> 0.4.12 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/circlize </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/circlize </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-08 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ckmeans.1d.dp </td>
   <td style="text-align:center;"> Ckmeans.1d.dp </td>
   <td style="text-align:center;"> 4.3.3 </td>
   <td style="text-align:center;"> 4.3.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Ckmeans.1d.dp </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Ckmeans.1d.dp </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-07-22 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cli </td>
   <td style="text-align:center;"> cli </td>
   <td style="text-align:center;"> 2.2.0 </td>
   <td style="text-align:center;"> 2.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/cli </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/cli </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-20 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> clue </td>
   <td style="text-align:center;"> clue </td>
   <td style="text-align:center;"> 0.3.58 </td>
   <td style="text-align:center;"> 0.3-58 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/clue </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/clue </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-03 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cluster </td>
   <td style="text-align:center;"> cluster </td>
   <td style="text-align:center;"> 2.1.0 </td>
   <td style="text-align:center;"> 2.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/cluster </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/cluster </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-06-19 </td>
   <td style="text-align:center;"> CRAN (R 4.0.3) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> codetools </td>
   <td style="text-align:center;"> codetools </td>
   <td style="text-align:center;"> 0.2.18 </td>
   <td style="text-align:center;"> 0.2-18 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/codetools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/codetools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-04 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> colorspace </td>
   <td style="text-align:center;"> colorspace </td>
   <td style="text-align:center;"> 2.0.0 </td>
   <td style="text-align:center;"> 2.0-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/colorspace </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/colorspace </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-11 </td>
   <td style="text-align:center;"> CRAN (R 4.0.3) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ComplexHeatmap </td>
   <td style="text-align:center;"> ComplexHeatmap </td>
   <td style="text-align:center;"> 2.6.2 </td>
   <td style="text-align:center;"> 2.6.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ComplexHeatmap </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ComplexHeatmap </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-12 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> corrplot </td>
   <td style="text-align:center;"> corrplot </td>
   <td style="text-align:center;"> 0.84 </td>
   <td style="text-align:center;"> 0.84 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/corrplot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/corrplot </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2017-10-16 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cowplot </td>
   <td style="text-align:center;"> cowplot </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/cowplot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/cowplot </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-30 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> crayon </td>
   <td style="text-align:center;"> crayon </td>
   <td style="text-align:center;"> 1.3.4 </td>
   <td style="text-align:center;"> 1.3.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/crayon </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/crayon </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2017-09-16 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> crosstalk </td>
   <td style="text-align:center;"> crosstalk </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/crosstalk </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/crosstalk </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-12 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> curl </td>
   <td style="text-align:center;"> curl </td>
   <td style="text-align:center;"> 4.3 </td>
   <td style="text-align:center;"> 4.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/curl </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/curl </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-12-02 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:center;"> data.table </td>
   <td style="text-align:center;"> 1.13.6 </td>
   <td style="text-align:center;"> 1.13.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/data.table </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/data.table </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-30 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DBI </td>
   <td style="text-align:center;"> DBI </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/DBI </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/DBI </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-15 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dbplyr </td>
   <td style="text-align:center;"> dbplyr </td>
   <td style="text-align:center;"> 2.0.0 </td>
   <td style="text-align:center;"> 2.0.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/dbplyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/dbplyr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-03 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DelayedArray </td>
   <td style="text-align:center;"> DelayedArray </td>
   <td style="text-align:center;"> 0.16.0 </td>
   <td style="text-align:center;"> 0.16.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/DelayedArray </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/DelayedArray </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DelayedMatrixStats </td>
   <td style="text-align:center;"> DelayedMatrixStats </td>
   <td style="text-align:center;"> 1.12.1 </td>
   <td style="text-align:center;"> 1.12.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/DelayedMatrixStats </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/DelayedMatrixStats </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-24 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dendextend </td>
   <td style="text-align:center;"> dendextend </td>
   <td style="text-align:center;"> 1.14.0 </td>
   <td style="text-align:center;"> 1.14.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/dendextend </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/dendextend </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-08-26 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> desc </td>
   <td style="text-align:center;"> desc </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/desc </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/desc </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2018-05-01 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DESeq2 </td>
   <td style="text-align:center;"> DESeq2 </td>
   <td style="text-align:center;"> 1.30.0 </td>
   <td style="text-align:center;"> 1.30.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/DESeq2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/DESeq2 </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> devtools </td>
   <td style="text-align:center;"> devtools </td>
   <td style="text-align:center;"> 2.3.2 </td>
   <td style="text-align:center;"> 2.3.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/devtools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/devtools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-09-18 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dichromat </td>
   <td style="text-align:center;"> dichromat </td>
   <td style="text-align:center;"> 2.0.0 </td>
   <td style="text-align:center;"> 2.0-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/dichromat </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/dichromat </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2013-01-24 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> digest </td>
   <td style="text-align:center;"> digest </td>
   <td style="text-align:center;"> 0.6.27 </td>
   <td style="text-align:center;"> 0.6.27 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/digest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/digest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-24 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DMRcate </td>
   <td style="text-align:center;"> DMRcate </td>
   <td style="text-align:center;"> 2.4.0 </td>
   <td style="text-align:center;"> 2.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/DMRcate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/DMRcate </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DMRcatedata </td>
   <td style="text-align:center;"> DMRcatedata </td>
   <td style="text-align:center;"> 2.8.0 </td>
   <td style="text-align:center;"> 2.8.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/DMRcatedata </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/DMRcatedata </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-29 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> doRNG </td>
   <td style="text-align:center;"> doRNG </td>
   <td style="text-align:center;"> 1.8.2 </td>
   <td style="text-align:center;"> 1.8.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/doRNG </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/doRNG </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-01-27 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dplyr </td>
   <td style="text-align:center;"> dplyr </td>
   <td style="text-align:center;"> 1.0.3 </td>
   <td style="text-align:center;"> 1.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/dplyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/dplyr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-15 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DSS </td>
   <td style="text-align:center;"> DSS </td>
   <td style="text-align:center;"> 2.38.0 </td>
   <td style="text-align:center;"> 2.38.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/DSS </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/DSS </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> edgeR </td>
   <td style="text-align:center;"> edgeR </td>
   <td style="text-align:center;"> 3.32.0 </td>
   <td style="text-align:center;"> 3.32.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/edgeR </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/edgeR </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ellipsis </td>
   <td style="text-align:center;"> ellipsis </td>
   <td style="text-align:center;"> 0.3.1 </td>
   <td style="text-align:center;"> 0.3.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ellipsis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ellipsis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-05-15 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ensembldb </td>
   <td style="text-align:center;"> ensembldb </td>
   <td style="text-align:center;"> 2.14.0 </td>
   <td style="text-align:center;"> 2.14.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ensembldb </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ensembldb </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> evaluate </td>
   <td style="text-align:center;"> evaluate </td>
   <td style="text-align:center;"> 0.14 </td>
   <td style="text-align:center;"> 0.14 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/evaluate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/evaluate </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-05-28 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ExperimentHub </td>
   <td style="text-align:center;"> ExperimentHub </td>
   <td style="text-align:center;"> 1.16.0 </td>
   <td style="text-align:center;"> 1.16.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ExperimentHub </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ExperimentHub </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> factoextra </td>
   <td style="text-align:center;"> factoextra </td>
   <td style="text-align:center;"> 1.0.7 </td>
   <td style="text-align:center;"> 1.0.7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/factoextra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/factoextra </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-04-01 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fansi </td>
   <td style="text-align:center;"> fansi </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/fansi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/fansi </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-15 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> farver </td>
   <td style="text-align:center;"> farver </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/farver </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/farver </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-01-16 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fastmap </td>
   <td style="text-align:center;"> fastmap </td>
   <td style="text-align:center;"> 1.0.1 </td>
   <td style="text-align:center;"> 1.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/fastmap </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/fastmap </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-10-08 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ff </td>
   <td style="text-align:center;"> ff </td>
   <td style="text-align:center;"> 4.0.4 </td>
   <td style="text-align:center;"> 4.0.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ff </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ff </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-13 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> forcats </td>
   <td style="text-align:center;"> forcats </td>
   <td style="text-align:center;"> 0.5.0 </td>
   <td style="text-align:center;"> 0.5.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/forcats </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/forcats </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-03-01 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> foreach </td>
   <td style="text-align:center;"> foreach </td>
   <td style="text-align:center;"> 1.5.1 </td>
   <td style="text-align:center;"> 1.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/foreach </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/foreach </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-15 </td>
   <td style="text-align:center;"> CRAN (R 4.0.3) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> foreign </td>
   <td style="text-align:center;"> foreign </td>
   <td style="text-align:center;"> 0.8.81 </td>
   <td style="text-align:center;"> 0.8-81 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/foreign </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/foreign </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-22 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Formula </td>
   <td style="text-align:center;"> Formula </td>
   <td style="text-align:center;"> 1.2.4 </td>
   <td style="text-align:center;"> 1.2-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Formula </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Formula </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-16 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fs </td>
   <td style="text-align:center;"> fs </td>
   <td style="text-align:center;"> 1.5.0 </td>
   <td style="text-align:center;"> 1.5.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/fs </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/fs </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-07-31 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gbRd </td>
   <td style="text-align:center;"> gbRd </td>
   <td style="text-align:center;"> 0.4.11 </td>
   <td style="text-align:center;"> 0.4-11 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/gbRd </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/gbRd </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2012-10-01 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> genefilter </td>
   <td style="text-align:center;"> genefilter </td>
   <td style="text-align:center;"> 1.72.0 </td>
   <td style="text-align:center;"> 1.72.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/genefilter </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/genefilter </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> geneplotter </td>
   <td style="text-align:center;"> geneplotter </td>
   <td style="text-align:center;"> 1.68.0 </td>
   <td style="text-align:center;"> 1.68.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/geneplotter </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/geneplotter </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> generics </td>
   <td style="text-align:center;"> generics </td>
   <td style="text-align:center;"> 0.1.0 </td>
   <td style="text-align:center;"> 0.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/generics </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/generics </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-31 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GenomeInfoDb </td>
   <td style="text-align:center;"> GenomeInfoDb </td>
   <td style="text-align:center;"> 1.26.2 </td>
   <td style="text-align:center;"> 1.26.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GenomeInfoDb </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GenomeInfoDb </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-08 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GenomeInfoDbData </td>
   <td style="text-align:center;"> GenomeInfoDbData </td>
   <td style="text-align:center;"> 1.2.4 </td>
   <td style="text-align:center;"> 1.2.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GenomeInfoDbData </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GenomeInfoDbData </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-16 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GenomicAlignments </td>
   <td style="text-align:center;"> GenomicAlignments </td>
   <td style="text-align:center;"> 1.26.0 </td>
   <td style="text-align:center;"> 1.26.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GenomicAlignments </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GenomicAlignments </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GenomicFeatures </td>
   <td style="text-align:center;"> GenomicFeatures </td>
   <td style="text-align:center;"> 1.42.1 </td>
   <td style="text-align:center;"> 1.42.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GenomicFeatures </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GenomicFeatures </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-11 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GenomicRanges </td>
   <td style="text-align:center;"> GenomicRanges </td>
   <td style="text-align:center;"> 1.42.0 </td>
   <td style="text-align:center;"> 1.42.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GenomicRanges </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GenomicRanges </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GEOquery </td>
   <td style="text-align:center;"> GEOquery </td>
   <td style="text-align:center;"> 2.58.0 </td>
   <td style="text-align:center;"> 2.58.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GEOquery </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GEOquery </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GetoptLong </td>
   <td style="text-align:center;"> GetoptLong </td>
   <td style="text-align:center;"> 1.0.5 </td>
   <td style="text-align:center;"> 1.0.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GetoptLong </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GetoptLong </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-15 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggplot2 </td>
   <td style="text-align:center;"> ggplot2 </td>
   <td style="text-align:center;"> 3.3.3 </td>
   <td style="text-align:center;"> 3.3.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ggplot2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ggplot2 </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-30 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggpubr </td>
   <td style="text-align:center;"> ggpubr </td>
   <td style="text-align:center;"> 0.4.0 </td>
   <td style="text-align:center;"> 0.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ggpubr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ggpubr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-06-27 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggrepel </td>
   <td style="text-align:center;"> ggrepel </td>
   <td style="text-align:center;"> 0.9.1 </td>
   <td style="text-align:center;"> 0.9.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ggrepel </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ggrepel </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-15 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggridges </td>
   <td style="text-align:center;"> ggridges </td>
   <td style="text-align:center;"> 0.5.3 </td>
   <td style="text-align:center;"> 0.5.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ggridges </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ggridges </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-08 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggsignif </td>
   <td style="text-align:center;"> ggsignif </td>
   <td style="text-align:center;"> 0.6.0 </td>
   <td style="text-align:center;"> 0.6.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ggsignif </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ggsignif </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-08-08 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GlobalOptions </td>
   <td style="text-align:center;"> GlobalOptions </td>
   <td style="text-align:center;"> 0.1.2 </td>
   <td style="text-align:center;"> 0.1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GlobalOptions </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GlobalOptions </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-06-10 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> glue </td>
   <td style="text-align:center;"> glue </td>
   <td style="text-align:center;"> 1.4.2 </td>
   <td style="text-align:center;"> 1.4.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/glue </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/glue </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-08-27 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO.db </td>
   <td style="text-align:center;"> GO.db </td>
   <td style="text-align:center;"> 3.12.1 </td>
   <td style="text-align:center;"> 3.12.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GO.db </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/GO.db </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-16 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gridExtra </td>
   <td style="text-align:center;"> gridExtra </td>
   <td style="text-align:center;"> 2.3 </td>
   <td style="text-align:center;"> 2.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/gridExtra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/gridExtra </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2017-09-09 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gtable </td>
   <td style="text-align:center;"> gtable </td>
   <td style="text-align:center;"> 0.3.0 </td>
   <td style="text-align:center;"> 0.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/gtable </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/gtable </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-03-25 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gtools </td>
   <td style="text-align:center;"> gtools </td>
   <td style="text-align:center;"> 3.8.2 </td>
   <td style="text-align:center;"> 3.8.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/gtools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/gtools </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-03-31 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gviz </td>
   <td style="text-align:center;"> Gviz </td>
   <td style="text-align:center;"> 1.34.0 </td>
   <td style="text-align:center;"> 1.34.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Gviz </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Gviz </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Haplin </td>
   <td style="text-align:center;"> Haplin </td>
   <td style="text-align:center;"> 7.2.3 </td>
   <td style="text-align:center;"> 7.2.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Haplin </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Haplin </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-09-07 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> haven </td>
   <td style="text-align:center;"> haven </td>
   <td style="text-align:center;"> 2.3.1 </td>
   <td style="text-align:center;"> 2.3.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/haven </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/haven </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-06-01 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HDF5Array </td>
   <td style="text-align:center;"> HDF5Array </td>
   <td style="text-align:center;"> 1.18.0 </td>
   <td style="text-align:center;"> 1.18.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/HDF5Array </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/HDF5Array </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> highr </td>
   <td style="text-align:center;"> highr </td>
   <td style="text-align:center;"> 0.8 </td>
   <td style="text-align:center;"> 0.8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/highr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/highr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-03-20 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hmisc </td>
   <td style="text-align:center;"> Hmisc </td>
   <td style="text-align:center;"> 4.4.2 </td>
   <td style="text-align:center;"> 4.4-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Hmisc </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Hmisc </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-29 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hms </td>
   <td style="text-align:center;"> hms </td>
   <td style="text-align:center;"> 1.0.0 </td>
   <td style="text-align:center;"> 1.0.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/hms </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/hms </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-13 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> htmlTable </td>
   <td style="text-align:center;"> htmlTable </td>
   <td style="text-align:center;"> 2.1.0 </td>
   <td style="text-align:center;"> 2.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/htmlTable </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/htmlTable </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-09-16 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> htmltools </td>
   <td style="text-align:center;"> htmltools </td>
   <td style="text-align:center;"> 0.5.1 </td>
   <td style="text-align:center;"> 0.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/htmltools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/htmltools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-12 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> htmlwidgets </td>
   <td style="text-align:center;"> htmlwidgets </td>
   <td style="text-align:center;"> 1.5.3 </td>
   <td style="text-align:center;"> 1.5.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/htmlwidgets </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/htmlwidgets </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-10 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> httpuv </td>
   <td style="text-align:center;"> httpuv </td>
   <td style="text-align:center;"> 1.5.5 </td>
   <td style="text-align:center;"> 1.5.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/httpuv </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/httpuv </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-13 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> httr </td>
   <td style="text-align:center;"> httr </td>
   <td style="text-align:center;"> 1.4.2 </td>
   <td style="text-align:center;"> 1.4.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/httr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/httr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-07-20 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> igraph </td>
   <td style="text-align:center;"> igraph </td>
   <td style="text-align:center;"> 1.2.6 </td>
   <td style="text-align:center;"> 1.2.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/igraph </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/igraph </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-06 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IlluminaHumanMethylation450kanno.ilmn12.hg19 </td>
   <td style="text-align:center;"> IlluminaHumanMethylation450kanno.ilmn12.hg19 </td>
   <td style="text-align:center;"> 0.6.0 </td>
   <td style="text-align:center;"> 0.6.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/IlluminaHumanMethylation450kanno.ilmn12.hg19 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/IlluminaHumanMethylation450kanno.ilmn12.hg19 </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-10 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IlluminaHumanMethylation450kmanifest </td>
   <td style="text-align:center;"> IlluminaHumanMethylation450kmanifest </td>
   <td style="text-align:center;"> 0.4.0 </td>
   <td style="text-align:center;"> 0.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/IlluminaHumanMethylation450kmanifest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/IlluminaHumanMethylation450kmanifest </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-05-19 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IlluminaHumanMethylationEPICanno.ilm10b4.hg19 </td>
   <td style="text-align:center;"> IlluminaHumanMethylationEPICanno.ilm10b4.hg19 </td>
   <td style="text-align:center;"> 0.6.0 </td>
   <td style="text-align:center;"> 0.6.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/IlluminaHumanMethylationEPICanno.ilm10b4.hg19 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/IlluminaHumanMethylationEPICanno.ilm10b4.hg19 </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-05-19 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IlluminaHumanMethylationEPICanno.ilm10b5.hg38 </td>
   <td style="text-align:center;"> IlluminaHumanMethylationEPICanno.ilm10b5.hg38 </td>
   <td style="text-align:center;"> 0.0.1 </td>
   <td style="text-align:center;"> 0.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/IlluminaHumanMethylationEPICanno.ilm10b5.hg38 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/IlluminaHumanMethylationEPICanno.ilm10b5.hg38 </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-06 </td>
   <td style="text-align:center;"> Github (achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38@3db0691) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IlluminaHumanMethylationEPICmanifest </td>
   <td style="text-align:center;"> IlluminaHumanMethylationEPICmanifest </td>
   <td style="text-align:center;"> 0.3.0 </td>
   <td style="text-align:center;"> 0.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/IlluminaHumanMethylationEPICmanifest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/IlluminaHumanMethylationEPICmanifest </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-16 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> illuminaio </td>
   <td style="text-align:center;"> illuminaio </td>
   <td style="text-align:center;"> 0.32.0 </td>
   <td style="text-align:center;"> 0.32.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/illuminaio </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/illuminaio </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> interactiveDisplayBase </td>
   <td style="text-align:center;"> interactiveDisplayBase </td>
   <td style="text-align:center;"> 1.28.0 </td>
   <td style="text-align:center;"> 1.28.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/interactiveDisplayBase </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/interactiveDisplayBase </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IRanges </td>
   <td style="text-align:center;"> IRanges </td>
   <td style="text-align:center;"> 2.24.1 </td>
   <td style="text-align:center;"> 2.24.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/IRanges </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/IRanges </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-12 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> iterators </td>
   <td style="text-align:center;"> iterators </td>
   <td style="text-align:center;"> 1.0.13 </td>
   <td style="text-align:center;"> 1.0.13 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/iterators </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/iterators </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-15 </td>
   <td style="text-align:center;"> CRAN (R 4.0.3) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> jpeg </td>
   <td style="text-align:center;"> jpeg </td>
   <td style="text-align:center;"> 0.1.8.1 </td>
   <td style="text-align:center;"> 0.1-8.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/jpeg </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/jpeg </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-10-24 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> jsonlite </td>
   <td style="text-align:center;"> jsonlite </td>
   <td style="text-align:center;"> 1.7.2 </td>
   <td style="text-align:center;"> 1.7.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/jsonlite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/jsonlite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-09 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> kableExtra </td>
   <td style="text-align:center;"> kableExtra </td>
   <td style="text-align:center;"> 1.3.1 </td>
   <td style="text-align:center;"> 1.3.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/kableExtra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/kableExtra </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-22 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> knitr </td>
   <td style="text-align:center;"> knitr </td>
   <td style="text-align:center;"> 1.30 </td>
   <td style="text-align:center;"> 1.30 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/knitr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/knitr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-09-22 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> labeling </td>
   <td style="text-align:center;"> labeling </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/labeling </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/labeling </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-20 </td>
   <td style="text-align:center;"> CRAN (R 4.0.3) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> later </td>
   <td style="text-align:center;"> later </td>
   <td style="text-align:center;"> 1.1.0.1 </td>
   <td style="text-align:center;"> 1.1.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/later </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/later </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-06-05 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lattice </td>
   <td style="text-align:center;"> lattice </td>
   <td style="text-align:center;"> 0.20.41 </td>
   <td style="text-align:center;"> 0.20-41 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/lattice </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/lattice </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-04-02 </td>
   <td style="text-align:center;"> CRAN (R 4.0.3) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> latticeExtra </td>
   <td style="text-align:center;"> latticeExtra </td>
   <td style="text-align:center;"> 0.6.29 </td>
   <td style="text-align:center;"> 0.6-29 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/latticeExtra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/latticeExtra </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-12-19 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lazyeval </td>
   <td style="text-align:center;"> lazyeval </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/lazyeval </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/lazyeval </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-03-15 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lifecycle </td>
   <td style="text-align:center;"> lifecycle </td>
   <td style="text-align:center;"> 0.2.0 </td>
   <td style="text-align:center;"> 0.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/lifecycle </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/lifecycle </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-03-06 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> limma </td>
   <td style="text-align:center;"> limma </td>
   <td style="text-align:center;"> 3.46.0 </td>
   <td style="text-align:center;"> 3.46.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/limma </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/limma </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> locfit </td>
   <td style="text-align:center;"> locfit </td>
   <td style="text-align:center;"> 1.5.9.4 </td>
   <td style="text-align:center;"> 1.5-9.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/locfit </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/locfit </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-03-25 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lubridate </td>
   <td style="text-align:center;"> lubridate </td>
   <td style="text-align:center;"> 1.7.9.2 </td>
   <td style="text-align:center;"> 1.7.9.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/lubridate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/lubridate </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-13 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> magrittr </td>
   <td style="text-align:center;"> magrittr </td>
   <td style="text-align:center;"> 2.0.1 </td>
   <td style="text-align:center;"> 2.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/magrittr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/magrittr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-17 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MASS </td>
   <td style="text-align:center;"> MASS </td>
   <td style="text-align:center;"> 7.3.53 </td>
   <td style="text-align:center;"> 7.3-53 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/MASS </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/MASS </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-09-09 </td>
   <td style="text-align:center;"> CRAN (R 4.0.3) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Matrix </td>
   <td style="text-align:center;"> Matrix </td>
   <td style="text-align:center;"> 1.3.2 </td>
   <td style="text-align:center;"> 1.3-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Matrix </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Matrix </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-06 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MatrixGenerics </td>
   <td style="text-align:center;"> MatrixGenerics </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/MatrixGenerics </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/MatrixGenerics </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> matrixStats </td>
   <td style="text-align:center;"> matrixStats </td>
   <td style="text-align:center;"> 0.57.0 </td>
   <td style="text-align:center;"> 0.57.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/matrixStats </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/matrixStats </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-09-25 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mclust </td>
   <td style="text-align:center;"> mclust </td>
   <td style="text-align:center;"> 5.4.7 </td>
   <td style="text-align:center;"> 5.4.7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/mclust </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/mclust </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-20 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> memoise </td>
   <td style="text-align:center;"> memoise </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/memoise </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/memoise </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2017-04-21 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MethylToSNP </td>
   <td style="text-align:center;"> MethylToSNP </td>
   <td style="text-align:center;"> 0.99.0 </td>
   <td style="text-align:center;"> 0.99.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/MethylToSNP </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/MethylToSNP </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-06 </td>
   <td style="text-align:center;"> Github (elnitskilab/MethylToSNP@92a4282) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mgcv </td>
   <td style="text-align:center;"> mgcv </td>
   <td style="text-align:center;"> 1.8.33 </td>
   <td style="text-align:center;"> 1.8-33 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/mgcv </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/mgcv </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-08-27 </td>
   <td style="text-align:center;"> CRAN (R 4.0.3) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mime </td>
   <td style="text-align:center;"> mime </td>
   <td style="text-align:center;"> 0.9 </td>
   <td style="text-align:center;"> 0.9 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/mime </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/mime </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-02-04 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> minfi </td>
   <td style="text-align:center;"> minfi </td>
   <td style="text-align:center;"> 1.36.0 </td>
   <td style="text-align:center;"> 1.36.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/minfi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/minfi </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> missMethyl </td>
   <td style="text-align:center;"> missMethyl </td>
   <td style="text-align:center;"> 1.24.0 </td>
   <td style="text-align:center;"> 1.24.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/missMethyl </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/missMethyl </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> modelr </td>
   <td style="text-align:center;"> modelr </td>
   <td style="text-align:center;"> 0.1.8 </td>
   <td style="text-align:center;"> 0.1.8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/modelr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/modelr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-05-19 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> multtest </td>
   <td style="text-align:center;"> multtest </td>
   <td style="text-align:center;"> 2.46.0 </td>
   <td style="text-align:center;"> 2.46.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/multtest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/multtest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> munsell </td>
   <td style="text-align:center;"> munsell </td>
   <td style="text-align:center;"> 0.5.0 </td>
   <td style="text-align:center;"> 0.5.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/munsell </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/munsell </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2018-06-12 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nlme </td>
   <td style="text-align:center;"> nlme </td>
   <td style="text-align:center;"> 3.1.151 </td>
   <td style="text-align:center;"> 3.1-151 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/nlme </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/nlme </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-10 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nnet </td>
   <td style="text-align:center;"> nnet </td>
   <td style="text-align:center;"> 7.3.14 </td>
   <td style="text-align:center;"> 7.3-14 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/nnet </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/nnet </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-04-26 </td>
   <td style="text-align:center;"> CRAN (R 4.0.3) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nor1mix </td>
   <td style="text-align:center;"> nor1mix </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> 1.3-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/nor1mix </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/nor1mix </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-06-13 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> openssl </td>
   <td style="text-align:center;"> openssl </td>
   <td style="text-align:center;"> 1.4.3 </td>
   <td style="text-align:center;"> 1.4.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/openssl </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/openssl </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-09-18 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> openxlsx </td>
   <td style="text-align:center;"> openxlsx </td>
   <td style="text-align:center;"> 4.2.3 </td>
   <td style="text-align:center;"> 4.2.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/openxlsx </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/openxlsx </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> CRAN (R 4.0.3) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> org.Hs.eg.db </td>
   <td style="text-align:center;"> org.Hs.eg.db </td>
   <td style="text-align:center;"> 3.12.0 </td>
   <td style="text-align:center;"> 3.12.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/org.Hs.eg.db </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/org.Hs.eg.db </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-16 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pacman </td>
   <td style="text-align:center;"> pacman </td>
   <td style="text-align:center;"> 0.5.1 </td>
   <td style="text-align:center;"> 0.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/pacman </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/pacman </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-03-11 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> permute </td>
   <td style="text-align:center;"> permute </td>
   <td style="text-align:center;"> 0.9.5 </td>
   <td style="text-align:center;"> 0.9-5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/permute </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/permute </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-03-12 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pheatmap </td>
   <td style="text-align:center;"> pheatmap </td>
   <td style="text-align:center;"> 1.0.12 </td>
   <td style="text-align:center;"> 1.0.12 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/pheatmap </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/pheatmap </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-01-04 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pillar </td>
   <td style="text-align:center;"> pillar </td>
   <td style="text-align:center;"> 1.4.7 </td>
   <td style="text-align:center;"> 1.4.7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/pillar </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/pillar </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-20 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgbuild </td>
   <td style="text-align:center;"> pkgbuild </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/pkgbuild </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/pkgbuild </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-15 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgconfig </td>
   <td style="text-align:center;"> pkgconfig </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/pkgconfig </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/pkgconfig </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-09-22 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgload </td>
   <td style="text-align:center;"> pkgload </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/pkgload </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/pkgload </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-05-29 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> plotly </td>
   <td style="text-align:center;"> plotly </td>
   <td style="text-align:center;"> 4.9.3 </td>
   <td style="text-align:center;"> 4.9.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/plotly </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/plotly </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-10 </td>
   <td style="text-align:center;"> CRAN (R 4.0.3) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> plyr </td>
   <td style="text-align:center;"> plyr </td>
   <td style="text-align:center;"> 1.8.6 </td>
   <td style="text-align:center;"> 1.8.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/plyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/plyr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-03-03 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> png </td>
   <td style="text-align:center;"> png </td>
   <td style="text-align:center;"> 0.1.7 </td>
   <td style="text-align:center;"> 0.1-7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/png </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/png </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2013-12-03 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> preprocessCore </td>
   <td style="text-align:center;"> preprocessCore </td>
   <td style="text-align:center;"> 1.52.1 </td>
   <td style="text-align:center;"> 1.52.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/preprocessCore </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/preprocessCore </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-08 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> prettyunits </td>
   <td style="text-align:center;"> prettyunits </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/prettyunits </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/prettyunits </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-01-24 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> processx </td>
   <td style="text-align:center;"> processx </td>
   <td style="text-align:center;"> 3.4.5 </td>
   <td style="text-align:center;"> 3.4.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/processx </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/processx </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-30 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> progress </td>
   <td style="text-align:center;"> progress </td>
   <td style="text-align:center;"> 1.2.2 </td>
   <td style="text-align:center;"> 1.2.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/progress </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/progress </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-05-16 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> promises </td>
   <td style="text-align:center;"> promises </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/promises </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/promises </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-06-09 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ProtGenerics </td>
   <td style="text-align:center;"> ProtGenerics </td>
   <td style="text-align:center;"> 1.22.0 </td>
   <td style="text-align:center;"> 1.22.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ProtGenerics </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ProtGenerics </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ps </td>
   <td style="text-align:center;"> ps </td>
   <td style="text-align:center;"> 1.5.0 </td>
   <td style="text-align:center;"> 1.5.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ps </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/ps </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-05 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> purrr </td>
   <td style="text-align:center;"> purrr </td>
   <td style="text-align:center;"> 0.3.4 </td>
   <td style="text-align:center;"> 0.3.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/purrr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/purrr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-04-17 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> quadprog </td>
   <td style="text-align:center;"> quadprog </td>
   <td style="text-align:center;"> 1.5.8 </td>
   <td style="text-align:center;"> 1.5-8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/quadprog </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/quadprog </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-11-20 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> R.methodsS3 </td>
   <td style="text-align:center;"> R.methodsS3 </td>
   <td style="text-align:center;"> 1.8.1 </td>
   <td style="text-align:center;"> 1.8.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/R.methodsS3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/R.methodsS3 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-08-26 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> R.oo </td>
   <td style="text-align:center;"> R.oo </td>
   <td style="text-align:center;"> 1.24.0 </td>
   <td style="text-align:center;"> 1.24.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/R.oo </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/R.oo </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-08-26 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> R.utils </td>
   <td style="text-align:center;"> R.utils </td>
   <td style="text-align:center;"> 2.10.1 </td>
   <td style="text-align:center;"> 2.10.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/R.utils </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/R.utils </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-08-26 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> R6 </td>
   <td style="text-align:center;"> R6 </td>
   <td style="text-align:center;"> 2.5.0 </td>
   <td style="text-align:center;"> 2.5.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/R6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/R6 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-28 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rappdirs </td>
   <td style="text-align:center;"> rappdirs </td>
   <td style="text-align:center;"> 0.3.1 </td>
   <td style="text-align:center;"> 0.3.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rappdirs </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rappdirs </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2016-03-28 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rbibutils </td>
   <td style="text-align:center;"> rbibutils </td>
   <td style="text-align:center;"> 2.0 </td>
   <td style="text-align:center;"> 2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rbibutils </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rbibutils </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-18 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RColorBrewer </td>
   <td style="text-align:center;"> RColorBrewer </td>
   <td style="text-align:center;"> 1.1.2 </td>
   <td style="text-align:center;"> 1.1-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/RColorBrewer </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/RColorBrewer </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2014-12-07 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rcpp </td>
   <td style="text-align:center;"> Rcpp </td>
   <td style="text-align:center;"> 1.0.6 </td>
   <td style="text-align:center;"> 1.0.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Rcpp </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Rcpp </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-15 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RCurl </td>
   <td style="text-align:center;"> RCurl </td>
   <td style="text-align:center;"> 1.98.1.2 </td>
   <td style="text-align:center;"> 1.98-1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/RCurl </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/RCurl </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-04-18 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rdpack </td>
   <td style="text-align:center;"> Rdpack </td>
   <td style="text-align:center;"> 2.1 </td>
   <td style="text-align:center;"> 2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Rdpack </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Rdpack </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-09 </td>
   <td style="text-align:center;"> CRAN (R 4.0.3) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> readr </td>
   <td style="text-align:center;"> readr </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/readr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/readr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-05 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> readxl </td>
   <td style="text-align:center;"> readxl </td>
   <td style="text-align:center;"> 1.3.1 </td>
   <td style="text-align:center;"> 1.3.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/readxl </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/readxl </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-03-13 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> remotes </td>
   <td style="text-align:center;"> remotes </td>
   <td style="text-align:center;"> 2.2.0 </td>
   <td style="text-align:center;"> 2.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/remotes </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/remotes </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-07-21 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> reprex </td>
   <td style="text-align:center;"> reprex </td>
   <td style="text-align:center;"> 0.3.0 </td>
   <td style="text-align:center;"> 0.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/reprex </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/reprex </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-05-16 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> reshape </td>
   <td style="text-align:center;"> reshape </td>
   <td style="text-align:center;"> 0.8.8 </td>
   <td style="text-align:center;"> 0.8.8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/reshape </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/reshape </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2018-10-23 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> reshape2 </td>
   <td style="text-align:center;"> reshape2 </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/reshape2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/reshape2 </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-04-09 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rhdf5 </td>
   <td style="text-align:center;"> rhdf5 </td>
   <td style="text-align:center;"> 2.34.0 </td>
   <td style="text-align:center;"> 2.34.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rhdf5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rhdf5 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rhdf5filters </td>
   <td style="text-align:center;"> rhdf5filters </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rhdf5filters </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rhdf5filters </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rhdf5lib </td>
   <td style="text-align:center;"> Rhdf5lib </td>
   <td style="text-align:center;"> 1.12.0 </td>
   <td style="text-align:center;"> 1.12.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Rhdf5lib </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Rhdf5lib </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rio </td>
   <td style="text-align:center;"> rio </td>
   <td style="text-align:center;"> 0.5.16 </td>
   <td style="text-align:center;"> 0.5.16 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rio </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rio </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2018-11-26 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rjson </td>
   <td style="text-align:center;"> rjson </td>
   <td style="text-align:center;"> 0.2.20 </td>
   <td style="text-align:center;"> 0.2.20 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rjson </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rjson </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2018-06-08 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rlang </td>
   <td style="text-align:center;"> rlang </td>
   <td style="text-align:center;"> 0.4.10 </td>
   <td style="text-align:center;"> 0.4.10 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rlang </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rlang </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-30 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rmarkdown </td>
   <td style="text-align:center;"> rmarkdown </td>
   <td style="text-align:center;"> 2.6 </td>
   <td style="text-align:center;"> 2.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rmarkdown </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rmarkdown </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-14 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rngtools </td>
   <td style="text-align:center;"> rngtools </td>
   <td style="text-align:center;"> 1.5 </td>
   <td style="text-align:center;"> 1.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rngtools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rngtools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-01-23 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rpart </td>
   <td style="text-align:center;"> rpart </td>
   <td style="text-align:center;"> 4.1.15 </td>
   <td style="text-align:center;"> 4.1-15 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rpart </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rpart </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-04-12 </td>
   <td style="text-align:center;"> CRAN (R 4.0.3) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rprojroot </td>
   <td style="text-align:center;"> rprojroot </td>
   <td style="text-align:center;"> 2.0.2 </td>
   <td style="text-align:center;"> 2.0.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rprojroot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rprojroot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-15 </td>
   <td style="text-align:center;"> CRAN (R 4.0.3) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rsamtools </td>
   <td style="text-align:center;"> Rsamtools </td>
   <td style="text-align:center;"> 2.6.0 </td>
   <td style="text-align:center;"> 2.6.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Rsamtools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/Rsamtools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RSQLite </td>
   <td style="text-align:center;"> RSQLite </td>
   <td style="text-align:center;"> 2.2.2 </td>
   <td style="text-align:center;"> 2.2.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/RSQLite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/RSQLite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-08 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstatix </td>
   <td style="text-align:center;"> rstatix </td>
   <td style="text-align:center;"> 0.6.0 </td>
   <td style="text-align:center;"> 0.6.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rstatix </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rstatix </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-06-18 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstudioapi </td>
   <td style="text-align:center;"> rstudioapi </td>
   <td style="text-align:center;"> 0.13 </td>
   <td style="text-align:center;"> 0.13 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rstudioapi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rstudioapi </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-12 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rtracklayer </td>
   <td style="text-align:center;"> rtracklayer </td>
   <td style="text-align:center;"> 1.50.0 </td>
   <td style="text-align:center;"> 1.50.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rtracklayer </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rtracklayer </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rvest </td>
   <td style="text-align:center;"> rvest </td>
   <td style="text-align:center;"> 0.3.6 </td>
   <td style="text-align:center;"> 0.3.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rvest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rvest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-07-25 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S4Vectors </td>
   <td style="text-align:center;"> S4Vectors </td>
   <td style="text-align:center;"> 0.28.1 </td>
   <td style="text-align:center;"> 0.28.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/S4Vectors </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/S4Vectors </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-09 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> scales </td>
   <td style="text-align:center;"> scales </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/scales </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/scales </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-05-11 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> scrime </td>
   <td style="text-align:center;"> scrime </td>
   <td style="text-align:center;"> 1.3.5 </td>
   <td style="text-align:center;"> 1.3.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/scrime </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/scrime </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2018-12-01 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sessioninfo </td>
   <td style="text-align:center;"> sessioninfo </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/sessioninfo </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/sessioninfo </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2018-11-05 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> shape </td>
   <td style="text-align:center;"> shape </td>
   <td style="text-align:center;"> 1.4.5 </td>
   <td style="text-align:center;"> 1.4.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/shape </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/shape </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-09-13 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> shiny </td>
   <td style="text-align:center;"> shiny </td>
   <td style="text-align:center;"> 1.5.0 </td>
   <td style="text-align:center;"> 1.5.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/shiny </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/shiny </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-06-23 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> shinyMethyl </td>
   <td style="text-align:center;"> shinyMethyl </td>
   <td style="text-align:center;"> 1.26.0 </td>
   <td style="text-align:center;"> 1.26.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/shinyMethyl </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/shinyMethyl </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> siggenes </td>
   <td style="text-align:center;"> siggenes </td>
   <td style="text-align:center;"> 1.64.0 </td>
   <td style="text-align:center;"> 1.64.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/siggenes </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/siggenes </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sparseMatrixStats </td>
   <td style="text-align:center;"> sparseMatrixStats </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/sparseMatrixStats </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/sparseMatrixStats </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> statmod </td>
   <td style="text-align:center;"> statmod </td>
   <td style="text-align:center;"> 1.4.35 </td>
   <td style="text-align:center;"> 1.4.35 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/statmod </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/statmod </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-19 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> stringi </td>
   <td style="text-align:center;"> stringi </td>
   <td style="text-align:center;"> 1.5.3 </td>
   <td style="text-align:center;"> 1.5.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/stringi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/stringi </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-09-09 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> stringr </td>
   <td style="text-align:center;"> stringr </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/stringr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/stringr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-02-10 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SummarizedExperiment </td>
   <td style="text-align:center;"> SummarizedExperiment </td>
   <td style="text-align:center;"> 1.20.0 </td>
   <td style="text-align:center;"> 1.20.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/SummarizedExperiment </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/SummarizedExperiment </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-27 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> survival </td>
   <td style="text-align:center;"> survival </td>
   <td style="text-align:center;"> 3.2.7 </td>
   <td style="text-align:center;"> 3.2-7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/survival </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/survival </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-09-28 </td>
   <td style="text-align:center;"> CRAN (R 4.0.3) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> testthat </td>
   <td style="text-align:center;"> testthat </td>
   <td style="text-align:center;"> 3.0.1 </td>
   <td style="text-align:center;"> 3.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/testthat </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/testthat </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-17 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tibble </td>
   <td style="text-align:center;"> tibble </td>
   <td style="text-align:center;"> 3.0.5 </td>
   <td style="text-align:center;"> 3.0.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/tibble </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/tibble </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-15 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidyr </td>
   <td style="text-align:center;"> tidyr </td>
   <td style="text-align:center;"> 1.1.2 </td>
   <td style="text-align:center;"> 1.1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/tidyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/tidyr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-08-27 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidyselect </td>
   <td style="text-align:center;"> tidyselect </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/tidyselect </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/tidyselect </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-05-11 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidyverse </td>
   <td style="text-align:center;"> tidyverse </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/tidyverse </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/tidyverse </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-11-21 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> usethis </td>
   <td style="text-align:center;"> usethis </td>
   <td style="text-align:center;"> 2.0.0 </td>
   <td style="text-align:center;"> 2.0.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/usethis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/usethis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-10 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> VariantAnnotation </td>
   <td style="text-align:center;"> VariantAnnotation </td>
   <td style="text-align:center;"> 1.36.0 </td>
   <td style="text-align:center;"> 1.36.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/VariantAnnotation </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/VariantAnnotation </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-28 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> vctrs </td>
   <td style="text-align:center;"> vctrs </td>
   <td style="text-align:center;"> 0.3.6 </td>
   <td style="text-align:center;"> 0.3.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/vctrs </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/vctrs </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-17 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> viridis </td>
   <td style="text-align:center;"> viridis </td>
   <td style="text-align:center;"> 0.5.1 </td>
   <td style="text-align:center;"> 0.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/viridis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/viridis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2018-03-29 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> viridisLite </td>
   <td style="text-align:center;"> viridisLite </td>
   <td style="text-align:center;"> 0.3.0 </td>
   <td style="text-align:center;"> 0.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/viridisLite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/viridisLite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2018-02-01 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> webshot </td>
   <td style="text-align:center;"> webshot </td>
   <td style="text-align:center;"> 0.5.2 </td>
   <td style="text-align:center;"> 0.5.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/webshot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/webshot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-11-22 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> wesanderson </td>
   <td style="text-align:center;"> wesanderson </td>
   <td style="text-align:center;"> 0.3.6 </td>
   <td style="text-align:center;"> 0.3.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/wesanderson </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/wesanderson </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2018-04-20 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> withr </td>
   <td style="text-align:center;"> withr </td>
   <td style="text-align:center;"> 2.4.0 </td>
   <td style="text-align:center;"> 2.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/withr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/withr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-16 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xfun </td>
   <td style="text-align:center;"> xfun </td>
   <td style="text-align:center;"> 0.20 </td>
   <td style="text-align:center;"> 0.20 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/xfun </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/xfun </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-06 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> XML </td>
   <td style="text-align:center;"> XML </td>
   <td style="text-align:center;"> 3.99.0.5 </td>
   <td style="text-align:center;"> 3.99-0.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/XML </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/XML </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-07-23 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xml2 </td>
   <td style="text-align:center;"> xml2 </td>
   <td style="text-align:center;"> 1.3.2 </td>
   <td style="text-align:center;"> 1.3.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/xml2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/xml2 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-04-23 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xtable </td>
   <td style="text-align:center;"> xtable </td>
   <td style="text-align:center;"> 1.8.4 </td>
   <td style="text-align:center;"> 1.8-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/xtable </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/xtable </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-04-21 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> XVector </td>
   <td style="text-align:center;"> XVector </td>
   <td style="text-align:center;"> 0.30.0 </td>
   <td style="text-align:center;"> 0.30.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/XVector </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/XVector </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-28 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> yaml </td>
   <td style="text-align:center;"> yaml </td>
   <td style="text-align:center;"> 2.2.1 </td>
   <td style="text-align:center;"> 2.2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/yaml </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/yaml </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-02-01 </td>
   <td style="text-align:center;"> CRAN (R 4.0.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> zip </td>
   <td style="text-align:center;"> zip </td>
   <td style="text-align:center;"> 2.1.1 </td>
   <td style="text-align:center;"> 2.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/zip </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/zip </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-08-27 </td>
   <td style="text-align:center;"> CRAN (R 4.0.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> zlibbioc </td>
   <td style="text-align:center;"> zlibbioc </td>
   <td style="text-align:center;"> 1.36.0 </td>
   <td style="text-align:center;"> 1.36.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/zlibbioc </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library/zlibbioc </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-28 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.0/Resources/library </td>
  </tr>
</tbody>
</table>

