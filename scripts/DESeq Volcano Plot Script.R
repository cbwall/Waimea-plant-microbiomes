#Code by Melissa Atkins
#Purpose: Differential abundance of ESVs by gradient in fungi and bacteria by class
#DESeq Script

#load libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(Viridis)

#import abundance data for fungi and bacteria
f_allCounts <- read.csv("~/BOT_662_VolcanoPlot/abundance_table_esv_volcanoPlot.csv", row.names=1, quote="")
View(f_allCounts)
b_allCounts <- read.csv("~/BOT_662_VolcanoPlot/abundance_ESV_bacteria.csv", row.names=1, quote="")
View(b_allCounts)
#import metadata for fungi and bacteria
f_conditions <- as.matrix(read.csv("~/BOT_662_VolcanoPlot/VP_Fungi_Metadata.csv", sep = ",", header = TRUE))
View(f_conditions)
b_conditions <- as.matrix(read.csv("~/BOT_662_VolcanoPlot/VP_Bacteria_Metadata.csv", sep = ",", header = TRUE))
View(b_conditions)

#construct DESeqDataSet Object
f_dds <- DESeqDataSetFromMatrix(countData = f_allCounts, colData = f_conditions, design = ~FieldSite)
b_dds <- DESeqDataSetFromMatrix(countData = b_allCounts, colData = b_conditions, design = ~FieldSite)
#check dds object
f_dds
b_dds
#set the factor levels for a control condition
f_dds$FieldSite <- relevel(f_dds$FieldSite, "s1","s2","s3")
b_dds$FieldSite <- relevel(b_dds$FieldSite, "s1","s2","s3")
#transform data to deal with 0's using geoMeans in DESeq
f_counts <- counts(f_dds)
geoMeans = apply(f_counts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
f_ddsGeo <- estimateSizeFactors(f_dds, geoMeans=geoMeans)
b_counts <- counts(b_dds)
geoMeans = apply(b_counts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
b_ddsGeo <- estimateSizeFactors(b_dds, geoMeans=geoMeans)
#Begin analysis
f_dds <- DESeq(f_ddsGeo)
b_dds <- DESeq(b_ddsGeo)
#get results
f_res <- results(f_dds)
b_res <- results(b_dds)
#check results
f_res
b_res

#match results with corresponding ID's for volcano plot
#import taxonomy information
fungi_taxonomy <- read.csv("~/BOT_662_VolcanoPlot/consensus_taxonomy.csv")
bacteria_taxonomy <- read.csv("~/BOT_662_VolcanoPlot/bacteria_parametricDESeq_AllOTU.csv") 
#check that data imported properly
View(fungi_taxonomy)
View(bacteria_taxonomy)
#create data frame to match results with OTU IDs
f_res.df <- as.data.frame(f_res)
f_res.df$OTU_ID <- rownames(f_res.df)
b_res.df <- as.data.frame(b_res)
b_res.df$OTU_ID <- rownames(b_res.df)
#check dataframe
head(f_res.df)
head(b_res.df)
#combine results with OTU IDs
f_res.df.tax <- merge(f_res.df, fungi_taxonomy, by="OTU_ID")
b_res.df.tax <- merge(b_res.df, bacteria_taxonomy, by="OTU_ID")

#create fungi volcano plot of significant classes using ggplot
#create top significant classes from fungi res data
f_significant_classes <- c("Agaricomycetes","Dothideomycetes","Eurotiomycetes",
                           "Lecanoromycetes","Mortierellomycetes","Sordariomycetes","Tremellomycetes")
#subset significant classes from merged data frames
f_res.df.tax.significant <- f_res.df.tax %>%
  subset(Class %in% f_significant_classes)
#drop unused levels from factors
f_res.df.tax.significant$Class <- droplevels(factor(f_res.df.tax.significant$Class))
#check
summary(f_res.df.tax.significant$Class)
#ggplot with cut off lines at p-value of .05 and log2 FC using viridis color scheme
ggplot(data=f_res.df.tax.significant) + 
    geom_point(aes(x = f_res.df.tax.significant$log2FoldChange, y = -log10(f_res.df.tax.significant$pvalue)), color = ifelse(-log10(f_res.df.tax.significant$pvalue) < -log10(0.05),"#39568CFF","#55C667FF"), alpha=0.5) +
    geom_hline(yintercept = -log10(0.05), color="#39568CFF") +
    geom_vline(xintercept = 2, color="#39568CFF") +
    geom_vline(xintercept = -2, color="#39568CFF") +
    labs(x="Log2(FoldChange)", y="-Log10(Pvalue)", color="Class") +
    theme_classic() + ggtitle("Differentially Abundant Fungal ESV's by Class") +
    facet_wrap(~Class)
  
#create bacteria volcano plot of significant classes using ggplot
#create top significant classes from bacteria res data
b_significant_classes <- c("Acidimicrobiia","Acidobacteriia","Actinobacteria","Alphaproteobacteria", 
"Anaerolineae","Bacteroidia","Blastocatellia","Chloroflexia","Deltaproteobacteria","Gammaproteobacteria", 
"Gemmatimondadetes","Latescibacteria_cl","NC10", "OM190","Oxyphotobacteria","Phycisphaerae","Planctomycetacia",
"Thermoleophilia","Verrucomicrobiae")
#subset significant classes from merged data frames
b_res.df.tax.significant <- b_res.df.tax %>%
  subset(Class %in% b_significant_classes)
#drop unused levels from factors
b_res.df.tax.significant$Class <- droplevels(factor(b_res.df.tax.significant$Class))
#check
summary(b_res.df.tax.significant$Class)
#ggplot with cut off lines at p-value of .05 and log2 FC using viridis color scheme
ggplot(data=b_res.df.tax.significant) + 
    geom_point(aes(x = b_res.df.tax.significant$log2FoldChange, y = -log10(b_res.df.tax.significant$pvalue)), color = ifelse(-log10(b_res.df.tax.significant$pvalue) < -log10(0.05),"#39568CFF","#55C667FF"), alpha=0.5) +
    geom_hline(yintercept = -log10(0.05), color="#39568CFF") +
    geom_vline(xintercept = 2, color="#39568CFF") +
    geom_vline(xintercept = -2, color="#39568CFF") +
    labs(x="Log2(FoldChange)", y="-Log10(Pvalue)", color="Class") +
    theme_classic() + ggtitle("Differentially Abundant Bacterial ESV's by Class") +
    facet_wrap(~Class)

#Export plots
tiff(file="VP.bacteria.classFacet.tiff",width=2100,height=1890,res=300)
plot(ggplot(data=b_res.df.tax.significant) + 
    geom_point(aes(x = b_res.df.tax.significant$log2FoldChange, y = -log10(b_res.df.tax.significant$pvalue)), color = ifelse(-log10(b_res.df.tax.significant$pvalue) < -log10(0.05),"#39568CFF","#55C667FF"), alpha=0.5) +
    geom_hline(yintercept = -log10(0.05), color="#39568CFF") +
    geom_vline(xintercept = 2, color="#39568CFF") +
    geom_vline(xintercept = -2, color="#39568CFF") +
    labs(x="Log2(FoldChange)", y="-Log10(Pvalue)", color="Class") +
    theme_classic() + ggtitle("Differentially Abundant Bacterial ESV's by Class") +
    facet_wrap(~Class))
dev.off()

tiff(file="VP.fungi.classFacet.tiff",width=2100,height=1890,res=300)
plot(ggplot(data=f_res.df.tax.significant) + 
    geom_point(aes(x = f_res.df.tax.significant$log2FoldChange, y = -log10(f_res.df.tax.significant$pvalue)), color = ifelse(-log10(f_res.df.tax.significant$pvalue) < -log10(0.05),"#39568CFF","#55C667FF"), alpha=0.5) +
    geom_hline(yintercept = -log10(0.05), color="#39568CFF") +
    geom_vline(xintercept = 2, color="#39568CFF") +
    geom_vline(xintercept = -2, color="#39568CFF") +
    labs(x="Log2(FoldChange)", y="-Log10(Pvalue)", color="Class") +
    theme_classic() + ggtitle("Differentially Abundant Fungal ESV's by Class") +
    facet_wrap(~Class))
dev.off()

#Export data
write.table(b_res.df.tax.significant, file="bacteria_parametricDESeq_significant.tsv", sep="\t",eol="\n")
write.table(f_res.df.tax.significant, file="fungi_parametricDESeq_significant.tsv", sep="\t",eol="\n")