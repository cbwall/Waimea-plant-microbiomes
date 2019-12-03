
##Load Libraries
if(!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")	#Load packages from Bioconductor
if(!require(phyloseq)){
	BiocManager::install("phyloseq")
	library(phyloseq)}
if(!require("pacman")){
	install.packages("pacman")
	library(pacman)}
p_load("plyr", "vegan")

##PERMANOVA Bacteria
PS_normal = readRDS("PS_normal")
PS_normal = transform_sample_counts(PS_normal, sqrt)
PS_normal = subset_taxa(PS_normal, Order!="Chloroplast")
PS_normal = subset_taxa(PS_normal, Kingdom!="Archaea")
PS_normal = subset_taxa(PS_normal, Family!="Mitochondria")
PS_normal = subset_taxa(PS_normal, Kingdom!="Unknown")
PS_normal <- subset_taxa(PS_normal, Kingdom != "unknown")	#Remove unclassified kingdoms
PS_normal <- subset_taxa(PS_normal, Kingdom != "Eukaryota")	#Remove eukaryotes

#Bacterial ALL
Site <- as.factor(PS_normal@sam_data$FieldSite)
Sample_Type <- PS_normal@sam_data$SampleType
bact_ps <- phyloseq::distance(PS_normal, method="bray") #Normalized bact data - phyloseq-class
bact_meta_all = as(sample_data(PS_normal), "data.frame") #Extract sample data and coerce to df
bact_all = adonis(bact_ps ~ Site + Sample_Type, by = "margin", data = bact_meta_all, permutations = 10000)
bact_all


##PERMANOVA Fungi
FPS_normal = readRDS("fungal_PS_normal")
fun_ps = transform_sample_counts(FPS_normal, sqrt)
site <- as.factor(FPS_normal@sam_data$FieldSite)
SampleType <- FPS_normal@sam_data$SampleType
fun_ps <- phyloseq::distance(fun_ps, method="bray") #Normalized bact data - phyloseq-class
fun_meta_all = as(sample_data(FPS_normal), "data.frame") # Extract sample data and coerce to df
fun_all = adonis(fun_ps ~ site + SampleType, by = "margin", data = fun_meta_all, permutations = 10000)
fun_all