
#Load Libraries
library("ggplot2")
library("plyr")
library("vegan")

##PERMANOVA Bacteria
PS_normal = readRDS("PS_normal")
PS_normal=transform_sample_counts(PS_normal, sqrt)
PS_normal = subset_taxa(PS_normal, Order!="Chloroplast")
PS_normal = subset_taxa(PS_normal, Kingdom!="Archaea")
PS_normal = subset_taxa(PS_normal, Family!="Mitochondria")
PS_normal = subset_taxa(PS_normal, Kingdom!="Unknown")
#Use base boxplot functin to look at sequencing depth
boxplot(sample_sums(PS_normal)~sample_data(PS_normal)$FieldSite)
boxplot(sample_sums(PS_normal)~sample_data(PS_normal)$SampleType)


#Bacterial ALL
Gradient <- PS_normal@sam_data$rain
Sample_Type <- PS_normal@sam_data$SampleType
bact_ps = PS_normal #Normalized bact data - phyloseq-class
bact_meta_all = as(sample_data(bact_ps), "data.frame") # Extract sample data and coerce to df
bact_all = adonis(distance(bact_ps, method="bray") ~ Gradient*Sample_Type, by = "margin", data = bact_meta_all, permutations = 10000)
bact_all


##PERMANOVA Fungi
FPS_normal = readRDS("fungal_PS_normal")
gradient <- FPS_normal@sam_data$rain
SampleType <- FPS_normal@sam_data$SampleType
fun_ps = transform_sample_counts(FPS_normal, sqrt)
fun_meta_all = as(sample_data(fun_ps), "data.frame") # Extract sample data and coerce to df
fun_all = adonis(distance(fun_ps, method="bray") ~ gradient*SampleType, by = "margin", data = bact_meta_all, permutations = 10000)

fun_all









