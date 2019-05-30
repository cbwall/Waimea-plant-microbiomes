##Helen Sung
##Creating Heatmaps for ITS Fungal Data 

###############################################################################################
####HEAT MAP ITS DATA####
###############################################################################################
##set working directory##
setwd("C:/Users/helen/Desktop/BOT662/Analysis/ITS_Data")

###############################################################################################
###Load required packages###
###############################################################################################

##Install 'phyloseq' package
source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')
library(phyloseq)

##Packages to make heat maps##
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
if (!require("plyr")) {
  install.packages("plyr", dependencies = TRUE)
  library(plyr)
}
###############################################################################################
###Load R phyloseq files###
###############################################################################################

PS_normal <- readRDS("fungal_PS_normal") ## normalized phyloseq data
physeq1 <- readRDS("fungal_physeq1")     ## phyloseq data
rarPhySeq <- readRDS("fungal_rarPhySeq") ## rarified phyloseq data 

###############################################################################################
###Fix Data to properly format for Heatmaps###
###############################################################################################

##Substitute all of the FieldSite 10 into 9 in all datasets 
a <- gsub("10", "9", (sample_data(PS_normal)$FieldSite))
b <- gsub("10", "9", (sample_data(physeq1)$FieldSite))
c <- gsub("10", "9", (sample_data(rarPhySeq)$FieldSite))

##Replace old fieldsite names with new fieldsite names in all datasets
a = sample_data(PS_normal)$FieldSite 
b = sample_data(physeq1)$FieldSite 
c = sample_data(rarPhySeq)$FieldSite 

##Organize sample type by levels of ground to air in all datasets
sample_data(PS_normal)$SampleType = factor(sample_data(PS_normal)$SampleType, levels = c('Soil', 'Root', 'Stem', 'Axil', 'Petiole', 'Leaf','Litter','Air'))
sample_data(physeq1)$SampleType = factor(sample_data(physeq1)$SampleType, levels = c('Soil', 'Root', 'Stem', 'Axil', 'Petiole', 'Leaf','Litter','Air'))
sample_data(rarPhySeq)$SampleType = factor(sample_data(rarPhySeq)$SampleType, levels = c('Soil', 'Root', 'Stem', 'Axil', 'Petiole', 'Leaf','Litter','Air'))

###############################################################################################
###Hellinger Transform Data###
###############################################################################################

##Hellinger transform data with Square root for any phyloseq dataset
phyobject.hell <- transform_sample_counts(rarPhySeq, sqrt)  ## I used rarPhyseq because it looked the best for ITS data

###############################################################################################
###Subset Data###
###############################################################################################

## NOTE: Can put any physeq dataset (i.e. PS_normal, physeq1, rarPhySeq, phyobject.hell) into tax_glom(physeqobject, taxrank = "Rank you want to subset to") to subset
## tax_glom() gloms things together at the specified rank and it will get rid of anything that's not classified at that level

##Pool all ESVs by Phylum 
phy_dat <- tax_glom(phyobject.hell, taxrank="Phylum")   ## can replace "phyobject.hell" with any physeq object
phy_dat1 <- tax_glom(rarPhySeq, taxrank="Phylum")

##Pool all ESVs by Order 
ord_dat <- tax_glom(phyobject.hell, taxrank="Order")
ord_dat1 <- tax_glom(rarPhySeq, taxrank="Order")

##Pool all ESVs by Class 
class_dat <- tax_glom(phyobject.hell, taxrank="Class")

##Pool all ESVs by Family
fam_dat <- tax_glom(phyobject.hell, taxrank="Family")

##Pool all ESVs by Genus 
gen_dat <- tax_glom(phyobject.hell, taxrank="Genus")

###############################################################################################
###Plot Heatmap###
###############################################################################################

## NOTE: Can put any Subsetted physeq dataset (i.e. phy_dat, ord_dat, etc.) or non-subsetted physeq dataset (i.e. physeq1, PS_normal, etc.) into plot_heatmap(physeqobject, ...)
## NOTE: If using a subsetted physeq dataset, make sure set the taxa.label = taxrank for whichever taxonomic rank you make heatmap for

##Ploting heatmap for Order in rarPhySeq subsetted data###

##Plot the heatmap clustering by NMDS and Bray curtis distance, labeling sample type and supressing Taxon labels
p1 <- plot_heatmap(ord_dat1, "NMDS", "bray", sample.label = "SampleType", sample.order= "FieldSite", taxa.label="Order", low = "yellow", high = "red", na.value= "white")
##Plot Heatmap organized by Fieldsites (1-9) and grouped by SampleType
q1 <- p1+ facet_grid(~SampleType, scales= "free_x", switch = "x")
##Plot Heatmap with tick marks and FieldSite labels removed
q2 <- q1 + theme(
  axis.text.x = element_blank(),
  axis.ticks = element_blank()
  ) 
plot(q2)

##Save Heatmap as .tiff file 
tiff("ITS_order_heatmap.tiff", units="in", width=12, height=12, res=300) ## Change "ITS_order_heatmap.tiff" to whatever title for the file you choose
plot(q2)
dev.off()





