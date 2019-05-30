#Code by Maria and Helen
#For NMDS plots and heat maps by Sample Type
library(ggplot2)
library(fossil)
library("phyloseq")
library(vegan)
library(RColorBrewer)
library("plyr")
library(devtools)
library(ggpubr)
library(gridExtra)
devtools::install_bitbucket("graumannlabtools/multipanelfigure")
#read in normalized bacterial data
PS_normal<-readRDS("PS_normal")
#read in bacteria phyloseq data
physeq1 <- readRDS("physeq1")  
#read in rarified phyloseq data
rarPhySeq <- readRDS("rarPhySeq") 
#Remove site 9
PS_normal=subset_samples(PS_normal, FieldSite!="s9")
#read in normalized fungal taxa
fungal_norm<-readRDS("fungal_PS_normal")
#read in fungal phyloseq data
funphys <- readRDS("fungal_physeq1") 
#read in rarified phyloseq data
rarfun <- readRDS("fungal_rarPhySeq")

####################################
##########NMDS PLOTS##############
#hellinger transform bacteria
phyobj.hell=transform_sample_counts(PS_normal, sqrt)
#hellinger transform fungus
fungal.hell<-transform_sample_counts(fungal_norm, sqrt)
##Organize sample type by levels of ground to air in all datasets
sample_data(PS_normal)$SampleType = factor(sample_data(PS_normal)$SampleType, levels = c('Soil', 'Root', 'Stem', 'Axil', 'Petiole', 'Leaf','Litter','Air'))
sample_data(physeq1)$SampleType = factor(sample_data(physeq1)$SampleType, levels = c('Soil', 'Root', 'Stem', 'Axil', 'Petiole', 'Leaf','Litter','Air'))
sample_data(rarPhySeq)$SampleType = factor(sample_data(rarPhySeq)$SampleType, levels = c('Soil', 'Root', 'Stem', 'Axil', 'Petiole', 'Leaf','Litter','Air'))
##Organize sample type by levels of ground to air in all datasets for fungus
sample_data(fungal_norm)$SampleType = factor(sample_data(fungal_norm)$SampleType, levels = c('Soil', 'Root', 'Stem', 'Axil', 'Petiole', 'Leaf','Litter','Air'))
sample_data(funphys)$SampleType = factor(sample_data(funphys)$SampleType, levels = c('Soil', 'Root', 'Stem', 'Axil', 'Petiole', 'Leaf','Litter','Air'))
sample_data(rarfun)$SampleType = factor(sample_data(rarfun)$SampleType, levels = c('Soil', 'Root', 'Stem', 'Axil', 'Petiole', 'Leaf','Litter','Air'))

#bacteria ordination
phyobj.ord<-ordinate(phyobj.hell, "NMDS", "bray")
#fungus ordination
funhell.ord<-ordinate(fungal.hell, "NMDS", "bray")
#define colorblind friendly color palette 
cbPalette1<- c("#F0E442","#D55E00","#E69F00","#009E73","#999999","#CC79A7","#000000","#56B4E9","#0072B2")
#bacteria NMDS by sample type
plot1<-plot_ordination(phyobj.hell, phyobj.ord, type="samples", color="SampleType")+ theme(legend.position="none")  #coord_equal(ratio=2) + theme_bw()
plot1+
  scale_fill_manual(values=cbPalette1) +  
  scale_colour_manual(values=cbPalette1)+
  stat_ellipse(geom="polygon",type="t", alpha=0.2, aes(fill=SampleType)) +
  theme(legend.position="none") + coord_equal(ratio=2)+
  annotate_figure(plot1,
                  top = text_grob("Bacteria", color = "black", face = "bold", size = 14))+
  theme_bw()
#fungal NMDS by sample type
plot2<-plot_ordination(fungal.hell, funhell.ord,type="samples", color="SampleType")
plot2 +
  scale_fill_manual(values=cbPalette1) +  
  scale_colour_manual(values=cbPalette1)+
  stat_ellipse(geom="polygon",type="t", alpha=0.2, aes(fill=SampleType)) +
  theme(legend.position="none")  coord_equal(ratio=2)+
  annotate_figure(plot2ell,op = text_grob("Fungi", color = "black", face = "bold", size = 14))
  theme_bw()


###################################
########HEATMAP PLOTS#############

####BACTERIA 16S###############

##Substitute all of the FieldSite 10 into 9 in all datasets 
a <- gsub("s10", "s9", (sample_data(PS_normal)$FieldSite))
b <- gsub("s10", "s9", (sample_data(physeq1)$FieldSite))
c <- gsub("s10", "s9", (sample_data(rarPhySeq)$FieldSite))
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
phyobject.hellrar <- transform_sample_counts(rarPhySeq, sqrt)  ## I used rarPhyseq because it looked the best for ITS data
###############################################################################################
###Subset Data###
###############################################################################################

## NOTE: Can put any physeq dataset (i.e. PS_normal, physeq1, rarPhySeq, phyobject.hell) into tax_glom(physeqobject, taxrank = "Rank you want to subset to") to subset
## tax_glom() gloms things together at the specified rank and it will get rid of anything that's not classified at that level

##Pool all ESVs by Phylum 
phy_dat <- tax_glom(phyobject.hell, taxrank="Phylum")   ## subsetting for Hellinger transformed normalized data
phy_dat1 <- tax_glom(phyobject.hellrar, taxrank="Phylum")   ## subsetting for Hellinger transformed rarified data
phy_dat2 <- tax_glom(rarPhySeq, taxrank="Phylum") ## subsetting for rarPhySeq data
phy_dat3 <- tax_glom(PS_normal, taxrank="Phylum") ## subsetting for normalized data

##Pool all ESVs by Order 
ord_dat <- tax_glom(phyobject.hell, taxrank="Order")  
ord_dat1 <- tax_glom(phyobject.hell1, taxrank="Order")   
ord_dat2 <- tax_glom(rarPhySeq, taxrank="Order")
ord_dat3 <- tax_glom(PS_normal, taxrank="Order")


## NOTE: Can put any Subsetted physeq dataset (i.e. phy_dat, ord_dat, etc.) or non-subsetted physeq dataset (i.e. physeq1, PS_normal, etc.) into plot_heatmap(physeqobject, ...)
## NOTE: If using a subsetted physeq dataset, make sure set the taxa.label = taxrank for whichever taxonomic rank you make heatmap for

##Ploting heatmap for Order in rarPhySeq subsetted data###

##Plot the heatmap clustering by NMDS and Bray curtis distance, labeling sample type and supressing Taxon labels
p1 <- plot_heatmap(phy_dat3, "NMDS", "bray", sample.label = "SampleType", sample.order= "FieldSite", taxa.label="Phylum", low = "yellow", high = "red", na.value= "white")
##Plot Heatmap organized by Fieldsites (1-9) and grouped by SampleType
q1 <- p1+ facet_grid(~SampleType, scales= "free_x", switch = "x")
##Plot Heatmap with tick marks and FieldSite labels removed
q2 <- q1 + theme(
  axis.text.x = element_blank(),
  axis.ticks = element_blank()
) 
plot(q2)


#########FUNGAL ITS##################


##Substitute all of the FieldSite 10 into 9 in all datasets 
d <- gsub("s10", "s9", (sample_data(fungal_norm)$FieldSite))
e <- gsub("s10", "s9", (sample_data(funphys)$FieldSite))
f <- gsub("s10", "s9", (sample_data(rarfun)$FieldSite))

##Replace old fieldsite names with new fieldsite names in all datasets
d = sample_data(fungal_norm)$FieldSite 
e = sample_data(funphys)$FieldSite 
f = sample_data(rarfun)$FieldSite 

##Organize sample type by levels of ground to air in all datasets
sample_data(fungal_norm)$SampleType = factor(sample_data(fungal_norm)$SampleType, levels = c('Soil', 'Root', 'Stem', 'Axil', 'Petiole', 'Leaf','Litter','Air'))
sample_data(funphys)$SampleType = factor(sample_data(funphys)$SampleType, levels = c('Soil', 'Root', 'Stem', 'Axil', 'Petiole', 'Leaf','Litter','Air'))
sample_data(rarfun)$SampleType = factor(sample_data(rarfun)$SampleType, levels = c('Soil', 'Root', 'Stem', 'Axil', 'Petiole', 'Leaf','Litter','Air'))

###############################################################################################
###Hellinger Transform Data###
###############################################################################################

##Hellinger transform data with Square root for any phyloseq dataset
rarfun.hell <- transform_sample_counts(rarfun, sqrt)  ## I used rarfun because it looked the best for ITS data

###############################################################################################
###Subset Data###
###############################################################################################

## NOTE: Can put any physeq dataset (i.e. fungal_norm, funphys, rarfun, phyobject.hell) into tax_glom(physeqobject, taxrank = "Rank you want to subset to") to subset
## tax_glom() gloms things together at the specified rank and it will get rid of anything that's not classified at that level

##Pool all ESVs by Phylum 
funphy_dat <- tax_glom(rarfun.hell, taxrank="Phylum")   ## can replace "phyobject.hell" with any physeq object
funphy_dat1 <- tax_glom(rarfun, taxrank="Phylum")
funphy_dat2<- tax_glom(fungal.hell, taxrank= "Phylum")

##Pool all ESVs by Order 
funord_dat <- tax_glom(rarfun.hell, taxrank="Order")
funord_dat1 <- tax_glom(rarfun, taxrank="Order")
funord_dat2 <- tax_glom(fungal.hell, taxrank="Order")

##Pool all ESVs by Class 
funclass_dat <- tax_glom(rarfun.hell, taxrank="Class")
funclass_dat1 <- tax_glom(rarfun, taxrank="Class")
##Pool all ESVs by Family
funfam_dat <- tax_glom(rarfun.hell, taxrank="Family")


##Pool all ESVs by Genus 
fungen_dat <- tax_glom(rarfun.hell, taxrank="Genus")

###############################################################################################
###Plot Heatmap###
###############################################################################################

## NOTE: Can put any Subsetted physeq dataset (i.e. phy_dat, ord_dat, etc.) or non-subsetted physeq dataset (i.e. funphys, fungal_norm, etc.) into plot_heatmap(physeqobject, ...)
## NOTE: If using a subsetted physeq dataset, make sure set the taxa.label = taxrank for whichever taxonomic rank you make heatmap for

##Ploting heatmap for Order in rarfun subsetted data###

##Plot the heatmap clustering by NMDS and Bray curtis distance, labeling sample type and supressing Taxon labels
fp1 <- plot_heatmap(funord_dat1, "NMDS", "bray", sample.label = "SampleType", sample.order= "FieldSite", taxa.label="Order", low = "yellow", high = "red", na.value= "white")
##Plot Heatmap organized by Fieldsites (1-9) and grouped by SampleType
fq1 <- fp1+ facet_grid(~SampleType,scales= "free_x", switch = "x")
##Plot Heatmap with tick marks and FieldSite labels removed
fq2 <- fq1 + theme(
  axis.text.x = element_blank(),
  axis.ticks = element_blank()
) 
plot(fq2)

###################################################################################################################
##################### BELOW HERE I AM LOST ON HOW TO MAKE PRETTY FIGURE############################################

####MAKING A COMPOSITE FIGURE#########
nmds<- ggarrange(qplot1,qplot2 + rremove("x.text"), common.legend = TRUE, legend="bottom", ncol=1, nrow=2)
heatmaps<-ggarrange(q2+rremove("xlab"),fq2 + rremove("xlab"), common.legend = FALSE, legend="right", ncol=1, nrow=2)
fullfig<-ggarrange(nmds, heatmaps + rremove("x.text"))
ggdraw() +
  draw_plot(nmds, x = 0, y = .5, width = .4, height = .4) +
  draw_plot(heatmaps, x = 0, y = 0, width = .6, height = 0.6) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.5))
plot(fullfig)         



######THIS IS ANTHONY'S CODE FOR EPS FILE##################################

#create EPS file
setEPS()
#using a postscript device
postscript("nestedplots.eps")
#make a two row composite figure, reguce spacing between inner and outer margins
par(mfrow=c(2,1), mai=c(.1,1.,1,1), oma=c(.1,1,.1,1))
#plot the fungal nested figure, suppress taxon names
plot(funnested, kind="incid", names=c(TRUE,FALSE), col=c("white", "green"), lwd=3, main="Fungi Nested Plot")
#add stats as margin text
mtext("Nested Temp=42.5, P=0.001")
#plot the Bacteria nested figure, suppress taxon names
plot(nested, kind="incid", names=c(TRUE,FALSE), col=c("white", "purple"), lwd=3, main="Bacteria Nested Plot")
#add stats are margin text
mtext("Nested Temp=39.8, P=0.001")
dev.off()
