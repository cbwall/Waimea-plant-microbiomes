# Code by Anthony; modified by Feresa Corazon
# Purpose: Generate abundance occupancy graphs with bacteria and fungal data
# and site locations which are scaled the mean-annual precipitation 

# Load in packages
require(phyloseq)
require(bipartite)
require(raster)
library("ggplot2")
library("plotrix")
library("viridis")
library("lattice")


#Section 1A. Calling all data and values for bacteria
##################################
## Bacterial Abundance Occupancy##
##################################

rarPhySeq=readRDS("rarPhySeq")
#Calculate distance to shore
sample_data(rarPhySeq)$Shore_dist=pointDistance(cbind(sample_data(rarPhySeq)$Long,sample_data(rarPhySeq)$Lat), c(-158.062848, 21.640741),lonlat=TRUE)
# Merge by sample type
agg=merge_samples(rarPhySeq, "SampleType")
#Load in standarization code
phyloseq_standardize_otu_abundance <- function(physeq, method="total", ...){
  
  ## Check the orientation of the OTU table
  trows <- phyloseq::taxa_are_rows(physeq)
  if(trows == TRUE){ marg <- 2 } else { marg <- 1 }
  
  ## Extact OTU table
  comm <- as(object = phyloseq::otu_table(physeq), Class = "matrix")
  
  ## Standardize community table
  comm_std <- vegan::decostand(comm, method, MARGIN = marg, ...)
  
  ## Replace old otu_table with the new one
  phyloseq::otu_table(physeq) <- phyloseq::otu_table(comm_std, taxa_are_rows = trows)
  
  return(physeq)
}
#How many habitats is each ESV found? Standardize to convert to presence absence
habitats=colSums(otu_table(phyloseq_standardize_otu_abundance(agg, method = "pa")))
habitats
#Merge by site location
sites=merge_samples(rarPhySeq, "FieldSite")
#Standardize to convert to presence absence and calculate site total
site=colSums(otu_table(phyloseq_standardize_otu_abundance(sites, method = "pa")))
site
#Convert the otu table to presence absence in order to calculate range size
binarysite=phyloseq_standardize_otu_abundance(sites, method = "pa")
#Multiply each OTU by it's distance to shore
range=otu_table(binarysite)*sample_data(binarysite)$Shore_dist
#Convert 0 to NA
is.na(range) <- range==0
#Calculate min and max distance to shore
rangespan=apply(range,2,range, na.rm=TRUE)
#Subtract min distance from max distance
rangespan=rangespan[2,]-rangespan[1,]
#Do same thing but log transform range
nozerorangespan=rangespan
is.na(nozerorangespan) <- nozerorangespan==0
#Calculate total abundance
abund=colSums(otu_table(rarPhySeq))     
#Calculate how many samples an OTU is present
occupancy=colSums(otu_table(phyloseq_standardize_otu_abundance(rarPhySeq, method = "pa")))
#Calculate mean abundance per sample (where present)
meanabund=otu_table(rarPhySeq)
is.na(meanabund) <- meanabund==0
meanabund=colMeans(meanabund, na.rm=TRUE)
#Calculate per site habitat diversity (where present)
habpersite=cbind(colSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeq, FieldSite=="1"), method = "pa"))),
      colSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeq, FieldSite=="2"), method = "pa"))),
      colSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeq, FieldSite=="3"), method = "pa"))),
      colSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeq, FieldSite=="4"), method = "pa"))),
      colSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeq, FieldSite=="5"), method = "pa"))),
      colSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeq, FieldSite=="6"), method = "pa"))),
      colSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeq, FieldSite=="7"), method = "pa"))),
      colSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeq, FieldSite=="8"), method = "pa"))),
      colSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeq, FieldSite=="10"), method = "pa"))))
is.na(habpersite) <- habpersite==0
habpersite=rowMeans(habpersite, na.rm=TRUE)


#Section 1B. Calling all data and values for fungi
###############################
## Fungal Abundance Occupancy##
###############################
physeqF=readRDS("fungal_physeq1")
rarPhySeqF=rarefy_even_depth(physeqF, sample.size = 2000)

sample_data(rarPhySeqF)$Shore_dist=pointDistance(cbind(sample_data(rarPhySeqF)$Long,sample_data(rarPhySeqF)$Lat), c(-158.063625, 21.640741),lonlat=TRUE)
aggF=merge_samples(rarPhySeqF, "SampleType")
phyloseq_standardize_otu_abundance <- function(physeq, method="total", ...){
  
  ## Check the orientation of the OTU table
  trows <- phyloseq::taxa_are_rows(physeq)
  if(trows == TRUE){ marg <- 2 } else { marg <- 1 }
  
  ## Extact OTU table
  comm <- as(object = phyloseq::otu_table(physeq), Class = "matrix")
  
  ## Standardize community table
  comm_std <- vegan::decostand(comm, method, MARGIN = marg, ...)
  
  ## Replace old otu_table with the new one
  phyloseq::otu_table(physeq) <- phyloseq::otu_table(comm_std, taxa_are_rows = trows)
  
  return(physeq)
}

#How many habitats is each ESV found? Standardize to convert to presence absence
habitatsF=colSums(otu_table(phyloseq_standardize_otu_abundance(aggF, method = "pa")))
habitatsF
#Merge by site location
sitesF=merge_samples(rarPhySeqF, "FieldSite")
#Standardize to convert to presence absence and calculate site total
siteF=colSums(otu_table(phyloseq_standardize_otu_abundance(sitesF, method = "pa")))
siteF
#Convert the otu table to presence absence in order to calculate range size
binarysiteF=phyloseq_standardize_otu_abundance(sitesF, method = "pa")
#Multiply each OTU by it's distance to shore
rangeF=otu_table(binarysiteF)*sample_data(binarysiteF)$Shore_dist
#Convert 0 to NA
is.na(rangeF) <- rangeF==0
#Calculate min and max distance to shore
rangespanF=apply(rangeF,2,range.default, na.rm=TRUE)
#Subtract min distance from max distance
rangespanF=rangespanF[2,]-rangespanF[1,]
#Do same thing but log transform range
nozerorangespanF=rangespanF
is.na(nozerorangespanF) <- nozerorangespanF==0
#Calculate total abundance
abundF=colSums(otu_table(rarPhySeqF)) 
#Calculate how many samples an OTU is present
occupancyF=colSums(otu_table(phyloseq_standardize_otu_abundance(rarPhySeqF, method = "pa")))
#Calculate mean abundance per sample (where present)
meanabundF=otu_table(rarPhySeqF)
is.na(meanabundF) <- meanabundF==0
meanabundF=rowMeans(meanabundF, na.rm=TRUE)
#Calculate per site habitat diversity (where present)
habpersiteF=cbind(rowSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeqF, FieldSite=="1"), method = "pa"))),
                 rowSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeqF, FieldSite=="2"), method = "pa"))),
                 rowSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeqF, FieldSite=="3"), method = "pa"))),
                 rowSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeqF, FieldSite=="4"), method = "pa"))),
                 rowSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeqF, FieldSite=="5"), method = "pa"))),
                 rowSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeqF, FieldSite=="6"), method = "pa"))),
                 rowSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeqF, FieldSite=="7"), method = "pa"))),
                 rowSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeqF, FieldSite=="8"), method = "pa"))),
                 rowSums(otu_table(phyloseq_standardize_otu_abundance(subset_samples(rarPhySeqF, FieldSite=="10"), method = "pa"))))
is.na(habpersiteF) <- habpersiteF==0
habpersiteF=rowMeans(habpersiteF, na.rm=TRUE)

#Section 2. Making all plots 
######################################################
########## ALL PLOTS AND STATISTICAL TESTS ###########
######################################################

#FIGURE 1. ESV RANGE RELATIONSHIPS
# Create a 3 x 3 plotting matrix
# The next  plots created will be plotted next to each other
par(mfrow = c(3, 2))

# First row: plotting # OF HABITATS VS ESV RANGE
##Plot Bacteria: # of Habitats vs. ESV Range (PLOT 1)
plot(jitter(rangespan,20),jitter(habitats,3),col="#481567FF", cex=.2, xlab="ESV range (m)", ylab="Number of Habitats ESV is Present", main="Bacteria: # of Habitats vs. ESV Range") + abline(lsfit(rangespan,habitats,), col="black", lwd=3)
#When plot panel says locator is active,
#click a vacant area of the plot where you want to annotate the correlation and p-value
temp <- locator(1)
#Add cut-off lines
ablineclip(lm(rangespan~habitats),lty = 2)
#Calculate correlation
cor.test(rangespan, habitats)
text(temp, expression(paste( italic("cor = 0.512 
p-value= 0.00000000000000022"))))
#Plot Fungi: # of Habitats vs. ESV Range (PLOT 2)
plot(jitter(rangespanF,20),jitter(habitatsF,3),col="#95D840FF", cex=.2, xlab="ESV range (M)", ylab="Number of Habitats ESV is Present", main="Fungi: # of Habitats vs. ESV Range") + abline(lsfit(rangespanF,habitatsF,), col="black", lwd=3)
#When plot panel says locator is active,
#click a vacant area of the plot where you want to annotate the correlation and p-value
temp <- locator(1)
#Add cut-off lines
ablineclip(lm(rangespanF~habitatsF), lty = 2)
#Calculate correlation
cor.test(rangespanF, habitatsF)
text(temp, expression(paste( italic("cor = 0.625
p-value= 0.00000000000000022"))))

#Second row: plotting MEAN ABUNDANCE VS ESV RANGE
##Plot Bacteria: Mean Abundance vs. ESV Range (PLOT 3)
plot(jitter(rangespan,20),jitter(log(meanabund),3),col="#481567FF", cex=.2, xlab="ESV range (m)", ylab="Mean Abundance Per Sample (Where Present)", main="Bacteria: Mean Abundance vs. ESV Range")
#When plot panel says locator is active,
#click a vacant area of the plot where you want to annotate the correlation and p-value
temp <- locator(1)
#Calculate correlation
cor.test(rangespan, meanabund)
text(temp,expression(paste( italic("cor = -0.003
p-value = 0.6783"))))
##Plot Fungi: Mean abundance vs.  ESV Range (PLOT 4 )
plot(jitter(rangespanF,20),jitter(log(meanabundF),3),col="#95D840FF", cex=.2, xlab="ESV range (m)", ylab="Mean Abundance Per Sample (Where Present)", main=" Fungi: Mean Abundance vs. ESV Range")+ abline(lsfit(rangespanF,log(meanabundF)), col="black", lwd=3)
#Add cut-off lines
ablineclip(lm(rangespanF~meanabundF), lty = 2)
#When plot panel says locator is active,
#click a vacant area of the plot where you want to annotate the correlation and p-value
temp <- locator(1)
#Calculate correlation
cor.test(rangespanF, meanabundF)
text(temp,expression(paste( italic("cor = 0.105
p-value = 0.000000000000021"))))

#Third row: plotting MEAN HABITAT OCCURENCE VS ESV RANGE
##Plot Bacteria: Mean Habitat Occurence vs. ESV Range (PLOT 5)
plot(jitter(rangespan,20),jitter(habpersite,3),col="#481567FF", cex=.2, xlab="ESV range (m)", ylab="Mean Habitat Occurence (Per Site)", main="Bacteria: Mean Habitat Occurence vs. ESV Range") + abline(lsfit(rangespan,log(meanabund)), col="black", lwd=3) 
#Add cut-off lines
ablineclip(lm(rangespan~habpersite), lty = 2)
#When plot panel says locator is active,
#click a vacant area of the plot where you want to annotate the correlation and p-value
temp <- locator(1)
#Calculate correlation
cor.test(rangespan, habpersite)
text(temp, expression(paste( italic("cor = 0.235
p-value= 0.00000000000000022"))))
##Plot Fungi: Mean Habitat Occurence vs. ESV Range (PLOT 6)
plot(jitter(rangespanF,20),jitter(habpersiteF,3),col="#95D840FF", cex=.2, xlab="ESV range (m)", ylab="Mean Habitat Occurence (Per Site)", main="Fungi: Mean Habitat Occurence vs. ESV Range") + abline(lsfit(rangespanF,log(meanabundF)), col="black", lwd=3)
#Add cut-off lines
ablineclip(lm(rangespanF~habpersiteF), lty = 2)
#When plot panel says locator is active,
#click a vacant area of the plot where you want to annotate the correlation and p-value
temp <- locator(1)
#Calculate correlation
cor.test(rangespanF, habpersiteF)
text(temp, expression(paste( italic("cor = 0.333
p-value= 0.00000000000000022"))))

##FIGURE 2. RANGESPAN OF PHYTOBIOME COMMUNITIES 
#Set layout to combine boxplot and histogram
#Next plots will be created next to each other 
layout(mat = matrix(c(1,2,3,4),2,2),  height = c(1,8))

#First plot: Bacterial rangespan
#Set margin size for boxplot
par(mar=c(0, 4, 1.1, 2.1)) 
boxplot(rangespan, horizontal=TRUE, vertical=TRUE, xaxt="n" , col="#481567FF", main="Bacteria", frame=F) 
#Set margin size for histogram
par(mar=c(4, 4, 1.1, 2.1)) 
hist(rangespan, breaks=40 , col="#481567FF", main="", ylab = "Frequency", xlab="ESV Range (m)", ylim=c(0,5000), xlim=c(0,6000))

#Second Plot: Fungal rangespan
#Set margin size for boxplot
par(mar=c(0, 4, 1.1, 2.1))
boxplot(rangespanF, horizontal=TRUE , xaxt="n" , col="#95D840FF", main="Fungi", frame=F)
#Set margin size for histogram
par(mar=c(4, 4, 1.1, 2.1))
hist(rangespanF, breaks=40 , col="#95D840FF", main="", ylab = "Frequency", xlab="ESV Range (m)", ylim=c(0,3500), xlim=c(0,6000))

#Add title to figure!
title(main = "Rangespan of Phytobiome Communities", outer = TRUE, line = -0.5)

##END


