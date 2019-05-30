#Code by Anthony
#To produce nestedness plots
require(phyloseq)
require(bipartite)
#read in rarefied Bacteria data
rarPhySeq=readRDS("rarPhySeq")
#Merge samples by sample type
agg=merge_samples(rarPhySeq, "SampleType")
#convert to dataframe
aggtab=as.data.frame(otu_table(agg))
#calculate nested temperature
nested=nestedtemp(aggtab)
#read in Fungal data (not rarefied)
fungalphyseq1=readRDS("ITS_Data/fungal_physeq1")
#rarefy to 20000 sequences (drops 3 samples :())
funrarPhySeq=rarefy_even_depth(fungalphyseq1, sample.size=20000)
#merge by sample type
funagg=merge_samples(funrarPhySeq, "SampleType")
#convert to dataframe
funaggtab=as.data.frame(otu_table(funagg))
#calculate nested temp
funnested=nestedtemp(funaggtab)
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
