#Code by Anthony
#To plot mantel correlation between fungi and Bacteri
require(vegan)
require(ggplot2)
#read normalized Bacteria Phylloseq object
PS_normal=readRDS("PS_normal")
#read normalized Fungi Phylloseq object
fungal_PS_normal=readRDS("./ITS_Data/fungal_PS_normal")
#sqrttransform
sqrtfun=transform_sample_counts(fungal_PS_normal, sqrt)
#remove the sample missing from Bacteria
sqrtfun=subset_samples(fungal_PS_normal, CollectionID !="WMEA_01223_Pl_1")
#calculate bray-curtis distance 
funhel.dist=phyloseq::distance(sqrtfun, "bray") 
#remove this sample wich was low abundance in fungi
PS_normal=subset_samples(PS_normal, CollectionID !="WMEA_01223_Pl_1")
#sqrttransform
sqrt=transform_sample_counts(PS_normal, sqrt)
#calculate bray curtis distance
hel.dist=phyloseq::distance(sqrt, "bray") 
#Calculate mantel
m=mantel(hel.dist,funhel.dist)
#make a dataframe for ggplot
mantel=as.data.frame(cbind(c(hel.dist), c(funhel.dist)))
#make EPS file
setEPS()
#create postscript device
postscript("community.eps", width=80, height=80, pointsize=300)
ggplot(mantel,aes(mantel$V1,mantel$V2))+
  geom_point()+
  geom_smooth()+
  theme_bw()+
  ggtitle("Fungi and Bacteria Community Similarity")+
  xlab("Fungal Dissimilarity")+
  ylab("Bacteria Dissimilarity")+
  annotate("text", x=.5, y=.3,label="r=0.434, P=0.001")+
  theme(plot.title = element_text(hjust = 1))
#shutoff device
dev.off()
#png(filename="community.mantel.png",
    units="mm",
    width=80,
    height=80, res=300)  
