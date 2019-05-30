#######################################################################################
######################## Mantel Tests and Slope Plots By Plant Part ###################
################################ Chad Wilhite 04/25/19 ################################
#######################################################################################

# Do you need to install any of the packages?
# install.packages(c("fossil"))

# Load required packages
require(fossil)
require(vegan)
require(phyloseq)
require(ggplot2)

# Set your working directory!
setwd("C:\\Users\\cwilh\\Documents\\01_School\\01_Graduate School\\01_UH Manoa\\Coursework\\BOT 662 High Thruput Sequencing\\Corrected_Data")


####################################
#### Mantel Tests by SampleType ####
####################################


#
##
###
#### Bacteria Analysis


# Load RDS file of normalized bacterial reads
PS_normal = readRDS("PS_normal") #"PS_normal"


# Hellinger (square root) transform the data!
PS_normal.hell = transform_sample_counts(PS_normal, sqrt)

# View SampleType names
unique( get_variable(PS_normal.hell, sample_variables(PS_normal.hell)[17]) )
#
##
### Subset out bacterial data by SampleType

# Subset bacterial data by plant part - the hard way... xD
stem_otu = subset_samples(PS_normal.hell, SampleType=="Stem")
root_otu = subset_samples(PS_normal.hell, SampleType=="Root")
air_otu  = subset_samples(PS_normal.hell, SampleType=="Air")
leaf_otu = subset_samples(PS_normal.hell, SampleType=="Leaf")
soil_otu = subset_samples(PS_normal.hell, SampleType=="Soil")
litt_otu = subset_samples(PS_normal.hell, SampleType=="Litter")
axil_otu = subset_samples(PS_normal.hell, SampleType=="Axil")
peti_otu = subset_samples(PS_normal.hell, SampleType=="Petiole")


#
##
### Generate Dissimilarity Matrices for bacterial OTU data

# Generate a dissimilarity matrix for grouped and individual plant part bacterial OTU data 
all_otu.dist  = distance(PS_normal.hell, "bray")
stem_otu.dist = distance(stem_otu, "bray") 
root_otu.dist = distance(root_otu, "bray") 
air_otu.dist  = distance(air_otu,  "bray") 
leaf_otu.dist = distance(leaf_otu, "bray") 
soil_otu.dist = distance(soil_otu, "bray")
litt_otu.dist = distance(litt_otu, "bray")
axil_otu.dist = distance(axil_otu, "bray") 
peti_otu.dist = distance(peti_otu, "bray")



#
##
### Create correctly ordered sized geographic distance dissimilarity matrix


#Gather bacterial spatial data from the phyloseq object
all.geo  = cbind(sample_data(PS_normal.hell)$Lon, sample_data(PS_normal.hell)$Lat)
stem.geo = cbind(sample_data(stem_otu)$Lon, sample_data(stem_otu)$Lat)
root.geo = cbind(sample_data(root_otu)$Lon, sample_data(root_otu)$Lat)
air.geo  = cbind(sample_data(air_otu)$Lon,  sample_data(air_otu)$Lat)
leaf.geo = cbind(sample_data(leaf_otu)$Lon, sample_data(leaf_otu)$Lat)
soil.geo = cbind(sample_data(soil_otu)$Lon, sample_data(soil_otu)$Lat)
litt.geo = cbind(sample_data(litt_otu)$Lon, sample_data(litt_otu)$Lat)
axil.geo = cbind(sample_data(axil_otu)$Lon, sample_data(axil_otu)$Lat) 
peti.geo = cbind(sample_data(peti_otu)$Lon, sample_data(peti_otu)$Lat)



# Generate a dissimilarity matrix for grouped (all plant parts) and by plant part 
# 	bacterial geographic distance data 
all.geodist  = earth.dist(all.geo)
stem.geodist = earth.dist(stem.geo)
root.geodist = earth.dist(root.geo)
air.geodist  = earth.dist(air.geo)
leaf.geodist = earth.dist(leaf.geo)
soil.geodist = earth.dist(soil.geo)
litt.geodist = earth.dist(litt.geo)
axil.geodist = earth.dist(axil.geo)
peti.geodist = earth.dist(peti.geo)



#Mantel test by SampleType
all_bacteria_mantel  = mantel(log( all.geodist+1),log( all_otu.dist), permutations= 999)
stem_bacteria_mantel = mantel(log(stem.geodist+1),log(stem_otu.dist), permutations= 999)
root_bacteria_mantel = mantel(log(root.geodist+1),log(root_otu.dist), permutations= 999)
air_bacteria_mantel  = mantel(log( air.geodist+1),log( air_otu.dist), permutations= 999)
leaf_bacteria_mantel = mantel(log(leaf.geodist+1),log(leaf_otu.dist), permutations= 999)
soil_bacteria_mantel = mantel(log(soil.geodist+1),log(soil_otu.dist), permutations= 999)
litt_bacteria_mantel = mantel(log(litt.geodist+1),log(litt_otu.dist), permutations= 999)
axil_bacteria_mantel = mantel(log(axil.geodist+1),log(axil_otu.dist), permutations= 999)
peti_bacteria_mantel = mantel(log(peti.geodist+1),log(peti_otu.dist), permutations= 999)




#Call the bacterial results
all_bacteria_mantel
stem_bacteria_mantel
root_bacteria_mantel
air_bacteria_mantel
leaf_bacteria_mantel
soil_bacteria_mantel
litt_bacteria_mantel
axil_bacteria_mantel
peti_bacteria_mantel



#############################################
################ Fungal #####################
#############################################

#Load RDS file of fungal reads
FPS_normal = readRDS("fungal_PS_normal") #"fungal_PS_normal"



#square root transform (Hellinger)
FPS_normal.hell = transform_sample_counts(FPS_normal, sqrt)




#Subset Fungal Data by SampleType
fstem_otu = subset_samples(FPS_normal.hell, SampleType=="Stem")
froot_otu = subset_samples(FPS_normal.hell, SampleType=="Root")
fair_otu  = subset_samples(FPS_normal.hell, SampleType=="Air")
fleaf_otu = subset_samples(FPS_normal.hell, SampleType=="Leaf")
fsoil_otu = subset_samples(FPS_normal.hell, SampleType=="Soil")
flitt_otu = subset_samples(FPS_normal.hell, SampleType=="Litter")
faxil_otu = subset_samples(FPS_normal.hell, SampleType=="Axil")
fpeti_otu = subset_samples(FPS_normal.hell, SampleType=="Petiole")


#Generate Distance Matrix
fall_otu.dist  = distance(FPS_normal.hell, "bray")
fstem_otu.dist = distance(fstem_otu, "bray")
froot_otu.dist = distance(froot_otu, "bray")
fair_otu.dist  = distance( fair_otu, "bray")
fleaf_otu.dist = distance(fleaf_otu, "bray")
fsoil_otu.dist = distance(fsoil_otu, "bray")
flitt_otu.dist = distance(flitt_otu, "bray")
faxil_otu.dist = distance(faxil_otu, "bray")
fpeti_otu.dist = distance(fpeti_otu, "bray")




#Gather fungal spatial data from the phyloseq object
fall.geo  = cbind(sample_data(FPS_normal.hell)$Lon, sample_data(FPS_normal.hell)$Lat)
fstem.geo = cbind(sample_data(fstem_otu)$Lon, sample_data(fstem_otu)$Lat)
froot.geo = cbind(sample_data(froot_otu)$Lon, sample_data(froot_otu)$Lat)
fair.geo  = cbind(sample_data( fair_otu)$Lon, sample_data( fair_otu)$Lat)
fleaf.geo = cbind(sample_data(fleaf_otu)$Lon, sample_data(fleaf_otu)$Lat)
fsoil.geo = cbind(sample_data(fsoil_otu)$Lon, sample_data(fsoil_otu)$Lat)
flitt.geo = cbind(sample_data(flitt_otu)$Lon, sample_data(flitt_otu)$Lat)
faxil.geo = cbind(sample_data(faxil_otu)$Lon, sample_data(faxil_otu)$Lat)
fpeti.geo = cbind(sample_data(fpeti_otu)$Lon, sample_data(fpeti_otu)$Lat)


#Create correct sized fungal geographic distance dissimilarity matrix
fall.geodist  = earth.dist( fall.geo)
fstem.geodist = earth.dist(fstem.geo)
froot.geodist = earth.dist(froot.geo)
fair.geodist  = earth.dist( fair.geo)
fleaf.geodist = earth.dist(fleaf.geo)
fsoil.geodist = earth.dist(fsoil.geo)
flitt.geodist = earth.dist(flitt.geo)
faxil.geodist = earth.dist(faxil.geo)
fpeti.geodist = earth.dist(fpeti.geo)





#Mantel test by SampleType
fall_mantel  = mantel(log( fall.geodist+1),log( fall_otu.dist), permutations= 999)
fstem_mantel = mantel(log(fstem.geodist+1),log(fstem_otu.dist), permutations= 999)
froot_mantel = mantel(log(froot.geodist+1),log(froot_otu.dist), permutations= 999)
fair_mantel  = mantel(log( fair.geodist+1),log( fair_otu.dist), permutations= 999)
fleaf_mantel = mantel(log(fleaf.geodist+1),log(fleaf_otu.dist), permutations= 999)
fsoil_mantel = mantel(log(fsoil.geodist+1),log(fsoil_otu.dist), permutations= 999)
flitt_mantel = mantel(log(flitt.geodist+1),log(flitt_otu.dist), permutations= 999)
faxil_mantel = mantel(log(faxil.geodist+1),log(faxil_otu.dist), permutations= 999)
fpeti_mantel = mantel(log(fpeti.geodist+1),log(fpeti_otu.dist), permutations= 999)




#call the fungal results
fall_mantel
fstem_mantel
froot_mantel
fair_mantel
fleaf_mantel
fsoil_mantel
flitt_mantel
faxil_mantel
fpeti_mantel







###################################################
# Save Mantel Correlation Statistics and P Values #
###################################################

#Call plant part results
mant_cor = c(stem_bacteria_mantel$statistic, root_bacteria_mantel$statistic, 
	air_bacteria_mantel$statistic, leaf_bacteria_mantel$statistic, 
	soil_bacteria_mantel$statistic, litt_bacteria_mantel$statistic, 
	axil_bacteria_mantel$statistic, peti_bacteria_mantel$statistic, 
	fstem_mantel$statistic, froot_mantel$statistic, fair_mantel$statistic, 
	fleaf_mantel$statistic, fsoil_mantel$statistic, flitt_mantel$statistic, 
	faxil_mantel$statistic, fpeti_mantel$statistic)

#Call plant part p-values
p_val = c(stem_bacteria_mantel$signif, root_bacteria_mantel$signif, 
	air_bacteria_mantel$signif, leaf_bacteria_mantel$signif, 
	soil_bacteria_mantel$signif, litt_bacteria_mantel$signif, 
	axil_bacteria_mantel$signif, peti_bacteria_mantel$signif, 
	fstem_mantel$signif, froot_mantel$signif, fair_mantel$signif, 
	fleaf_mantel$signif, fsoil_mantel$signif, flitt_mantel$signif, 
	faxil_mantel$signif, fpeti_mantel$signif)

#Create results data frame with mantel R and p-values for each plant part
results = data.frame( c( rep('Bacterial',8), rep('Fungal',8) ), 
	c( rep( c('Stem', 'Root', 'Air', 'Leaf', 'Soil', 'Litter', 'Axil', 'Petiole'),2 ) ),
	mant_cor, p_val
	)

#Name the columns intelligible names
colnames(results) = c('type', 'Part', 'mant_cor', 'p_val')

#Adjust p-values for repeated tests on the within plant part groups for bacteria
bact.cor.p = p.adjust(results[results$type == 'Bacterial',4], method = 'bonferroni')
bact.cor.p = data.frame(cor.p = bact.cor.p, type = rep('Bacterial',8), 
	Part = as.character(unique(results$Part)) )

#Adjust p-values for repeated tests on the within plant part groups for fungi
fung.cor.p = p.adjust(results[results$type == 'Fungal',4], method = 'bonferroni')
fung.cor.p = data.frame(cor.p = fung.cor.p, type = rep('Fungal',8), 
	Part = as.character(unique(results$Part)) )

#Join corrected p-values to for fungus and bacteria together
cor.p = rbind(bact.cor.p, fung.cor.p)


#Join corrected p-values to results data frame
results = merge(results, cor.p)



#Add overall test (all parts together) to data frame
results2 = rbind( results, data.frame( type = c('Bacterial', 'Fungal'), Part = rep('All',2),
	mant_cor = c(all_bacteria_mantel$statistic, fall_mantel$statistic), 
	p_val = c(all_bacteria_mantel$signif, fall_mantel$signif),
	cor.p = rep('NA', 2))
	)



#Save all result data as a csv
write.csv(results2, file = 'Mantel_Plant_Part_Results_Hell_Bray.csv')















####################################################
################### Slopes plot ####################
####################################################


# The colorblind friendly palette with grey:
cbPalette = c('#999999', '#E69F00', '#56B4E9', 
	'#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7')



###################
#### Bacterial ####
###################

#Combine all plant part bacterial OTU dissimilarity matrices
ball.otudist = c(stem_otu.dist, root_otu.dist, air_otu.dist,
	leaf_otu.dist, soil_otu.dist, litt_otu.dist, axil_otu.dist, peti_otu.dist)


#Combine all plant part bacterial geographic dissimilarity matrices
ball.geodist = c(stem.geodist, root.geodist, air.geodist,
	leaf.geodist, soil.geodist, litt.geodist, axil.geodist, peti.geodist)

#Set up properly ordered names by plant part
bplpt = c( rep('Stem', 36), rep('Root', 36), rep('Air', 36), rep('Leaf', 36),
	 rep('Soil', 36), rep('Litter', 36), rep('Axil', 36), rep('Petiole', 36) ) 


#Combine the OTU and geographic dissimilarity vectors with their plant part names
bac.all.dist = data.frame(b.geo.dist = ball.geodist, b.otu.dist = ball.otudist, b.pl.part = bplpt)


#Give each plant part our particular colorblind friendly color
b.part.col = data.frame(cbPalette, levels(bac.all.dist$b.pl.part))
colnames(b.part.col) = c('b.col', 'b.pl.part')
b.part.col$b.col = as.character(b.part.col$b.col)

#Merge assigned plant part colors in the main data frame
bac.all.dist = merge(bac.all.dist, b.part.col)


#Create and save a plot of slopes for each plant part!
bac.slope.plot = ggplot(bac.all.dist, aes(x=b.geo.dist, y=b.otu.dist, color = b.pl.part)) +
	scale_color_manual( values = unique(bac.all.dist$b.col), name = "Sample Type") +
	geom_smooth(method = lm, se = FALSE) + geom_point() +
	theme(axis.line = element_line(color = "black"), legend.background = element_rect(),
		panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), legend.key=element_blank(),
		legend.justification=c(1,0), legend.position=c(1,0)) +
	scale_x_continuous(name = "Pairwise Geographic Dissimilarity", breaks = 0:5) + 
	ylab("Pairwise Community Dissimilarity")

#Show plot
bac.slope.plot

#Save plot
ggsave(filename = "Bacterial_Slopes_by_Plant_Part.eps", width = 10, height = 10)

dev.off()













################
#### Fungal ####
################

#Combine all fungal plant part OTU dissimilarity matrices
fall.otudist = c(fstem_otu.dist, froot_otu.dist, fair_otu.dist,
	fleaf_otu.dist, fsoil_otu.dist, flitt_otu.dist, faxil_otu.dist, fpeti_otu.dist)


#Combine all fungal plant part geographic dissimilarity matrices
fall.geodist = c(fstem.geodist, froot.geodist, fair.geodist,
	fleaf.geodist, fsoil.geodist, flitt.geodist, faxil.geodist, fpeti.geodist)

#Set up properly ordered names by plant part
plpt = c( rep('Stem', 36), rep('Root', 36), rep('Air', 36), rep('Leaf', 36),
	 rep('Soil', 36), rep('Litter', 36), rep('Axil', 36), rep('Petiole', 36) ) 


#Combine the fungal OTU and geographic dissimilarity vectors with their plant part names
fun.all.dist = data.frame(geo.dist = fall.geodist, otu.dist = fall.otudist, pl.part = plpt)


#Give each plant part our particular colorblind friendly color
part.col = data.frame(cbPalette, levels(fun.all.dist$pl.part))
colnames(part.col) = c('col', 'pl.part')
part.col$col = as.character(part.col$col)

fun.all.dist = merge(fun.all.dist, part.col)


#Make and save a plot of slopes for each fungal plant part!
fun.slope.plot = ggplot(fun.all.dist, aes(x=geo.dist, y=otu.dist, color = pl.part)) +
	scale_color_manual( values = unique(fun.all.dist$col), name = "Sample Type") +
	geom_smooth(method = lm, se = FALSE) + geom_point() +
	theme(axis.line = element_line(color = "black"), legend.background = element_rect(),
		panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), legend.key=element_blank(),
		legend.justification=c(1,0), legend.position=c(1,0)) +
	 scale_x_continuous(name = "Pairwise Geographic Dissimilarity", breaks = 0:5) + 
	ylab("Pairwise Community Dissimilarity")

#Show plot
fun.slope.plot

#Save plot
ggsave(filename = "Fungal_Slopes_by_Plant_Part.eps", width = 10, height = 10)

dev.off()




