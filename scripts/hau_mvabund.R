####			Multi- and Univariate Analyses of Microbial ESVs
####				Associated with Hibiscus tiliaceus
####				Along an Environmental Gradient

###						Jared Bernard

##Load packages (automatically install if needed)
if(!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")	#Load packages from Bioconductor
if(!require(phyloseq)){
	BiocManager::install("phyloseq")
	library(phyloseq)}
if(!require(qvalue)){
	BiocManager::install("qvalue")
	library(qvalue)}
if(!require("pacman")){
	install.packages("pacman")
	library(pacman)}
p_load("plyr", "reshape2", "mvabund", "graphics", "MASS", "effects",
	"BaylorEdPsych", "ResourceSelection", "dplyr", "stringr", "viridis",
	"multipanelfigure", "magrittr")


memory.limit(100000)	#May be necessary for faster manyglm runs


###BACTERIAL DATA

rarPhySeq = readRDS("../data/rarPhySeq")	#Read in rarefied count data
bac_count_df <- psmelt(rarPhySeq)	#Make it an ordinary dataframe, not phyloseq
bac_count_df$taxon <- paste(bac_count_df$Class, bac_count_df$OTU)	#Create field combining OTU and Class
str(bac_count_df)	#Check structure
write.csv(bac_count_df, file = "../data/bac_count_df.csv", row.names = F)	#Save
#bac_count_df <- read.csv("../data/bac_count_df.csv", header = T)	#Read in if saved


##Subsample bacterial dataframe
newdata <- subset(bac_count_df, Abundance >= (0.05 * max(bac_count_df$Abundance)),
	select = c(rain, SampleType, taxon, Abundance))	#Remove lowest 5% and subset
tmp <- ddply(newdata, .(SampleType, taxon), transform,
	newid = paste(SampleType, seq_along(taxon)))	#Prepare df for transposing taxon column
sub_bacteria <- dcast(tmp, rain + SampleType + newid ~ taxon, na.rm = T,
	value.var = "Abundance")[,-3]	#Transforms long df to wide df, with OTUs as columns
sub_bacteria[is.na(sub_bacteria)] <- 0	#Replace NA with 0
colnames(sub_bacteria) <- gsub(x = colnames(sub_bacteria),
	pattern = "Otu", replacement = "ESV")	#Rename OTUs to ESVs
colnames(sub_bacteria) <- make.names(colnames(sub_bacteria), unique = TRUE)	#For ease of matching mvabund objects
str(sub_bacteria)	#Check structure
write.csv(sub_bacteria, file = "../data/sub_bacteria.csv", row.names = F)	#Save df
#sub_bacteria <- read.csv("../data/sub_bacteria.csv")	#Read in if saved


###FUNGAL DATA

fun_rarPhySeq = readRDS("../data/fungal_rarPhySeq")	#Read in rarefied count data
fun_count_df <- psmelt(fun_rarPhySeq)	#Make it an ordinary dataframe, not phyloseq
fun_count_df$taxon <- paste(fun_count_df$Class, fun_count_df$OTU)	#Create field combining OTU and Class
str(fun_count_df)	#Check structure
write.csv(fun_count_df, file = "../data/fun_count_df.csv", row.names = F)	#Save df
#fun_count_df <- read.csv("../data/fun_count_df.csv", header = T)	#Read in df if saved


##Subsample fungal dataframe
newdata1 <- subset(fun_count_df, Abundance >= (0.05 * max(bac_count_df$Abundance)),
	select = c(rain, SampleType, taxon, Abundance))	#Remove lowest 5% and subset
tmp1 <- ddply(newdata1, .(SampleType, taxon), transform,
	newid = paste(SampleType, seq_along(taxon)))	#Prepare df for transposing taxon column
sub_fungi <- dcast(tmp1, rain + SampleType + newid ~ taxon, na.rm = T,
	value.var = "Abundance")[,-3]	#Transforms long df to wide df, with OTUs as columns
sub_fungi[is.na(sub_fungi)] <- 0	#Replace NA with 0
colnames(sub_fungi) <- gsub(x = colnames(sub_fungi),
	pattern = "esv_", replacement = "ESV")	#Rename esvs to ESVs for consistency
colnames(sub_fungi) <- make.names(colnames(sub_fungi), unique = TRUE)	#For ease of matching mvabund objects
str(sub_fungi)	#Check structure
write.csv(sub_fungi, file = "../data/sub_fungi.csv", row.names = F)	#Save df
#sub_fungi <- read.csv("../data/sub_fungi.csv")	#Read in if saved



###BACTERIAL ANALYSIS


##Convert abundance data to mvabund object.
sub_bac <- mvabund(sub_bacteria[,3:length(sub_bacteria)])

##Bacterial ESV distributions explained by environmental gradient (rain)
bac_gradient <- manyglm(sub_bac ~ sub_bacteria$rain,
	family = "negative.binomial")	#Perform manyglm
pdf(file = "../figures/bac_gradient_residuals.pdf")	#Prepare to save image
plot(bac_gradient)	#Plot of residuals
dev.off()	#Deactivate pdf device
bac_gradient.anova <- anova(bac_gradient, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
bac_gradient.df <- as.data.frame(bac_gradient.anova$uni.p)	#Access p-values
write.csv(bac_gradient.df, file = "../data/bac_gradient.df.csv", row.names = F)	#Save p-values
bac_gradient_spp = colnames(bac_gradient.df)[bac_gradient.df[2,]<=0.05]	#Select significant p-values
write.csv(bac_gradient_spp, file = "../data/bac_gradient_spp.csv", row.names = F)	#Save list of significant ESVs
bac_gradient.pi <- pi0est(p = bac_gradient.anova$uni.p, 
	lambda = seq(0.1, 0.9, 0.1), pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_bac_gradient <- qvalue(p = bac_gradient.anova$uni.p,
	pi0 = bac_gradient.pi$pi0)	#Calculate q-values
q_bac_gradient.df <- as.data.frame(q_bac_gradient$qvalues)	#Access q-values
write.csv(q_bac_gradient.df, file = "../data/q_bac_gradient.df.csv", row.names = F)	#Save q-values
bac_gradient.q <- colnames(q_bac_gradient.df)[q_bac_gradient.df[2,] <= 0.05]	#Select significant q-values
sigbac_gradient <- subset(bac_gradient_spp, x %in% bac_gradient.q$x)
write.csv(sigbac_gradient, file = "../data/sigbac_gradient.csv", row.names = F)	#Save list of significant ESVs

##Represent strength of relationship between environment gradient
##and significantly associated bacterial ESVs
sigbac_gradient <- as.data.frame(sigbac_gradient)	#Change to df
sigbacgrad.df <- sub_bacteria[,c(1:2, na.omit(match(sigbac_gradient$x, colnames(sub_bacteria))))]	#Subset dataframe by significant ESVs
write.csv(sigbacgrad.df, file = "../data/sigbacgrad.df.csv", row.names = F)	#Save df
#sigbacgrad.df <- read.csv("../data/sigbacgrad.df.csv")	#Read in significant ESV df
Gradient <- sigbacgrad.df[,1]	#Rename variable for stats
sigbacgrad.abund <- mvabund(sigbacgrad.df[,3:length(sigbacgrad.df)])	#mvabund
sigbacgrad.glm <- manyglm(sigbacgrad.abund ~ Gradient,
	family = "negative.binomial")	#Assess GLMs for significant ESVs
pdf(file = "../figures/sigbacgrad_coefplot.pdf")	#Prepare to save plot
par(mar=c(5,8,3,1))	#Change margins of imaging device
coefplot.manyglm(sigbacgrad.glm, which.Xcoef=2, mfrow = NULL)	#Plot GLM coefficients
dev.off()	#Deactivate plotting device

##Assess those with significant coefficient CIs
bacgrad_glm <- c(grep("00020", names(sub_bacteria), value=TRUE))	#Select significant ESV
bacgrad_glm.df <- sub_bacteria[,c(1:2, na.omit(match(bacgrad_glm, colnames(sub_bacteria))))]	#Subset dataframe by significant ESVs
Gradient <- bacgrad_glm.df[,1]	#Rename variable for stats
sigbac <- bac_count_df %>%
	select(Family, Genus) %>% filter(bac_count_df$OTU == "Otu00020")	#Fetch finer taxonomic names
unique(sigbac)	#Call finer taxa
colnames(bacgrad_glm.df)[3] <- "Amnibacterium.ESV00020"	#Replace class with finer taxon
sigbacgrad <- glm.nb(Amnibacterium.ESV00020 ~ Gradient, data = bacgrad_glm.df)	#Examine only significant GLMs
pdf(file = "../figures/sigbacgrad.pdf")	#Prepare to save plot
plot(predictorEffect("Gradient", sigbacgrad), main="",
	axes = list(x = list(Gradient = list(lab = "Environmental Gradient (mm rain/year)")),
		y = list(type = "response", lab = "Amnibacterium ESV 20 (abundance)")),
	lines = list(col = "#481567FF"))	#Plot GLM
dev.off()	#Deactivate plotting device
bacgrad_hl <- hoslem.test(bacgrad_glm.df$Amnibacterium.ESV00020,
	fitted(sigbacgrad), g=10)	#Fetch Hosmer-Lemeshow goodness-of-fit
bacgrad_stats <- PseudoR2(sigbacgrad)	#Assess McFadden's pseudo-R-squared


##Bacterial ESV distributions explained by plant part, air, or soil
bac_air <- manyglm(sub_bac ~ sub_bacteria$SampleType=="Air",
	family = "negative.binomial")	#Perform manyglm
bac_air.anova <- anova(bac_air, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
bac_air.df <- as.data.frame(bac_air.anova$uni.p)	#Access p-values
write.csv(bac_air.df, file = "../data/bac_air.df.csv", row.names = F)	#Save p-values
bac_air_spp = colnames(bac_air.df)[bac_air.df[2,]<=0.05]	#Select significant p-values
write.csv(bac_air_spp, file = "../data/bac_air_spp.csv", row.names = F)	#Save list of significant ESVs
bac_air.pi <- pi0est(p = bac_air.anova$uni.p, 
	lambda = seq(0.1, 0.9, 0.1), pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_bac_air <- qvalue(p = bac_air.anova$uni.p, pi0 = bac_air.pi$pi0)	#Calculate q-values
q_bac_air.df <- as.data.frame(q_bac_air$qvalues)	#Access q-values
write.csv(q_bac_air.df, file = "../data/q_bac_air.df.csv", row.names = F)	#Save q-values
bac_air.q <- colnames(q_bac_air.df)[q_bac_air.df[2,] <= 0.05]	#Select significant q-values
sigbac_air <- subset(bac_air_spp, x %in% bac_air.q$x)
write.csv(sigbac_air, file = "../data/sigbac_air.csv", row.names = F)	#Save list of significant ESVs

bac_soil <- manyglm(sub_bac ~ sub_bacteria$SampleType=="Soil",
	family = "negative.binomial")	#Perform manyglm
bac_soil.anova <- anova(bac_soil, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
bac_soil.df <- as.data.frame(bac_soil.anova$uni.p)	#Access p-values
write.csv(bac_soil.df, file = "../data/bac_soil.df.csv", row.names = F)	#Save p-values
bac_soil_spp = colnames(bac_soil.df)[bac_soil.df[2,]<=0.05]	#Select significant p-values
write.csv(bac_soil_spp, file = "../data/bac_soil_spp.csv", row.names = F)	#Save list of significant ESVs
bac_soil.pi <- pi0est(p = bac_soil.anova$uni.p, 
	lambda = seq(0.1, 0.9, 0.1), pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_bac_soil <- qvalue(p = bac_soil.anova$uni.p, pi0 = bac_soil.pi$pi0)	#Calculate q-values
q_bac_soil.df <- as.data.frame(q_bac_soil$qvalues)	#Access q-values
write.csv(q_bac_soil.df, file = "../data/q_bac_soil.df.csv", row.names = F)	#Save q-values
bac_soil.q <- colnames(q_bac_soil.df)[q_bac_soil.df[2,] <= 0.05]	#Select significant q-values
sigbac_soil <- subset(bac_soil_spp, x %in% bac_soil.q$x)
write.csv(sigbac_soil, file = "../data/sigbac_soil.csv", row.names = F)	#Save list of significant ESVs

bac_root <- manyglm(sub_bac ~ sub_bacteria$SampleType=="Root",
	family = "negative.binomial")	#Perform manyglm
bac_root.anova <- anova(bac_root, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
bac_root.df <- as.data.frame(bac_root.anova$uni.p)	#Access p-values
write.csv(bac_root.df, file = "../data/bac_root.df.csv", row.names = F)	#Save p-values
bac_root_spp = colnames(bac_root.df)[bac_root.df[2,]<=0.05]	#Select significant p-values
write.csv(bac_root_spp, file = "../data/bac_root_spp.csv", row.names = F)	#Save list of significant ESVs
bac_root.pi <- pi0est(p = bac_root.anova$uni.p, 
	lambda = seq(0.1, 0.9, 0.1), pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_bac_root <- qvalue(p = bac_root.anova$uni.p, pi0 = bac_root.pi$pi0)	#Calculate q-values
q_bac_root.df <- as.data.frame(q_bac_root$qvalues)	#Access q-values
write.csv(q_bac_root.df, file = "../data/q_bac_root.df.csv", row.names = F)	#Save q-values
bac_root.q <- colnames(q_bac_root.df)[q_bac_root.df[2,] <= 0.05]	#Select significant q-values
sigbac_root <- subset(bac_root_spp, x %in% bac_root.q$x)
write.csv(sigbac_root, file = "../data/sigbac_root.csv", row.names = F)	#Save list of significant ESVs

bac_stem <- manyglm(sub_bac ~ sub_bacteria$SampleType=="Stem",
	family = "negative.binomial")	#Perform manyglm
bac_stem.anova <- anova(bac_stem, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
bac_stem.df <- as.data.frame(bac_stem.anova$uni.p)	#Access p-values
write.csv(bac_stem.df, file = "../data/bac_stem.df", row.names = F)	#Save p-values
bac_stem_spp = colnames(bac_stem.df)[bac_stem.df[2,]<=0.05]	#Select significant p-values
write.csv(bac_stem_spp, file = "../data/bac_stem_spp.csv", row.names = F)	#Save list of significant ESVs
bac_stem.pi <- pi0est(p = bac_stem.anova$uni.p,
	lambda = 0.5, pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_bac_stem <- qvalue(p = bac_stem.anova$uni.p, pi0 = bac_stem.pi$pi0)	#Calculate q-values
q_bac_stem.df <- as.data.frame(q_bac_stem$qvalues)	#Access q-values
write.csv(q_bac_stem.df, file = "../data/q_bac_stem.df.csv", row.names = F)	#Save q-values
bac_stem.q <- colnames(q_bac_stem.df)[q_bac_stem.df[2,] <= 0.05]	#Select significant q-values
sigbac_stem <- subset(bac_stem_spp, x %in% bac_stem.q$x)
write.csv(sigbac_stem, file = "../data/sigbac_stem.csv", row.names = F)	#Save list of significant ESVs

bac_axil <- manyglm(sub_bac ~ sub_bacteria$SampleType=="Axil",
	family = "negative.binomial")	#Perform manyglm
bac_axil.anova <- anova(bac_axil, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
bac_axil.df <- as.data.frame(bac_axil.anova$uni.p)	#Access p-values
write.csv(bac_axil.df, file = "../data/bac_axil.df", row.names = F)	#Save p-values
bac_axil_spp = colnames(bac_axil.df)[bac_axil.df[2,]<=0.05]	#Select significant p-values
write.csv(bac_axil_spp, file = "../data/bac_axil_spp", row.names = F)	#Save list of significant ESVs
bac_axil.pi <- pi0est(p = bac_axil.anova$uni.p,
	lambda = 0.5, pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_bac_axil <- qvalue(p = bac_axil.anova$uni.p, pi0 = bac_axil.pi$pi0)	#Calculate q-values
q_bac_axil.df <- as.data.frame(q_bac_axil$qvalues)	#Access q-values
write.csv(q_bac_axil.df, file = "../data/q_bac_axil.df.csv", row.names = F)	#Save q-values
bac_axil.q <- colnames(q_bac_axil.df)[q_bac_axil.df[2,] <= 0.05]	#Select significant q-values
sigbac_axil <- subset(bac_axil_spp, x %in% bac_axil.q$x)
write.csv(sigbac_axil, file = "../data/sigbac_axil.csv", row.names = F)	#Save list of significant ESVs

bac_petiole <- manyglm(sub_bac ~ sub_bacteria$SampleType=="Petiole",
	family = "negative.binomial")	#Perform manyglm
bac_petiole.anova <- anova(bac_petiole, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
bac_petiole.df <- as.data.frame(bac_petiole.anova$uni.p)	#Access p-values
write.csv(bac_petiole.df, file = "../data/bac_petiole.df.csv", row.names = F)	#Save p-values
bac_petiole_spp = colnames(bac_petiole.df)[bac_petiole.df[2,]<=0.05]	#Select significant p-values
write.csv(bac_petiole_spp, file = "../data/bac_petiole_spp.csv", row.names = F)	#Save list of significant ESVs
bac_petiole.pi <- pi0est(p = bac_petiole.anova$uni.p, 
	lambda = seq(0.1, 0.9, 0.1), pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_bac_petiole <- qvalue(p = bac_petiole.anova$uni.p, pi0 = bac_petiole.pi$pi0)	#Calculate q-values
q_bac_petiole.df <- as.data.frame(q_bac_petiole$qvalues)	#Access q-values
write.csv(q_bac_petiole.df, file = "../data/q_bac_petiole.df.csv", row.names = F)	#Save q-values
bac_petiole.q <- colnames(q_bac_petiole.df)[q_bac_petiole.df[2,] <= 0.05]	#Select significant q-values
sigbac_petiole <- subset(bac_petiole_spp, x %in% bac_petiole.q$x)
write.csv(sigbac_petiole, file = "../data/sigbac_petiole.csv", row.names = F)	#Save list of significant ESVs

bac_leaf <- manyglm(sub_bac ~ sub_bacteria$SampleType=="Leaf",
	family = "negative.binomial")	#Perform manyglm
bac_leaf.anova <- anova(bac_leaf, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
bac_leaf.df <- as.data.frame(bac_leaf.anova$uni.p)	#Access p-values
write.csv(bac_leaf.df, file = "../data/bac_leaf.df.csv", row.names = F)	#Save p-values
bac_leaf_spp = colnames(bac_leaf.df)[bac_leaf.df[2,]<=0.05]	#Select significant p-values
write.csv(bac_leaf_spp, file = "../data/bac_leaf_spp.csv", row.names = F)	#Save list of significant ESVs
bac_leaf.pi <- pi0est(p = bac_leaf.anova$uni.p, lambda = 0.5,
	pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_bac_leaf <- qvalue(p = bac_leaf.anova$uni.p, pi0 = bac_leaf.pi$pi0)	#Calculate q-values
q_bac_leaf.df <- as.data.frame(q_bac_leaf$qvalues)	#Access q-values
write.csv(q_bac_leaf.df, file = "../data/q_bac_leaf.df.csv", row.names = F)	#Save q-values
bac_leaf.q <- colnames(q_bac_leaf.df)[q_bac_leaf.df[2,] <= 0.05]	#Select significant q-values
sigbac_leaf <- subset(bac_leaf_spp, x %in% bac_leaf.q$x)
write.csv(sigbac_leaf, file = "../data/sigbac_leaf.csv", row.names = F)	#Save list of significant ESVs

bac_litter <- manyglm(sub_bac ~ sub_bacteria$SampleType=="Litter",
	family = "negative.binomial")	#Perform manyglm
bac_litter.anova <- anova(bac_litter, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
bac_litter.df <- as.data.frame(bac_litter.anova$uni.p)	#Access p-values
write.csv(bac_litter.df, file = "../data/bac_litter.df.csv", row.names = F)	#save p-values
bac_litter_spp = colnames(bac_litter.df)[bac_litter.df[2,]<=0.05]	#Select significant p-values
write.csv(bac_litter_spp, file = "../data/bac_litter_spp.csv", row.names = F)	#Save list of significant ESVs
bac_litter.pi <- pi0est(p = bac_litter.anova$uni.p, lambda = seq(0.1, 0.9, 0.1),
	pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_bac_litter <- qvalue(p = bac_litter.anova$uni.p, pi0 = bac_litter.pi$pi0)	#Calculate q-values
q_bac_litter.df <- as.data.frame(q_bac_litter$qvalues)	#Access q-values
write.csv(q_bac_litter.df, file = "../data/q_bac_litter.df.csv", row.names = F)	#Save q-values
bac_litter.q <- colnames(q_bac_litter.df)[q_bac_litter.df[2,] <= 0.05]	#Select significant q-values
sigbac_litter <- subset(bac_litter_spp, x %in% bac_litter.q$x)
write.csv(sigbac_litter, file = "../data/sigbac_litter.csv", row.names = F)	#Save list of significant ESVs


pdf(file = "../figures/bac_SampleType_residuals.pdf")	#Prepare pdf file
par(mfrow = c(2,4), bty = "l", pty = "m")	#Set plotting parameters
plot(bac_air)	#Plot of air residuals
plot(bac_leaf)	#Plot of leaf residuals
plot(bac_petiole)	#Plot of petiole residuals
plot(bac_axil)	#Plot of axil residuals
plot(bac_stem)	#Plot of stem residuals
plot(bac_root)	#Plot of root residuals
plot(bac_litter)	#Plot of litter residuals
plot(bac_soil)	#Plot of soil residuals
dev.off()	#Deactivate pdf device


##Graphically represent number of bacterial ESVs significantly associated with
##each sample type
bacair_class <- as.data.frame(str_split_fixed(sigbac_air$x, ".ESV", 2))	#Split ESV class and number
bacair_class <- as.data.frame(bacair_class %>%
	group_by(V1) %>%
	summarise(n()))	#Summarise ESVs by class
colnames(bacair_class)[2] <- "Air"	#Rename column

bacleaf_class <- as.data.frame(str_split_fixed(sigbac_leaf$x, ".ESV", 2))	#Split ESV class and number
bacleaf_class <- as.data.frame(bacleaf_class %>%
	group_by(V1) %>%
	summarise(n()))	#Summarise ESVs by class
colnames(bacleaf_class)[2] <- "Leaf"	#Rename column

bacpetiole_class <- as.data.frame(str_split_fixed(sigbac_petiole$x, ".ESV", 2))	#Split ESV class and number
bacpetiole_class <- as.data.frame(bacpetiole_class %>%
	group_by(V1) %>%
	summarise(n()))	#Summarise ESVs by class
colnames(bacpetiole_class)[2] <- "Petiole"	#Rename column

bacaxil_class <- as.data.frame(str_split_fixed(sigbac_axil$x, ".ESV", 2))	#Split ESV class and number
bacaxil_class <- as.data.frame(bacaxil_class %>%
	group_by(V1) %>%
	summarise(n()))	#Summarise ESVs by class
colnames(bacaxil_class)[2] <- "Axil"	#Rename column

bacstem_class <- as.data.frame(str_split_fixed(sigbac_stem$x, ".ESV", 2))	#Split ESV class and number
bacstem_class <- as.data.frame(bacstem_class %>%
	group_by(V1) %>%
	summarise(n()))	#Summarise ESVs by class
colnames(bacstem_class)[2] <- "Stem"	#Rename column

bacroot_class <- as.data.frame(str_split_fixed(sigbac_root$x, ".ESV", 2))	#Split ESV class and number
bacroot_class <- as.data.frame(bacroot_class %>%
	group_by(V1) %>%
	summarise(n()))	#Summarise ESVs by class
colnames(bacroot_class)[2] <- "Root"	#Rename column

baclitter_class <- as.data.frame(str_split_fixed(sigbac_litter$x, ".ESV", 2))	#Split ESV class and number
baclitter_class <- as.data.frame(baclitter_class %>%
	group_by(V1) %>%
	summarise(n()))	#Summarise ESVs by class
colnames(baclitter_class)[2] <- "Litter"	#Rename column

bacsoil_class <- as.data.frame(str_split_fixed(sigbac_soil$x, ".ESV", 2))	#Split ESV class and number
bacsoil_class <- as.data.frame(bacsoil_class %>%
	group_by(V1) %>%
	summarise(n()))	#Summarise ESVs by class
colnames(bacsoil_class)[2] <- "Soil"	#Rename column

##Merge significant ESVs by class and sample type
bacsamp <- merge(merge(merge(merge(merge(merge(merge(
	bacair_class, bacleaf_class, by = "V1", all = T),
	bacpetiole_class, by = "V1", all = T), bacaxil_class, by = "V1", all = T),
	bacstem_class, by = "V1", all = T), bacroot_class, by = "V1", all = T),
	baclitter_class, by = "V1", all = T), bacsoil_class, by = "V1", all = T)
bacsamp[is.na(bacsamp)] <- 0	#Replace NA with 0
bacsamp <- bacsamp[-4,]	#Remove unnecessary NA row
rownames(bacsamp) <- bacsamp[,1]	#Make the classes row names
bacsamp <- bacsamp[,-1]	#Remove redundant column
write.csv(bacsamp, file = "../data/bacsamp.csv", row.names = T)	#Save df of significant ESVs

##Bar graph of ESV numbers per class for each sample type
pdf(file = "../figures/bacsamp.pdf")	#Prepare pdf file
par(mfrow=c(1, 1), mar = c(5, 5, 4, 9))	#Set margins
barplot(as.matrix(bacsamp), col = c("#440154FF","#404788FF", "#287D8EFF"),
	border = "white", space = 0.04, ylim = c(0,35), font.axis = 2,
	legend.text = rownames(bacsamp), args.legend =
	list(x = "topright", bty = "n", inset = c(-0.4, 0), xpd = TRUE),
	xlab = "Sample Type", ylab = "Number of Associated Bacterial ESVs")
dev.off()	#Deactivate pdf device


##Assess interactive effects of gradient and sample type on bacterial ESVs.
bac_rainsamp <- manyglm(sub_bac ~ sub_bacteria$rain * sub_bacteria$SampleType,
	family="negative.binomial")	#Perform manyglm
bac_rainsamp.anova = anova(bac_rainsamp, nBoot=1000, test="LR",
	p.uni="unadjusted", resamp="montecarlo")	#Perform ANOVA on manyglm
bac_rainsamp.anova.df <- as.data.frame(bac_rainsamp.anova$uni.p)	#Access p-values
write.csv(bac_rainsamp.anova.df, file = "../data/bac_rainsamp.summary.df.csv", row.names = F)
bac_rainsamp.p = bac_rainsamp.anova$uni.p	#Access p-values
bac_rainsamp.p.df <- as.data.frame(bac_rainsamp.p)	#Change to df
write.csv(bac_rainsamp.p.df, file = "../data/bac_rainsamp.p.df.csv")	#Save p-values
sigbacrainsamp_spp = colnames(bac_rainsamp.p.df)[bac_rainsamp.p.df[4,]<=0.05]	#Select significant p-values
write.csv(sigbacrainsamp_spp, file = "../data/sigbacrainsamp_spp.csv")	#Save list of significant ESVs
sigbacrainsamp.pi <- pi0est(p = bac_rainsamp.anova$uni.p, 
	lambda = seq(0.1, 0.9, 0.1), pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_sigbacrainsamp <- qvalue(p = bac_rainsamp.anova$uni.p,
	pi0 = sigbacrainsamp.pi$pi0)	#Calculate q-values
q_sigbacrainsamp.df <- as.data.frame(q_sigbacrainsamp$qvalues)	#Access q-values
write.csv(q_sigbacrainsamp.df, file = "../data/q_sigbacrainsamp.df.csv", row.names = F)	#Save q-values
sigbacrainsamp.q <- colnames(q_sigbacrainsamp.df)[q_sigbacrainsamp.df[4,] <= 0.05]	#Select significant q-values
sigbacrainsamp <- subset(sigbacrainsamp_spp, x %in% sigbacrainsamp.q$x)
write.csv(sigbacrainsamp, file = "../data/sigbacrainsamp.csv", row.names = F)	#Save list of significant ESVs

##Represent strength of interactive effects between explanatory
##variables and significantly associated bacterial ESVs.
sigbacrainsamp.df <- sub_bacteria[,c(1:2, na.omit(match(sigbacrainsamp$x,
	colnames(sub_bacteria))))]	#Subset dataframe by significant ESVs
write.csv(sigbacrainsamp.df, file = "../data/sigbacrainsamp.df.csv", row.names = F)	#Save df
sigbacrainsamp.abund <- mvabund(sigbacrainsamp.df[,3:length(sigbacrainsamp.df)]) #mvabund
Gradient <- sigbacrainsamp.df[,1]	#Rename variable for x-axis
Sample_type <- sigbacrainsamp.df$SampleType	#Rename variable for x-axis
sigbacrainsamp.glm <- manyglm(sigbacrainsamp.abund ~ Gradient * Sample_type,
	family = "negative.binomial")	#Assess GLMs for significant ESVs
pdf(file = "../figures/sigbacrainsamp_coefplot1.pdf")	#Prepare to save plot
par(mar=c(5,8,3,1))	#Change margins of imaging device
coefplot.manyglm(sigbacrainsamp.glm, which.Xcoef=10, mfrow = NULL)	#Plot GLM coefficients
dev.off()	#Deactivate plotting device



###FUNGAL ANALYSIS


##Convert abundance data to mvabund object.
sub_fun <- mvabund(sub_fungi[,3:length(sub_fungi)])

##Fungal ESV distributions explained by environmental gradient (rain)
fun_gradient <- manyglm(sub_fun ~ sub_fungi$rain,
	family = "negative.binomial")	#Perform manyglm
pdf(../figures/file = "fun_gradient_residuals.pdf")	#Prepare to save image
plot(fun_gradient)	#Plot of residuals
dev.off()	#Deactivate pdf device
fun_gradient.anova <- anova(fun_gradient, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
fun_gradient.df <- as.data.frame(fun_gradient.anova$uni.p)	#Access p-values
write.csv(fun_gradient.df, file = "../data/fun_gradient.df.csv", row.names = F)	#Save p-values
fun_gradient_spp = colnames(fun_gradient.df)[fun_gradient.df[2,]<=0.05]	#Select significant p-values
write.csv(fun_gradient_spp, file = "../data/fun_gradient_spp.csv", row.names = F)	#Save list of significant ESVs
fun_gradient.pi <- pi0est(p = fun_gradient.anova$uni.p, 
	lambda = seq(0.1, 0.9, 0.1), pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_fun_gradient <- qvalue(p = fun_gradient.anova$uni.p,
	pi0 = fun_gradient.pi$pi0)	#Calculate q-values
q_fun_gradient.df <- as.data.frame(q_fun_gradient$qvalues)	#Access q-values
write.csv(q_fun_gradient.df, file = "../data/q_fun_gradient.df.csv", row.names = F)	#Save q-values
fun_gradient.q <- colnames(q_fun_gradient.df)[q_fun_gradient.df[2,] <= 0.05]	#Select significant q-values
sigfun_gradient <- subset(fun_gradient_spp, x %in% fun_gradient.q$x)
write.csv(sigfun_gradient, file = "../data/sigfun_gradient.csv", row.names = F)	#Save list of significant ESVs


##Represent strength of relationship between environment gradient
##and significantly associated bacterial ESVs.
sigfun_gradient <- as.data.frame(sigfun_gradient)	#Change to df
sigfungrad.df <- sub_fungi[,c(1:2, na.omit(match(sigfun_gradient$x, colnames(sub_fungi))))]	#Subset dataframe by significant ESVs
write.csv(sigfungrad.df, file = "../data/sigfungrad.df.csv", row.names = F)	#Save df
sigfungrad.abund <- mvabund(sigfungrad.df[,3:length(sigfungrad.df)])	#mvabund
Gradient <- sigfungrad.df[,1]	#Rename variable for coefficient assessment
sigfungrad.glm <- manyglm(sigfungrad.abund ~ Gradient,
	family = "negative.binomial")	#Assess GLMs for significant ESVs
pdf(file = "../figures/sigfungrad_coefplot.pdf")	#Prepare to save plot
par(mar=c(5,7,3,1))	#Change margins of imaging device
coefplot.manyglm(sigfungrad.glm, which.Xcoef = 2, mfrow = NULL)	#Plot GLM coefficients
dev.off()	#Deactivate plotting device

##Assess those with significant coefficient CIs
fungrad_glm <- c("Dothideomycetes.ESV5")	#Select significant ESVs
fungrad_glm.df <- sub_fungi[,c(1:2, na.omit(match(fungrad_glm, colnames(sub_fungi))))]	#Subset dataframe by significant ESVs
Gradient <- fungrad_glm.df[,1]	#Rename variable for stats
sigfun <- fun_count_df %>%
	select(Family, Genus) %>% filter(fun_count_df$OTU == "esv_5")	#Fetch finer taxonomic names
unique(sigfun)	#Call finer taxa
colnames(fungrad_glm.df)[3] <- c("Mycosphaerellaceae.ESV5")	#Replace classes with finer taxa
sigfungrad <- glm.nb(Mycosphaerellaceae.ESV5 ~ Gradient, data = fungrad_glm.df)	#Examine only significant GLMs
pdf(file = "../figures/sigfungrad.pdf")	#Prepare to save plot
plot(predictorEffect("Gradient", sigfungrad), main="",
	axes = list(x = list(Gradient = list(lab = "Environmental Gradient (mm rain/year)")),
		y = list(type = "response", lab = "Mycosphaerellaceae ESV 5 (abundance)")),
	lines = list(col = "#95D840FF"))	#Plot GLM
dev.off()	#Deactivate plotting device
fungrad_hl <- hoslem.test(fungrad_glm.df$Mycosphaerellaceae.ESV5,
	fitted(sigfungrad), g=10)	#Fetch Hosmer-Lemeshow goodness-of-fit
fungrad_stats <- PseudoR2(sigfungrad)	#Assess McFadden's pseudo-R-squared



##Fungal ESV distributions explained by plant part, air, or soil
fun_air <- manyglm(sub_fun ~ sub_fungi$SampleType=="Air",
	family = "negative.binomial")	#Perform manyglm
fun_air.anova <- anova(fun_air, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
fun_air.df <- as.data.frame(fun_air.anova$uni.p)	#Access p-values
write.csv(fun_air.df, file = "../data/fun_air.df.csv", row.names = F)	#Save p-values
fun_air_spp = colnames(fun_air.df)[fun_air.df[2,]<=0.05]	#Select significant p-values
write.csv(fun_air_spp, file = "../data/fun_air_spp.csv", row.names = F)	#Save list of significant ESVs
fun_air.pi <- pi0est(p = fun_air.anova$uni.p, 
	lambda = seq(0.1, 0.9, 0.1), pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_fun_air <- qvalue(p = fun_air.anova$uni.p, pi0 = fun_air.pi$pi0)	#Calculate q-values
q_fun_air.df <- as.data.frame(q_fun_air$qvalues)	#Access q-values
write.csv(q_fun_air.df, file = "../data/q_fun_air.df.csv", row.names = F)	#Save q-values
fun_air.q <- colnames(q_fun_air.df)[q_fun_air.df[2,] <= 0.05]	#Select significant q-values
sigfun_air <- subset(fun_air_spp, x %in% fun_air.q$x)
write.csv(sigfun_air, file = "../data/sigfun_air.csv", row.names = F)	#Save list of significant ESVs

fun_soil <- manyglm(sub_fun ~ sub_fungi$SampleType=="Soil",
	family = "negative.binomial")	#Perform manyglm
fun_soil.anova <- anova(fun_soil, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
fun_soil.df <- as.data.frame(fun_soil.anova$uni.p)	#Access p-values
write.csv(fun_soil.df, file = "../data/fun_soil.df.csv", row.names = F)	#Save p-values
fun_soil_spp = colnames(fun_soil.df)[fun_soil.df[2,]<=0.05]	#Select significant p-values
write.csv(fun_soil_spp, file = "../data/fun_soil_spp.csv", row.names = F)	#Save list of significant ESVs
fun_soil.pi <- pi0est(p = fun_soil.anova$uni.p, 
	lambda = seq(0.1, 0.9, 0.1), pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_fun_soil <- qvalue(p = fun_soil.anova$uni.p, pi0 = fun_soil.pi$pi0)	#Calculate q-values
q_fun_soil.df <- as.data.frame(q_fun_soil$qvalues)	#Access q-values
write.csv(q_fun_soil.df, file = "../data/q_fun_soil.df.csv", row.names = F)	#Save q-values
fun_soil.q <- colnames(q_fun_soil.df)[q_fun_soil.df[2,] <= 0.05]	#Select significant q-values
sigfun_soil <- subset(fun_soil_spp, x %in% fun_soil.q$x)
write.csv(sigfun_soil, file = "../data/sigfun_soil.csv", row.names = F)	#Save list of significant ESVs

fun_root <- manyglm(sub_fun ~ sub_fungi$SampleType=="Root",
	family = "negative.binomial")	#Perform manyglm
fun_root.anova <- anova(fun_root, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
fun_root.df <- as.data.frame(fun_root.anova$uni.p)	#Access p-values
write.csv(fun_root.df, file = "../data/fun_root.df.csv", row.names = F)	#Save p-values
fun_root_spp = colnames(fun_root.df)[fun_root.df[2,]<=0.05]	#Select significant p-values
write.csv(fun_root_spp, file = "../data/fun_root_spp.csv", row.names = F)	#Save list of significant ESVs
fun_root.pi <- pi0est(p = fun_root.anova$uni.p, 
	lambda = seq(0.1, 0.9, 0.1), pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_fun_root <- qvalue(p = fun_root.anova$uni.p, pi0 = fun_root.pi$pi0)	#Calculate q-values
q_fun_root.df <- as.data.frame(q_fun_root$qvalues) #Access q-values
write.csv(q_fun_root.df, file = "../data/q_fun_root.df.csv", row.names = F)	#Save q-values
fun_root.q <- colnames(q_fun_root.df)[q_fun_root.df[2,] <= 0.05]	#Select significant q-values
sigfun_root <- subset(fun_root_spp, x %in% fun_root.q$x)
write.csv(sigfun_root, file = "../data/sigfun_root.csv", row.names = F)	#Save list of significant ESVs

fun_stem <- manyglm(sub_fun ~ sub_fungi$SampleType=="Stem",
	family = "negative.binomial")	#Perform manyglm
fun_stem.anova <- anova(fun_stem, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
fun_stem.df <- as.data.frame(fun_stem.anova$uni.p)	#Access p-values
write.csv(fun_stem.df, file = "../data/fun_stem.df.csv", row.names = F)	#Save p-values
fun_stem_spp = colnames(fun_stem.df)[fun_stem.df[2,]<=0.05]	#Select significant p-values
write.csv(fun_stem_spp, file = "../data/fun_stem_spp.csv", row.names = F)	#Save list of significant ESVs
fun_stem.pi <- pi0est(p = fun_stem.anova$uni.p, lambda = seq(0.1, 0.9, 0.1),
	pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_fun_stem <- qvalue(p = fun_stem.anova$uni.p, pi0 = fun_stem.pi$pi0)	#Calculate q-values
q_fun_stem.df <- as.data.frame(q_fun_stem$qvalues)	#Access q-values
write.csv(q_fun_stem.df, file = "../data/q_fun_stem.df.csv", row.names = F)	#Save q-values
fun_stem.q <- colnames(q_fun_stem.df)[q_fun_stem.df[2,] <= 0.05]	#Select significant q-values
sigfun_stem <- subset(fun_stem_spp, x %in% fun_stem.q$x)
write.csv(sigfun_stem, file = "../data/sigfun_stem.csv", row.names = F)	#Save list of significant ESVs

fun_axil <- manyglm(sub_fun ~ sub_fungi$SampleType=="Axil",
	family = "negative.binomial")	#Perform manyglm
fun_axil.anova <- anova(fun_axil, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
fun_axil.df <- as.data.frame(fun_axil.anova$uni.p)	#Access p-values
write.csv(fun_axil.df, file = "../data/fun_axil.df.csv", row.names = F)	#Save p-values
fun_axil_spp = colnames(fun_axil.df)[fun_axil.df[2,]<=0.05]	#Select significant p-values
write.csv(fun_axil_spp, file = "../data/fun_axil_spp.csv", row.names = F)	#Save list of significant ESVs
fun_axil.pi <- pi0est(p = fun_axil.anova$uni.p, lambda = seq(0.1, 0.9, 0.1),
	pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_fun_axil <- qvalue(p = fun_axil.anova$uni.p, pi0 = fun_axil.pi$pi0)	#Calculate q-values
q_fun_axil.df <- as.data.frame(q_fun_axil$qvalues)	#Access q-values
write.csv(q_fun_axil.df, file = "../data/q_fun_axil.df.csv", row.names = F)	#Save q-values
fun_axil.q <- colnames(q_fun_axil.df)[q_fun_axil.df[2,] <= 0.05]	#Select significant q-values
sigfun_axil <- subset(fun_axil_spp, x %in% fun_axil.q$x)
write.csv(sigfun_axil, file = "../data/sigfun_axil.csv", row.names = F)	#Save list of significant ESVs

fun_petiole <- manyglm(sub_fun ~ sub_fungi$SampleType=="Petiole",
	family = "negative.binomial")	#Perform manyglm
fun_petiole.anova <- anova(fun_petiole, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
fun_petiole.df <- as.data.frame(fun_petiole.anova$uni.p)	#Access p-values
write.csv(fun_petiole.df, file = "../data/fun_petiole.df.csv", row.names = F)	#Save p-values
fun_petiole_spp = colnames(fun_petiole.df)[fun_petiole.df[2,]<=0.05]	#Select significant p-values
write.csv(fun_petiole_spp, file = "../data/fun_petiole_spp.csv", row.names = F)	#Save list of significant ESVs
fun_petiole.pi <- pi0est(p = fun_petiole.anova$uni.p, 
	lambda = seq(0.1, 0.9, 0.1), pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_fun_petiole <- qvalue(p = fun_petiole.anova$uni.p, pi0 = fun_petiole.pi$pi0)	#Calculate q-values
q_fun_petiole.df <- as.data.frame(q_fun_petiole$qvalues)	#Access q-values
write.csv(q_fun_petiole.df, file = "../data/q_fun_petiole.df.csv", row.names = F)	#Save q-values
fun_petiole.q <- colnames(q_fun_petiole.df)[q_fun_petiole.df[2,] <= 0.05]	#Select significant q-values
sigfun_petiole <- subset(fun_petiole_spp, x %in% fun_petiole.q$x)
write.csv(sigfun_petiole, file = "../data/sigfun_petiole.csv", row.names = F)	#Save list of significant ESVs

fun_leaf <- manyglm(sub_fun ~ sub_fungi$SampleType=="Leaf",
	family = "negative.binomial")	#Perform manyglm
fun_leaf.anova <- anova(fun_leaf, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
fun_leaf.df <- as.data.frame(fun_leaf.anova$uni.p)	#Access p-values
write.csv(fun_leaf.df, file = "../data/fun_leaf.df.csv", row.names = F)	#Save p-values
fun_leaf_spp = colnames(fun_leaf.df)[fun_leaf.df[2,]<=0.05]	#Select significant p-values
write.csv(fun_leaf_spp, file = "../data/fun_leaf_spp.csv", row.names = F)	#Save list of significant ESVs
fun_leaf.pi <- pi0est(p = fun_leaf.anova$uni.p, lambda = seq(0.1, 0.9, 0.1),
	pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_fun_leaf <- qvalue(p = fun_leaf.anova$uni.p, pi0 = fun_leaf.pi$pi0)	#Calculate q-values
q_fun_leaf.df <- as.data.frame(q_fun_leaf$qvalues)	#Access q-values
write.csv(q_fun_leaf.df, file = "../data/q_fun_leaf.df.csv", row.names = F)	#Save q-values
fun_leaf.q <- colnames(q_fun_leaf.df)[q_fun_leaf.df[2,] <= 0.05]	#Select significant q-values
sigfun_leaf <- subset(fun_leaf_spp, x %in% fun_leaf.q$x)
write.csv(sigfun_leaf, file = "../data/sigfun_leaf.csv", row.names = F)	#Save list of significant ESVs

fun_litter <- manyglm(sub_fun ~ sub_fungi$SampleType=="Litter",
	family = "negative.binomial")	#Perform manyglm
fun_litter.anova <- anova(fun_litter, nBoot = 1000, p.uni = "unadjusted",
	test = "LR", resamp = "montecarlo")	#Perform ANOVA on manyglm
fun_litter.df <- as.data.frame(fun_litter.anova$uni.p)	#Access p-values
write.csv(fun_litter.df, file = "../data/fun_litter.df.csv", row.names = F)	#save p-values
fun_litter_spp = colnames(fun_litter.df)[fun_litter.df[2,]<=0.05]	#Select significant p-values
write.csv(fun_litter_spp, file = "../data/fun_litter_spp.csv", row.names = F)	#Save list of significant ESVs
fun_litter.pi <- pi0est(p = fun_litter.anova$uni.p, lambda = seq(0.1, 0.9, 0.1),
	pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_fun_litter <- qvalue(p = fun_litter.anova$uni.p, pi0 = fun_litter.pi$pi0)	#Calculate q-values
q_fun_litter.df <- as.data.frame(q_fun_litter$qvalues)	#Access q-values
write.csv(q_fun_litter.df, file = "../data/q_fun_litter.df.csv", row.names = F)	#Save q-values
fun_litter.q <- colnames(q_fun_litter.df)[q_fun_litter.df[2,] <= 0.05]	#Select significant q-values
sigfun_litter <- subset(fun_litter_spp, x %in% fun_litter.q$x)
write.csv(sigfun_litter, file = "../data/sigfun_litter.csv", row.names = F)	#Save list of significant ESVs


pdf(file = "../figures/fun_SampleType_residuals.pdf")	#Prepare pdf file
par(mfrow = c(2,4), bty = "l", pty = "m")	#Set plotting parameters
plot(fun_air)	#Plot of air residuals
plot(fun_leaf)	#Plot of leaf residuals
plot(fun_petiole)	#Plot of petiole residuals
plot(fun_axil)	#Plot of axil residuals
plot(fun_stem)	#Plot of stem residuals
plot(fun_root)	#Plot of root residuals
plot(fun_litter)	#Plot of litter residuals
plot(fun_soil)	#Plot of soil residuals
dev.off()	#Deactivate pdf device


##Graphically represent number of bacterial ESVs significantly associated with
##each sample type
funair_class <- as.data.frame(str_split_fixed(sigfun_air$x, ".ESV", 2))	#Split ESV class and number
funair_class <- as.data.frame(funair_class %>%
	group_by(V1) %>%
	summarise(n()))	#Summarise ESVs by class
colnames(funair_class)[2] <- "Air"	#Rename column

funleaf_class <- as.data.frame(str_split_fixed(sigfun_leaf$x, ".ESV", 2))	#Split ESV class and number
funleaf_class <- as.data.frame(funleaf_class %>%
	group_by(V1) %>%
	summarise(n()))	#Summarise ESVs by class
colnames(funleaf_class)[2] <- "Leaf"	#Rename column

funpetiole_class <- as.data.frame(str_split_fixed(sigfun_petiole$x, ".ESV", 2))	#Split ESV class and number
funpetiole_class <- as.data.frame(funpetiole_class %>%
	group_by(V1) %>%
	summarise(n()))	#Summarise ESVs by class
colnames(funpetiole_class)[2] <- "Petiole"	#Rename column

funaxil_class <- as.data.frame(str_split_fixed(sigfun_axil$x, ".ESV", 2))	#Split ESV class and number
funaxil_class <- as.data.frame(funaxil_class %>%
	group_by(V1) %>%
	summarise(n()))	#Summarise ESVs by class
colnames(funaxil_class)[2] <- "Axil"	#Rename column

funstem_class <- as.data.frame(str_split_fixed(sigfun_stem$x, ".ESV", 2))	#Split ESV class and number
funstem_class <- as.data.frame(funstem_class %>%
	group_by(V1) %>%
	summarise(n()))	#Summarise ESVs by class
colnames(funstem_class)[2] <- "Stem"	#Rename column

funroot_class <- as.data.frame(str_split_fixed(sigfun_root$x, ".ESV", 2))	#Split ESV class and number
funroot_class <- as.data.frame(funroot_class %>%
	group_by(V1) %>%
	summarise(n()))	#Summarise ESVs by class
colnames(funroot_class)[2] <- "Root"	#Rename column

funlitter_class <- as.data.frame(str_split_fixed(sigfun_litter$x, ".ESV", 2))	#Split ESV class and number
funlitter_class <- as.data.frame(funlitter_class %>%
	group_by(V1) %>%
	summarise(n()))	#Summarise ESVs by class
colnames(funlitter_class)[2] <- "Litter"	#Rename column

funsoil_class <- as.data.frame(str_split_fixed(sigfun_soil$x, ".ESV", 2))	#Split ESV class and number
funsoil_class <- as.data.frame(funsoil_class %>%
	group_by(V1) %>%
	summarise(n()))	#Summarise ESVs by class
colnames(funsoil_class)[2] <- "Soil"	#Rename column

##Merge significant ESVs by class and sample type
funsamp <- merge(merge(merge(merge(merge(merge(merge(
	funair_class, funleaf_class, by = "V1", all = T),
	funpetiole_class, by = "V1", all = T), funaxil_class, by = "V1", all = T),
	funstem_class, by = "V1", all = T), funroot_class, by = "V1", all = T),
	funlitter_class, by = "V1", all = T), funsoil_class, by = "V1", all = T)
funsamp[is.na(funsamp)] <- 0	#Replace NA with 0
funsamp <- funsamp[-9,]	#Remove unnecessary NA row
rownames(funsamp) <- funsamp[,1]	#Make the classes row names
funsamp <- funsamp[,-1]	#Remove redundant column
write.csv(funsamp, file = "../data/funsamp.csv", row.names = T)	#Save df of significant ESVs

##Bar graph of ESV numbers per class for each sample type
pdf(file = "../figures/funsamp.pdf")	#Prepare pdf file
par(mfrow=c(1, 1), mar = c(5, 5, 4, 9))	#Set margins
barplot(as.matrix(funsamp), col = c("#29AF7FFF","#3CBB75FF", "#55C667FF",
	"#73D055FF", "#95D840FF", "#B8DE29FF", "#DCE319FF", "#FDE725FF"),
	border = "white", space = 0.04, ylim = c(0,35), font.axis = 2,
	legend.text = rownames(funsamp), args.legend =
	list(x = "topright", bty = "n", inset = c(-0.4, -0.1), xpd = TRUE),
	xlab = "Sample Type", ylab = "Number of Associated Fungal ESVs")
dev.off()	#Deactivate pdf device


##Assess interactive effects of gradient and sample type on fungal ESVs.
fun_rainsamp <- manyglm(sub_fun ~ sub_fungi$rain * sub_fungi$SampleType,
	family="negative.binomial")	#Perform manyglm
fun_rainsamp.anova = anova(fun_rainsamp, nBoot=1000, test="LR",
	p.uni="unadjusted", resamp="montecarlo")	#Perform ANOVA on manyglm
fun_rainsamp.anova.df <- as.data.frame(fun_rainsamp.anova$uni.p)	#Access p-values
write.csv(fun_rainsamp.anova.df, file = "../data/fun_rainsamp.summary.df.csv", row.names = F)
fun_rainsamp.p = fun_rainsamp.anova$uni.p	#Access p-values
fun_rainsamp.p.df <- as.data.frame(fun_rainsamp.p)	#Change to df
fun_rainsamp.p.df <- read.csv("fun_rainsamp.p.df.csv")	#Read in significant ESVs.
sigfunrainsamp_spp = colnames(fun_rainsamp.p.df)[fun_rainsamp.p.df[4,]<=0.05]	#Select significant p-values
write.csv(sigfunrainsamp_spp, file = "../data/sigfunrainsamp_spp.csv")	#Save list of significant ESVs
sigfunrainsamp.pi <- pi0est(p = fun_rainsamp.anova$uni.p, 
	lambda = seq(0.1, 0.9, 0.1), pi0.method = "smoother", smooth.log.pi0="TRUE")	#Find optimal pi0 value
q_sigfunrainsamp <- qvalue(p = fun_rainsamp.anova$uni.p,
	pi0 = sigfunrainsamp.pi$pi0)	#Calculate q-values
q_sigfunrainsamp.df <- as.data.frame(q_sigfunrainsamp$qvalues)	#Access q-values
write.csv(q_sigfunrainsamp.df, file = "../data/q_sigfunrainsamp.df.csv", row.names = F)	#Save q-values
sigfunrainsamp.q <- colnames(q_sigfunrainsamp.df)[q_sigfunrainsamp.df[4,] <= 0.05]	#Select significant q-values
sigfunrainsamp <- subset(sigfunrainsamp_spp, x %in% sigfunrainsamp.q$x)
write.csv(sigfunrainsamp, file = "../data/sigfunrainsamp.csv", row.names = F)	#Save list of significant ESVs

##Graphically represent strength of interactive effects between explanatory
##variables and significantly associated fungal ESVs.
sigfunrainsamp.df <- sub_fungi[,c(1:2, na.omit(match(sigfunrainsamp$x,
	colnames(sub_fungi))))]	#Subset dataframe by significant ESVs
write.csv(sigfunrainsamp.df, file = "../data/sigfunrainsamp.df.csv", row.names = F)	#Save df
#sigfunrainsamp.df <- read.csv("sigfunrainsamp.df.csv")	#Read in significant ESV df
sigfunrainsamp.abund <- mvabund(sigfunrainsamp.df[,3:length(sigfunrainsamp.df)])	#mvabund
Gradient <- sigfunrainsamp.df[,1]	#Rename variable for x-axis
Sample_type <- sigfunrainsamp.df$SampleType	#Rename variable for x-axis
sigfunrainsamp.glm <- manyglm(sigfunrainsamp.abund ~ Gradient * Sample_type,
	family = "negative.binomial")	#Assess GLMs for significant ESVs
pdf(file = "../figures/sigfunrainsamp_coefplot1.pdf")	#Prepare to save plot
par(mar=c(5,7,3,1))	#Change margins of imaging device
coefplot.manyglm(sigfunrainsamp.glm, which.Xcoef=10, mfrow = NULL)	#Plot GLM coefficients
dev.off()	#Deactivate plotting device



###Plot numbers of significant microbes by explanatory variable

significant_microbes <- multi_panel_figure(width = c(178, 356),
	height = c(178, 178)) %>%
fill_panel("../figures/sigbacgrad.pdf") %>%
fill_panel("../figures/bacsamp.pdf") %>%
fill_panel("../figures/sigfungrad.pdf") %>%
fill_panel("../figures/funsamp.pdf")

significant_microbes %>% save_multi_panel_figure(filename = "../figures/significant_microbes.tiff")

rm(list = ls())

##End (*not run*)