# Code by Austin Greene
# Purpose: Generate a map of Waimea Valley, Oahu with elevation data, 
# and site locations which are scaled the mean-annual precipitation 

# Load in packages
library(FedData)
library(latticeExtra)
library(scales)
library(viridis)
library(ggmap)
library(phyloseq)
library(raster)
library(rasterVis)
library(rgdal)
library(RColorBrewer)
library(ggplot2)

# Set wd to where physeq object lives
setwd("/Users/austingreene/Documents/UH Manoa/BOT662/BOT662_Analysis/Microbial_Data_Anthony")

# Load in physeq object
physeq1=readRDS("physeq1")

# Can exclude these lines. Changing location to a place where elevation data will be stored after download. 
# Also where final plots are saved. 
setwd("/Users/austingreene/Documents/UH Manoa/BOT662/BOT662_Analysis/ElevData")

# Getting base map via Stamen of location within bounding box. 
# Using Stamen maps which do not require a Google API key
map3 <- get_map(location = c(left = -158.070631, bottom = 21.598356, right = -157.998464, top = 21.655101), zoom=13, source="stamen", maptype = c("terrain-background"))
Waimea_map_3 <- ggmap(map3) # Turn into ggmap graphical object
Waimea_map_3 # Check that it looks good 

# Generating extents of map to pull elevation
bb <- attr(map3, "bb") # Mapp attribute object, from which we extract map extents 
extentB <- polygon_from_extent(raster::extent(bb$ll.lon, bb$ur.lon, bb$ll.lat, bb$ur.lat),
                               proj4string = "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs")

# Confirm bounding box for our elevation data
ggmap(map3) + geom_polygon(data=extentB, fill=NA, aes(x=long,y=lat), color="red", size=3)

# Download elevation data tiles
elev_waimea <- get_ned(template = extentB, label = "ned_waimea", res="1", force.redo = F)

# Make df of elevation data with lat and long included
elev_raster_df <- raster::as.data.frame(elev_waimea, xy=TRUE) #xy True 

# Make a dataframe of geographic location and mean annual precip at each site (from existing metadata)
geo_p=cbind(sample_data(physeq1)$Long, sample_data(physeq1)$Lat, sample_data(physeq1)$rain, sample_data(physeq1)$FieldSite)
geo3 = data.frame(geo_p) #we convert it to a dataframe
colnames(geo3)=c("Long", "Lat","Rain", "Site") # Repair column names
geo3 <- subset(geo3, Site != 9) # Remove Site 9 
geo3$Site <- factor(geo3$Site) # Convert sites to factors, not critical

# Plot elevational data and site locations scaled by precipitation on top of existing ggmap of Waimea Valley
Waimea_map_6_sizescaled <- Waimea_map_3 + # Existing map
  geom_raster(data = elev_raster_df, aes(x=x, y=y, fill=elev_raster_df$ned_waimea_NED_1)) + # Raster elevation data
  geom_point(data = geo3, aes(x=Long, y=Lat, size=Rain), fill="white", color="black", pch=21, alpha=0.5) + # Points for sites, scaled by mean-annual precipitation
  coord_cartesian() +  # Set to cartesian coordinates
  scale_fill_viridis(option = "plasma", alpha = 0.7) + # Set fill colors to be colorblind-friendly
  coord_fixed(1.3) + # Aspect ratio of cartesian coordinates is 1.3
  labs(x="Longitude", y="Latitude", size="Precipitation (mm/year)", fill="Elevation (m)") + # Labels
  theme_classic() # Classic minimalist theme

Waimea_map_6_sizescaled # Print the plot

# Save the plot as a tiff image with 300 dpi. 
# Setting to 6x6 inches of the plot with the expectation that this is 2x larger than publication size 
# per publication guidelines (expected print size 3x3 in)
ggsave("Map6_SizeScaled.tiff", Waimea_map_6_sizescaled, height = 6, width = 6, units = "in", device = "tiff", dpi=300)
