##### Script to create Figures 3, 5, and 5 and conduct the SHAP analysis for van Galen et. al. (2025) "The biogeography and conservation of Earth’s ‘dark’ ectomycorrhizal fungi"

library(terra) # v1.7-78
library(raster) # v3.6-26
library(berryFunctions) # v1.22.5
library(viridisLite) # v0.4.2
library(scales) # v1.3.0
library(ggplot2) # v3.5.1
library(h2o) # v3.44.0.3
library(tidyverse) # v2.0.0
library(gridExtra) # v2.3
library(ggpubr) # v0.6.0
library(geodata) # v0.6-2
library(webr) # v0.1.5
library(sf) # v1.0-19


##################################################################################################################################################
# SHAP analysis
##################################################################################################################################################

#### SHAP model ####################################################################################################################
sample1=read.csv("./Supps_data/Dark EcM taxa richness maps/SHAP_training_data.csv") # values from 1000 random grid cells selected across the globe
family.list=readRDS("./host_sp_download_files/data_families/family.list.rds") # list of host plant families
family.list2=paste0(family.list," host presence") # edit host plant names
names(family.list2)<-family.list

# Initialize H2O
h2o.init()

# Subset only the columns of interest for making a random forest model
training_data=sample1[,c("dark_taxa_hotspot_index","CHELSA_BIO_Annual_Mean_Temperature","CHELSA_BIO_Annual_Precipitation","SG_Soil_pH_H2O_005cm",
                         "SG_SOC_Content_005cm","GRIP4_DistanceToAllRoads","GHS_Population_Density","ecm_sampleintensity_5degrees",
                         "veg_prop_ECM_Jinsu","host_family_richness",family.list)] 

training_data <- training_data %>% mutate_at(vars(c(11:35)), as.factor) # change host vars to factor
training_data=training_data[ , -which(names(training_data) %in% c("Ticodendraceae","Sarcolaenaceae","Asteropeiaceae"))] # drop hosts with too small distributions

regressionMatrix <- as.h2o(training_data) # Convert to H2O regression matrix

# Train a random forest model
rfModel <- h2o.randomForest(training_frame = regressionMatrix,
                            y = "dark_taxa_hotspot_index",
                            ntrees = 999, 
                            seed = 42)

h2o.performance(rfModel)

# SHAP summary plot (takes time, 5 or 10 mins)
p1=h2o.shap_summary_plot(
  model = rfModel,
  top_n_features = 40,
  newdata = regressionMatrix
)
# p1

saveRDS(rfModel,"output_shap_model.rds")
saveRDS(p1,"output_shap_plot.rds")

rfModel=readRDS("output_shap_model.rds")
h2o.performance(rfModel)
h2o.r2(rfModel, xval = FALSE)  # model r2 performance

############## Donut plot for Figure 3E #####################################################################################
shap_plot=readRDS("output_shap_plot.rds")
shap_data=shap_plot$data

shap_data$category=as.factor(ifelse(shap_data$feature=="Dipterocarpaceae"|shap_data$feature=="Juglandaceae"|shap_data$feature=="Rosaceae"|shap_data$feature=="Myrtaceae"|shap_data$feature=="Malvaceae"
                                    |shap_data$feature=="Polygonaceae"|shap_data$feature=="Fagaceae"|shap_data$feature=="Fabaceae"|shap_data$feature=="Gnetaceae"
                                    |shap_data$feature=="host_family_richness"|shap_data$feature=="veg_prop_ECM_Jinsu"
                                    |shap_data$feature=="Cistaceae"|shap_data$feature=="Asteraceae"|shap_data$feature=="Betulaceae"|shap_data$feature=="Pinaceae"
                                    |shap_data$feature=="Ticodendraceae"|shap_data$feature=="Sarcolaenaceae"|shap_data$feature=="Salicaceae"|shap_data$feature=="Nyctaginaceae"
                                    |shap_data$feature=="Apiaceae"|shap_data$feature=="Asteropeiaceae"|shap_data$feature=="Goodeniaceae"|shap_data$feature=="Cyperaceae"
                                    |shap_data$feature=="Phyllanthaceae"|shap_data$feature=="Sarcolaenaceae"|shap_data$feature=="Nothofagaceae"|shap_data$feature=="Casuarinaceae"
                                    |shap_data$feature=="Ericaceae"|shap_data$feature=="Rhamnaceae"|shap_data$feature=="Achatocarpaceae","Vegetation",
                                    ifelse(shap_data$feature=="ecm_sampleintensity_5degrees"|shap_data$feature=="GRIP4_DistanceToAllRoads"|shap_data$feature=="ConsensusLandCover_Human_Development_Percentage"
                                           |shap_data$feature=="GHS_Population_Density","Human",
                                           ifelse(shap_data$feature=="CHELSA_BIO_Annual_Precipitation"|shap_data$feature=="CHELSA_BIO_Annual_Mean_Temperature"
                                                  |shap_data$feature=="SG_SOC_Content_005cm"|shap_data$feature=="SG_Soil_pH_H2O_005cm","Environment","!!!!!"))))

# calculate the overall pie segment percentages for each category
sum=aggregate(contribution~category+id.x,data=shap_data,FUN="sum") # add up the shap values in each group
sum_final=aggregate(abs(contribution)~category, data=sum, FUN=mean) # calculate group means of absolute values
names(sum_final)<-c("category","contribution")
PieDonut(sum_final, aes(pie=category,count=contribution),showRatioThreshold=0.01)
sum_final$perc=(sum_final$contribution/sum(sum_final$contribution))*100

# apply the adjustment for the values of individual variables
shap_data3=aggregate(abs(contribution)~feature+category, data=shap_data, FUN=mean) # calculate mean shap for individual variables
names(shap_data3)<-c("feature","category","contribution")
sum_orig=aggregate(contribution~category, data=shap_data3, FUN="sum") # add them up to find how much each needs to be adjusted
shap_data3$contribution_perc <- ifelse(shap_data3$category=="Environment",shap_data3$contribution/sum_orig$contribution[1],
                                       ifelse(shap_data3$category=="Human",shap_data3$contribution/sum_orig$contribution[2],shap_data3$contribution/sum_orig$contribution[3])) # the % of each cat that each var makes up (divide each by the original sum - sum of mean abs shap)
shap_data3$contribution_adjusted <- ifelse(shap_data3$category=="Environment",shap_data3$contribution_perc*sum_final$contribution[1],
                                           ifelse(shap_data3$category=="Human",shap_data3$contribution_perc*sum_final$contribution[2],shap_data3$contribution_perc*sum_final$contribution[3])) # multiply by the new sum (mean for each cat after summing the raw shap values)

shap_data3 <- shap_data3 %>% mutate(feature = fct_reorder(feature, contribution_adjusted, .desc = F)) # order segments by size

PieDonut(shap_data3, aes(pie=category,donut=feature,count=contribution_adjusted),showRatioThreshold=0.01)


pdf("./figure_pieDonut_shap.pdf",width=4,height=4)
PieDonut(shap_data3, aes(pie=category,donut=feature,count=contribution_adjusted),
         r0 = 0.45, r1 = 0.9, showRatioThreshold = 1, pieLabelSize = 0, # Note, switch showRatioThreshold = 1 to remove labels
         showPieName = F, showRatioPie = FALSE, addDonutLabel = F)
dev.off()


############## Figure S3 #####################################################################################
shap_plot=readRDS("output_shap_plot.rds")

p_supps2=ggplot(shap_plot$data, aes(.data$feature, .data$contribution, color = .data$normalized_value, text = .data$row)) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = h2o:::position_jitter_density(), alpha = 0.5, size = 0.5) +
  scale_color_gradient(low = "#00AAEE", high = "#FF1166") +
  coord_flip() +
  labs(y = "SHAP value", x = NULL, color = "Normalised\ncovariate\nvalues",title="B)") +
  scale_x_discrete(labels=labels) +
  theme(
    axis.title.y = element_blank(),  # Remove redundant y-axis label
    axis.text.y = element_blank(),   # Remove redundant y-axis text
    axis.ticks.y = element_blank()  # Remove redundant y-axis ticks
  )+
  theme(legend.position = c(0.8, 0.17),legend.title=element_text(size=8),legend.text=element_text(size=7),plot.title = element_text(size=12)) 

### calculate mean importance ###
shap_data=shap_plot$data
shap_ag=aggregate(abs(contribution)~feature, data=shap_data, FUN=mean) # calculate group means of absolute values
names(shap_ag)<-c("feature","contribution")

family.list=readRDS("./host_sp_download_files/data_families/family.list.rds")
family.list2=paste0(family.list," host presence")
names(family.list2)<-family.list
labels=c("ectomycorrhizal_richness_V1"="Total EcM fungal richness","CHELSA_BIO_Annual_Mean_Temperature"="Annual mean temperature","CHELSA_BIO_Annual_Precipitation"="Annual mean precipitation","SG_Soil_pH_H2O_005cm"="Soil pH",
         "SG_SOC_Content_005cm"="Soil organic carbon","GRIP4_DistanceToAllRoads"="Distance to roads","GHS_Population_Density"="Population density","ConsensusLandCover_Human_Development_Percentage"="Human development percentage",
         "ecm_sampleintensity_5degrees"="Training data sampling density",
         "veg_prop_ECM_Jinsu"="% vegetation biomass comprising EcM hosts","ECM_rwr_sampledensity"="EcM fungal rarity index","host_family_richness"="Host plant family richness",family.list2)


p_supps3=ggplot(shap_ag, aes(feature, contribution)) +
  theme_bw() +
  coord_flip()+
  geom_bar(stat="identity")+
  scale_x_discrete(labels=labels) +
  labs(y = "Mean absolute SHAP value", x = NULL,title="A)") +
  theme(plot.title = element_text(size=12))

png("./figure_shap_supps.png",width=8,height=6.5,res=600,units="in")
ggarrange(p_supps3, p_supps2, ncol = 2, widths = c(1.8, 1), align="h")
dev.off()



###################################################################################################################################################
# Make maps for Figure 3
###################################################################################################################################################
# V1 refers to the total EcM richness map, U1 to the dark taxa EcM richness map

maps <- rast("./original_maps/analysis_layers_masked_95_all_10km_robinson.tif") # the layers reprojected to Robinson projection. 
background=rast("./Background/Background_robinson_10km.tif")
er <- rast(ext(-17500000,17500000,-6100000,8565000), resolution=res(background),crs=crs(background)) # empty raster for outline, set extent manually to chop off bottom
values(er) <- 0


######### colour settings #########
n_viridis <- 30 
viridis_colors <- viridisLite::viridis(n_viridis)
repetition_pattern <- rep(seq_len(ceiling(length(viridis_colors)/2)), each = 2)[1:length(viridis_colors)]
skewed_viridis <- rep(viridis_colors, times = repetition_pattern)
skewed_viridis <- c(skewed_viridis,rep(skewed_viridis[length(skewed_viridis)],50))

# Viridis palette settings
min2 <- min(minmax(maps$ectomycorrhizal_richness_V1)[[1]], minmax(maps$ectomycorrhizal_richness_U1)[[1]])
max2 <- max(minmax(maps$ectomycorrhizal_richness_V1)[[2]], minmax(maps$ectomycorrhizal_richness_U1)[[2]])
color_sequence2 <- seq(min2, max2, length.out = length(skewed_viridis))

lim_num=minmax(maps$ectomycorrhizal_richness_U1)[[2]]/minmax(maps$ectomycorrhizal_richness_V1)[[2]] # what proportion is U1 max of V1 max
lim_num2=length(color_sequence2)*lim_num # now use this to adjust the seqeunce2 value of the U1 plot so the histogram only goes as far as the max

graticule <- st_graticule(lat = seq(-60, 90, 30), lon = seq(-180, 180, 30))
graticule <- st_transform(graticule, crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

########## Plots for Figure 3A, B, C ##########################################################################################################################################
# png("./figure_base_maps.png",width=12,height=2,units="in",res=800)
# pdf("./figure_base_maps.pdf",width=12,height=2)

cex1=0.8
cex2 <- 1
lwd1=0.5
mar1=c(0,0.1,0.5,0.1)
maxcell1=10^7 # 10^7, set lower when experimenting
par(mfrow=c(1,3),mar=c(0.5,0.5,0.5,0.5))

plot(er,legend=F,col=NA,maxcell=maxcell1,mar=mar1,box=F,pax=list(col=NA,col.axis="grey30"),axes=F) # plot empty raster for setting the margins
plot(graticule$geometry,add=T, col = "grey80", lty = 1,lwd=lwd1) # add a grid
plot(background,add=T,legend=F,axes=F,col=c(NA,"grey85"),maxcell=maxcell1) # col NA makes it transparent
plot(maps$ectomycorrhizal_richness_V1, add=T, maxcell=maxcell1, col = skewed_viridis, breaks = color_sequence2, legend = FALSE)
smallPlot(hist(maps$ectomycorrhizal_richness_V1, breaks = color_sequence2, col = skewed_viridis, border = skewed_viridis, yaxt = 'n', ylab = "", xlab = "", main = "", cex.lab = cex2, cex.axis = cex2, mgp = c(1, 0.2, 0), tck = -0.05),
          x1 = 0.05, x2 = 0.27, y1 = 0.27, y2 = 0.47, mar = c(0, 0, 0, 0), mgp = c(1.3, 0.5, 0), border = FALSE, bg = NA)

plot(er,legend=F,col=NA,maxcell=maxcell1,mar=mar1,box=F,pax=list(col=NA,col.axis="grey30"),axes=F) # plot empty raster for setting the margins
plot(graticule$geometry,add=T, col = "grey80", lty = 1,lwd=lwd1) # add a grid
plot(background,add=T,legend=F,axes=F,col=c(NA,"grey85"),maxcell=maxcell1) # col NA makes it transparent
plot(maps$ectomycorrhizal_richness_U1, add=T, maxcell=maxcell1, col = skewed_viridis, breaks = color_sequence2, legend = FALSE)
smallPlot(hist(maps$ectomycorrhizal_richness_U1, breaks = color_sequence2[c(1:lim_num2)], col = skewed_viridis, border = skewed_viridis, yaxt = 'n', ylab = "", xlab = "", 
               main = "", cex.lab = cex2, cex.axis = cex2, mgp = c(1, 0.2, 0), tck = -0.05, xlim=c(0,max2)),
          x1 = 0.05, x2 = 0.27, y1 = 0.27, y2 = 0.47, mar = c(0, 0, 0, 0), mgp = c(1.3, 0.5, 0), border = FALSE, bg = NA)

plot(er,legend=F,col=NA,maxcell=maxcell1,mar=mar1,box=F,pax=list(col=NA,col.axis="grey30"),axes=F) # plot empty raster for setting the margins
plot(graticule$geometry,add=T, col = "grey80", lty = 1,lwd=lwd1) # add a grid
plot(background,add=T,legend=F,axes=F,col=c(NA,"grey85"),maxcell=maxcell1) # col NA makes it transparent
plot(maps$percentage_undescribed_clamp, add=T, maxcell=maxcell1, col = viridis(100), legend = FALSE)
smallPlot(hist(maps$percentage_undescribed_clamp, breaks = seq(minmax(maps$percentage_undescribed_clamp)[[1]], minmax(maps$percentage_undescribed_clamp)[[2]], length.out = 100), col = viridis(100), border = viridis(100), yaxt = 'n', ylab = "", xlab = "", main = "", cex.lab = cex2, cex.axis = cex2, mgp = c(1, 0.2, 0), tck = -0.05),
          x1 = 0.05, x2 = 0.27, y1 = 0.27, y2 = 0.47, mar = c(0, 0, 0, 0), mgp = c(1.3, 0.5, 0), border = FALSE, bg = NA)

# dev.off()


########## Plot for Figure 3D ##########################################################################################################################################

# Skewed magma palette with more yellow
n=100
magma_colors <- viridisLite::inferno(n)
values <- seq(0, 1, length.out = n)
skewed_values <- scales::rescale(values, to = c(0, 1), from = c(0, 0.5)) # adjust last number from to skew
skewed_magma <- scales::gradient_n_pal(magma_colors)(skewed_values)
skewed_magma[is.na(skewed_magma)] <- magma_colors[length(magma_colors)]

# png("./figure_metric_maps.png",width=8,height=3.8,units="in",res=800)
# pdf("./figure_metric_maps.pdf",width=8,height=3.8)

cex1=0.8
cex2 <- 0.7
lwd1=0.5
mar1=c(0,0.1,0.5,0.1)
maxcell1=10^7 # 10^7
par(mfrow=c(1,1),mar=c(0.5,0.5,0.5,0.5))

plot(er,legend=F,col=NA,maxcell=maxcell1,mar=mar1,box=F,pax=list(col=NA,col.axis="grey30"),axes=F) # plot empty raster for setting the margins
plot(graticule$geometry,add=T, col = "grey80", lty = 1,lwd=lwd1) # add a grid
plot(background,add=T,legend=F,axes=F,col=c(NA,"grey85"),maxcell=maxcell1) # col NA makes it transparent
plot(maps$dark_taxa_hotspot_index, add=T, maxcell=maxcell1, col = skewed_magma, breaks = seq(minmax(maps$dark_taxa_hotspot_index)[[1]], minmax(maps$dark_taxa_hotspot_index)[[2]], length.out = n), legend = FALSE)
smallPlot(hist(maps$dark_taxa_hotspot_index, breaks = seq(minmax(maps$dark_taxa_hotspot_index)[[1]], minmax(maps$dark_taxa_hotspot_index)[[2]], length.out = n),maxcell=10^7,
               col = skewed_magma, border = skewed_magma, yaxt = 'n', ylab = "", xlab = "", main = "", cex.lab = cex2, cex.axis = cex2, mgp = c(1, 0.2, 0), tck = -0.05),
          x1 = 0.05, x2 = 0.27, y1 = 0.27, y2 = 0.47, mar = c(0, 0, 0, 0), mgp = c(1.3, 0.5, 0), border = FALSE, bg = NA)
# dev.off()





###################################################################################################################################################
# Figure 4 - 95th percentile priority area 
###################################################################################################################################################
# First, transform the equal-earth projection so that grid cells are all equal size
# Then subset the research priority area metric map to select grid cells in the upper 95th percentile
# Then calculate the area (km2) of the subsetted region that falls within each biome and host family range. Do this at both global level and within each continent

# # reproject to EE projection
# # maps <- rast("./original_maps/analysis_layers_masked_95_all.tif") 
# # template=rast("./host_sp_download_files/data_families/family_rasters/raster_Achatocarpaceae_thresh200km_buffer150km.tif") # EE template
# # mapsEE=project(maps$dark_taxa_hotspot_index,template) # change to EE
# # templateMask=rast("./original_maps/ectomycorrhizal_richness_Classified_MultibandImage_final_EE.tif")["ectomycorrhizal_richness_Ensemble_mean"] # EE template
# # mapsEE2=mask(mapsEE,templateMask) # mask to remove excess area (can't use this for orig template because crs not showing)
# # writeRaster(mapsEE2, filename="./original_maps/analysis_layers_masked_95_metric_EE.tif",overwrite=TRUE) # save as tif
# 
# biome=rast("Resolve_Biome_combined_1km_EE.tif")
# maps=rast("./original_maps/analysis_layers_masked_95_metric_EE.tif")
# 
# # select top 95th percentile
# coVar_quantile <- quantile(raster(maps), c(0.95), na.rm=T) # 95th percentile
# metric_hotspot=maps>coVar_quantile
# plot(metric_hotspot,nr=1,col=c("grey","red"))
# metric_hotspot2=subst(metric_hotspot, from=0, to=NA) # make zeros NA
# writeRaster(metric_hotspot2, filename="./original_maps/darkspot_map_EE.tif",overwrite=TRUE) # save as tif
# 
# 
# ## extract continents #############
# countries <- world(resolution = 5, path = "maps")  # you may choose a smaller (more detailed) resolution for the polygon borders, and a different folder path to save the imported map
# cntry_codes <- country_codes() # import a table with country codes and continents
# countries <- merge(countries, cntry_codes, by.x = "GID_0", by.y = "ISO3", all.x = TRUE) # add this table to the countries map attributes
# countries.EE=project(countries,template) # reproject to EE
# continent.list=c(unique(countries$continent))
# continent.list=continent.list[continent.list != "Antarctica"] # remove antarctica
# 
# #### biomes ###################################################################################################################################
# area_output=NULL
# biome_df=data.frame(value=c(1:14))
# for(j in 1:length(continent.list)){ # 5 mins
#   print(continent.list[j])
#   print(Sys.time())
#   continent1=subset(countries.EE, countries.EE$continent == continent.list[j])
#   biome_continent_masked=mask(biome,continent1) # takes time
#   biome_continent_masked2=mask(biome_continent_masked,metric_hotspot2) # mask by hotspot raster - select hotspot regions for that continent
# 
#   biome_vals=freq(biome_continent_masked2) # extract count of 1km2 cells for each biome category
#   area=merge(biome_df,biome_vals,by="value",all.x=T) # merge with base df to add missing biomes
#   area$count[is.na(area$count)] <- 0 # replace NAs with zero
#   area$labels=rev(c("Mangroves","Deserts","Mediterranean Forests","Tundra","Montane grasslands","Flooded Grasslands","Temperate Grasslands","Tropical Grasslands","Boreal Forests","Temperate Conifer Forests","Temperate Broadleaf Forests","Tropical Conifer Forests","Tropical Dry Forests","Tropical Moist Forests"))
#   area$continent=continent.list[j] # ad continent name
#   
#   area_output=rbind(area_output,area)
# }
# 
# # add global
# biome_mask=mask(biome,metric_hotspot2)
# biome_vals=freq(biome_mask)
# area=merge(biome_df,biome_vals,by="value",all.x=T) # merge with base df to add missing biomes
# area$count[is.na(area$count)] <- 0 # replace NAs with zero
# area$labels=rev(c("Mangroves","Deserts","Mediterranean Forests","Tundra","Montane grasslands","Flooded Grasslands","Temperate Grasslands","Tropical Grasslands","Boreal Forests","Temperate Conifer Forests","Temperate Broadleaf Forests","Tropical Conifer Forests","Tropical Dry Forests","Tropical Moist Forests"))
# area$continent="Global"
# area_output2=rbind(area_output,area)
# area_output2$continent=factor(area_output2$continent,ordered=T,levels=c("Global","Africa","Asia","Europe","North America","Oceania","South America"))
# area_output2$value=as.factor(area_output2$value)
# saveRDS(area_output2,"darkspot_area_biome.rds")
# 
# #### host plants ###################################################################################################################################
# # read in host rasters
# family.list=readRDS("./host_sp_download_files/data_families/family.list.rds")
# family.list=family.list[ ! family.list %in% c("Cyperaceae","Apiaceae","Asteropeiaceae","Sarcolaenaceae","Ticodendraceae") ] # drop families that we are no longer using
# # family.list=family.list[c(1:2)]
# 
# raster_list = list() # using EE projection
# for (i in 1:length(family.list)) {
#   raster_path = paste0("./host_sp_download_files/data_families/family_rasters/raster_",gsub(" ", "_", family.list[i]),"_thresh200km_buffer150km.tif")
#   raster_list[[i]] = rast(raster_path)
#   names(raster_list[[i]])<-family.list[i]
# }
# family_stack = do.call(c, raster_list)
# 
# area_output_host=NULL
# family_df=data.frame(layer=c(1:length(family.list)))
# for(j in 1:length(continent.list)){ # 10 mins
#   print(continent.list[j])
#   print(Sys.time())
#   continent1=subset(countries.EE, countries.EE$continent == continent.list[j])
#   fam_continent_masked=mask(family_stack,continent1) # takes time
#   fam_continent_masked2=mask(fam_continent_masked,metric_hotspot2) # mask by hotspot raster - select hotspot regions for that continent
#   
#   fam_vals=freq(fam_continent_masked2) # extract count of 1km2 cells for each biome category
#   area=merge(family_df,fam_vals,by="layer",all.x=T) # merge with base df to add missing fams
#   area$count[is.na(area$count)] <- 0 # replace NAs with zero
#   area$labels=family.list
#   area$continent=continent.list[j] # ad continent name
#   
#   area_output_host=rbind(area_output_host,area)
# }
# 
# # add global
# full_mask=mask(family_stack,metric_hotspot2)
# full_vals=freq(full_mask)
# area=merge(family_df,full_vals,by="layer",all.x=T) # merge with base df to add missing fams
# area$count[is.na(area$count)] <- 0 # replace NAs with zero
# area$labels=family.list
# area$continent="Global"
# area_output_host2=rbind(area_output_host,area)
# area_output_host2$continent=factor(area_output_host2$continent,ordered=T,levels=c("Global","Africa","Asia","Europe","North America","Oceania","South America"))
# area_output_host2$value=as.factor(area_output_host2$value)
# saveRDS(area_output_host2,"darkspot_area_host.rds")


## Make Figure 4 ###############################
host=readRDS("darkspot_area_host.rds")
bio=readRDS("darkspot_area_biome.rds")

labels1=rev(c("Mangroves","Deserts","Mediterranean forests","Tundra","Montane grasslands","Flooded grasslands","Temperate grasslands","Tropical grasslands","Boreal forests","Temperate conifer forests","Temperate broadleaf forests","Tropical conifer forests","Tropical dry forests","Tropical moist forests"))
bio <- bio %>% mutate(labels = labels1[value]) # add labels, so they are lower case

bio$labels <- gsub("Temperate", "Temp.", bio$labels)
bio$labels <- gsub("Tropical", "Trop.", bio$labels)
bio$labels <- gsub("broadleaf", "broad.", bio$labels)


# Select top 5 categories per continent and order within each group
df_top_bio <- bio %>%
  group_by(continent) %>%
  arrange(desc(count)) %>%
  slice_head(n = 5) %>%
  mutate(labels = paste0(labels, " (", continent, ")")) %>%
  mutate(labels = factor(labels, levels = labels[order(count)])) %>%
  mutate(y_indent = 0.03 * max(count / 100000)) %>%  # Small indent relative to max value
  ungroup()

# Create bar plot
plotb3=ggplot(df_top_bio, aes(x = labels, y = count / 100000)) +
  geom_col(fill=rgb(255/255, 140/255, 0, 0.5)) + # darkorange with 50% opacity for fill
  geom_text(aes(label = sub(" \\(.*\\)", "", labels), y = y_indent), color = "black", size = 2.6,hjust=0) +  # Labels inside bars
  coord_flip() +
  facet_wrap(~continent, scales = "free", nrow = 1) +
  scale_x_discrete(labels = function(x) NULL) +  # Remove axis labels
  theme_minimal() +
  theme(axis.text.y = element_blank(),  # Remove text labels on y-axis
        legend.position = "none",plot.title = element_text(size=11)) +
  labs(x = NULL, y = NULL, title="A) Biomes")

# Select top 5 categories per continent and order within each group
df_top <- host %>%
  group_by(continent) %>%
  arrange(desc(count)) %>%
  slice_head(n = 10) %>%
  mutate(labels = paste0(labels, " (", continent, ")")) %>%
  mutate(labels = factor(labels, levels = labels[order(count)])) %>%
  mutate(y_indent = 0.03 * max(count / 100000)) %>%  # Small indent relative to max value
  ungroup()

# Create bar plot
ploth3=ggplot(df_top, aes(x = labels, y = count / 100000)) +
  geom_col(fill=rgb(34/255, 139/255, 34/255, 0.5)) + # Forest green with 50% opacity for fill
  geom_text(aes(label = sub(" \\(.*\\)", "", labels), y = y_indent), color = "black", size = 2.6,hjust=0) +  # Labels inside bars
  coord_flip() +
  facet_wrap(~continent, scales = "free", nrow = 1) +
  scale_x_discrete(labels = function(x) NULL) +  # Remove axis labels
  theme_minimal() +
  theme(axis.text.y = element_blank(),  # Remove text labels on y-axis
        legend.position = "none") +
  theme( strip.text.x = element_blank(),plot.title = element_text(size=11)) +# remove facet labels
  labs(x = NULL, y = "Priority 'darkspot' area (x 100,000 km2)",title="B) Host plant families")

png("./figure_darkspot_biome_host_panel_top.png",width=8.6,height=3.5,res=600,units="in")
pdf("./figure_darkspot_biome_host_panel_top.pdf",width=8.6,height=4)
ggarrange(plotb3, ploth3, ncol = 1, heights = c(0.67, 1), align="v")
dev.off()


#### Make Figure S5 - darkspot map for supps ################################################################################################
# # change to robinson for plotting
# template=rast("original_maps/AM_richness_robinson.tif") # units in km
# crs(template) <-  gsub("LENGTHUNIT\\[\"kilometre\",1000", "LENGTHUNIT[\"metre\",1", crs(template)) # update units to meters
# ext(template) <- ext(template) * 1000
# spot=rast("./original_maps/darkspot_map_EE.tif") # save as tif
# spot2=project(spot,template,method="near",res=1000)
# spot3=mask(spot2,template) # remove corner bits
# spot3[is.na(spot3)] <- 0 # replace nas with 0 for aggregating
# spot4=aggregate(spot3,fact=10,fun="modal")
# spot4[spot4 == 0] <- NA # put NAs back
# plot(spot4,col="red")
# writeRaster(spot4, filename="./original_maps/darkspot_map_robinson_10km.tif",overwrite=TRUE) # save as tif

spot=rast("./original_maps/darkspot_map_robinson_10km.tif") # save as tif

background=rast("./Background/Background_robinson_10km.tif")
er <- rast(ext(-17500000,17500000,-6100000,8565000), resolution=res(background),crs=crs(background)) # empty raster for outline, set extent manually to chop off bottom
values(er) <- 0
graticule <- st_graticule(lat = seq(-60, 90, 30), lon = seq(-180, 180, 30))
graticule <- st_transform(graticule, crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

png("./figure_darkspot_map_supps.png",width=8,height=4.2,res=800,units="in")
cex1=0.8
cex2 <- 1
lwd1=0.5
mar1=c(0,0.1,0.5,0.1)
maxcell1=10^7 # 10^7
par(mfrow=c(1,1),mar=c(0.5,0.5,0.5,0.5))

plot(er,legend=F,col=NA,maxcell=maxcell1,mar=mar1,box=F,pax=list(col=NA,col.axis="grey30"),axes=F) # plot empty raster for setting the margins
plot(graticule$geometry,add=T, col = "grey80", lty = 1,lwd=lwd1) # add a grid
plot(background,add=T,legend=F,axes=F,col=c(NA,"grey85"),maxcell=maxcell1) # col NA makes it transparent
plot(spot$dark_taxa_hotspot_index, add=T, maxcell=maxcell1, col = "red", legend = FALSE)
# mtext(expression("Upper 95"^"th"*" percentile of dark taxa research priority metric values"),line=-1,adj=0,cex=1.3)

dev.off()





##################################################################################################################################################
# Figure 5 - mean metric scores for each biome and host family range
##################################################################################################################################################
# first, extract 1000 grid cells from priority metric map with each host family range and biome. Do this at the global level, but also for each continent.

# maps <- rast("./original_maps/analysis_layers_masked_95_all.tif")["dark_taxa_hotspot_index"] # load the research priority metric raster

# ###### random sample of families #############################################################################################
# family.list=readRDS("./host_sp_download_files/data_families/family.list.rds") # list of host families
# 
# # take random sample within each family
# family.coords=NULL
# for (i in 1:length(family.list)) {
#   print(family.list[i])
#   family_rast = rast(paste0("./host_sp_download_files/data_families/family_rasters/raster_",gsub(" ", "_", family.list[i]),"_thresh200km_buffer150km_equirec.tif"))
#   maps2=mask(maps,family_rast)
#   n=1000 # set sample size
#   non_na_cells <- which(!is.na(values(maps2)))
#   set.seed(4534)
#   cell.num.sample=sample(non_na_cells,n)
#   coords=terra::extract(maps2,cell.num.sample,xy=T)
#   coords$Family=family.list[i]
#   family.coords=rbind(family.coords,coords)
# }
# saveRDS(family.coords,"random.sample.values.families.rds")
# 
# ## take random sample within each family but split by continent #############
# countries <- world(resolution = 5, path = "maps")  # you may choose a smaller (more detailed) resolution for the polygon borders, and a different folder path to save the imported map
# cntry_codes <- country_codes() # import a table with country codes and continents
# countries <- merge(countries, cntry_codes, by.x = "GID_0", by.y = "ISO3", all.x = TRUE) # add this table to the countries map attributes
# continent.list=c(unique(countries$continent))
# continent.list=continent.list[continent.list != "Antarctica"] # remove antarctica
# 
# output=NULL
# output_area=NULL
# for(j in 1:length(continent.list)){ # ~6 hours
#   print(continent.list[j])
#   print(Sys.time())
#   continent1=subset(countries, countries$continent == continent.list[j])
#   maps_continent_masked=mask(maps,continent1)
# 
#   results=c()
#   host_area=NULL
#   for(i in 1:length(family.list)){
#     print(paste(i,family.list[i]))
#     family_rast = rast(paste0("./host_sp_download_files/data_families/family_rasters/raster_",gsub(" ", "_", family.list[i]),"_thresh200km_buffer150km_equirec.tif"))
#     maps2=mask(maps_continent_masked,family_rast)
# 
#     # calculate area of each family
#     area=cellSize(maps2, unit = "km", mask=T)/1000 # surface area size of each cell. Divide by 1000 so units are 1000 km and don't get too big
#     sum1=global(area, sum, na.rm = TRUE)
#     sum1$Family=family.list[i]
#     sum1$continent=continent.list[j]
#     host_area=rbind(host_area,sum1)
# 
#     n=1000 # set sample size
#     non_na_cells <- which(!is.na(values(maps2)))
# 
#     if(length(non_na_cells)>0){ # only move on if more than zero cells
#       if(length(non_na_cells)>n){ # if more than n, subsample them
#         set.seed(4534)
#         cell.num.sample=sample(non_na_cells,n)
#         coords=terra::extract(maps2,cell.num.sample,xy=T)
#         coords$Family=family.list[i]
#         coords$continent=continent.list[j]
#         results=rbind(results,coords)
#       } else{ # otherwise use all
#         cell.num.sample=non_na_cells
#         coords=terra::extract(maps2,cell.num.sample,xy=T)
#         coords$Family=family.list[i]
#         coords$continent=continent.list[j]
#         results=rbind(results,coords)
#       }
#     }
# 
#   }
#   output=rbind(output,results)
#   output_area=rbind(output_area,host_area)
# }
# saveRDS(output,"./random.sample.values.families.continent.rds")
# saveRDS(output_area,"./random.sample.families.continent.area.1000km2.rds")
# 
# ###### random sample of biomes #############################################################################################
# # first selecting cells, then extract values after
# biome=rast("Resolve_Biome_combined_1km.tif") # resolve biome raster https://ecoregions.appspot.com/
# maps <- rast("./original_maps/analysis_layers_masked_95_all.tif") # the dark taxa maps
# biome_masked=mask(biome,maps[[1]]) # mask the biome in the same way as the maps
# 
# ## biomes at the global level ##############
# results=c()
# for(i in 1:14){ # took 3 minutes
#   print(i)
#   n=1000 # set sample size
#   cell.num=cells(biome_masked,i) # get cells that match value i
#   set.seed(534)
#   cell.num.sample=sample(cell.num$Resolve_Biome,n)
#   results=c(results,cell.num.sample)
# }
# saveRDS(results,"./random.sample.cells.biomes.global.rds")
# Sys.time()
# 
# ## biomes at the continent level #############
# countries <- world(resolution = 5, path = "maps")  # you may choose a smaller (more detailed) resolution for the polygon borders, and a different folder path to save the imported map
# cntry_codes <- country_codes() # import a table with country codes and continents
# countries <- merge(countries, cntry_codes, by.x = "GID_0", by.y = "ISO3", all.x = TRUE) # add this table to the countries map attributes
# continent.list=c(unique(countries$continent))
# continent.list=continent.list[continent.list != "Antarctica"] # remove antarctica
# 
# output=NULL
# biome_area=NULL
# for(j in 1:length(continent.list)){ # 20 mins
#   print(continent.list[j])
#   print(Sys.time())
#   continent1=subset(countries, countries$continent == continent.list[j])
#   biome_continent_masked=mask(biome_masked,continent1)
# 
#   # calculate area of each biome
#   area=cellSize(biome_continent_masked, unit = "km", mask=T)/1000 # surface area size of each cell. Divide by 1000 so units are 1000 km and don't get too big
#   zonal_sum <- zonal(area, biome_continent_masked, fun = "sum") # add up area for each biome
#   zonal_sum$continent=continent.list[j]
#   biome_area=rbind(biome_area,zonal_sum)
# 
#   results=c()
#   for(i in 1:14){
#     print(i)
#     n=1000 # set sample size
#     cell.num=cells(biome_continent_masked,i) # get cells that match value i
#     if(length(cell.num$Resolve_Biome)>0){ # only move on if more than zero cells
#       if(length(cell.num$Resolve_Biome)>n){ # if more than n, subsample them
#         set.seed(534)
#         cell.num.sample=sample(cell.num$Resolve_Biome,n)
#         results=c(results,cell.num.sample)
#       } else{
#         results=c(results,cell.num$Resolve_Biome) # otherwise keep all of them
#       }
#     }
#   }
#   results2=data.frame(cellID=results,continent=continent.list[j]) # save the continent ID
#   output=rbind(output,results2)
# }
# 
# saveRDS(output,"./random.sample.cells.biomes.continent.rds")
# saveRDS(biome_area,"./random.sample.biome.continent.area.1000km2.rds")
# Sys.time()
#
# ### extract values #########################################################################
# biome=rast("Resolve_Biome_combined_1km.tif")
# maps <- rast("./original_maps/analysis_layers_masked_95_all.tif")
# stack2=c(biome,maps)
# cells_bio_global=readRDS("./random.sample.cells.biomes.global.rds")
# cells_bio_cont=readRDS("./random.sample.cells.biomes.continent.rds")
#
# # global biomes (1 min)
# sample1=terra::extract(stack2, cells_bio_global, xy=T)
# saveRDS(sample1,"./random.sample.values.biomes.global.rds")
#
# # continent biomes (3 mins)
# sample2=terra::extract(stack2, cells_bio_cont$cellID, xy=T) # 3 mins
# sample2$continent=cells_bio_cont$continent # add on continent
# saveRDS(sample2,"./random.sample.values.biomes.continent.rds")
#


# Biome codes
# 1 Tropical & Subtropical Moist Broadleaf Forests
# 2 Tropical & Subtropical Dry Broadleaf Forests
# 3 Tropical & Subtropical Coniferous Forests
# 4 Temperate Broadleaf & Mixed Forests
# 5 Temperate Conifer Forests 
# 6 Boreal Forests/Taiga
# 7 Tropical & Subtropical Grasslands, Savannas & Shrublands
# 8 Temperate Grasslands, Savannas & Shrublands
# 9 Flooded Grasslands & Savannas
# 10 Montane Grasslands & Shrublands
# 11 Tundra
# 12 Mediterranean Forests, Woodlands & Scrub
# 13 Deserts & Xeric Shrublands
# 14 Mangroves



##### Make Figure 5 - mean metric values ##############################################################################################################################################

######## biomes ###############################################################
biome.coords.cont=readRDS("random.sample.values.biomes.continent.rds")
biome.all=readRDS("random.sample.values.biomes.global.rds")
biome_area=readRDS("random.sample.biome.continent.area.1000km2.rds")

to_remove=biome_area[which(biome_area$area<150),] # remove those with less than 150,000 km2 of land area.

# drop all the combinations in the "to_remove" df
biome.coords.cont3=biome.coords.cont
for(i in 1:nrow(to_remove)){
  biome.coords.cont3 <- biome.coords.cont3[!(biome.coords.cont3$continent == to_remove$continent[i] & biome.coords.cont3$Resolve_Biome == to_remove$Resolve_Biome[i]),] 
}

# add global
biome.all$continent="Global"
biome.coords.cont4=rbind(biome.coords.cont3,biome.all)

# get mean and SE
biome.coords.cont5 <- biome.coords.cont4 %>%
  group_by(continent,Resolve_Biome) %>%
  summarise(
    mean_value = mean(dark_taxa_hotspot_index, na.rm = TRUE),
    se_value = sd(dark_taxa_hotspot_index, na.rm = TRUE) / sqrt(n())
  )

labels=rev(c("Mangroves","Deserts","Mediterranean forests","Tundra","Montane grasslands","Flooded grasslands","Temperate grasslands","Tropical grasslands","Boreal forests","Temperate conifer forests","Temperate broadleaf forests","Tropical conifer forests","Tropical dry forests","Tropical moist forests"))
bio <- biome.coords.cont5 %>% mutate(labels = labels[Resolve_Biome]) # add labels
bio$continent=factor(bio$continent,ordered=T,levels=c("Global","Africa","Asia","Europe","North America","Oceania","South America"))

bio$labels <- gsub("Temperate", "Temp.", bio$labels)
bio$labels <- gsub("Tropical", "Trop.", bio$labels)
bio$labels <- gsub("broadleaf", "broad.", bio$labels)

# Select top 5 categories per continent and order within each group
df_top_bio <- bio %>%
  group_by(continent) %>%
  arrange(desc(mean_value)) %>%
  slice_head(n = 5) %>%
  mutate(labels = paste0(labels, " (", continent, ")")) %>%
  mutate(labels = factor(labels, levels = labels[order(mean_value)])) %>%
  mutate(y_indent = 0.03 * max(mean_value / 10)) %>%  # Small indent relative to max value
  ungroup()

# Create bar plot
plotb3=ggplot(df_top_bio, aes(x = labels, y = mean_value)) +
  geom_col(fill=rgb(255/255, 140/255, 0, 0.5)) + # darkorange with 50% opacity for fill
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value),  width = 0.2,size=0.5, color = "darkorange1") +  # Add error bars
  geom_text(aes(label = sub(" \\(.*\\)", "", labels), y = y_indent), color = "black", size = 2.6,hjust=0) +  # Labels inside bars
  coord_flip() +
  facet_wrap(~continent, scales = "free_y", nrow = 1) +
  scale_x_discrete(labels = function(x) NULL) +  # Remove axis labels
  theme_minimal() +
  theme(axis.text.y = element_blank(), legend.position = "none",plot.title = element_text(size=11)) +  # Remove text labels on y-axis
  labs(x = NULL, y = NULL,title="A) Biomes")

######## host families ###############################################################
host=readRDS("random.sample.values.families.continent.rds")
host.all=readRDS("random.sample.values.families.rds")
host.area=readRDS("random.sample.families.continent.area.1000km2.rds")

to_remove=host.area[which(host.area$sum <150),] # remove those with less than 150,000 km2 of land area.

# drop all the combinations in the "to_remove" df
host2=host
for(i in 1:nrow(to_remove)){
  host2 <- host2[!(host2$continent == to_remove$continent[i] & host2$Family == to_remove$Family[i]),] 
}

# add global
host.all$continent="Global"
host3=rbind(host2,host.all)

host4 <- host3 %>%
  group_by(continent,Family) %>%
  summarise(
    mean_value = mean(dark_taxa_hotspot_index, na.rm = TRUE),
    se_value = sd(dark_taxa_hotspot_index, na.rm = TRUE) / sqrt(n())
  )
host4$continent=factor(host4$continent,ordered=T,levels=c("Global","Africa","Asia","Europe","North America","Oceania","South America"))


# Select top 5 categories per continent and order within each group
df_top <- host4 %>%
  group_by(continent) %>%
  arrange(desc(mean_value)) %>%
  slice_head(n = 10) %>%
  mutate(Family = paste0(Family, " (", continent, ")")) %>%
  mutate(Family = factor(Family, levels = Family[order(mean_value)])) %>%
  mutate(y_indent = 0.03 * max(mean_value / 10)) %>%  # Small indent relative to max value
  ungroup()

# Create bar plot
ploth3=ggplot(df_top, aes(x = Family, y = mean_value)) +
  geom_col(fill=rgb(34/255, 139/255, 34/255, 0.5)) + # Forest green with 50% opacity for fill
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value),  width = 0.2, size=0.5, color = "forestgreen") +  # Add error bars
  geom_text(aes(label = sub(" \\(.*\\)", "", Family), y = y_indent), color = "black", size = 2.6,hjust=0) +  # Labels inside bars
  coord_flip() +
  facet_wrap(~continent, scales = "free_y", nrow = 1) +
  scale_x_discrete(labels = function(x) NULL) +  # Remove axis labels
  theme_minimal() +
  theme(axis.text.y = element_blank(),  legend.position = "none",plot.title = element_text(size=11)) + # Remove text labels on y-axis
  theme( strip.text.x = element_blank() ) + # remove facet labels
  labs(x = NULL, y = "Mean research priority metric score",title="B) Host plant families")

png("./figure_mean_metric_biome_host_panel_top.png",width=8.6,height=4,res=600,units="in")
pdf("./figure_mean_metric_biome_host_panel_top.pdf",width=8.6,height=4)
ggarrange(plotb3, ploth3, ncol = 1, heights = c(0.67, 1), align="v")
dev.off()

