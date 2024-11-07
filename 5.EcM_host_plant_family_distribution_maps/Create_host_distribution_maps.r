# code used to create host plant family rasters for van Galen et. al. (2025) "The biogeography and conservation of Earth’s ‘dark’ ectomycorrhizal fungi"

library(terra) # v1.7-78
library(rgbif) # v3.8.0
library(CoordinateCleaner) #v3.0.1
library(sf) # v1.0-17
library(raster) # v3.6-26
library(dplyr) # v1.1.4
library(nngeo) # v0.4.8 (for nearest distances)

#################################################################################################################################################################
## download gbif data
#################################################################################################################################################################
# https://docs.ropensci.org/rgbif/articles/gbif_credentials.html used this to set username and password of gbif account

###### to get ECM host genera ########################
# used supplementary file from Soudzilovskaia et. al. (2020) fungalRoot database containing mycorrhizal assignments for plant genera, and filtered those recorded as ECM https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.16569
# extract GBIF records for each of the genera
#######################################################

ecto.host.list=read.csv("fungalRoot_ECM.csv") # character vector of genus names from fungalRoot

failed.genera=NULL # things that fail because not in backbone (no usageKey found)
download.details=list()
citation.df=NULL
too.big=NULL # things that fail because GBIF tries to download all records due to a glitch (>50gb)

for(i in c(1:length(ecto.host.list))){
  print(paste(i, ecto.host.list[i]))
  setwd("./data_downloaded") # change directory to where the download files should be stored
  
  # Reset variables to avoid carry-over from previous iterations
  usage_key <- NULL
  num_records <- NULL
  occ_result <- NULL
  
  # Attempt to get the taxonKey
  usage_key <- name_backbone(ecto.host.list[i])$usageKey
  
  # Check if usage_key is NULL
  if (is.null(usage_key)) {
    print(paste("No usageKey found for", ecto.host.list[i]))
    failed.genera <- c(failed.genera, ecto.host.list[i])
    next # Skip remaining code, move to next iteration
  }
  
  # Estimate size based on the number of records
  occ_result <- occ_search(taxonKey = usage_key, limit = 0, facet = "none") # limit = 0 to get only the count
  num_records <- occ_result$meta$count
  
  if (num_records > 90000000) {  # sometimes rgbif tries to download 9 digits worth of records, which seems to be everything in GBIF. So skip if that's the case
    print(paste("Skipping", ecto.host.list[i], "due to estimated file size too large"))
    too.big <- c(too.big, ecto.host.list[i])
    next # Skip remaining code, move to next iteration
  }
  
  # Proceed to download the data
  gbif_download <- occ_download(pred("taxonKey", usage_key), format = "SIMPLE_CSV")
  download1 <- occ_download_wait(gbif_download)
  download1[length(download1) + 1] <- ecto.host.list[i] # Add genus name
  download.details[[i]] <- download1 # Save the full download details
  
  # Save citation details
  cite1 <- gbif_citation(download1$key) # Using the download key
  cite2 <- data.frame(Genus = ecto.host.list[i], Citation = cite1$download)
  citation.df <- rbind(citation.df, cite2)
  
  # Get the downloaded file path and import the data
  download_file <- occ_download_get(gbif_download)
  d <- occ_download_import(download_file) # Read the data from the file
  saveRDS(d, paste0("./data_raw/data_", gsub(" ", "_", ecto.host.list[i]), ".rds")) # Save as .rds
  
  # Delete the zip file to save space
  file.remove(download_file)
}

saveRDS(failed.genera,"./data_downloaded/failed.genera.rds")
saveRDS(download.details,"./data_downloaded/download.details.rds")
saveRDS(citation.df,"./data_downloaded/citation.df.rds")
saveRDS(too.big,"./data_downloaded/too.big.rds")

setwd("origingal_directory_file_path") # change wd back


## repeat for problem genera ###########################################################################################################################################################################
# some genera don't work if there is a duplicate with the same name (an animal or something). Therefore need to specify the taxonomic authority.
# went through this list of failed genera by hand and updated them after searching in GBIF https://www.gbif.org/

failed.genera=readRDS("./data_downloaded/failed.genera.rds") # character vector of the genera that failed because no usageKey could be found
too.big=readRDS("./data_downloaded/too.big.rds") # character vector of the genera that failed because rgbif tried to download the entirety of GBIF

failed.genera2=c("Acanthocladium F.Muell.","Afzelia Sm.","Agonis (DC.) Sweet","Aldina Endl.","Allotropa Torr. & A.Gray ex A.Gray","Aotus Sm.","Astus Trudgen & Rye","Beaufortia R.Br.","Betula L.",
                 "Brachysema R.Br.","Chrysolepis Hjelmq.","Corymbia K.D.Hill & L.A.S.Johnson","Cowania D.Don",     
                 "Dampiera R.Br.","Dryas L.","Eremaea Lindl.","Eutaxia R.Br.","Gastrolobium R.Br.","Haeckeria F.Muell.","Hudsonia L.","Jacksonia R.Br. ex Sm.",
                 "Jansonia Kippist","Kunzea Rchb.","Leptorhynchos Less.","Lithocarpus Blume","Marquesia Gilg",    
                 "Melaleuca L.","Mirbelia Sm.","Pinus L.","Platylobium Sm.","Pterochaeta Steetz","Pterospora Nutt.","Raoulia Hook.f.","Regelia S.Schauer","Salix L.","Seorsus Rye & Trudgen",
                 "Tristania R.Br.","Urodon Turcz.","Vatica L.","Verreauxia Benth.","Xylococcus Nutt.","Callistachys Vent.","Cedrus Trew","Chrysocoryne Endl.","Lechea Kalm ex L.","Poranthera Rudge","Waitzia J.C.Wendl.") # updated list of failed genera with their taxonomic authorities 
saveRDS(failed.genera2,"./data_downloaded/failed.genera.updated.names.rds") # save the list

# re-run download code for the failed genera
failed.genera2=readRDS("./data_downloaded/failed.genera.updated.names.rds")
download.details=readRDS("./data_downloaded/download.details.rds")
citation.df=readRDS("./data_downloaded/citation.df.rds")
for(i in c(1:length(failed.genera2))){
  print(paste(i, failed.genera2[i]))
  setwd("./data_downloaded") # change directory
  
  # Reset variables to avoid carry-over from previous iterations
  usage_key <- NULL
  num_records <- NULL
  occ_result <- NULL
  
  # Attempt to get the taxonKey
  usage_key <- name_backbone(failed.genera2[i])$usageKey
  
  # Proceed to download the data
  gbif_download <- occ_download(pred("taxonKey", usage_key), format = "SIMPLE_CSV")
  download1 <- occ_download_wait(gbif_download)
  download1[length(download1) + 1] <- failed.genera2[i] # Add genus name
  download.details[[length(download.details)+1]] <- download1 # Save the full download details (adding on to end)
  
  # Save citation details
  cite1 <- gbif_citation(download1$key) # Using the download key
  cite2 <- data.frame(Genus = failed.genera2[i], Citation = cite1$download)
  citation.df <- rbind(citation.df, cite2)
  
  # Get the downloaded file path and import the data
  download_file <- occ_download_get(gbif_download)
  d <- occ_download_import(download_file) # Read the data from the file
  saveRDS(d, paste0("./data_raw/data_", gsub(" ", "_", failed.genera2[i]), ".rds")) # Save as .rds
  
  # Delete the zip file to save space
  file.remove(download_file)
}

saveRDS(download.details,"./data_downloaded/download.details.rds")
saveRDS(citation.df,"./data_downloaded/citation.df.rds")
setwd("origingal_directory_file_path") # change wd back

#################################################################################################################################################################
## clean the gbif records
#################################################################################################################################################################
#loads a file that shows which parts of of the world are seas (available at https://github.com/nvkelso/natural-earth-vector/blob/master/50m_physical/ne_50m_land.shp)
seadata <- shapefile("./ne_50m_land/ne_50m_land.shp")

#This function removes NA values in lat and/or longitude,
#Inputs: a data frame with the values and the name of the column you want to remove NAs from 
completeFun <- function(data, desiredColsName) {
  colNum <-which(colnames(data)==desiredColsName)
  completeVec <- complete.cases(data[,colNum])
  return(data[completeVec, ])
}

########################################################################

citation.df=readRDS("./data_downloaded/citation.df.rds")
ecto.host.list=citation.df$Genus # using the updated names for those that have taxonomic authorities when needed

no.records.left=NULL # open a file to save names of genera that have no records left after cleaning
for(i in c(1:length(ecto.host.list))){
  print(i)
  data <-readRDS(paste0("./data_raw/data_",gsub(" ", "_", ecto.host.list[i]),".rds")) # occurrence records directly downloaded above
  
  #removes NAs from the read-in file
  data1 <- completeFun(data, "decimalLongitude")
  data1 <- completeFun(data1, "decimalLatitude")
  
  # removes duplicates from the file using the CoordinateCleaner cc_dupl function
  data2<-cc_dupl(data1, lon = "decimalLongitude", lat="decimalLatitude", species="genus") # removing anything with duplicate records. Using "genus" because don't care if multiple species of the same genus are present
  
  # remove pre-1970 records
  if(min(data2$year,na.rm=T)<1970){ # doesn't work if there are no rows to delete
    data3=data2[-which(data2$year<1970),] # this keeps NAs
  } else {
    data3=data2
  }
  
  # remove those with spatial uncertainty > 1km
  if(max(data3$coordinateUncertaintyInMeters,na.rm=T)>1000){ # doesn't work if there are no rows to delete
    data4=data3[-which(data3$coordinateUncertaintyInMeters>1000),] # this keeps NAs
  } else {
    data4=data3
  }
  
 if(nrow(data4)>0){
   # uses Coordinate Cleaner clean_coordinates function, which removes a number of suspicious coordinates--those from the capitals, centres of cities
   # located in the GBIF headquarters, etc--all the things that could have been likely to been put in by default and not where the species was.
   dataClean<-clean_coordinates(data4, lon = "decimalLongitude", lat="decimalLatitude", species="species", tests = c("capitals", "centroids", "equal", "gbif", "institutions", "seas", "zeros"), value="clean", seas_ref = seadata, verbose = F)
   if(nrow(dataClean)>0){
     saveRDS(dataClean,paste0("./data_clean/clean_",gsub(" ", "_", ecto.host.list[i]),".rds")) # saves as .rds
   } else {no.records.left=c(no.records.left,ecto.host.list[i])}
 } else {no.records.left=c(no.records.left,ecto.host.list[i])}
}

saveRDS(no.records.left,"./data_clean/no.records.left.rds")

# # have a look at them
# for(i in c(1:length(ecto.host.list[1:2]))){
#   assign(paste0("clean_",gsub(" ", "_", ecto.host.list[i])), readRDS(paste0("./data_clean/clean_",gsub(" ", "_", ecto.host.list[i]),".rds")))
# 
# }

#################################################################################################################################################################
## reproject points (using equal earth)
#################################################################################################################################################################
##### transform #####
citation.df=readRDS("./data_downloaded/citation.df.rds")
no.records.left=readRDS("./data_clean/no.records.left.rds")
ecto.host.list=citation.df$Genus # using just the ones that worked
ecto.host.list <- setdiff(ecto.host.list, no.records.left) # remove those that were lost last time

template=rast("original_maps/ectomycorrhizal_richness_Classified_MultibandImage_final_EE.tif")["ectomycorrhizal_richness_Ensemble_mean"] # use a template raster with equal-earth projection

for(i in c(1:length(ecto.host.list))){
  print(i)
  data <- readRDS(paste0("./data_clean/clean_",gsub(" ", "_", ecto.host.list[i]),".rds")) # occurrence records
  
  # transform coordinates
  data2 <- data %>% st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs =  4326) # change to sf object
  data3 <- data2 %>% st_transform(crs = crs(template)) #transform to the same projection as the template layer
  data_reprojected <- data.frame(st_coordinates(data3)) # extract coordinates
  data4=cbind(data,data_reprojected) # add x y coords to end
  saveRDS(data4,paste0("./data_repro/clean_",gsub(" ", "_", ecto.host.list[i]),"_reprojected.rds")) # saves as .rds
}

#################################################################################################################################################################
## remove duplicates (per genera)
#################################################################################################################################################################
citation.df=readRDS("./data_downloaded/citation.df.rds")
no.records.left=readRDS("./data_clean/no.records.left.rds")
ecto.host.list=citation.df$Genus # using just the ones that worked
ecto.host.list <- setdiff(ecto.host.list, no.records.left) # remove those that were lost last time

template=rast("original_maps/ectomycorrhizal_richness_Classified_MultibandImage_final_EE.tif")["ectomycorrhizal_richness_Ensemble_mean"] # use a template raster with equal-earth projection and resolution of 1km2

for(i in c(1:length(ecto.host.list))){
  print(i)
  data <-readRDS(paste0("./data_repro/clean_",gsub(" ", "_", ecto.host.list[i]),"_reprojected.rds")) # occurrence records
  
  ####### remove cell duplicates ########
  occurrence.cells <- terra::extract(template, data[c("X","Y")], cells = TRUE)
  occurrence.cellDups <- duplicated(occurrence.cells$cell)
  occurrence_noDup <- data[!occurrence.cellDups,]
  saveRDS(occurrence_noDup,paste0("./data_noDup/clean_",gsub(" ", "_", ecto.host.list[i]),"_reprojected_noDup.rds"))
}

#################################################################################################################################################################
## merge families
#################################################################################################################################################################
citation.df=readRDS("./data_downloaded/citation.df.rds")
no.records.left=readRDS("./data_clean/no.records.left.rds")
ecto.host.list=citation.df$Genus # using just the ones that worked
ecto.host.list <- setdiff(ecto.host.list, no.records.left) # remove those that were lost last time

# merge the dataframes
data.all=NULL
for(i in c(1:length(ecto.host.list))){
  if(i%%10==0){print(i)}
  data <-readRDS(paste0("./data_noDup/clean_",gsub(" ", "_", ecto.host.list[i]),"_reprojected_noDup.rds")) # occurrence records
  data.all=rbind(data.all,data)
}

# save list of families
families=unique(data.all$family)
saveRDS(families,paste0("./data_families/family.list.rds"))
family.count=data.frame(table(data.all$family))
names(family.count)<-c("family","num_records")
genera_count <- data.all %>%
  group_by(family) %>%
  summarize(num_genera = n_distinct(genus), .groups = 'drop')
family.summary=merge(family.count,genera_count,by="family")
write.csv(family.summary,paste0("./data_families/family.summary.csv"))

# split and save each family
split_data <- split(data.all, data.all$family) # split each family into a list
for (family_name in names(split_data)) {
  print(family_name)
  saveRDS(split_data[[family_name]], paste0("./data_families/data_final_", family_name, ".rds"))
}

#################################################################################################################################################################
## remove duplicates again (per family)
#################################################################################################################################################################
family.list=readRDS("./data_families/family.list.rds")
template=rast("original_maps/ectomycorrhizal_richness_Classified_MultibandImage_final_EE.tif")["ectomycorrhizal_richness_Ensemble_mean"] # use a template raster with equal-earth projection and resolution of 1km2

for(i in c(1:length(family.list))){
  print(family.list[i])
  data <-readRDS(paste0("./data_families/data_final_", family.list[i], ".rds")) # occurrence records
  
  ####### remove cell duplicates ########
  occurrence.cells <- terra::extract(template, data[c("X","Y")], cells = TRUE)
  occurrence.cellDups <- duplicated(occurrence.cells$cell)
  occurrence_noDup <- data[!occurrence.cellDups,]
  saveRDS(occurrence_noDup,paste0("./data_families/data_final_", family.list[i], "_noDup.rds"))
}

#################################################################################################################################################################
## calculate minimum distance between points
#################################################################################################################################################################
family.list=readRDS("./data_families/family.list.rds")
template=rast("original_maps/ectomycorrhizal_richness_Classified_MultibandImage_final_EE.tif")["ectomycorrhizal_richness_Ensemble_mean"]

for(i in c(1:length(family.list))){
  print(family.list[i])
  print(Sys.time())
  coords <-readRDS(paste0("./data_families/data_final_",gsub(" ", "_", family.list[i]),"_noDup.rds"))[c("X","Y")] # occurrence records
  sf_coords <- st_as_sf(coords, coords = c("X", "Y"), crs = crs(template)) # Convert the dataframe to an sf object with UTM CRS

  nearest_neighbors <- st_nn(sf_coords, sf_coords, k = 2, returnDist = TRUE) # Find the nearest neighbor for each point using spatial indexing, k = 2 to skip the point itself and find the nearest other point
  nearest_distances <- sapply(nearest_neighbors$dist, function(x) x[2]) / 1000 #  Extract the distances (in meters) and convert to kilometers
  
  saveRDS(nearest_distances,paste0("./data_families/min_distance_data/min_distance_",gsub(" ", "_", family.list[i]),"",".rds"))
}

#################################################################################################################################################################
## spatRaster layers of coordinates with buffer
#################################################################################################################################################################
family.list=readRDS("./data_families/family.list.rds")
template=rast("original_maps/ectomycorrhizal_richness_Classified_MultibandImage_final_EE.tif")["ectomycorrhizal_richness_Ensemble_mean"]

# ran this three times, with 50000, 100000, and 150000 buffer sizes
for(i in c(1:length(family.list))){
  print(family.list[i])
  for(radius in c(50000, 100000, 150000)){ # units m
    for(dist.thresh in c(Inf, 500, 200)){ # units km
      print(paste(radius,dist.thresh))
      print(Sys.time())
      coords <-readRDS(paste0("./data_families/data_final_",gsub(" ", "_", family.list[i]),"_noDup.rds"))[c("X","Y")] # occurrence records
      sf_coords <- st_as_sf(coords, coords = c("X", "Y"), crs = crs(template)) # Convert the dataframe to an sf object with UTM CRS
      min_distances <- readRDS(paste0("./data_families/min_distance_data/min_distance_",gsub(" ", "_", family.list[i]),"",".rds")) # Units are in km
      
      # delete points where the closest point is > xxx km away
      dist.thresh=dist.thresh # units km
      sf_coords <- sf_coords[min_distances <= dist.thresh, ]
      
      # define buffers
      radius <- radius # Define the buffer radius in meters (100km = 100000 meters)
      buffers <- st_buffer(sf_coords, dist = radius) # Create buffers (circles) around each point
      buffers <- st_transform(buffers, crs = crs(template)) # Ensure buffers are in the same CRS as the template raster
      
      # Rasterize the buffer zones using the template raster parameters
      raster_circles <- rasterize(buffers, template, field = 1, background = NA)
      raster_circles2=mask(raster_circles,template) # clip to land area
      writeRaster(raster_circles2, filename=paste0("./data_families/family_rasters/raster_",gsub(" ", "_", family.list[i]),"_thresh",dist.thresh,"km_buffer",radius/1000,"km.tif"),overwrite=TRUE) # save as tif
    }
  }
}

# make plots for each family to compare parameters
for(i in c(1:length(family.list))){
  print(family.list[i])
  png(paste0("./data_families/figure_points_comparison_",family.list[i],".png"),width=14,height=8,res=300,units="in")
  par(mfrow=c(3,3))
  for(radius in c(50000, 100000, 150000)){
    for(dist.thresh in c(Inf, 500, 200)){
      data=rast(paste0("./data_families/family_rasters/raster_",gsub(" ", "_", family.list[i]),"_thresh",dist.thresh,"km_buffer",radius/1000,"km.tif"))
      plot(template,legend=F,main=paste0("Family = ",family.list[i],", buffer = ",radius/1000,", thresh = ",dist.thresh))
      plot(data,add=T,legend=F)
    }
  }
  dev.off()
}

#################################################################################################################################################################
## calculate family statistics for supplementary table, and prepare file of GBIF download citations
#################################################################################################################################################################

####### GBIF citation file ########################################################################
citation.df=readRDS("./data_downloaded/citation.df.rds")
ecto.host.list=citation.df$Genus # using just the ones that worked
ecto.host.list.end.removed <- sub(" .*", "", ecto.host.list) # remove taxonomic authorities from those that have them

citation.df$family=""
citation.df$num_records=NA

for(i in c(1:length(ecto.host.list))){
  if(i%%10==0){print(i)}
  orig_data <- readRDS(paste0("./data_raw/data_",gsub(" ", "_", ecto.host.list[i]),".rds")) # occurrence records
  if(nrow(orig_data)>0){citation.df$family[i]=unique(orig_data$family)}# add family if there are any rows
  
  file_list <- list.files("./data_noDup") # those that have cleaned records
  if (any(grepl(ecto.host.list.end.removed[i], file_list))) { # only do this if clean records exist
    data <-readRDS(paste0("./data_noDup/clean_",gsub(" ", "_", ecto.host.list[i]),"_reprojected_noDup.rds")) # occurrence records
    citation.df$num_records[i]=nrow(data) # add number of records
  } else {
    citation.df$num_records=ifelse(citation.df$Genus==ecto.host.list[i],0, citation.df$num_records) # add number of records
  }
}

citation.df2=citation.df[-which(citation.df$family==""),] # drop rows that don't have a family, i.e. that there was no data downloaded for
citation.df3=data.frame(cbind(Family=citation.df2$family,Genus=citation.df2$Genus,Clean_records=citation.df2$num_records,Citation=citation.df2$Citation)) # reorder and rename
citation.df4=citation.df3[-which(citation.df3$Family=="Cyperaceae"|citation.df3$Family=="Apiaceae"),] # remove unwanted families
write.csv(citation.df4,"./file_genus_family_citation.csv",row.names=F)

####### Supplementary table ########################################################################
family.list=readRDS("./data_families/family.list.rds")
template=rast("original_maps/ectomycorrhizal_richness_Classified_MultibandImage_final_EE.tif")["ectomycorrhizal_richness_Ensemble_mean"]
s1=read.csv("./file_genus_family_citation.csv") # for genera

family.list2=family.list[ ! family.list %in% c("Cyperaceae","Apiaceae") ] # drop families that we are no longer using (see manuscript Methods)

output=NULL
for(i in c(1:length(family.list2))){
  print(family.list2[i])
  coords <-readRDS(paste0("./data_families/data_final_",gsub(" ", "_", family.list2[i]),"_noDup.rds")) # occurrence records
  sf_coords <- st_as_sf(coords, coords = c("X", "Y"), crs = crs(template)) # Convert the dataframe to an sf object with UTM CRS
  min_distances <- readRDS(paste0("./data_families/min_distance_data/min_distance_",gsub(" ", "_", family.list2[i]),"",".rds")) # Units are in km
  
  # delete points where the closest point is > xxx km away
  dist.thresh=200 # units km
  sf_coords <- sf_coords[min_distances <= dist.thresh, ]
  
  # genera (can't use genus col in coords because some are lost when points are deleted)
  s2=s1[which(s1$Family==family.list2[i]),]
  
  result=data.frame(Family=family.list2[i],Num_genera=nrow(s2),Cleaned_records=nrow(sf_coords))
  output=rbind(result,output)
}
output2 <- output[order(output$Family),]
write.csv(output2,"./family_table_supps.csv",row.names=F)

#################################################################################################################################################################
## transform to Robinson projection for use in analysis 
#################################################################################################################################################################
family.list=readRDS("./data_families/family.list.rds")
template=rast("original_maps/AM_richness_robinson.tif") # template raster for Robinson (units in km)
crs(template) <-  gsub("LENGTHUNIT\\[\"kilometre\",1000", "LENGTHUNIT[\"metre\",1", crs(template)) # update units to meters
ext(template) <- ext(template) * 1000

# just doing 200 and 150, because that's what we decided to use
for(i in c(2:length(family.list))){
  print(family.list[i])
  for(radius in c(150000)){ 
    for(dist.thresh in c(200)){
      print(paste(radius,dist.thresh))
      print(Sys.time())
      data=rast(paste0("./data_families/family_rasters/raster_",gsub(" ", "_", family.list[i]),"_thresh",dist.thresh,"km_buffer",radius/1000,"km.tif"))
      data_repro=terra::project(data,template,method="near") # "near", categorical
      data_repro=subst(data_repro, from=0, to=NA)
      writeRaster(data_repro, filename=paste0("./data_families/family_rasters/raster_",gsub(" ", "_", family.list[i]),"_thresh",dist.thresh,"km_buffer",radius/1000,"km_robinson.tif"),overwrite=TRUE) # save as tif
    }
  }
}

#################################################################################################################################################################
## make plot of host distribution maps for Supplementary Figures
#################################################################################################################################################################
family.list=readRDS("./data_families/family.list.rds")
family.list2=sort(family.list)
family.list2=family.list2[ ! family.list2 %in% c("Cyperaceae","Apiaceae") ] # drop families that we are no longer using

# # make Robinson background raster
# template=rast("original_maps/AM_richness_robinson.tif") # units in km
# crs(template) <-  gsub("LENGTHUNIT\\[\"kilometre\",1000", "LENGTHUNIT[\"metre\",1", crs(template)) # update units to meters
# ext(template) <- ext(template) * 1000
# template_10=aggregate(template,fact=10) # doesn't matter about NAs, purely just for the grid cell size
# background=rast("./Background/Background_equirectangular_10km.tif") # raster of land area to use as map background
# background_rob=terra::project(background,template_10,method="near")
# background_rob=subst(background_rob, from=NA, to=0)
# plot(background_rob)
# writeRaster(background_rob, filename=paste0("./Background/Background_robinson_10km.tif"),overwrite=TRUE) # save as tif

background=rast("./Background/Background_robinson_10km.tif")
coltab(background) <- data.frame(value=c(0,1),col=c("NA","grey85"))
er <- rast(ext(-17500000,17500000,-6100000,8565000), resolution=res(background),crs=crs(background)) # empty raster for outline, set extent manually to exclude Antarctica
values(er) <- 0

graticule <- st_graticule(lat = seq(-60, 90, 30), lon = seq(-180, 180, 30))
graticule <- st_transform(graticule, crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

png("figure_host_supps.png",width=8,height=7.6,res=800,units="in")
par(mfrow=c(7,4),mar=c(0,0,0,0),oma=c(0,0,0,0))

cex1=0.8
cex2 <- 1
lwd1=0.5
mar1=c(0,0.1,1.5,0.1)

for(i in c(1:length(family.list2))){
  plot(er,legend=F,col=NA,mar=mar1,box=F,pax=list(col=NA,col.axis="grey30"),axes=F) # plot empty raster for setting the margins
  plot(graticule$geometry,add=T, col = "grey80", lty = 1,lwd=lwd1) # add a grid
  plot(background,add=T,legend=F,axes=F) # col NA makes it transparent

  data=rast(paste0("./data_families/family_rasters/raster_",gsub(" ", "_", family.list2[i]),"_thresh200km_buffer150km_robinson.tif"))
  plot(data,add=T,col="forestgreen",legend=F,)
  mtext(text=family.list2[[i]],line=-1.5,adj=0,cex=0.8)
}
dev.off()


# Some host families have tiny distributions. Calculate land area of those:
Asteropeiaceae=rast(paste0("./data_families/family_rasters/raster_Asteropeiaceae_thresh200km_buffer150km.tif")) # EE projection
Asteropeiaceae_freq=freq(Asteropeiaceae) # units km
Asteropeiaceae_freq

Sarcolaenaceae=rast(paste0("./data_families/family_rasters/raster_Sarcolaenaceae_thresh200km_buffer150km.tif")) # EE projection
Sarcolaenaceae_freq=freq(Sarcolaenaceae) # units km
Sarcolaenaceae_freq

Ticodendraceae=rast(paste0("./data_families/family_rasters/raster_Ticodendraceae_thresh200km_buffer150km.tif")) # EE projection
Ticodendraceae_freq=freq(Ticodendraceae) # units km
Ticodendraceae_freq

# all cover less than 600,000 km2

#################################################################################################################################################################
## transform back to equirectangular for use in analysis
#################################################################################################################################################################
family.list=readRDS("./data_families/family.list.rds")
template=rast("original_maps/ectomycorrhizal_richness_Classified_MultibandImage_final.tif")["ectomycorrhizal_richness_Ensemble_mean"] # equirectangular template

for(i in c(2:length(family.list))){
  print(family.list[i])
  for(radius in c(100000, 150000)){ # not doing the 50km buffer because decided not to use it
    for(dist.thresh in c(Inf, 500, 200)){
      print(paste(radius,dist.thresh))
      print(Sys.time())
      data=rast(paste0("./data_families/family_rasters/raster_",gsub(" ", "_", family.list[i]),"_thresh",dist.thresh,"km_buffer",radius/1000,"km.tif"))
      data_repro=terra::project(data,template)
      writeRaster(data_repro, filename=paste0("./data_families/family_rasters/raster_",gsub(" ", "_", family.list[i]),"_thresh",dist.thresh,"km_buffer",radius/1000,"km_equiRec.tif"),overwrite=TRUE) # save as tif
    }
  }
}

#################################################################################################################################################################
## make host diversity layer
#################################################################################################################################################################
family.list=readRDS("./data_families/family.list.rds")
family.list=family.list[ ! family.list %in% c("Cyperaceae","Apiaceae") ] # drop families that we are no longer using

raster_list = list()
for (i in 1:length(family.list)) {
  raster_path = paste0("./data_families/family_rasters/raster_",gsub(" ", "_", family.list[i]),"_thresh200km_buffer150km_equirec.tif")
  raster_list[[i]] = rast(raster_path)
  names(raster_list[[i]])<-family.list[i]
}
family_stack = do.call(c, raster_list)

family_sum=sum(family_stack,na.rm=T)
names(family_sum)<-"host_family_richness"
plot(family_sum)
writeRaster(family_sum, filename=paste0("./data_families/family_rasters/raster_host_family_richness_equirectangular.tiff"),overwrite=TRUE) # save as tif
