##### Script to make the list of Australian ECM species names based on NSL + GBIF + GlobalFungi + GSMc for van Galen et. al. (2025) "The biogeography and conservation of Earth’s ‘dark’ ectomycorrhizal fungi"

library(rgbif) # v3.8.0
library(dplyr) # v1.1.4
library(ggplot2) # v3.5.1
library(stringr) # v1.5.1
library(VennDiagram) # v1.7.3
library(eulerr) # v7.0.2
library(seqinr) # v4.2-36

############################################################################################################################################
# Australian ECM species in the National Species List (NSL)
############################################################################################################################################
# download NSL data
nsl = data.frame(
  read.csv("https://fungi.biodiversity.org.au/nsl/services/export/taxonCsv",encoding = "UTF-8"))

# Select species names (subspecies/var not taken into account here)
unique(nsl$taxonRank)
nsl_sp = subset(nsl, taxonRank %in% c("Species"))
unique(nsl_sp$taxonRank)

# Filter out names as "sp."
nsl_sp = nsl_sp[!grepl("sp. ", nsl_sp$scientificName),]
nsl_sp = with(nsl_sp,  nsl_sp[order(scientificName) , ])   # reorder alphabetically (scientific name)

# Create new column with genus name
nsl_sp$genus = word(nsl_sp$acceptedNameUsage, 1)

# Upload FungalTraits to filter out EcM species
FungalTraits = data.frame(read.csv("13225_2020_466_MOESM4_ESM.csv",encoding = "UTF-8")) # FungalTraits database from https://link.springer.com/article/10.1007/s13225-020-00466-2
colnames(FungalTraits)[colnames(FungalTraits) == 'GENUS'] <- 'genus' 

# Append guild info to NSL species list
nsl_guild = left_join(
  x = nsl_sp,
  y = FungalTraits,
  by = "genus",
  copy = FALSE,
  suffix = c(".x", ".y"),
  keep = NULL
)

# Select ECM entries
nsl_ecm = nsl_guild[nsl_guild$primary_lifestyle %in% "ectomycorrhizal",]

nsl_ecm = nsl_ecm %>% select(scientificName, canonicalName, acceptedNameUsage,taxonomicStatus)
nsl_ecm$NSL = c("yes")
colnames(nsl_ecm)[colnames(nsl_ecm) == 'scientificName'] <- 'scientificName_NSL'
View(nsl_ecm)

# Select accepted ECM taxa = 1270
nsl_ecm_valid = nsl_ecm[nsl_ecm$taxonomicStatus %in% "accepted",]
View(nsl_ecm_valid)

############################################################################################################################################
# Australian ECM species in GBIF
############################################################################################################################################
?download_predicate_dsl   # see list of search fields
name_backbone("Fungi")   # get taxonKey for species or group of interest

# Accessed on 3 October 2024
gbif_download = occ_download(
    pred("country", "AU"),
    pred("HAS_GEOSPATIAL_ISSUE",FALSE),
    pred("taxonKey",5),
    pred_not(pred_in("BASIS_OF_RECORD",
                     c("MATERIAL_SAMPLE"))),
    format = "SIMPLE_CSV",
    user = Sys.getenv("gbif_user"),
    pwd = Sys.getenv("gbif_pwd"),
    email = Sys.getenv("gbif_email")
    )

# Check download status
occ_download_wait(gbif_download) 

# citation details 
my_download<-occ_download_get("0038895-240906103802322",overwrite=TRUE)
gbif_citation(my_download)

# Upload table
gbif = data.frame(
        occ_download_get(gbif_download) %>%   
        occ_download_import(),
        row.names = 1,
        encoding = "UTF-8")   # convert integer64 to numericals

# Select species rank
unique(gbif$taxonRank)
gbif_sp = subset(gbif, taxonRank %in% "SPECIES")
unique(gbif_sp$taxonRank)

# Upload FungalTraits to filter out EcM species
FungalTraits = data.frame(read.csv("13225_2020_466_MOESM4_ESM.csv",encoding = "UTF-8")) # FungalTraits database from https://link.springer.com/article/10.1007/s13225-020-00466-2
colnames(FungalTraits)[colnames(FungalTraits) == 'GENUS'] <- 'genus' 

# Append guild info to  species list
gbif_guild = left_join(
  x = gbif_sp,
  y = FungalTraits,
  by = "genus",
  copy = FALSE,
  suffix = c(".x", ".y"),
  keep = NULL)

# Select ECM entries
gbif_ecm = gbif_guild[gbif_guild$primary_lifestyle %in% "ectomycorrhizal",]
View(gbif_ecm)

# Count the number of GBIF records in each category
unique(gbif_ecm$basisOfRecord)

gbif_ecm_record = gbif_ecm %>% 
  count(species, basisOfRecord) %>%
  tidyr::spread(key = basisOfRecord, value = n)

# Merge categories by observation and specimen: 1767 ECM species in gbif
gbif_ecm_record = gbif_ecm_record %>%
  mutate(OBSERVATION_TOT = rowSums(across(c(HUMAN_OBSERVATION,MACHINE_OBSERVATION,OBSERVATION, OCCURRENCE)),
                                   na.rm=TRUE)) %>%
  mutate(SPECIMEN_TOT = rowSums(across(c(LIVING_SPECIMEN, PRESERVED_SPECIMEN)),
                                   na.rm=TRUE))

gbif_ecm_record = gbif_ecm_record %>% select(species, OBSERVATION_TOT, SPECIMEN_TOT)
colnames(gbif_ecm_record)[colnames(gbif_ecm_record) == 'species'] <- 'canonicalName'
gbif_ecm_record$GBIF = c("yes")
View(gbif_ecm_record)

gbif_ecm_record[, c(2:3)] <- sapply(gbif_ecm_record[, c(2:3)], as.numeric)
gbif_ecm_record[is.na(gbif_ecm_record)]<-0
colSums(Filter(is.numeric, gbif_ecm_record))

############################################################################################################################################
# Australian EcM species in GlobalFungi version 5
############################################################################################################################################
# GlobalFungi version 5 EcM OTUs were obtained from the authors of the GlobalFungi database (https://globalfungi.com/)

gf_new=readRDS("GFv5_ECM_OTU_TAXA_METADAT.rds") # Metadata of GlobalFungi version 5
gf_new_aus=gf_new[which(gf_new$country=="Australia"&gf_new$continent=="Australia"),] # filter to keep samples from Australia
gf_new_aus$SIMILARITY=as.numeric(gf_new_aus$SIMILARITY)

# Only including species names if similarity (match between OTU representative sequence and taxonomy reference database) is > 97% 
gf_new_aus2 <- gf_new_aus %>%
  mutate(Species = ifelse(SIMILARITY > 97,Species,paste0(sub("_.+", "", Species), "_sp")))

df_new_aus_named=gf_new_aus2[!grepl("_sp$", gf_new_aus2$Species), ] # drop unnamed
df_new_aus_named_unique=data.frame(canonicalName=gsub("_", " ", unique(df_new_aus_named$Species)),GlobalFungi="yes") # unique names

############################################################################################################################################
# Australian EcM species in GSMc
############################################################################################################################################
# GSMc data downloaded from here: https://doi.plutof.ut.ee/doi/10.15156/BIO/2263453
meta=read.csv("Tedersoo L, Mikryukov V, Anslan S et al. Fungi_GSMc_sample_metadata.csv") # note: manually copied txt file into excel and save as .csv because otherwise was an issue with read.table
meta_aus <- meta[grepl("Australia", meta$area), ] # select Australian samples
meta_aus=meta_aus[-which(meta_aus$plot=="G4905"),] # drop outlier from Cuba

# #subset OTU to Aus samples (takes a long time to read in otu table). Have saved output
# otu=read.table("Fungi_GSMc_OTU_Table.txt")
# colnames(otu)<-otu[1,]
# otu_aus <- otu[, names(otu)[(names(otu) %in% c("OTU",meta_aus$plot))]] # extract Aus plots
# otu_aus=otu_aus[-1,] # drop first row
# rownames(otu_aus)=otu_aus$OTU
# otu_aus=otu_aus[,-1] # drop first col
# otu_aus2 <- data.frame(sapply( otu_aus, as.numeric )) # convert cols to numeric
# rownames(otu_aus2)<-rownames(otu_aus)
# otu_aus3=otu_aus2[rowSums(otu_aus2[])>0,] # keep only rows with non zero values
# otu_aus3$OTU=rownames(otu_aus3)
# write.csv(otu_aus3,"GSMc_Aus_subset.csv",row.names=F)

otu_aus=read.csv("GSMc_Aus_subset.csv")
tax=read.csv("Tedersoo L, Mikryukov V, Anslan S et al. Fungi_GSMc_taxonomy-function table_final.csv")
tax_aus=tax[tax$OTU %in% otu_aus$OTU,] # extract Aus
tax_aus_named=tax_aus[-which(tax_aus$species=="."),] # those named to species level
tax_aus_named_ECM = tax_aus_named[grepl("ectomycorrhizal", tax_aus_named$primary_lifestyle), ] # ECM
tax_aus_named_ECM$canonicalName = paste(tax_aus_named_ECM$genus,tax_aus_named_ECM$species)
GSMc=data.frame(canonicalName=unique(tax_aus_named_ECM$canonicalName),GSMc="yes")

############################################################################################################################################
# merge data and remove synonyms
############################################################################################################################################
nsl_gbif <- merge(nsl_ecm,gbif_ecm_record,by=c("canonicalName"), all.x=TRUE, all.y=TRUE) # merge NSL and GBIF
nsl_gbif_GF <- merge(nsl_gbif,df_new_aus_named_unique,by=c("canonicalName"), all=TRUE) # add GlobalFungi
nsl_gbif_GF_GSMc <- merge(nsl_gbif_GF,GSMc,by=c("canonicalName"), all=TRUE) # add GSMc

##### Remove synonyms ####################
# acceptedNameUsage column has NAs, so fill these with the canonicalName values
nsl_gbif_GF_GSMc$acceptedNameUsage <- ifelse(is.na(nsl_gbif_GF_GSMc$acceptedNameUsage), nsl_gbif_GF_GSMc$canonicalName, nsl_gbif_GF_GSMc$acceptedNameUsage)

df_summary <- nsl_gbif_GF_GSMc %>%
  group_by(acceptedNameUsage) %>% 
  summarise(
    NSL = ifelse(any(NSL == "yes"), "yes", NA), # if yes occurs in any row, keep yes, otherwise NA
    GBIF = ifelse(any(GBIF == "yes"), "yes", NA), 
    OBSERVATION_TOT = sum(OBSERVATION_TOT, na.rm = TRUE), # sum the values of the rows, excluding NAs
    SPECIMEN_TOT  = sum(SPECIMEN_TOT , na.rm = TRUE), 
    GlobalFungi = ifelse(any(GlobalFungi == "yes"), "yes", NA), # if yes occurs in any row, keep yes, otherwise NA
    GSMc = ifelse(any(GSMc == "yes"), "yes", NA), # if yes occurs in any row, keep yes, otherwise NA
    synonym_canonicalNames = paste(unique(canonicalName,), collapse = ", ") # so we can double check it's worked, keep a list of the synonym names
  )
df_summary$canonicalName = word(df_summary$acceptedNameUsage, 1,2)   # create a variable with species name without authority 

# write.csv(df_summary,"ECM_species_list_AUS_noSynonyms.csv",row.names=F,fileEncoding = "UTF-8")

############################################################################################################################################
# merge with species names in UNITE version 10
############################################################################################################################################
# downloaded from the UNITE database https://unite.ut.ee/repository.php
fasta_sequences <- read.fasta(file = "sh_general_release_dynamic_s_04.04.2024_dev.fasta", seqtype = "DNA")
headers <- sapply(getAnnot(fasta_sequences), function(x) x)
split_taxonomy <- strsplit(headers, "\\|")# Split the taxonomy information
genus_species <- sapply(split_taxonomy, function(x) head(x, 1)) # Extract the last two elements (genus and species)
unite <- as.data.frame((genus_species))# Convert to a data frame for easy viewing
colnames(unite) <- c("canonicalName")
unite$canonicalName <- gsub(">", "", unite$canonicalName) # remove >
unite <- data.frame(canonicalName=unite[!grepl("_sp$", unite$canonicalName), ]) # remove _sp
unite$canonicalName <- gsub("_", " ", unite$canonicalName)
unite=data.frame(canonicalName=unique(unite$canonicalName))
unite$UNITE = c("yes")

# merge with dataframe
df_summary2=read.csv("ECM_species_list_AUS_noSynonyms.csv")
data <- merge(df_summary2,unite,by=c("canonicalName"), all.x=TRUE) # keep only those from unite that were detected in the other Aus records

# write.csv(data,"ECM_final_species_list_w_UNITE.csv",row.names=F,fileEncoding = "UTF-8")

############################################################################################################################################
# make euler diagrams
############################################################################################################################################

####### Main text figure  #################################################################################################
data=read.csv("ECM_final_species_list_w_UNITE.csv")
data$NSL_GBIF=ifelse(data$NSL=="yes"|data$GBIF=="yes","yes",NA) # merge NSL and GBIF
data$GlobalFungi_GSMc=ifelse(data$GlobalFungi=="yes"|data$GSMc=="yes","yes",NA) # merge GlobalFungi and GSMc

# create sets
set_NSL_GBIF <- data[which(data$NSL_GBIF=="yes"),"canonicalName"]
set_UNITE <- data[which(data$UNITE=="yes"),"canonicalName"]
set_GlobalFungi_GSMc <- data[which(data$GlobalFungi_GSMc=="yes"),"canonicalName"]

# Euler diagram
fit <- euler(list(NSL_GBIF=set_NSL_GBIF, GlobalFungi_GSMc=set_GlobalFungi_GSMc, UNITE=set_UNITE))
plot=plot(fit, 
          quantities = list(type = c("counts", "percent"), font=1, round=2, cex=0.75),
          shape = "ellipse",
          labels=list(font=2, cex=1))
plot

pdf("./figure_euler_mainText.pdf",width=5,height=5)
plot
dev.off()

# # use this to work out what's missing from the euler
# venn.diagram(
#   x = list(set_NSL_GBIF, set_GlobalFungi_GSMc, set_UNITE),
#   category.names = c("Physical","GlobalFungi","UNITE"),
#   print.mode=c("raw","percent"),
#   filename = 'venn_diagramm.png',
#   output=TRUE
# )


####### Make figure with split categories for supps #################################################################################################
data=read.csv("ECM_final_species_list_w_UNITE.csv")

# create sets
set_NSL <- data[which(data$NSL=="yes"),"canonicalName"]
set_GBIF <- data[which(data$GBIF=="yes"),"canonicalName"]
set_GlobalFungi <- data[which(data$GlobalFungi=="yes"),"canonicalName"]
set_GSMc <- data[which(data$GSMc=="yes"),"canonicalName"]
set_UNITE <- data[which(data$UNITE=="yes"),"canonicalName"]


# Euler diagram
fit <- euler(list(NSL=set_NSL, GBIF=set_GBIF, GlobalFungi=set_GlobalFungi, GSMc=set_GSMc, UNITE=set_UNITE))
plot_supps=plot(fit,
     quantities = list(type = c("counts", "percent"), font=1, round=2, cex=0.75),
     shape = "ellipse",
     labels=list(font=2, cex=1))
plot_supps

pdf("./figure_euler_supps.pdf",width=5,height=5)
plot_supps
dev.off()

# # use this to check what's missing from euler
# venn.diagram(
#   x = list(set_NSL, set_GBIF, set_GlobalFungi, set_GSMc,set_UNITE),
#   category.names = c("Checklist" , "GBIF","GlobalFungi", "GSMc","UNITE"),
#   print.mode=c("raw"),
#   filename = 'venn_diagramm_supps.png',
#   output=TRUE
# )


