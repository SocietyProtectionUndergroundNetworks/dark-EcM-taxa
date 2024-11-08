## Data and code used to conduct the analyses and create the figures for van Galen et. al. (2025) “The biogeography and conservation of Earth’s ‘dark’ ectomycorrhizal fungi”

### 1.Proportion_of_EcM_fungal_OTUs_unassigned (will become availabe shortly, in time for publication)
Data used to create Figure 1.

Files:
-	“GFv5_EcM_unassigned_per_sample.csv”: Metadata of the 39,495 samples in the GlobalFungi v.5 dataset, including the proportion of EcM fungal OTUs that were assigned and unassigned to species-level identities.
-	“GFv5_EcM_order_unassigned_continents.csv”: breakdown of the proportion of EcM fungal OTUs assigned and unassigned to species-level by fungal orders and continents (donut plots in Figure 1).

### 2.EcM_fungal_species_in_the_Catalogue_of_Life
Code and output used to estimate the number of EcM fungal species in the Catalogue of Life.

Files:
- "README.md": Full instrcutions for dowloading and filtering all data
- "ECM_genera2.csv": The list of 327 EcM fungal genera from FungalTraits v.1.1
- "NameUsage.tsv": List of all fungi names (accepted and synonyms) downloaded from the Catalogue of Life
- "metadata.yaml": Metadata from Catalogue of Life download
- "ECM_accepted_match.tsv": Output file containing the final list of species used

### 3.EcM_fungal_species_from_Australia
Code and data used to download/extract the number of named EcM fungal species that have been recorded in Australia from the National Species List, GBIF, GlobalFungi v.5, and the Global Soil Mycobiome consortium (GSMc). 

Files:
-	“ECM_species_names_Australia.r”: R code used to download and summarise EcM fungal species recorded in Australia, and to create Figure 2.
-	“ECM_final_species_list_w_UNITE.csv”: The output file from the R code of EcM fungal species names present in the different databases used to create Figure 2.

### 4.Dark_taxa_richness_maps
Raster layers of total EcM fungal richness, dark taxa EcM richness, percentage dark taxa, and the dark taxa research priority metric plotted in Figure 3.

Files:
-	“Figure_3_rasters.tif”: raster layers (equirectangular projection, 30 arc-second resolution) of the maps presented in Figure 3.

### 5.EcM_host_plant_family_distribution_maps
Code and data used to create the EcM host plant family distribution maps.

Files:
-	“GBIF_record_details_and_download_DOIs.xlx”: List of the EcM host plant genera, and their corresponding families, extracted from the FungalRoot database. Occurrence records for each genus were downloaded from GBIF (see citation DOIs), and rigorously cleaned (see manuscript Methods). The “Clean_records” column shows the number of occurrence records remaining after cleaning, and therefore used to create the host plant distribution maps.
-	“Create_host_distribution_maps.r”: R code used to download and clean the GBIF records, and create the raster layers of host plant family distributions.


