### Cleaning MOOREA CORAL REEF (MCR) data

# --------------------------------------------------------------------------------------------------------------------------------

# Max Castorani
# 8 June 2017

## Data manipulation packages
library(dplyr)
library(tidyr)

#source("Group2-explore-data/format_data/pull_data_gdrive_fun.R")

# --------------------------------------------------------------------------------------------------------------------------------

### MCR Algae Data ###

# ---------------------------------------------------------------------------------------------
## Read in the data from EDI Data Portal

# Package ID: knb-lter-mcr.8.28 Cataloging System:https://pasta.lternet.edu.
# Data set title: MCR LTER: Coral Reef: Long-term Population and Community Dynamics: Benthic Algae and Other Community Components, ongoing since 2005.
# Data set creator:    - Moorea Coral Reef LTER 
# Data set creator:  Robert Carpenter - Moorea Coral Reef LTER 
# Contact:    - Information Manager LTER Network Office  - tech-support@lternet.edu
# Contact:    - Information Manager Moorea Coral Reef LTER  - mcrlter@msi.ucsb.edu
# Metadata Link: https://portal.lternet.edu/nis/metadataviewer?packageid=knb-lter-mcr.8.28
# Stylesheet for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

infile1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-mcr/8/28/54d54c25616a48b9ec684118df9d6fca"
infile1 <- sub("^https","http",infile1)
mcr.algae <-read.csv(infile1,header=F
              ,skip=1
              ,sep=","
              , col.names=c(
                "Year",
                "Date",
                "Location",
                "Site",
                "Habitat",
                "Transect",
               "Quadrat",
                "Taxonomy_Substrate_Functional_Group",
                "Percent_Cover"), check.names=TRUE, stringsAsFactors = FALSE)

# Replace underscores with dots for convenience. Also convert to lowercase.
colnames(mcr.algae) <- tolower(gsub("_", ".", colnames(mcr.algae)))

# Fix species 
mcr.algae <- mcr.algae %>%
  dplyr::rename(species = taxonomy.substrate.functional.group) %>%
  dplyr::filter(species != "Coral",  # Remove the non-focal taxa (non-algae)
                species != "Ascidians",
                species != "Sponge",
                species != "Shell Debris",
                species != "Coral Rubble",
                species != "No data",
                species != "Bare Space",
                species != "Soft Coral",
                species != "Sand",
                species != "",
                !is.na(species)
  ) %>%
  droplevels()

# Clean the dataset
mcr.algae$year <- as.numeric(as.character(mcr.algae$year))
mcr.algae$year[is.na(mcr.algae$year)] <- 2006

# Aggregate data by year, taxon, habitat, transect and site (6 sites total)
mcr.algae_clean <- mcr.algae %>%
  group_by(year, site, habitat, transect, species) %>%
  dplyr::summarise(percent.cover = mean(percent.cover, na.rm = TRUE)) %>%
  ungroup() %>%
  droplevels() 

# Convert from long to wide and back to long to be sure that we have fully propagated taxa
mcr.algae_clean_wide <- spread(mcr.algae_clean, key = species, value = percent.cover, fill = 0)
mcr.algae_clean_long <- gather(mcr.algae_clean_wide, key = species, value = percent.cover, -year, -site, -transect, -habitat)

# Finish cleaning data by renaming and adding columns
mcr.algae_clean <- mcr.algae_clean_long %>%
  mutate(project = "algae",  # rename what they called site to what we call project
         plot = site) %>%
  select(-site) %>%
  mutate(site = "mcr",   # format column names
         plot = paste0("location_", sapply(strsplit(as.character(plot), " "), tail, 1)),
         subplot = paste0("transect_", transect),
         abundance = percent.cover,  
         unitAbund = "mean.percent.cover",
         scaleAbund = "0.25_m2",
         growth = NA,
         uniqueID = paste(site, project, plot, subplot, sep = "_"),
         guild = "algae") %>%
  select(year, site, habitat, project, plot, subplot, uniqueID, guild, species, abundance, unitAbund, scaleAbund) #, growth)

# Remove unneeded files
rm(mcr.algae, mcr.algae_clean_long, mcr.algae_clean_wide)

# --------------------------------------------------------------------------------------------------------------------------------

# Reformat column names
mcr.algae_reformat <- mcr.algae_clean %>%
  dplyr::mutate(OBSERVATION_TYPE = "TAXON_COUNT",
                VARIABLE_UNITS = paste0(unitAbund, ".per.", scaleAbund),
                UNIQUE_SPATIAL_ID = paste(plot, habitat, subplot, sep = "_")) %>%
  dplyr::rename(VALUE = abundance,
                VARIABLE_NAME = species,
                DATE = year,
                SITE_ID = plot,
                HABITAT = habitat,
                SUB_SITE_ID = subplot) %>%
  dplyr::select(OBSERVATION_TYPE,
                SITE_ID, 
                HABITAT, 
                SUB_SITE_ID, 
                UNIQUE_SPATIAL_ID,
                DATE,
                VARIABLE_NAME,
                VARIABLE_UNITS,
                VALUE)

# --------------------------------------------------------------------------------------------------------------------------------
# Aggregate by site, then add spatial information
mcr.algae_L3 <- mcr.algae_reformat %>%
  dplyr::group_by(SITE_ID, DATE, VARIABLE_NAME) %>%
  dplyr::summarise(OBSERVATION_TYPE = unique(OBSERVATION_TYPE),
                   VARIABLE_UNITS = unique(VARIABLE_UNITS),
                   VALUE = mean(VALUE, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(OBSERVATION_TYPE, SITE_ID, DATE, VARIABLE_NAME, VARIABLE_UNITS, VALUE) %>%
  as.data.frame(.)

# Replace underscores with dots in location IDs for future plotting. 
mcr.algae_L3$SITE_ID <- gsub("_", "", mcr.algae_L3$SITE_ID)


mcr.algae_L3_final <- rbind(mcr.algae_L3)

# Write CSV file for cleaned data (L3)
write.csv(mcr.algae_L3_final, file = "data/L3-mcr-algae-castorani.csv", row.names = F)
