### Cleaning MOOREA CORAL REEF (MCR) data

# --------------------------------------------------------------------------------------------------------------------------------

# Max Castorani
# 8 June 2017

## Data manipulation packages
library(dplyr)
library(tidyr)

#source("Group2-explore-data/format_data/pull_data_gdrive_fun.R")

# --------------------------------------------------------------------------------------------------------------------------------

### MCR Coral Data ###
# ---------------------------------------------------------------------------------------------
## Read in the data from EDI Data Portal
# Package ID: knb-lter-mcr.4.35 Cataloging System:https://pasta.edirepository.org.
# Data set title: MCR LTER: Coral Reef: Long-term Population and Community Dynamics: Corals, ongoing since 2005.
# Data set creator:    - Moorea Coral Reef LTER 
# Data set creator:  Peter Edmunds - Moorea Coral Reef LTER 
# Contact:    - Information Manager Moorea Coral Reef LTER  - mcrlter@msi.ucsb.edu
# Stylesheet for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

infile1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-mcr/4/35/a9c9466ca7831118dc9f1eabe179e61e"
infile1 <- sub("^https","http",infile1)
mcr.coral <-read.csv(infile1,header=F
              ,skip=1
              ,sep=","
              , col.names=c(
                "Date",
                "Location",
                "Site",
                "Habitat",
                "Section_of_Transect",
                "Quadrat",
                "Taxonomy_Substrate_or_Functional_Group",
                "Percent_Cover"    ), check.names=TRUE, stringsAsFactors = FALSE)
rm(infile1)

# Replace underscores with dots for convenience. Also convert to lowercase.
colnames(mcr.coral) <- tolower(gsub("_", ".", colnames(mcr.coral)))

# Code species guild
mcr.coral <- mcr.coral %>%
  dplyr::rename(species = "taxonomy.substrate.or.functional.group",
                transect = "section.of.transect") %>%
  dplyr::filter(species != "Sand",              # Drop species that are non-coral
                species != "Unknown or Other",
                species != "Macroalgae",
                species != "Crustose Coralline Algae / Bare Space",
                species != "Turf",
                species != "Non-coralline Crustose Algae",
                !is.na(species),
                !is.na(percent.cover)) %>%
  droplevels() %>%
  # Convert date to year
  mutate(year = as.numeric(strtrim(as.character(date), 4)))

# For each species, average the abundance data by year, habitat, plot ('site'), and subplot ('transect')
mcr.coral_clean <- mcr.coral %>%
  group_by(year, site, habitat, transect, species) %>%
  dplyr::summarise(percent.cover = mean(percent.cover, na.rm = TRUE)) %>%
  ungroup() %>%
  droplevels() 

# Convert from long to wide and back to long to be sure that we have fully propagated taxa
mcr.coral_clean_wide <- spread(mcr.coral_clean, key = species, value = percent.cover, fill = 0)
mcr.coral_clean_long <- gather(mcr.coral_clean_wide, key = species, value = percent.cover, -year, -site, -transect, -habitat)

# Finish cleaning data by renaming and adding columns
mcr.coral_clean <- mcr.coral_clean_long %>%
  mutate(project = "coral",  # rename what they called site to what we call project
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
         guild = "coral") %>%
  select(year, site, habitat, project, plot, subplot, uniqueID, guild, species, abundance, unitAbund, scaleAbund) #, growth)

# Remove unneeded files
rm(mcr.coral, mcr.coral_clean_long, mcr.coral_clean_wide)

# --------------------------------------------------------------------------------------------------------------------------------

# Reformat column names
mcr.coral_reformat <- mcr.coral_clean %>%
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
mcr.coral_L3 <- mcr.coral_reformat %>%
  dplyr::group_by(SITE_ID, DATE, VARIABLE_NAME) %>%
  dplyr::summarise(OBSERVATION_TYPE = unique(OBSERVATION_TYPE),
                   VARIABLE_UNITS = unique(VARIABLE_UNITS),
                   VALUE = mean(VALUE, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(OBSERVATION_TYPE, SITE_ID, DATE, VARIABLE_NAME, VARIABLE_UNITS, VALUE) %>%
  as.data.frame(.)

# Replace underscores with dots in location IDs for future plotting. 
mcr.coral_L3$SITE_ID <- gsub("_", "", mcr.coral_L3$SITE_ID)


mcr.coral_L3_final <- rbind(mcr.coral_L3)

# Write CSV file for cleaned data (L3)
write.csv(mcr.coral_L3_final, file = "data/L3-mcr-coral-castorani.csv", row.names = F)
