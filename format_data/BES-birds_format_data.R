# --------------------------------------------------------- #
# Format raw data from L0 into L3 format                    #
# BES Birds Nilon                                           #
# Revised 8 Nov 2018 by Shawn Taylor                        #
# --------------------------------------------------------- #

# Contributors: Riley Andrade, Max Castorani, Nina Lany, Sydne Record, Nicole Voelker

#Here is some template or example code that may be helpful for putting datasets into the format for the LTER metacommunities working group synthesis project and performing a few data checks. 
#Make one script for each dataset that converts the data from the form we originally got it in to the one we need. This might involve aggregating subplots and/or multiple sampling occasins by year, subsetting out sites or species, etc. It may also be necessary to fill in an abundance of zero for any taxa that appear in the dataset but were not observed at a particular site-year.
#Save the resulting file in the directory "~/Google Drive/LTER Metacommunities/LTER-DATA/L3-aggregated-by-year-and-space". This file should be ready for preliminary analalysis.
#The next step in the pipeline after this one is to run each dataset through the script 3_explore_comm_dat.R. This script plots the sampling effort (to be sure overvations are balanced), the species accumulation curve, and a time series of species richness at each site as well as in aggragate. If these preliminiary visualizations turn up something problematic, go back to the script spcific to the dataset and add in annotated code to fix the problem.


#####################################################
##### Cleaning steps ################################
# I applied the following filters to the BES birds dataset
# 
# - Dropped all observations of distance 'FT", which stands for "fly thru".
# - removied uknown species
# - where there was multipel surveys for a site in a single year,
#   keep the maximum count for each species among all the surveys.
# - only kept sites with surveys every year from 2005 - 2009 (126 sites)
# 
# A note on spatial coordinates
# The only indicator for location is Park Names, or numbered locations small regions.
# (ie. WS 263 Plot 57 stands for Watershed 263). I tried pull lat/longs from a lookup 
# service, but it only worked on a small subset of sites. If lat/longs are required
# someone will likely have to go thru and locate them manually, or ask the data provider 
# directly. 

# Clear environment

require(tidyverse)
# ---------------------------------------------------------------------------------------------------
#IMPORT AND FORMAT DATA 

#IMPORT DATA FROM A DIFFERENT SOURCE THAT IS STORED ON GOOGLE DRIVE (you'll have to format it yourself)
#Here is an example of reading in a dataset using the Google Drive ID. This is for datasets that we obtained somehow and then put in the directory "~/Google Drive/LTER Metacommunities/LTER-DATA/L0-raw". 

# Assign data set of interest in the L0 data folder.
# NOTE: Google Drive file ID is different for each dataset

# Baltimore Urban LTER
data.set <- "BES-birds-Nilon"
observations <- read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-bes.543.170&entityid=0329310e827da4f12986224abad4d0a8")
sites <- read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-bes.543.170&entityid=7197e52c126c5881c69c194b0d4f5edd")
surveys <- read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-bes.543.170&entityid=88715faecbee458c7c16b7db47bd8a41")
species_list <- read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-bes.543.170&entityid=494044ed31ec80c9d751d491929b71fc")

surveys$survey_date <- as.Date(surveys$survey_date, format='%m/%d/%Y')

# For some reason there is a whole bunch of repeated rows in the surveys file
surveys <- surveys %>%
  distinct()



########################################
# Only keep species which were seen at < 40ft
# 'FT' indicates a "fly through", where a bird was seen flying thru
# the area (but most likely > 40 feet away)
observations <- observations %>%
  filter(!distance %in% c('FT'))

# Drop unknown or unspecified species codes
species_codes_to_drop = c('UNAC','UNSP','temp','UNKN')
observations <- observations %>%
  filter(!species_id %in% species_codes_to_drop)

########################################
# Initial filters, potentially dropping sites due to lack of spatial coordinates

# aggregate obserations to total counts per survey
# accounts for an observer counting 2 groups of the same species
observations_aggregated <- observations %>%
  group_by(survey_id, species_id) %>%
  summarise(count = sum(bird_count)) %>% 
  ungroup()

# fill in 0's for all species/sites/years
observations_aggregated <- observations_aggregated %>%
  complete(survey_id, species_id, fill = list(count = 0))

# Add in the dates and keep only the year from the survey info
observations_aggregated <- observations_aggregated %>%
  left_join(select(surveys, site_id, survey_id, survey_date), by=c('survey_id')) %>%
  mutate(year = lubridate::year(survey_date)) %>%
  select(-survey_id, -survey_date)


# For multiple surveys in a single year, take the maximium count
# for each species/site
observations_aggregated <- observations_aggregated %>%
  group_by(year, site_id, species_id) %>%
  top_n(1, -count)  %>%
  ungroup() %>%
  distinct()

# only keep sites with surveys 2005 - 2009
years_to_keep = 2005:2009
observations_aggregated <- observations_aggregated %>%
  filter(year %in% years_to_keep)

sites_to_keep <- observations_aggregated %>%
  group_by(site_id) %>%
  summarise(n_years = n_distinct(year)) %>%
  filter(n_years == length(years_to_keep)) %>%
  pull(site_id) %>%
  unique()

observations_aggregated <- observations_aggregated %>%
  filter(site_id %in% sites_to_keep)

# Drop any species which aren't ever actually observed.
observations_aggregated <- observations_aggregated %>%
  group_by(species_id) %>%
  filter(sum(count) > 0) %>%
  ungroup()

# Drop sites that had 0 species recorded in a year. These sites
# won't work with the meta analsysis. 
sites_with_no_species_in_atleast_one_year = observations_aggregated %>%
  group_by(site_id, year) %>%
  summarise(total_abund = sum(count)) %>%
  ungroup() %>%
  filter(total_abund == 0) %>%
  pull(site_id) %>%
  unique()

observations_aggregated <- observations_aggregated %>%
  filter(!site_id %in% sites_with_no_species_in_atleast_one_year)

################################################
# QA Check, the sites/years in the final version
# actually have recorded surveys
surveys <- surveys %>%
  mutate(year = lubridate::year(survey_date)) %>%
  select(site_id, year) %>%
  distinct() %>%
  mutate(was_surveyed = TRUE)

observations_aggregated <- observations_aggregated %>%
  left_join(surveys, by=c('site_id','year'))

if(!all(observations_aggregated$was_surveyed)) stop('Potentially unsurveyed year/sites in final dataset')
################################################


# Convert to final working group format
final_format <- observations_aggregated %>%
  mutate(OBSERVATION_TYPE = 'TAXON_COUNT',
         VARIABLE_UNITS = 'count') %>%
  select(OBSERVATION_TYPE,
         SITE_ID = site_id,
         DATE = year,
         VARIABLE_NAME = species_id,
         VARIABLE_UNITS,
         VALUE = count) %>%
  arrange(DATE)


write.csv(final_format, 
          file = here("data/L3-bes-birds-nilon.csv"),
          row.names = F)

