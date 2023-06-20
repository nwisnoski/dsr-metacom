# --------------------------------------------------------- #
# Format raw data                       #
# FCE diatoms                                                       #
# Revised 01 Jun 2017                                       #
# --------------------------------------------------------- #

# Contributors: Aldo Compagnoni, Chris Catano, Riley Andrade, Max Castorani, Nina Lany, Sydne Record, Nicole Voelker

# modified by N Wisnoski 20 Jun 2023

library(tidyverse)

# Clear environment
rm(list = ls())

options( stringsAsFactors = F )

diat <- read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-fce.1211.4&entityid=c55358e65d766fc8f2b2b3b28dec4600")
diat$Year <- year(diat$OBS_DATE)
diat <- diat %>% filter(
  WETLAND_BASIN %in% c("SRS", "TSL"),
  Year %in% 2005:2014) 


# format diatom data
vars <- c("OBSERVATION_TYPE","SITE_ID","DATE",
          "VARIABLE_NAME","VARIABLE_UNITS","VALUE")

# coordinates ----------------------------------------------------------------

# find a centroid for each PSU
centr   <- diat %>%
            mutate(PRIMARY_SAMPLING_UNIT = as.factor(PRIMARY_SAMPLING_UNIT) ) %>%
            select(PRIMARY_SAMPLING_UNIT, EASTING_UTM, NORTHING_UTM) %>%
            group_by(PRIMARY_SAMPLING_UNIT) %>%
            summarise( east_centroid = mean(EASTING_UTM),
                       north_centroid = mean(NORTHING_UTM) )

# check we didn't loose the PSU values (through "as.factor")
expect_equal(length(intersect(unique(as.character(diat$PRIMARY_SAMPLING_UNIT)),
                       as.character(centr$PRIMARY_SAMPLING_UNIT))),
             length(unique(diat$PRIMARY_SAMPLING_UNIT)))


# coordinates data frame
coord <- data.frame( OBSERVATION_TYPE = "SPATIAL_COORDINATE",
                      SITE_ID = as.character(centr$PRIMARY_SAMPLING_UNIT),
                      DATE = NA,
                      VARIABLE_NAME = "EASTING",
                      VARIABLE_UNITS = "METERS",
                      VALUE = centr$east_centroid) %>%
            rbind( data.frame( OBSERVATION_TYPE = "SPATIAL_COORDINATE",
                               SITE_ID = as.character(centr$PRIMARY_SAMPLING_UNIT),
                               DATE = NA,
                               VARIABLE_NAME = "NORTHING",
                               VARIABLE_UNITS = "METERS",
                               VALUE = centr$north_centroid)
            )



# taxon count 

# calculate mean of diatom species across samples
# (NOT NEEDED, number of rows does not change)
mean_diat <- diat %>% 
                # # modify the one instance of sampling in january
                # mutate( Date = replace(Date, 
                #          grepl('^1/[0-9]{1,2}/2011', Date),
                #          '12/31/2010') ) %>% 
                select(Year, PRIMARY_SAMPLING_UNIT, ACHADNADN:UNKWNVALV) %>%
                group_by(Year, PRIMARY_SAMPLING_UNIT) %>%
                summarise_all( funs(mean) )

# test equality of two data sets
expect_equal(nrow(mean_diat), nrow(diat) )

# long-form this data set
diat_long <- diat %>%
              select(Year, PRIMARY_SAMPLING_UNIT, ACHADNADN:UNKWNVALV) %>%
              gather(VARIABLE_NAME, VALUE,-Year,-PRIMARY_SAMPLING_UNIT) %>%
              mutate(VARIABLE_UNITS = "PERCENT",
                     OBSERVATION_TYPE = "TAXON_RELATIVE_ABUNDANCE") %>%
              rename(DATE = Year,
                     SITE_ID = PRIMARY_SAMPLING_UNIT) %>%
              select_(.dots = vars)
  


# put it all together
dat.long <- Reduce(function(...) rbind(...),
                   list(coord, diat_long)) %>%
                # convert SITE_ID to numeric
                mutate( SITE_ID = as.numeric(SITE_ID) )


# PART coded by Chris  ---------------------------------------------------------------------------------------------------

str(dat.long)
levels(dat.long$OBSERVATION_TYPE)

# Change 'TAXON_RELATIVE_ABUNDANCE' to 'TAXON_COUNT'
dat.long$OBSERVATION_TYPE <- gsub("TAXON_RELATIVE_ABUNDANCE", "TAXON_COUNT", dat.long$OBSERVATION_TYPE)

# Write this to the L3 folder in Google Drive 
#write.csv(dat.long, file = "~/Google Drive File Stream/My Drive/LTER Metacommunities/LTER-DATA/L3-aggregated_by_year_and_space/L3-fce-diatoms-marazzi.csv")


# MAKE DATA LIST
dat <- list()

# COMMUNITY DATA 
comm.long <- dat.long[dat.long$OBSERVATION_TYPE == "TAXON_COUNT", ] 
comm.long <- comm.long %>%
  droplevels()

str(comm.long)  # Inspect the structure of the community data

# Add number of unique taxa and number of years to data list:
dat$n.spp <- length(unique(comm.long$VARIABLE_NAME))
dat$n.years <- length(unique(comm.long$DATE))

# Ensure that community data VALUE and DATE are coded as numeric
comm.long <- comm.long %>%   # Recode if necessary
  mutate_at(vars(c(DATE, VALUE)), as.numeric)

# Ensure that community character columns coded as factors are re-coded as characters
comm.long <- comm.long %>%   # Recode if necessary
  mutate_if(is.factor, as.character)
  
# Ensure that SITE_ID is a character: recode numeric as character 
comm.long <- comm.long %>%   # Recode if necessary
  mutate_at(vars(SITE_ID), as.character)

# Double-check that all columns are coded properly
ifelse(FALSE %in% 
   c(
     class(comm.long$OBSERVATION_TYPE) == "character",
     class(comm.long$SITE_ID) == "character",
     class(comm.long$DATE) == "numeric",
     class(comm.long$VARIABLE_NAME) == "character",
     class(comm.long$VARIABLE_UNITS) == "character",
     class(comm.long$VALUE) == "numeric"
     #class(comm.long$TAXON_GROUP) == "character")
   ),
  "ERROR: Community columns incorrectly coded.", 
  "OK: Community columns correctly coded.")


# checking for weird taxa names and removing suspicious ones that are rare
comm.long_gamma_summary <- comm.long %>% group_by(VARIABLE_NAME) %>% 
  filter(!is.na(VARIABLE_NAME) & VALUE>0) %>%
  summarize(mean_value = mean(VALUE, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(RA = mean_value / sum(mean_value))

comm.long_gamma_summary_remove_unk <- comm.long_gamma_summary %>% 
  filter(!grepl('(?i)unk',VARIABLE_NAME))

if(!nrow(comm.long_gamma_summary) == nrow(comm.long_gamma_summary_remove_unk)){
  message('WARNING: suspicous taxa removed -- taxaID had "unk"')
}

# Subset data if necessary
keep <- comm.long_gamma_summary_remove_unk$VARIABLE_NAME
comm.long.ENP2 <- comm.long[comm.long$VARIABLE_NAME %in% keep, ] # records dropped
length(unique(comm.long.ENP2$VARIABLE_NAME)) # removed 2 species (UNKNGIRD, UNKNVALV)


# ---------------------------------------------------------------------------------------------------
# Check balanced sampling of species across space and time by inspecting table, and add to data list

if(length(unique(xtabs(~ SITE_ID + DATE, data = comm.long.ENP2))) > 1){
  comm.long2 <- comm.long.ENP2 %>%
    spread(VARIABLE_NAME, VALUE, fill = 0) %>%
    gather('VARIABLE_NAME','VALUE', -c(SITE_ID, DATE, VARIABLE_UNITS))
}

(cont.table <- xtabs(~ SITE_ID + DATE, data = comm.long2)) #number of taxa should be same in all
cont.table <- as.data.frame(cont.table)
hist(na.omit(comm.long2$DATE)) #frequency should be same (even distribution)


ifelse(length(unique(xtabs(~ SITE_ID + DATE, data = comm.long.ENP2))) == 1,
       "OK: Equal number of taxa recorded across space and time.", 
       "ERROR: Unequal numbers of observations across space and time, or taxa list not fully propagated across space and time. Inspect contingency table.")

# keep sites that have same sample coverage/spacing (remove those with gaps)
a <- cont.table[cont.table$Freq == 0, ]
irregular <- unique(a$SITE_ID)
comm.long.ENP3 <- comm.long.ENP2[!(comm.long.ENP2$SITE_ID %in% irregular), ]

xtabs(~ SITE_ID + DATE, data = comm.long.ENP3) #number of taxa should be same in all
hist(na.omit(comm.long.ENP3$DATE)) #frequency should be same (even distribution)



# ---------------------------------------------------------------------------------------------------
# Convert community data to wide form
comm.wide <- comm.long.ENP3 %>%
  select(-VARIABLE_UNITS) %>%
  spread(VARIABLE_NAME,  VALUE)

# check to make sure each species is sampled at least once
zero <- colSums(comm.wide[, c(4:285)]) == 0
zero <- as.data.frame(zero)
zero$species <- rownames(zero)

# Remove species with zero occurence (zero = True)
zero2 <- zero[!(zero$zero == TRUE), ]
species <- zero2[, -1]

# long form data
dat.long.ENP <- comm.long.ENP3[comm.long.ENP3$VARIABLE_NAME %in% species, ] # records dropped

# Convert community data to wide form
comm.wide.ENP <- dat.long.ENP %>%
  select(-VARIABLE_UNITS) %>%
  spread(VARIABLE_NAME,  VALUE)

length(unique(dat.long.ENP$SITE_ID)) # 30 sites 
unique(dat.long.ENP$DATE) # 7 years (2005-2011)
length(unique(dat.long.ENP$VARIABLE_NAME)) # 192 species


write.csv(dat.long.ENP, 
          file = here("data/L3-fce-diatoms-catano.csv"), 
          row.names = F)


