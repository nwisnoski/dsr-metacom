# --------------------------------------------------------- #
# Format raw data from L0 into L3 format                    #
# handy code TEMPLATE                                       #
# Revised 15 May 2018 by ERS                                #
# --------------------------------------------------------- #

options(stringsAsFactors = FALSE)

# Clear environment
rm(list = ls())


library(tidyverse)

# read data from EDI portal
# Package ID: knb-lter-sgs.527.1 Cataloging System:https://pasta.lternet.edu.
# Data set title: SGS-LTER Effects of grazing on ecosystem structure and function (GZTX): Vegetation basal cover on the Central Plains Experimental Range, Nunn, Colorado, USA 1992-2011, ARS Study Number 32.
# Data set creator:  Daniel Milchunas - Natural Resource Ecology Lab 
# Metadata Provider:    - Colorado State University 
# Contact:    - Information Manager LTER Network Office  - tech-support@lternet.edu
# Contact:  Daniel Milchunas -  Natural Resource Ecology Lab  - daniel.milchunas@colostate.edu
# Metadata Link: https://portal.lternet.edu/nis/metadataviewer?packageid=knb-lter-sgs.527.1
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

infile1 <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sgs/527/1/c3dc7ed391a6f6c6eca4512d27d0d091" 
infile1 <- sub("^https","http",infile1) 
dt1     <-read.csv(infile1,header=F 
          ,skip=1
            ,sep="\t"  
        , col.names=c(
                    "Year",     
                    "Site",     
                    "Treatment",     
                    "Plot",     
                    "Species",     
                    "Basal_Cover"    ), check.names=TRUE)
               
# format based on ltermetacomm data format
sgs  <- dt1 %>% 
            # only retain sites that WERE ALWAYS GRAZED
            subset( Treatment == 'GZGZ' ) %>% 
            # data gap before 1995
            subset( Year > 1994 ) %>% 
            # aggregated by treatment
            group_by( Year, Site, Species ) %>% 
            summarise( Basal_Cover = mean(Basal_Cover) ) %>% 
            ungroup() %>% 
            rename(SITE_ID = Site, 
                   DATE = Year,
                   VARIABLE_NAME = Species,
                   VALUE = Basal_Cover) %>%
            mutate(VARIABLE_UNITS = 'PERCENT_COVER',
                   OBSERVATION_TYPE = 'TAXON_COUNT') %>%
            select(OBSERVATION_TYPE,SITE_ID, DATE, VARIABLE_NAME, VARIABLE_UNITS, VALUE) %>% 
            # insert zeros
            spread(key = VARIABLE_NAME, value = VALUE, fill = 0) %>% 
            gather(key = VARIABLE_NAME, value = VALUE, -DATE, -SITE_ID, -OBSERVATION_TYPE, -VARIABLE_UNITS) 
          
sgs_out <- sgs %>% 
  filter(VARIABLE_NAME != "ASOX",
         VARIABLE_NAME != "BARE",
         VARIABLE_NAME != "LICH",
         VARIABLE_NAME != "LITT",
         VARIABLE_NAME != "UNFB")

# write it
write.csv(sgs_out, 
          "data/L3-sgs-plants-compagnoni.csv",
          row.names=F)
