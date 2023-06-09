# Aldo Compagnoni Nov 2018
#modified by NKL 03/25/2019
rm(list = ls())

library(dplyr)
library(tidyr)


#Data are from Ecologcial Archives
## http://esapubs.org/archive/ecol/E088/161/
##Peter B. Adler, William R. Tyburczy, and William K. Lauenroth. 2007. Long-term mapped quadrats from Kansas prairie: demographic information for herbaceaous plants. Ecology 88:2673.


hays <- read.csv('https://figshare.com/ndownloader/files/5599934', stringsAsFactors=F) 
spp_list <- read.csv('https://figshare.com/ndownloader/files/5599943', stringsAsFactors=F) 
quad_info <- read.csv('https://figshare.com/ndownloader/files/5599937', stringsAsFactors=F) 
quad_inv <- read.csv('https://figshare.com/ndownloader/files/5599940', stringsAsFactors=F)


# hays raw
hays_raw <- hays %>% 
              mutate( plotyear = as.character(plotyear) ) %>% 
              # rename plot that messes with pattern
              mutate( plotyear = replace(plotyear,
                                         grepl( 'e2qo-10',plotyear),
                                         'e2qo-X')) %>% 
              separate( plotyear, c('plot','year'), sep='-') %>% 
              mutate( plot_n = substr(year,1,1) ) %>% 
              mutate( year   = substr(year,2,3) ) %>% 
              # re-rename plot that messed with pattern
              mutate( plotyear = replace(plot,
                                         plot == 'e2qo-X',
                                         'e2qo-10')) %>% 
              mutate( plot = paste(plot,plot_n,sep='-') ) %>% 
              # count species by year/plot
              count(plot,year,species) %>% 
              # remove non-plants and unknowns
              subset( !(species %in% c('Unknown.',
                                       'Unknown',
                                       'Short grass',
                                       'Fragment',
                                       'Bare ground')) )

# visually examine items with unclear taxonomic ID
hays_raw %>% 
  group_by(species) %>% 
  summarise(all_n = sum(n) ) %>% 
  print(n=200)

# genera to remove (less than 5% of all counts)
r_gens <- c('Ambrosia spp.','Oxalis spp.', 'Solidago spp.')
# genera to lump (species identified lass than 95% of the time)
lump_gens <- c('Allium', 'Chamaesyce', 'Opuntia', 'Polygala')
         
# final hays data frame    
hays_out <- hays_raw  %>% 
              # remove a few individuals identified at the genera level
              subset( !(species %in% r_gens) ) %>%
              # "lump" gens. iddentified at spp. level less than 95% of counts
              mutate( species = replace(species,
                                        grepl('Allium',species),
                                        'Allium spp.') ) %>% 
              mutate( species = replace(species,
                                        grepl('Chamaesyce',species),
                                        'Chamaesyce spp.') ) %>%
              mutate( species = replace(species,
                                        grepl('Opuntia',species),
                                        'Opuntia spp.') ) %>% 
              mutate( species = replace(species,
                                        grepl('Opuntia',species),
                                        'Opuntia spp.') ) %>% 
              # re-count individuals 
              group_by(plot,year,species) %>% 
              summarise( count = sum(n) ) %>% 
              ungroup %>% 
              rename(SITE_ID = plot, 
                   DATE = year,
                   VARIABLE_NAME = species,
                   VALUE = count) %>%
              mutate(VARIABLE_UNITS = 'COUNT',
                     OBSERVATION_TYPE = 'TAXON_COUNT') %>%
              select(OBSERVATION_TYPE,SITE_ID, DATE, 
                     VARIABLE_NAME, VARIABLE_UNITS, VALUE) %>% 
              # insert zeros
              spread(key = VARIABLE_NAME, value = VALUE, fill = 0) %>% 
              gather(key = VARIABLE_NAME, value = VALUE, -DATE, -SITE_ID, -OBSERVATION_TYPE, -VARIABLE_UNITS) 
      
# keep plots with continuous rep
keep_plot <- hays_out %>% 
              mutate( DATE = as.numeric(DATE) ) %>% 
              subset( DATE > 37 & DATE < 73) %>% 
              select(SITE_ID,DATE) %>% 
              unique %>% 
              count(SITE_ID) %>% 
              subset( n == 35) %>% 
              .$SITE_ID

# write it out
hays_out %>% 
  subset( SITE_ID %in% keep_plot) %>% 
  subset( DATE >37 ) %>% 

write.csv(file = 'data/L3-hays-plants-compagnoni.csv',row.names=F)  

