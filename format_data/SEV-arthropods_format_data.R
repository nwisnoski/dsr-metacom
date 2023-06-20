# --------------------------------------------------------- #
# Format raw data from L0 into L3 format                    #
# handy code TEMPLATE                                       #
# Revised 15 May 2018 by ERS                                #
# --------------------------------------------------------- #


# Contributors: Riley Andrade, Max Castorani, Nina Lany, Sydne Record, Nicole Voelker
# revised by Eric Sokol
# revised by Aldo Compagnoni

options(stringsAsFactors = FALSE)

# Clear environment
rm(list = ls())

library(tidyverse)

# ---------------------------------------------------------------------------------------------------
#IMPORT AND FORMAT DATA FROM EDI
#1. read in flat tables from EDI. The code below was modified from the code that is available on the EDI page for each dataset (look for an 'import data from R' link). Once you obtain this code for the dataset of interest, the rest of the steps should run without modification. Be sure to add 'stringsAsFactors = F' to the code that imports each csv.

infile1   <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sev/29/175390/1a8ca1a97279d2e189452071665ae584" 
infile1   <- sub("^https","http",infile1) 
dt1       <-read.csv(infile1,header=T ,sep=",", 
                     col.names=c(
                    "Year",     
                    "Month",     
                    "Day",     
                    "Site",     
                    "Line.",     
                    "Trap",     
                    "Order",     
                    "Family",     
                    "Genus",     
                    "Species",     
                    "Count",     
                    "comments"    ), 
        check.names=TRUE,
        stringsAsFactors = F)

# function to remove NAs (-888)
remove_888 <- function(x) replace(x, x == -888, NA)
 
# format the SEV dataset
sev <- dt1 %>% 
          # remove missing site information (!!!)
          subset( !(Site %in% '') ) %>%
          # remove sites with less inconsistent sampling
          subset( !(Site %in% c('P','B')) ) %>% 
          # remove "pooled" line information
          subset( !(Line. %in% 'pooled') ) %>%
          # remove missing Trap information
          subset( !is.na(Trap) ) %>%
          # remove missing Trap information
          subset( !(Trap %in% c(2,4,6)) ) %>% 
          # trim white spaces
          mutate( Line.   = trimws(Line.),
                  Order   = trimws(Order),
                  Order   = trimws(Family),
                  Genus   = trimws(Genus),
                  Species = trimws(Species) ) %>% 
          # substitute -888 with NAs
          mutate( Order   = remove_888(Order) ) %>% 
          mutate( Family  = remove_888(Family) ) %>% 
          mutate( Genus   = remove_888(Genus) ) %>% 
          mutate( Species = remove_888(Species) ) %>% 
          mutate( Count   = remove_888(Count) ) %>% 
          # remove taxa not identified to the species level
          subset( !grepl('[0-9]',Order ) ) %>%        
          subset( !grepl('[0-9]',Family ) ) %>%
          subset( !grepl('[0-9]',Genus ) ) %>% 
          subset( !grepl('[0-9]',Species ) ) %>% 
          # "The traps within each trap line are subsamples, 
          # and data from those should be summed or averaged 
          # for a single value per line, per sample period
          # I get means because sampling is not 6 times/year as stated
          group_by( Year, Site, Line.,    
                    Order, Family, Genus, Species ) %>% 
          summarise( mean_count = mean(Count, na.rm=T) ) %>% 
          ungroup %>% 
          # collapse taxonomic information and spatial rep value
          mutate( taxa_id  = paste(Order, Family, 
                                   Genus, Species, sep='_') ) %>% 
          mutate( spat_rep = paste(Site,Line.,sep='_') ) %>% 
          select( -Site, -Line., -Order, -Family, -Genus, -Species) %>% 
          # ltermetacomm format
          rename( SITE_ID          = spat_rep, 
                  DATE             = Year,
                  VARIABLE_NAME    = taxa_id,
                  VALUE            = mean_count ) %>% 
          mutate( VARIABLE_UNITS   = 'Average_count',
                  OBSERVATION_TYPE = 'TAXON_COUNT') %>% 
          select(OBSERVATION_TYPE,SITE_ID, DATE, VARIABLE_NAME, VARIABLE_UNITS, VALUE) %>% 
          # introduce zeros (if need be!)
          spread(VARIABLE_NAME, VALUE, fill = 0) %>% 
          gather(VARIABLE_NAME, VALUE, -DATE, -SITE_ID, -OBSERVATION_TYPE, -VARIABLE_UNITS)
        

# check replication
ggplot2::ggplot(data = sev, 
                ggplot2::aes(x = DATE,
                             y = SITE_ID)) +
    ggplot2::geom_point(size = 5) +
    ggplot2::theme_bw() + 
    ggplot2::xlab("Year with available data") + 
    ggplot2::ylab("Site")

# write file out
write.csv(sev, 'data/L3-sev-arthropods-compagnoni.csv', row.names=F)