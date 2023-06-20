# --------------------------------------------------------- #
# Format raw data from L0 into L3 format                                                         #
# Revised 15 May 2018 by ERS                                #
# --------------------------------------------------------- #


# Contributors: Riley Andrade, Max Castorani, Nina Lany, Sydne Record, Nicole Voelker
# revised by Eric Sokol


options(stringsAsFactors = FALSE)
library(testthat)

# Clear environment
rm(list = ls())


# ---------------------------------------------------------------------------------------------------
#IMPORT AND FORMAT DATA FROM EDI

# Package ID: knb-lter-jrn.210351002.75 Cataloging System:https://pasta.edirepository.org.
# Data set title: Jornada Experimental Range permanent quadrat chart data beginning 1915 - plant density.
# Data set creator:  William Chapline -  
# Metadata Provider:    - USDA ARS Jornada Experimental Range (JER) 
# Contact:    - Data Manager Jornada Basin LTER  - datamanager@jornada-vmail.nmsu.edu
# Contact:  Darren James -  USDA ARS Jornada Experimental Range (JER)  - 
# Stylesheet for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

#read in observations data table
infile1 <- "https://pasta.lternet.edu/package/data/eml/knb-lter-jrn/210351002/75/10c1bf759b5700581368e64387c2a347" 
infile1 <- sub("^https","http",infile1) 
main_df <- read.csv( infile1,
                     header=F,
                     skip=1,
                     sep=",",
                     quot='"', 
                     col.names=c( "quadrat",     
                                  "year",     
                                  "month",     
                                  "USDA_code",     
                                  "scientific_name",     
                                  "common_name",     
                                  "duration",     
                                  "form",     
                                  "density"), 
                     check.names=TRUE)
 


main_tax <- main_df %>% 
              # fix mistake in the name of a species 
              mutate( scientific_name = gsub('Hoffmanseggia','Hoffmannseggia', scientific_name) ) %>% 
              mutate( scientific_name = replace(scientific_name, 
                                                scientific_name == 'Opuntia santa rita',
                                                'Opuntia santarita') ) %>%
              mutate( scientific_name = replace(scientific_name, 
                                                scientific_name == 'Ammocodon chenopodiodes',
                                                'Ammocodon chenopodioides') ) %>%
              # remove unknown species
              subset( !grepl('Unknown', scientific_name) ) %>% 
              # homogenize how genus-level identification is "flagged" ('spp.' instead of 'species')
              mutate( scientific_name = gsub("species","spp.",scientific_name) ) %>% 
              # fix mistakes in taxon_id
              mutate( USDA_code = replace(USDA_code, USDA_code == 'SPHA', 'SPHAE') ) %>% 
              # fix mistake with Spheraclea. I call 'hastulata' as. 'spp.' because in this genus
              # level ID is predominant (see below)
              mutate( scientific_name = replace( scientific_name, 
                                                 scientific_name == 'Sphaeralcea hastulata',
                                                 'Sphaeralcea spp.') )


# Take means of density by year
# select the columns of interest
jrn <- main_tax %>%
          group_by(quadrat, year, USDA_code, scientific_name) %>% 
          # take average of replicate observations in a given location and time
          summarize(value = mean(density, na.rm = TRUE)) %>%
          ungroup() %>%
          rename(SITE_ID          = quadrat, 
                 DATE             = year,
                 VARIABLE_NAME    = USDA_code,
                 VALUE            = value) %>% 
          mutate(OBSERVATION_TYPE = 'TAXON_COUNT',
                 VARIABLE_UNITS   = 'count (average)' ) %>%
          select(OBSERVATION_TYPE,
                 SITE_ID, 
                 DATE, 
                 VARIABLE_NAME, 
                 VARIABLE_UNITS, 
                 VALUE) 

# 3. resolve issues with genus-level IDs -------------------------------------------------------
  
# get genera where observations were identified at the genus level
genus_ids <- main_tax %>% 
                select( scientific_name ) %>% 
                mutate( scientific_name = gsub("species","spp.",scientific_name) ) %>% 
                subset( grepl('spp\\.|species', scientific_name) ) %>% 
                separate( scientific_name, into = c('genus', 'spp'), sep =" " ) %>% 
                unique

# Test: codes of taxa identified at the genus level correspond to a unique 'scientific_name'
main_tax %>% 
  mutate( scientific_name = gsub("species","spp.",scientific_name) ) %>% 
  subset( grepl('spp.', scientific_name) ) %>% 
  select(USDA_code, scientific_name) %>% 
  unique %>% 
  .$USDA_code %>% 
  table %>% 
  unique %>% 
  expect_equal(1)

# get the genera to lump (those that are identified at the species level AS WELL)
codes_df <- main_tax %>% 
              select( USDA_code, scientific_name ) %>% 
              unique %>% 
              subset( grepl(paste(genus_ids$genus,collapse="|"), 
                            scientific_name) ) %>% 
              separate( scientific_name, into = c('genus', 'species'), sep =" " ) %>% 
              group_by( genus ) %>% 
              summarise( rep=n() ) %>%
              # select IDs identified at genus AND species level
              subset( rep > 1 ) %>% 
              ungroup %>% 
              .$genus %>% 
              lapply(function(x) subset(main_tax, grepl(x,scientific_name) ) )

# calculate the proportion density measurements at the genus-level
prop_gen <- function(x){
  
  # transform taxonomic data in order to merge it to the (formatted) observation data frame
  x <- x %>% 
        rename(VARIABLE_NAME = USDA_code) %>% 
        select(VARIABLE_NAME, scientific_name) %>% 
        unique
  
  # calculate proportions of counts for each ID category
  # use this to calculate proportion of counts IDed to genus rather than species
  # dat1 %>% 
  jrn %>% 
    subset( grepl(paste(x$VARIABLE_NAME,collapse="|"), 
                        VARIABLE_NAME) ) %>% 
    group_by(VARIABLE_NAME) %>% 
    summarise( total = sum(VALUE) ) %>% 
    arrange( total ) %>% 
    mutate( prop = total/sum(total) ) %>% 
    arrange( desc(prop) ) %>% 
    as.data.frame %>% 
    left_join( x )
  
}

# proportion of species identified to genus level
prop_l <- lapply(codes_df, prop_gen) 

# remove genus level IDed that make up less than 5% of counts in a genus
r_gen  <- prop_l %>% 
            lapply( function(x){ 
                      out <- subset(x, grepl('spp\\.', scientific_name))
                      if(out$prop < 0.05) return(T) else return(F)} 
                  ) %>% 
            unlist %>% 
            which

# Lump entities to the genus level 
# because more than 5% of counts within a genus cannot not be identified at the species level
l_gen  <- prop_l %>% 
            lapply( function(x){ 
                      out <- subset(x, grepl('spp\\.', scientific_name))
                      if(out$prop >= 0.05) return(T) else return(F)} 
                  ) %>% 
            unlist %>% 
            which

# 3A. remove IDS of species identified to genus level only.
r_codes<- prop_l[r_gen] %>% 
            # select genus-level data (flagged by "spp.")
            lapply( function(x) subset(x, grepl('spp\\.', scientific_name)) ) %>% 
            # re-combine into a data frame
            Reduce( function(...) rbind(...), .) %>% 
            # get codes to remove
            .$VARIABLE_NAME

# test that you've not lost any identifier
expect_true( identical(1:length(prop_l), c(r_gen,l_gen) %>% sort ) )

# 3B. lump IDs to the genus level
lump_spp <- function(x){
  
  mutate(x, new_var = subset(x, grepl('spp\\.', scientific_name))$VARIABLE_NAME ) %>% 
    select(-total, -prop, -scientific_name)
  
}

# new variables (that lump individuals at the genus level)
newvar_df <- lapply(prop_l[l_gen], lump_spp) %>% 
                Reduce(function(...) rbind(...), .) %>% 
                subset( !(VARIABLE_NAME == new_var) )


# update taxonomic IDs
count_d   <- jrn %>% 
                # 3A. remove IDS of species identified to genus level only. 
                subset( !(VARIABLE_NAME %in% r_codes) ) %>% 
                # 3B. lump IDs to the genus level
                left_join( newvar_df ) %>% 
                # substitute VARIABLE_NAME with new_var only if new_var IS NOT an NA.
                mutate( VARIABLE_NAME = plyr::mapvalues(VARIABLE_NAME,
                                                        from = newvar_df$VARIABLE_NAME,
                                                        to   = newvar_df$new_var )

                        ) %>% 
                dplyr::select( -new_var ) %>% 
                
                # take means AGAIN: lumped taxonomy need be 
                group_by(OBSERVATION_TYPE, SITE_ID, DATE, VARIABLE_NAME, VARIABLE_UNITS) %>% 
                summarise( VALUE = mean(VALUE, na.rm=T) ) %>% 
                ungroup() %>% 
  
                # introduce zeros
                spread(key = VARIABLE_NAME, value = VALUE, fill = 0) %>% 
                gather(key = VARIABLE_NAME, value = VALUE, -DATE, -SITE_ID, -OBSERVATION_TYPE, -VARIABLE_UNITS) %>% 
                
                # 4. select sites and year BY HAND
                subset( SITE_ID %in% c( 'a1', 
                                        'a2',
                                        'a4',
                                        'a2',
                                        'b4',
                                        'b5',
                                        'g4',
                                        'g5',
                                        'i1',
                                        'i2',
                                        'i4') ) %>% 
                subset( DATE < 1939 ) %>% 
                subset( !(DATE %in% c(1918,1922,1926,1929,1930,1934)) ) 

# OVERKILL: test that you did substitude characters correctly
test_post  <- jrn %>% 
                # 3A. remove IDS of species identified to genus level only. 
                subset( !(VARIABLE_NAME %in% r_codes) ) %>% 
                # 3B. lump IDs to the genus level
                left_join( newvar_df )  %>% 
                # substitute VARIABLE_NAME with new_var only if new_var IS NOT an NA.
                mutate( VARIABLE_NAME = plyr::mapvalues(VARIABLE_NAME,
                                                        from = newvar_df$VARIABLE_NAME,
                                                        to   = newvar_df$new_var )

                        ) %>% 
                select( -new_var )
test_pre <- subset(jrn,!(VARIABLE_NAME %in% r_codes) )

# perform the actual tests
expect_equal( nrow(test_pre), nrow(test_post) )
for(ii in 1:nrow(newvar_df)){
  expect_true( (which(test_post$VARIABLE_NAME == newvar_df$VARIABLE_NAME[ii]) %in% 
                  which(test_pre$VARIABLE_NAME == newvar_df$new_var[ii]) ) %>% all )
}

# combine data into one dataframe
out <- Reduce(function(...) rbind(...), 
              list(count_d) )

write.csv(out, file = "data/L3-jrn-plants-compagnoni.csv", row.names = F)    
