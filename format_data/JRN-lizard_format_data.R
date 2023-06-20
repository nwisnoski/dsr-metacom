# --------------------------------------------------------- #
# Format raw data as a list of tables                       #
# JORNADA lizard data                                                      # Riley Andrade (main contact)
# Revised 09 Jun 2017  & 04 Nov 2018 by Nina Lany                                     #
# --------------------------------------------------------- #

# Contributors: Riley Andrade, Max Castorani, Nina Lany, Sydne Record, Nicole Voelker

# Clear environment
rm(list = ls())

library(tidyverse)

# ---------------------------------------------------------------------------------------------------


# Package ID: knb-lter-jrn.210007001.38 Cataloging System:https://pasta.edirepository.org.
# Data set title: Lizard pitfall trap data from 11 NPP study locations at the Jornada Basin LTER site, 1989-2006.
# Data set creator:  David C Lightfoot - Jornada Basin LTER (now at Univ. of NM) 
# Data set creator:  Walter G Whitford - Jornada Basin LTER/New Mexico State University 
# Metadata Provider:    - Jornada Basin LTER 
# Contact:    - Information Manager Jornada Basin LTER  - jornada.data@nmsu.edu
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-jrn/210007001/38/731f52d77045dfc5957589d35c2e6227" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "date",     
                 "zone",     
                 "site",     
                 "plot",     
                 "pit",     
                 "spp",     
                 "sex",     
                 "rcap",     
                 "toe_num",     
                 "SV_length",     
                 "total_length",     
                 "weight",     
                 "tail",     
                 "pc"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

# attempting to convert dt1$date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1date<-as.Date(dt1$date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1date) == length(tmp1date[!is.na(tmp1date)])){dt1$date <- tmp1date } else {print("Date conversion failed for dt1$date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1date) 
if (class(dt1$zone)!="factor") dt1$zone<- as.factor(dt1$zone)
if (class(dt1$site)!="factor") dt1$site<- as.factor(dt1$site)
if (class(dt1$plot)!="factor") dt1$plot<- as.factor(dt1$plot)
if (class(dt1$pit)=="factor") dt1$pit <-as.numeric(levels(dt1$pit))[as.integer(dt1$pit) ]               
if (class(dt1$pit)=="character") dt1$pit <-as.numeric(dt1$pit)
if (class(dt1$spp)!="factor") dt1$spp<- as.factor(dt1$spp)
if (class(dt1$sex)!="factor") dt1$sex<- as.factor(dt1$sex)
if (class(dt1$rcap)!="factor") dt1$rcap<- as.factor(dt1$rcap)
if (class(dt1$toe_num)=="factor") dt1$toe_num <-as.numeric(levels(dt1$toe_num))[as.integer(dt1$toe_num) ]               
if (class(dt1$toe_num)=="character") dt1$toe_num <-as.numeric(dt1$toe_num)
if (class(dt1$SV_length)=="factor") dt1$SV_length <-as.numeric(levels(dt1$SV_length))[as.integer(dt1$SV_length) ]               
if (class(dt1$SV_length)=="character") dt1$SV_length <-as.numeric(dt1$SV_length)
if (class(dt1$total_length)=="factor") dt1$total_length <-as.numeric(levels(dt1$total_length))[as.integer(dt1$total_length) ]               
if (class(dt1$total_length)=="character") dt1$total_length <-as.numeric(dt1$total_length)
if (class(dt1$weight)=="factor") dt1$weight <-as.numeric(levels(dt1$weight))[as.integer(dt1$weight) ]               
if (class(dt1$weight)=="character") dt1$weight <-as.numeric(dt1$weight)
if (class(dt1$tail)!="factor") dt1$tail<- as.factor(dt1$tail)
if (class(dt1$pc)!="factor") dt1$pc<- as.factor(dt1$pc)

# Convert Missing Values to NA for non-dates

dt1$plot <- as.factor(ifelse((trimws(as.character(dt1$plot))==trimws("NA")),NA,as.character(dt1$plot)))
dt1$pit <- ifelse((trimws(as.character(dt1$pit))==trimws("NA")),NA,dt1$pit)               
suppressWarnings(dt1$pit <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$pit))==as.character(as.numeric("NA"))),NA,dt1$pit))
dt1$sex <- as.factor(ifelse((trimws(as.character(dt1$sex))==trimws("NA")),NA,as.character(dt1$sex)))
dt1$rcap <- as.factor(ifelse((trimws(as.character(dt1$rcap))==trimws("NA")),NA,as.character(dt1$rcap)))
dt1$toe_num <- ifelse((trimws(as.character(dt1$toe_num))==trimws("NA")),NA,dt1$toe_num)               
suppressWarnings(dt1$toe_num <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$toe_num))==as.character(as.numeric("NA"))),NA,dt1$toe_num))
dt1$SV_length <- ifelse((trimws(as.character(dt1$SV_length))==trimws("NA")),NA,dt1$SV_length)               
suppressWarnings(dt1$SV_length <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$SV_length))==as.character(as.numeric("NA"))),NA,dt1$SV_length))
dt1$total_length <- ifelse((trimws(as.character(dt1$total_length))==trimws("NA")),NA,dt1$total_length)               
suppressWarnings(dt1$total_length <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$total_length))==as.character(as.numeric("NA"))),NA,dt1$total_length))
dt1$weight <- ifelse((trimws(as.character(dt1$weight))==trimws("NA")),NA,dt1$weight)               
suppressWarnings(dt1$weight <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$weight))==as.character(as.numeric("NA"))),NA,dt1$weight))
dt1$tail <- as.factor(ifelse((trimws(as.character(dt1$tail))==trimws("NA")),NA,as.character(dt1$tail)))


#in this table, each row represents a lizard. Remove columns with data measured on indivudal lizards (i.e. svl, sex, toe number)
data <- dt1 %>%
	select(date, zone, site, plot, spp,toe_num, pc)

data$datetime <- as.POSIXct(data$date, format = "%m/%d/%Y")
data$year <- as.numeric(format(data$datetime,"%Y"))
str(data)

############################
#re-do the data cleaning that Andrew Hope did (descibed in his metadata sheet)
#remove first and last years due to incomplete sampling
data <- subset(data, data$year > 1989 & data$year < 2006)


#Look at sampling:
ss <- function(x) {length(unique(x))}
tapply(data$site,data$year, ss) #two sites began in 1995
tapply(data$site,data$year, unique) #they were "NORT" and "SUMM"

#remove the two sites that began later, NORT and SUMM
data <- data %>%
	filter(site != "NORT" & site != "SUMM")

#Now remove the rows where no lizards were observed.
data <- data %>%
	filter(spp != "NONE")
	
#Now remove unknown lizards

data <- data %>%
	filter(spp != "UKLI" & spp != "UKCN" & spp != "UKPH")
	

	#NOTE: pc is a problem code. Look it up. 

#assume equal sampling effort (Pitfall traps were opened for 2 week sessions  4x per year (quarterly) at each site; trapping began at sites XXXX and XXXX in 1996) and calculate the number of unique individuals of each species captured per year: 
length(unique(data$spp))  

#make expanded grid of all spp-plot-year combinations
species <- unique(data$spp)
sites <- unique(data$site)
years <- unique(data$year)
grid <- as.data.frame(expand.grid(species, sites, years))
names(grid) <- c("VARIABLE_NAME", "SITE_ID", "DATE")

#sum number of unique individuals counted per spp-site-year:
sum_dat <- data %>%
  group_by(year, site, spp) %>%
  summarize(count = ss(toe_num))%>%
  as.data.frame()

names(sum_dat) <- c("DATE", "SITE_ID", "VARIABLE_NAME","VALUE")
str(sum_dat)

#merge the two together, replace NA with 0, and name according to Erics format:
comm.long <- merge(sum_dat, grid, by = c("VARIABLE_NAME", "SITE_ID", "DATE"), all = T)
comm.long$OBSERVATION_TYPE = rep("TAXON_COUNT", length(comm.long$DATE))
comm.long$VARIABLE_UNITS = rep("count", length(comm.long$DATE))
comm.long$VALUE <- as.numeric(ifelse(is.na(comm.long$VALUE), 0, comm.long$VALUE))

comm.long$VALUE  <- ifelse((comm.long$SITE_ID == "NORT" |comm.long$SITE_ID == "SUMM") & comm.long$DATE < 1995, "NA", comm.long$VALUE)

comm.long$VALUE <- as.numeric(comm.long$VALUE)
comm.long <- na.omit(comm.long)
str(comm.long)

#rbind the L3 dataset
dat.tog <- rbind(comm.long) %>%
  select(OBSERVATION_TYPE, SITE_ID, DATE, VARIABLE_NAME, VARIABLE_UNITS, VALUE)

#Write out the L3 dataset
write.csv(dat.tog, "data/L3-jrn-lizards-hope.csv", row.names=F)

