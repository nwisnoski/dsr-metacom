# ------------------------------------- #
# Make metadata table -2018-03-21       #
# edited by NKL 2023-04-19              #
# ------------------------------------- #

# Clear environment
rm(list = ls())

#load packages
library(xtable)
library(here)

#make a list of files to be appended to the table:
file.list <- list.files(here("supplement/metadata_tables"))

# first create the first row of the table:
metadata_table <- read.csv(here(paste("supplement/metadata_tables/", file.list[1], sep="")))

#then add following rows in a loop:
for (file in file.list[2:length(file.list)]){
	 a <- read.csv(here(paste("supplement/metadata_tables/", file, sep="")))
	 metadata_table <- rbind(metadata_table,a)
}

#save the table into the folder with the Supp Info TeX file:
write.csv(metadata_table, file = here("supplement/metadata_table.csv"), row.names=F)

#also generate LaTeX table that can be pasted in to the doc:
xtable(metadata_table)

