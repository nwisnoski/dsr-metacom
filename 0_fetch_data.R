# Fetch all data from EDI and create local instances of the data in the /data directory

## finds all .R and .r files within a folder and sources them
source_folder <- function(folder, recursive = FALSE, ...) 
{ 
  files <- list.files(folder, pattern = "[.][rR]$", 
                      full.names = TRUE, recursive = recursive)
  if (!length(files))
    stop(simpleError(sprintf('No R files in folder "%s"', folder)))
  src <- invisible(lapply(files, source, ...))
  message(sprintf('%s files sourced from folder "%s"', length(src), folder))
}


source_folder("format_data/")