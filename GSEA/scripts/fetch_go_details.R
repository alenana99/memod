#!/usr/bin/env Rscript

#import library
library(GO.db)

if (!require("optparse")) {
  install.packages("optparse", repos='http://cran.us.r-project.org')
  library(optparse)
}

# Define command line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "data/GO_ids",
              help = "Path to the input file containing GO IDs [default: %default]"),
  make_option(c("-o", "--output"), type = "character", default = "data/GO_description.csv",
              help = "Path to the output CSV file [default: %default]")
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Definition of the function to retrieve GO IDs information
fetch_go_details <- function(go_ids) {
  # Helper function to retrieve the information of a single GO ID
  get_go_info <- function(id) {
    if (!id %in% keys(GOTERM)) {
      return(list(GO_ID = id, Ontology = NA, Name = NA, Definition = NA))
    } else {
      ontology <- as.character(Ontology(GOTERM[[id]]))
      name <- as.character(Term(GOTERM[[id]]))
      definition <- as.character(Definition(GOTERM[[id]]))
      return(list(GO_ID = id, Ontology = ontology, Name = name, Definition = definition))
    }
  }
  
  # Apply the helper function to all GO IDs and collect the results
  go_info_list <- lapply(go_ids, get_go_info)
  
  # Convert the list of lists to a dataframe
  go_info_df <- do.call(rbind, lapply(go_info_list, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
  return(go_info_df)
}


# Load GO IDs from input file
GO_ids <- read.table(opt$input, header = FALSE)
GO_ids <- as.vector(as.matrix(GO_ids))

# Obtain details on GO terms
GO_description <- fetch_go_details(GO_ids)

# Save the result to output file
write.csv(GO_description, file = opt$output, row.names = FALSE)
