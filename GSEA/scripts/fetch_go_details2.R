#!/usr/bin/env Rscript

# Import library
library(GO.db)
library(optparse)

# Define command line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "GO_ids",
              help = "Path to the input file containing GO IDs [default: %default]"),
  make_option(c("-o", "--output"), type = "character", default = "GO_description.csv",
              help = "Path to the output CSV file [default: %default]")
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Function to retrieve GO IDs information
fetch_go_details <- function(go_ids) {
  get_go_info <- function(id) {
    if (!id %in% keys(GOTERM)) {
      return(list(GO_ID = id, Ontology = NA, Name = NA, Definition = NA))
    } else {
      ontology <- as.character(Ontology(GOTERM[[id]]))
      name <- as.character(Term(GOTERM[[id]]))
      definition <- as.character(Definition(GOTTERM[[id]]))
      return(list(GO_ID = id, Ontology = ontology, Name = name, Definition = definition))
    }
  }
  go_info_list <- lapply(go_ids, get_go_info)
  go_info_df <- do.call(rbind, lapply(go_info_list, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
  return(go_info_df)
}

# Load GO IDs
GO_ids <- read.table(opt$input, header = FALSE)
GO_ids <- as.vector(as.matrix(GO_ids))

# Obtain details on GO terms
GO_description <- fetch_go_details(GO_ids)

# Save the result
write.csv(GO_description, file = opt$output, row.names = FALSE)
