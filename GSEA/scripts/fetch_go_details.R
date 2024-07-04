#!/usr/bin/env Rscript

#import library
library(GO.db)
library(optparse)

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


# Load GO IDs
GO_ids <- read.table("data/GO_ids", header = FALSE)
GO_ids <- as.vector(as.matrix(GO_ids))

# Obtain details on GO terms
GO_description <- fetch_go_details(GO_ids)

# Save the result
write.csv(GO_description, file = "data/GO_description.csv", row.names = FALSE)
