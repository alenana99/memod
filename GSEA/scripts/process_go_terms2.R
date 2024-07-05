#!/usr/bin/env Rscript

# Import library
library(optparse)

# Define command line options
option_list <- list(
  make_option(c("-g", "--go_description"), type = "character", default = "GO_description.csv",
              help = "Path to the input file containing GO description [default: %default]"),
  make_option(c("-i", "--input"), type = "character", default = "GO.out",
              help = "Path to the input file containing GO terms [default: %default]"),
  make_option(c("-o", "--output"), type = "character", default = "data",
              help = "Directory for output files [default: %default]")
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Upload GO term information
GO_description <- read.csv(opt$go_description)

# Load gene-GO associations
GO.out <- read.table(opt$input, header = FALSE)
geneID_goID <- GO.out[, c(1, 3)]
names(geneID_goID) <- c("gene_ID", "GO_ID")

# Merge GO information with genes
geneID_goID_merged <- merge(geneID_goID, GO_description, by = "GO_ID")
geneID_goID_merged <- unique(geneID_goID_merged)

# Divide by category
BP <- subset(geneID_goID_merged, Ontology == "BP")
MF <- subset(geneID_goID_merged, Ontology == "MF")
CC <- subset(geneID_goID_merged, Ontology == "CC")

# Save the results
write.table(BP, file = file.path(opt$output, "BP"), row.names = FALSE)
write.table(MF, file = file.path(opt$output, "MF"), row.names = FALSE)
write.table(CC, file = file.path(opt$output, "CC"), row.names = FALSE)

# Create S for each category
S_BP <- data.frame(BP$gene_ID)
names(S_BP) <- c("g")

# Save S
write.table(S_BP, file = file.path(opt$output, "S_BP"), row.names = FALSE)
