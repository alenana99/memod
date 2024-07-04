#!/usr/bin/env Rscript

# Import library
library(optparse)

# Define command line options
option_list <- list(
  make_option(c("-p", "--pmsm_count_file"), type = "character", default = "pmsm_count_df",
              help = "Path to the pmsm count file [default: %default]"),
  make_option(c("-c", "--category_file"), type = "character", default = "BP",
              help = "Path to the BP category file [default: %default]"),
  make_option(c("-S", "--S_file"), type = "character", default = "S_BP",
              help = "Path to the S file [default: %default]"),
  make_option(c("-o", "--output"), type = "character", default = "list_S.rds",
              help = "Path to the output list S file [default: %default]")
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load pmsm count data
pmsm_count_df <- read.delim(opt$pmsm_count_file, header = FALSE)
names(pmsm_count_df) <- c("gene_ID", "r")

# Load S data
S_BP <- read.table(opt$S_file, header = TRUE)

#Load category data
BP <- read.table(opt$category_file, header = TRUE)

# Perform operations
S_BP$i <- match(S_BP$g, pmsm_count_df$g)
S_BP_GO <- data.frame(BP$gene_ID, BP$GO_ID)
names(S_BP_GO) <- c("g", "GO_ID")
S_BP_GO$i <- match(S_BP_GO$g, pmsm_count_df$g)
S_BP_GO <- merge(S_BP_GO, pmsm_count_df, by = "g")
unique_GO_BP <- unique(S_BP_GO$GO_ID)

# Create list_S
list_S <- list()
for (go_id in unique_GO_BP) {
  ss <- S_BP_GO[S_BP_GO$GO_ID == go_id, c("g", "i", "r")]
  ss <- ss[order(ss$r, decreasing = TRUE), ]
  ss <- ss[order(ss$i), ]
  ss <- ss[, c("g", "i")]
  list_S[[paste0("", gsub(":", "_", go_id))]] <- ss
}

# Save list_S
saveRDS(list_S, file = opt$output)
