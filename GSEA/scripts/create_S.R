#!/usr/bin/env Rscript

# create_S.R
# Load libraries
if (!require("optparse")) {
  install.packages("optparse", repos='http://cran.us.r-project.org')
  library(optparse)
}

option_list = list(
  make_option(c("-p", "--pmsm_count_file"), type="character", default=NULL, 
              help="Count methylation per gene", metavar="character"),
  make_option(c("-S", "--S_file"), type="character", default=NULL, 
              help="File BP", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default="list_S_BP.rds", 
              help="Output file [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Verifica che tutti i file di input siano specificati
if (is.null(opt$pmsm_count_file) | is.null(opt$bp_file)) {
  print_help(opt_parser)
  stop("Both input files must be specified", call.=FALSE)
}


pmsm_count_df <- read.table("data/pmsm_count_df", header=TRUE)
BP <- read.table("data/BP")
S_BP$i <- match(S_BP$g, pmsm_count_df$g)
S_BP_GO <- data.frame(BP$gene_ID, BP$GO_ID)
names(S_BP_GO) <- c("g", "GO_ID")
S_BP_GO$i <- match(S_BP_GO$g, pmsm_count_df$g)
S_BP_GO <- merge(S_BP_GO, pmsm_count_df, by = "g")
unique_GO_BP <- unique(S_BP_GO$GO_ID)

# Create list_S
list_S <- list()
for (go_id in unique_GO_BP) {
  ss <- S_BP_GO[S_BP_GO$GO_ID == go_id, c("g", "i", "n_meth")]
  ss <- ss[order(ss$n_meth, decreasing = TRUE), ]
  ss <- ss[order(ss$i), ]
  ss <- ss[, c("g", "i")]
  list_S[[paste0("", gsub(":", "_", go_id))]] <- ss
}

# Save list_S
saveRDS(list_S, file = "data/list_S_BP.rds")
