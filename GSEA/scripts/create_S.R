# create_S.R

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
saveRDS(list_S, file = "data/list_S_BP")
