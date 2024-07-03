# process_go_terms.R

# Upload GO term information
GO_description <- read.csv("data/GO_description.csv")

# Load gene-GO associations
geneID_goID <- read.csv("data/geneID_goID", sep="\t")
names(geneID_goID) <- c("gene_ID", "GO_ID")

# Merge GO information with genes
geneID_goID_merged <- merge(geneID_goID, GO_description, by = "GO_ID")
geneID_goID_merged <- unique(geneID_goID_merged)

# Divide by category
BP <- subset(geneID_goID_merged, Ontology == "BP")
MF <- subset(geneID_goID_merged, Ontology == "MF")
CC <- subset(geneID_goID_merged, Ontology == "CC")

# Save the results
write.csv(BP, file = "data/BP.csv", row.names = FALSE)
write.csv(MF, file = "data/MF.csv", row.names = FALSE)
write.csv(CC, file = "data/CC.csv", row.names = FALSE)
