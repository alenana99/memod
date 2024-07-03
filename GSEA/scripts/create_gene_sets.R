# create_gene_sets.R

# Carica i dati delle categorie
BP <- read.csv("data/BP.csv")
MF <- read.csv("data/MF.csv")
CC <- read.csv("data/CC.csv")

# Crea set di geni per ciascuna categoria
S_BP <- data.frame(BP$gene_ID)
names(S_BP) <- c("g")

# Salva i set di geni
write.csv(S_BP, file = "data/S_BP.csv", row.names = FALSE)
