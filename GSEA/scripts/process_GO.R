#!/usr/bin/env Rscript

# Definizione della funzione per recuperare le informazioni dei GO IDs
fetch_go_details <- function(go_ids) {
          # Funzione helper per recuperare le informazioni di un singolo GO ID
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

  # Applica la funzione helper a tutti i GO IDs e raccoglie i risultati
  go_info_list <- lapply(go_ids, get_go_info)

    # Converti l'elenco di liste in un dataframe
    go_info_df <- do.call(rbind, lapply(go_info_list, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
    return(go_info_df)
}

# Funzione per processare i dati GO
process_go_data <- function(anno_path, L_path, output_path) {
  library(AnnotationDbi)
  library(GO.db)
  
  # Read files
  anno.df <- read.table(anno_path, header = TRUE, sep = "\t", quote = "", colClasses = c(goid = "character"))
  anno.df$goid <- paste("GO:", anno.df$goid, sep = "")
  
  L <- read.table(L_path, header = FALSE, sep = "\t")
  names(L) <- c("g", "r")
  
  # Recupera i dettagli dei GO
  anno.df <- cbind(anno.df$qpid, fetch_go_details(anno.df$goid))
  names(anno.df)[names(anno.df) == "anno.df$qpid"] <- "gene_ID"
  
  # Crea una lista per memorizzare S per ogni categoria
  S_ontology <- list()
  for (ontology in unique(anno.df$Ontology)) {
    S_ontology[[ontology]] <- data.frame(
      g = subset(anno.df, Ontology == ontology)$gene_ID,
      GO_ID = subset(anno.df, Ontology == ontology)$GO_ID
    )
  }
  
  # Stampa i data frame per controllo
  #print(S_ontology)
  
  S_GO <- lapply(S_ontology, function(df) {
    df$i <- match(df$g, L$g)
    merge_df <- merge(df, L, by = "g")
    list(merge_df = merge_df)
  })
  
  # Visualizza i risultati
  #print(S_GO)
  
  unique_go_ids <- lapply(S_GO, function(item) {
    unique(item$merge_df$GO_ID)
  })
  
  list_S <- list()
  for (ontology in names(S_GO)) {
    for (go_id in unique_go_ids[[ontology]]) {
      ss <- S_GO[[ontology]]$merge_df[S_GO[[ontology]]$merge_df$GO_ID == go_id, c("g", "i", "r")]
      ss <- ss[order(ss$r, decreasing = TRUE), ]
      ss <- ss[order(ss$i), ]
      ss <- ss[, c("g", "i")]
      list_S[[ontology]][[paste0("", gsub(":", "_", go_id))]] <- ss
    }
  }

 #return(list_S)
    # Salva list_S in un file
    saveRDS(list_S, file = output_path)
}

args <- commandArgs(trailingOnly = TRUE)
anno_path <- args[1]
L_path <- args[2]
output_path <- args[3]

# Esegui la funzione con i percorsi specificati
list_S <- process_go_data(anno_path, L_path, output_path)
