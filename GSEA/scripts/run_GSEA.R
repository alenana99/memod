# run_gsea.R

#Load file
list_S <- readRDS("data/list_S_BP")
pmsm_count_df <- read.table("data/pmsm_count_df")

#ES_significance2 function
ES_significance2 <- function(S,L,nbootstrap=1000){
  ES_ <- ES_wrapper_fast(S,L)
  S_ <- S
  ES_pi <- sapply(1:nbootstrap,function(x){ x$i <- sort(sample(1:dim(L)[1])[1:dim(x)[1]])
  ES_wrapper_fast(x,L)})
  pos_av <- mean(ES_pi[ES_pi>0])
  neg_av <- mean(ES_pi[ES_pi<0])
  sig <- sum(abs(ES_pi[sign(ES_pi)==sign(ES_)])>abs(ES_)) / nbootstrap
  NES_pi <- c(
    ES_pi[ES_pi>0]/pos_av,
    ES_pi[ES_pi<0]/neg_av)
  if(ES_ < 0){NES <- ES_ / neg_av}
  else{NES <- ES_ / pos_av}
  return(list(ES_ = ES_, sig = sig, NES = NES, NES_pi = NES_pi))
}

#ES_wrapper_fast function -> calculate Enrichment Score
ES_wrapper_fast <- function(S,L){
  Pmiss <- 0
  latest_i <- 1
  latest <- S$i[latest_i] # S$i = match(S$g,L$g)
  memoise <- 0
  last_i <- S$i[length(S$i)]
  Nr <- sum(L[S$i,2]) 
  N <- dim(L)[1]
  Nh <- dim(S)[1]
  
  ES_final <- 0
  
  for (i in 1:length(L$r)){
    #print(c(i,latest))
    if( i < latest){
      Pmiss <- Pmiss + (1/(N - Nh))
      Phit <- memoise
    }
    else {
      #print(c(i,Phit,Pmiss,Nr))
      Phit <- memoise + (L[latest,2]/Nr)
      latest_i <- min(latest_i + 1, dim(S)[1])
      latest <- S$i[latest_i]
      memoise <- Phit 
      #print(c(i,Phit,Pmiss,Nr))
    }
    #print(c(i,Phit,Pmiss))
    ES <- Phit - Pmiss
    ES_final <- max(abs(ES_final),abs(ES))
    #if(i==last_i){print(i)}
     }
return(ES_final)
}

#GSEA_fdr function
GSEA_fdr <- function(x,NES,NES_pi)
  # x is the NES value for a given term
  # NES is the distribution of NES for all terms
  # NES_pi is the distribution of NES for all permutation of all terms
  NES_pi = NES_pi[sign(NES_pi) == sign(x)]
  NES <- NES[sign(NES) == sign(x)]
  FDR <- (sum(abs(x) >= abs(NES_pi))/(length(NES_pi)))/
          (sum(abs(x) >= abs(NES))/(length(NES)))
  return(FDR)

#FDR_wrapper function
FDR_wrapper <- function(NES,NES_pi){
  FDRs <- sapply(NES,function(x){GSEA_fdr(x,NES,NES_pi)})
  return(FDRs)
}

# Calculate ES for each gene set
res_ES <- lapply(list_S, function(x) ES_significance2(x, pmsm_count_df))

# Extract NES and NES_pi
NES_values <- sapply(res_ES, function(x) x$NES)
NES_pi_values <- do.call(c, lapply(res_ES, function(x) x$NES_pi))             

# Calculate FDR                                  
FDR_values <- FDR_wrapper(NES_values, NES_pi_values)

# Save the result
write.table(FDR_values, file = "data/FDR_values", row.names = FALSE)                                   
