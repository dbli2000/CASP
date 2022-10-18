# CASP-σ Motif Prediction

# Refactored version of Horia Todor's https://github.com/horiatodor/predictECF
# From the paper Todor et al (PNAS, 2020) paper
# Rewiring the specificity of extracytoplasmic function sigma factors

# Rationale:
# As the CASP-σs appear to be ECF sigma factors,
# we can leverage existing ECF sigma motif predictions.

#Load Data
setwd('/Users/davidli/Documents/GitHub/CASP/De Novo Sigma Motif Predictions')
load("applied_weights_3.RData")
# For sigma 2 domain
load("domain2.RData")
load("minus10.RData")
# For sigma 4 domain
load("domain4.RData")
load("minus35.RData")

# Define helper functions, largely untouched

loglikelihood <- function(actualresidue, protein, dna, row_weights){
  #first make sure that the number of residues = the size of the proteinmatrix
  if (length(actualresidue)!=dim(as.matrix(protein))[2]){return(c(NA,NA,NA,NA))}

  #this piece of code gets all of the things in protein which are the same as the thing of interest
  #it starts with everything in play and then narrows down only the things that match the sequence
  prots_in_play <- which(as.matrix(protein)[,1]==actualresidue[1])

  if (length(actualresidue) >1){
    for (i in 2:length(actualresidue)){
      prots_in_play <- intersect(prots_in_play, which(as.matrix(protein)[,i]==actualresidue[i]))
    }
  }

  #if the combination doesn't exist, then what?
  if (sum(row_weights[prots_in_play]) <= 100 & length(actualresidue)>1){
    return(loglikelihood_separate(actualresidue, protein, dna, row_weights))
  }

  #for log probability only we will optimize separately from the others
  #since this is all there is to this, and since this the the one I really use
  dna_of_interest <- dna[prots_in_play]
  row_weights_of_int <- row_weights[prots_in_play]
  table_dna <- rep(NA,4)
  table_dna[1] <- sum(row_weights_of_int[which(dna_of_interest == "a")])
  table_dna[2] <- sum(row_weights_of_int[which(dna_of_interest == "c")])
  table_dna[3] <- sum(row_weights_of_int[which(dna_of_interest == "g")])
  table_dna[4] <- sum(row_weights_of_int[which(dna_of_interest == "t")])

  result <- log(table_dna/sum(table_dna),2)

  for (i in 1:4){
    if (is.na(result[i])){result[i] <- log(0.01 ,2)}
    if (result[i] < -10000){result[i] <- log(0.01 ,2)}
  }

  return(result)
}

loglikelihood_separate <- function(actualresidue, protein, dna, row_weights = NA){

  temp_results <- NULL
  quantity <- c(0,0)

  for (i in 1:length(actualresidue)){
    res <- loglikelihood(actualresidue[i], protein[,i], dna, row_weights)
    temp_results <- cbind(temp_results, res)
  }

  return(log(rowMeans(2^(temp_results)),2))
}

predict10w <- function(sequence, dnamatrix, proteinmatrix, row_weights = NA){

  motif <- matrix(0, 4,7)
  rownames(motif) <- c('a','c','g','t')
  sequence <- tolower(sequence)

  #Work through bases from -13 to -7
  #Annotations from original sometimes use wrong AA.
  motif[,1] <- loglikelihood(sequence[53], proteinmatrix[,53], dnamatrix[,9],row_weights)
  motif[,2] <- loglikelihood(sequence[53],proteinmatrix[,53], dnamatrix[,10],row_weights)
  motif[,3] <- loglikelihood(c(sequence[12],sequence[48]),
                             cbind(proteinmatrix[,12],proteinmatrix[,48]), dnamatrix[,11],row_weights)
  motif[,4] <- loglikelihood(c(sequence[22],sequence[26]),
                             cbind(proteinmatrix[,22],proteinmatrix[,26]),
                             dnamatrix[,12],row_weights)
  motif[,5] <- loglikelihood(sequence[40],proteinmatrix[,40], dnamatrix[,13],row_weights)
  motif[,6] <- loglikelihood(sequence[45],proteinmatrix[,45], dnamatrix[,14],row_weights)
  motif[,7] <- loglikelihood(sequence[41],proteinmatrix[,41], dnamatrix[,15],row_weights)
  return(motif)
}

predict35w <- function(sequence, dnamatrix, proteinmatrix, row_weights=NA){

  motif <- matrix(0, 4,7)
  rownames(motif) <- c('a','c','g','t')
  sequence <- tolower(sequence)

  #Work through bases from -36 to -30
  #Annotations from original sometimes use wrong AA.
  motif[,1] <- loglikelihood(sequence[43],proteinmatrix[,43], dnamatrix[,6], row_weights)
  motif[,2] <- loglikelihood(sequence[46],proteinmatrix[,46], dnamatrix[,7], row_weights)
  motif[,3] <- loglikelihood(sequence[42:43],proteinmatrix[,42:43], dnamatrix[,8],row_weights)
  motif[,4] <- loglikelihood(sequence[c(38,42)], proteinmatrix[,c(38,42)], dnamatrix[,9],row_weights)
  motif[,5] <- loglikelihood(sequence[c(45,41)], proteinmatrix[,c(45,41)], dnamatrix[,10],row_weights)
  motif[,6] <- loglikelihood(sequence[41],proteinmatrix[,41], dnamatrix[,11],row_weights)
  motif[,7] <- loglikelihood(sequence[41],proteinmatrix[,41], dnamatrix[,12],row_weights)
  return(motif)
}

# De novo motif predictions
#Predict the -10 motif for Di CASP-σ
DiD2seq <- unlist(strsplit("dkedilsefvlqiq--sstssfkgqna----------kfetwayrifqnvwadffrkikn",""))
predict10w(DiD2seq, minus10, domain2, applied_weights_3)

# Predict the -35 motif for Di CASP-σ
DiD4seq <- unlist(strsplit("frkavkdpkkafcyhllvhytwglkqdemaaryemrhdtfrkrlsecrsrfkqyl",""))
predict35w(DiD4seq, minus35, domain4, applied_weights_3)

#Predict the -35 motif for Csb CASP-σ
CsbD4seq <- unlist(strsplit("gegilnaedialikalyegnengctqsemseemglkpdafkqrkcrlikklessgfgkdtlskilmsd",""))
predict35w(CsbD4seq, minus35, domain4, applied_weights_3)
