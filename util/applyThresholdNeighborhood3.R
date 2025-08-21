library(SubCellBarCode)
library(caret)
library(e1071)
library(stats)
library(dplyr)
library(RColorBrewer)


all.repA <- all3.A 
all.repB <- all3.B 
all.repC <- all3.C
threshold.df <- t.n.df3


applyThresholdNeighborhood3 <- function(all.repA, all.repB, all.repC, threshold.df){
  
  #upgrade compartment labels to neighborhood labels for prediction
  all.n.repA <- SubCellBarCode::replacePrediction(df = all.repA,
                                                  column = "svm.pred.all")
  all.n.repB <- SubCellBarCode::replacePrediction(df = all.repB,
                                                  column = "svm.pred.all")
  all.n.repC <- SubCellBarCode::replacePrediction(df = all.repC,
                                                  column = "svm.pred.all")
  #sum up compartment level predictions to neighborhood predictions
  m.all.repA <- SubCellBarCode::mergeProbability(all.n.repA)
  m.all.repB <- SubCellBarCode::mergeProbability(all.n.repB)
  m.all.repC <- SubCellBarCode::mergeProbability(all.n.repC)
  
  # same row order for every replicate 
  m.all.repB <- m.all.repB[rownames(m.all.repA), ]
  m.all.repC <- m.all.repC[rownames(m.all.repA), ]
  
  #merge the replicates by averaging the probabilities, the unclassified list
  all.repAB.mean <- rbind(m.all.repA, m.all.repB, m.all.repC)[,-2]
  
  # aggregate the measurement per how many protein (same gene) are repeated in the df
  all.repAB.mean <- aggregate(.~Proteins, data = all.repAB.mean, mean)
  
  # create an unclassified variable for all
  all.repAB <- data.frame(Proteins = all.repAB.mean$Proteins,
                          svm.pred.all = rep("Unclassified", nrow(all.repAB.mean)),
                          all.repAB.mean[,2:5])
  rownames(all.repAB) <- all.repAB$Proteins
  all.repAB <- all.repAB[rownames(m.all.repA), ]
  
  
  # common neighborhood prediction across all replicates - in boolean so i can subset 
  common_pred <-cbind(m.all.repA$svm.pred.all, m.all.repB$svm.pred.all, m.all.repC$svm.pred.all)
  common_bool <- apply(common_pred, 1, function(r){all(r==r[1])})
  
  # merge all replicates and average them using those that have same classification in each replicate
  n.repA.m <- m.all.repA[common_bool, ]
  n.repB.m <- m.all.repB[common_bool, ]
  n.repC.m <- m.all.repC[common_bool, ]
  
  
  # those with the same neighborhood classification level are avergaed across the 3 replicates
  # so the existence of more replicates will not decrese the unclassified since the aggregation
  # happens after comparing the classification neighborhood in each replicate.
  
  combined.reps <- rbind(n.repA.m, n.repB.m, n.repC.m)
  combined.df <- data.frame(Proteins = combined.reps$Proteins,
                            combined.reps[, 3:6])
  averaged.reps <- aggregate(.~Proteins, data = combined.df, mean)
  rownames(averaged.reps) <- averaged.reps$Proteins
  averaged.reps <- averaged.reps[rownames(n.repA.m), ]
  
  combined.rep.A.B <- data.frame(Proteins = averaged.reps$Proteins,
                                 svm.pred.all = n.repA.m$svm.pred.all, # predictions are the same with B and C here
                                 averaged.reps[, 2:5])
  
  #apply the thresholds for each neighborhood, generated from calculate threshold function 
  neighborhoods <- c("Secretory","Nuclear", "Cytosol", "Mitochondria")
  
  confident.classification <- lapply(neighborhoods, function(m){
    # temp precision
    t.p <- unname(unlist(threshold.df[threshold.df$Neighborhood == m, ][2]))
    #temp recall
    t.r <- unname(unlist(threshold.df[threshold.df$Neighborhood == m, ][3]))
    t.value <- max(t.p, t.r)
    if (is.numeric(t.p)){
      temp.df <- combined.rep.A.B[combined.rep.A.B$svm.pred.all == m, ]
      up.threshold.df <- temp.df[temp.df[m] >= t.value, ]
    }
  }) # adds the neighborhood if and onle if the averaged probability > threshold
  
  conf.df <- do.call("rbind", confident.classification)
  
  ##adding "unclassified proteins" those that had different neighborhood and those below the threshold
  no.class <- subset(all.repAB, ! rownames(all.repAB) %in% rownames(conf.df))
  
  n.cls.df <- rbind(conf.df, no.class)
  
}
