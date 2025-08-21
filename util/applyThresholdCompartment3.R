library(SubCellBarCode)
library(caret)
library(e1071)
library(stats)
library(dplyr)
library(RColorBrewer)


applyThresholdCompartment3 <- function(all.repA, all.repB, all.repC, threshold.df){
  
  if( ! identical(rownames(all.repA), rownames(all.repB)) )
    stop('Replicates have to possess same indexed rownames')
  
  if( ! identical(rownames(all.repA), rownames(all.repC)) )
    stop('Replicates have to possess same indexed rownames')
  
  if( ! identical(rownames(all.repB), rownames(all.repC)) )
    stop('Replicates have to possess same indexed rownames')
  
  all.repA$Proteins <- rownames(all.repA)
  all.repB$Proteins <- rownames(all.repB)
  all.repC$Proteins <- rownames(all.repC)
  
  #merge the replicates by averaging the probabilities
  
  # This part is going to be used for the unclassified proteins 
  repAB.mean <- rbind(all.repA, all.repB, all.repC)[,-1]
  repAB.mean <- aggregate(.~Proteins, data = repAB.mean, mean)
  repAB <- data.frame(Proteins = repAB.mean$Proteins,
                      svm.pred = rep("Unclassified", nrow(repAB.mean)),
                      repAB.mean[,seq_len(16)])
  rownames(repAB) <- repAB$Proteins
  repAB <- repAB[rownames(all.repA), ]
  
  common_prot <- cbind(all.repA$svm.pred.all, all.repB$svm.pred.all, all.repC$svm.pred.all)
  
  # This part is used for the classified that belong to all the replicates
  all.repA.match <- all.repA[apply(common_prot, 1, function(r){all(r == r[1])}),]
  all.repB.match <- all.repB[apply(common_prot, 1, function(r){all(r == r[1])}),]
  all.repC.match <- all.repC[apply(common_prot, 1, function(r){all(r == r[1])}),]
  
  combined.repAB <- rbind(all.repA.match, all.repB.match, all.repC.match)
  combined.df <- data.frame(Proteins = combined.repAB$Proteins,
                            combined.repAB[, 2:16])
  averaged.reps <- aggregate(.~Proteins, data = combined.df, mean)
  rownames(averaged.reps) <- averaged.reps$Proteins
  averaged.reps <- averaged.reps[rownames(all.repA.match), ]
  averaged.repAB <- data.frame(Proteins = averaged.reps$Proteins,
                               svm.pred = all.repA.match$svm.pred.all,
                               averaged.reps[, 2:16])
  
  
  #apply the theresholds for each compartment
  compartments <- c("S1", "S2", "S3", "S4", "N1", "N2", "N3", "N4",
                    "C1", "C2", "C3", "C4", "C5", "M1", "M2")
  
  confident.classification <- lapply(compartments, function(m){
    #temp precision
    t.p <- unname(unlist(threshold.df[threshold.df$Compartment == m, ][2]))
    #temp recall
    t.r <- unname(unlist(threshold.df[threshold.df$Compartment == m, ][3]))
    t.value <- max(t.p, t.r, na.rm = TRUE)
    if(is.numeric(t.value)){
      temp.df <- averaged.repAB[averaged.repAB$svm.pred == m, ]
      up.threshold.df <- temp.df[temp.df[m] >= t.value, ]
    }
  })
  
  confident.df <- do.call("rbind", confident.classification)
  confident.df <- confident.df <- confident.df[complete.cases(confident.df),]
  
  # adding "unclassified proteins"
  no.clss <- subset(repAB, ! rownames(repAB) %in% rownames(confident.df))[-3]
  
  c.cls.df <- rbind(confident.df, no.clss)
  
}
