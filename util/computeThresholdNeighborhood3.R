library(SubCellBarCode)
library(caret)
library(e1071)
library(stats)
library(dplyr)
library(RColorBrewer)



computeThresholdNeighborhood3 <- function(test.repA, test.repB, test.repC){
  
  #upgrade compartment labels to neighborhood labels for prediction
  replaceObervation <- SubCellBarCode::replacePrediction(df = test.repA,
                                                         column = "Observation")
  neighborhood.repA <- SubCellBarCode::replacePrediction(df=replaceObervation,
                                                         column = "svm.pred")
  
  replaceObervation <- SubCellBarCode::replacePrediction(df = test.repB,
                                                         column = "Observation")
  neighborhood.repB <- SubCellBarCode::replacePrediction(df=replaceObervation,
                                                         column = "svm.pred")
  
  replaceObervation <- SubCellBarCode::replacePrediction(df = test.repC,
                                                         column = "Observation")
  neighborhood.repC <- SubCellBarCode::replacePrediction(df = replaceObervation,
                                                         column = "svm.pred")
  
  # concatenate the probabilities for corresponding neighborhood
  
  sum.repA <- SubCellBarCode::sumProbability(neighborhood.repA)
  sum.repB <- SubCellBarCode::sumProbability(neighborhood.repB)
  sum.repC <- SubCellBarCode::sumProbability(neighborhood.repC)
  
  # order of rownames must be the same for all replicates 
  sum.repB <- sum.repB[rownames(sum.repA), ]
  sum.repC <- sum.repC[rownames(sum.repA), ]
  
  # # sanity check
  # identical(rownames(sum.repA), rownames(sum.repB))
  # identical(rownames(sum.repB), rownames(sum.repC))
  # identical(rownames(sum.repA), rownames(sum.repC))
  
  #common prediction across replicates
  common_pred <- cbind(sum.repA$svm.pred, sum.repB$svm.pred, sum.repC$svm.pred)
  common_bool <- apply(common_pred, 1, function(r){all(r==r[1])})
  
  #merge replicates and average them based on neighborhood classification
  n.test.repA.match <- sum.repA[common_bool, ]
  n.test.repA.match$Proteins <- rownames(n.test.repA.match)
  n.test.repB.match <- sum.repB[common_bool, ]
  n.test.repB.match$Proteins <- rownames(n.test.repB.match)
  n.test.repC.match <- sum.repC[common_bool, ]
  n.test.repC.match$Proteins <- rownames(n.test.repC.match)
  
  
  combined.reps <- rbind(n.test.repA.match, n.test.repB.match, n.test.repC.match)
  combined.df <- data.frame(Proteins = combined.reps$Proteins,
                            combined.reps[, 4:7])
  averaged.reps <- aggregate(.~Proteins, data = combined.df, mean)
  rownames(averaged.reps) <- averaged.reps$Proteins
  averaged.reps <- averaged.reps[rownames(n.test.repA.match), ]
  
  combined.rep.A.B <- data.frame(Proteins = averaged.reps$Proteins,
                                 Observation = n.test.repA.match$Observation,
                                 svm.pred = n.test.repA.match$svm.pred,
                                 averaged.reps[, 2:5])
  
  #estimate the neihborhood theresholds by calculating precision and recall
  cls.levels <- c("Secretory", "Nuclear", "Cytosol", "Mitochondria")
  
  results <- lapply(cls.levels, function(l){
    cls.df <- combined.rep.A.B[combined.rep.A.B$svm.pred == l, ]
    cls.obs.df <- neighborhood.repA[neighborhood.repA$Observation == l, ]
    parameters <- lapply(seq(0, 1, 0.005), function(t){
      u.df <- cls.df[cls.df[l] >= t, ]
      p.cls <- sum(u.df$Observation == u.df$svm.pred)/nrow(u.df)
      cls.down.df <- cls.df[cls.df[l] < t,]
      r.cls <- (sum(u.df$Observation == u.df$svm.pred))/nrow(cls.obs.df)
      f.score <- (2 * p.cls * r.cls) / (p.cls + r.cls)
      values <- list(Precision = p.cls,
                     Recall = r.cls,
                     fscore= f.score,
                     threshold = t,
                     Compartment = l)
      
    })
    result.df <- data.frame(do.call(rbind.data.frame, parameters))
    up.precision <- result.df[result.df$Precision >= 0.95, ]
    recall.threshold <- result.df[!duplicated(result.df[, c('Recall')]), ]
    threshold.df <- list(Neighborhood = l,
                         Precision = up.precision$threshold[1],
                         Recall = recall.threshold$threshold[2])
  })
  
  threshold.neighborhood.df <- data.frame(do.call(rbind, results))
  
  colnames(threshold.neighborhood.df)[2:3] <- c("PrecisionBasedThreshold",
                                                "RecallBasedThreshold")
  threshold.neighborhood.df$OptedThreshold <- apply(threshold.neighborhood.df,
                                                    1, function(x) max(unlist(x[2]), unlist(x[3])))
  return(threshold.neighborhood.df)
  
  
}
