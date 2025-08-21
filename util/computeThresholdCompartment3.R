library(SubCellBarCode)
library(caret)
library(e1071)
library(stats)
library(dplyr)
library(RColorBrewer)


computeThresholdCompartment3 <- function(test.repA, test.repB, test.repC){
  
  # The protein columns must be on the same order in all three replicates
  if( ! identical(rownames(test.repA), rownames(test.repB)) )
    stop('Replicates A - B have to possess same indexed rownames')
  
  if(! identical(rownames(test.repB), rownames(test.repC)))
    stop('Replicates B - C have to possess same indexed rownames')
  
  if(! identical(rownames(test.repA), rownames(test.repC)))
    stop('Replicates A - C have to possess same indexed rownames')
  
  # Create dataframe from all replicates - keeps the factor numbering and not the string literal 
  common <- cbind(test.repA$svm.pred, test.repB$svm.pred, test.repC$svm.pred)
  
  #Check for the common classifications and merge three replicates and average them
  
  #check if all the entries in the row are the same as the first one. 
  test.repA.match <- test.repA[apply(common, 1, function(row){all(row == row[1])}), ]
  test.repA.match$Proteins <- rownames(test.repA.match)
  
  test.repB.match <- test.repB[apply(common, 1, function(row){all(row == row[1])}), ]
  test.repB.match$Proteins <- rownames(test.repB.match)
  
  test.repC.match <- test.repC[apply(common, 1, function(row){all(row == row[1])}), ]
  test.repC.match$Proteins <- rownames(test.repC.match)
  
  # Combine and create the dataframe with same columns and entries as the original ones
  combined.reps <- rbind(test.repA.match, test.repB.match, test.repC.match)
  combined.df <- data.frame(Proteins = combined.reps$Proteins,
                            combined.reps[, 3:17])
  averaged.reps <- aggregate(.~Proteins, data = combined.df, mean)
  rownames(averaged.reps) <- averaged.reps$Proteins
  averaged.reps <- averaged.reps[rownames(test.repA.match),]
  
  combined.rep.A.B <- data.frame(Proteins = averaged.reps$Proteins,
                                 Observation = test.repA.match$Observation,
                                 svm.pred = test.repA.match$svm.pred,
                                 averaged.reps[, 2:16])
  
  
  
  #estimate the compartments thresholds  by calculating precision and recall
  cls.levels <- c("S1", "S2", "S3", "S4", "N1", "N2", "N3", "N4",
                  "C1", "C2", "C3", "C4", "C5", "M1", "M2")
  
  results <- lapply(cls.levels,function(l){
    
    # subset based on compartment eg. "S1", both the combined probabilities dataframe and observed classification 
    cls.df <- combined.rep.A.B[combined.rep.A.B$svm.pred==l, ]
    cls.obs.df <- test.repA[test.repA$Observation==l, ]
    
    # create the parameters list, where we test probability thresholds for each eg. "S1" compartment svm classifications from 0 to 1 in 0.005 incremental  
    if(nrow(cls.df) > 0 & nrow(cls.obs.df) > 0){
      parameters <- lapply(seq(0, 1, 0.005), function(t){
        u.df <- cls.df[cls.df[l] >= t, ]
        p.cls <- sum(u.df$Observation == u.df$svm.pred) / nrow(u.df)
        cls.down.df <- cls.df[cls.df[l] < t, ]
        r.cls <- (sum(u.df$Observation == u.df$svm.pred)) / nrow(cls.obs.df)
        f.score <- (2 * p.cls * r.cls) / (p.cls + r.cls)
        values <- list(Precision = p.cls,
                       Recall = r.cls,
                       fscore = f.score,
                       threshold = t,
                       Compartment = l)
      })
      
      result.df <- data.frame(do.call(rbind.data.frame, parameters))
      up.precision <- result.df[result.df$Precision>= 0.9, ]
      recall.threshold <- result.df[!duplicated(result.df[,c('Recall')]),]
      threshold.df <- list(Compartment = l,
                           Precision = up.precision$threshold[1],
                           Recall = recall.threshold$threshold[2])
    }})
  
  threshold.compartment.df <- data.frame(do.call(rbind, results))
  
  missing.comp <- setdiff(cls.levels, threshold.compartment.df$Compartment)
  
  if(length(missing.comp) > 1){
    miss.df <- data.frame(Compartment = missing.comp,
                          Precision = rep("NA", length(missing.comp)),
                          Recall = rep("NA", length(missing.comp)))
    threshold.compartment.df <- rbind(threshold.compartment.df, miss.df)
    comps <- threshold.compartment.df$Compartment
    rownames(threshold.compartment.df) <- comps
    threshold.compartment.df <- threshold.compartment.df[cls.levels,]
    
  }else if(length(missing.comp) == 1){
    miss.df <- data.frame(Compartment = missing.comp,
                          Precision = "NA",
                          Recall = "NA")
    threshold.compartment.df <- rbind(threshold.compartment.df, miss.df)
    comps <- threshold.compartment.df$Compartment
    rownames(threshold.compartment.df) <- comps
    threshold.compartment.df <- threshold.compartment.df[cls.levels,]
  }else{
    threshold.compartment.df <- threshold.compartment.df
  }
  colnames(threshold.compartment.df)[2:3] <- c("PrecisionBasedThreshold",
                                               "RecallBasedThreshold")
  threshold.compartment.df$OptedThreshold <- apply(threshold.compartment.df,
                                                   1, function(x) max(unlist(x[2]), unlist(x[3]), na.rm = TRUE))
  
  return(threshold.compartment.df)
}
