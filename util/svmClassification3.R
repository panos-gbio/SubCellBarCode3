library(SubCellBarCode)
library(caret)
library(e1071)
library(stats)
library(dplyr)
library(RColorBrewer)

# returns a list with 3 elements, one per replicate and each replicate contains: 
# df coontaining classification probabilities at compartment level for test data
# df coontaining classification probabilities at compertment level for whole data
# model parameters
# 
# 
# markerProteins <- r.markers3
# protein.data <- dfr3
# markerprot.df <- SubCellBarCode::markerProteins


svmClassification3 <- function(markerProteins, protein.data, markerprot.df){
  
  # subset marker proteins that passed QC and get MS signals 
  prot.df <- protein.data[markerProteins, ]
  
  # get the color for each marker-protein using the 5CL refrence df 
  marker.df <- markerprot.df[markerProteins, ] 
  
  # data from experiment and ground truth from marker proteins. 
  # Subset with those after quality control, r.markers, this is the training data 
  subcell.loc.data <- cbind(prot.df, marker.df[seq_len(2)])
  
  #splitting the data randomly using caret function 
  label.train <- caret::createDataPartition(subcell.loc.data$Compartments,
                                            times = 1,
                                            p=0.7,
                                            list = FALSE)
  
  
  df.train <- subcell.loc.data[label.train, ]
  df.test <- subcell.loc.data[-label.train, ]
  
  
  # Build classifier for  replicate A, B, and C respectively. 
  # Choose replicate name and use it as a placeholder for selecting columns in the second function
  # it loops three time, a replicate each time, generating classification probabilities per replicate 
  build.classifier <- lapply(c("A", "B", "C"), function(x) {
    rep.train.df <- df.train[, grepl(sprintf("\\.%s\\.", x), # regex
                                     colnames(df.train))]
    replicate.train.label <- df.train$Compartments
    
    # search for character %s, which is a placeholder for x 
    rep.test.df <- df.test[, grepl(sprintf("\\.%s\\.", x),
                                   colnames(df.test))]
    replicate.test.label <- df.test$Compartments
    
    # tune by grid-search the best parameters of RBF-SVM, the model uses factors as labels 
    svm.tune <- e1071::tune(svm,
                            train.x = rep.train.df,
                            train.y = factor(replicate.train.label),
                            kernel="radial",
                            ranges=list(cost = 10^(seq_len(4) - 2),
                                        gamma = c(.5, 1, 1.5, 2)))
    
    
    # here we fit the SVM model using the grid-search CV hyperparameters and using whole training set
    model.svm <- e1071::svm(factor(replicate.train.label) ~ .,
                            data = rep.train.df,
                            probability=TRUE,
                            kernel="radial",
                            cost=svm.tune$best.parameters[1],
                            gamma=svm.tune$best.parameters[2],
                            scale = FALSE, # data is scaled already 
                            cross=10)
    
    # Evaluate on the held-out test set  
    svm.pred <- predict(model.svm, rep.test.df, probability = TRUE)
    svm.pred.test.df <- data.frame(svm.pred)
    svm.test.prob <- data.frame(attr(svm.pred, "probabilities")) # predicted propabilities for each class
    test.observation <- data.frame(df.test$Compartments) # actual classes
    svm.test.prob.out <- cbind(test.observation,
                               svm.pred.test.df,
                               svm.test.prob)
    colnames(svm.test.prob.out)[1] <- 'Observation'
    svm.test.prob.out <- svm.test.prob.out[,
                                           c("Observation", "svm.pred", "S1", "S2",
                                             "S3", "S4", "N1", "N2", "N3", "N4", "C1",
                                             "C2", "C3", "C4", "C5", "M1", "M2")]
    # accuracy check 
    if(sum(svm.test.prob.out$Observation ==
           svm.test.prob.out$svm.pred) / nrow(df.test) < 0.50){
      
      cat("Overall prediction accuracy is < % 50.
                        Downstream analysis will not be accurate enough.
                        We highly recommend you to perform wet-lab analyis
                        again.")
    }
    
    # predict all proteins for this replicate
    all.prot.df <- protein.data[,grepl(sprintf("\\.%s\\.", x),
                                       colnames(df.train))]
    svm.pred.all <- predict(model.svm,
                            all.prot.df,
                            probability = TRUE)
    
    svm.all.label <- data.frame(svm.pred.all)
    svm.all.prob <- data.frame(attr(svm.pred.all, "probabilities"))
    all.prot.pred <- cbind(svm.all.label, svm.all.prob)
    all.prot.pred <- all.prot.pred[,c("svm.pred.all", "S1", "S2",
                                      "S3", "S4", "N1", "N2", "N3", "N4", "C1",
                                      "C2", "C3", "C4", "C5", "M1", "M2")]
    
    message("Extracting model ", paste0(x))
    
    # merge both test predictions and all protein predictions + model parameters 
    all.classifications <- list(svm.test.prob.out=svm.test.prob.out,
                                all.prot.pred = all.prot.pred,
                                model = model.svm)
    
    return(all.classifications)

    
  })
}
