# The MS data should be log transformed, mean subtracted across the TMT-channels 
# and median centered, as stated in nature protocols (Arslan et. al. 2022).

# We will use the same examples used from Bioconductor package in duplicates and
# triplicates. 

library(SubCellBarCode)
library(tidyverse)


# paths 
processed_path <- paste0(getwd(), "/data/processed/")
util_path <- paste0(getwd(), "/util/")
figs_path <- paste0(getwd(), "/figs/")
list.files(util_path)


# helper functions 
source(paste0(util_path,"markerQualityControl3.R")) # markerprotein analysis for triplicates
source(paste0(util_path,"tsneVisualization3.R")) # tsne for triplicates including pdf results 
source(paste0(util_path,"svmClassification3.R"))
source(paste0(util_path,"computeThresholdCompartment3.R"))
source(paste0(util_path,"applyThresholdCompartment3.R"))
source(paste0(util_path,"computeThresholdNeighborhood3.R"))
source(paste0(util_path,"applyThresholdNeighborhood3.R"))


# The data is hcc927ctrl from the SCBC data and triplicates generated from script #1 

df2 <- loadData(protein.data = hcc827Ctrl)
df3 <- read.csv(paste0(processed_path,"triplicates.csv"), header = TRUE) %>%
  column_to_rownames("X")
  

set.seed(2)
dfr2 <- df2[sample(nrow(df2), 6000),]
dfr3 <- df3[sample(nrow(df3), 6000),]



## ------- Step 1. Covered Marker Proteins ---------------

# A calculation of the percentage coverage of the Markerproteins in the specific experiment using
# the 5 CL as reference

c.prots2 <- SubCellBarCode::calculateCoveredProtein(rownames(dfr2), markerproteins = markerProteins[,1])
c.prots3 <- SubCellBarCode::calculateCoveredProtein(rownames(dfr3), markerproteins = markerProteins[,1])




## ------- Step 2. MarkerProtein QC ---------------

# For the triplicates, we can find the results of the analysis in the figure section of this repository 
# We need to provide a path to save the figures, it is a new argument

r.markers2 <- SubCellBarCode::markerQualityControl(coveredProteins = c.prots2, protein.data = dfr2)
r.markers3 <- markerQualityControl3(coveredProteins = c.prots3, protein.data = dfr3,
                                                   figure_path = figs_path)


## ------ Step 3 tSNE coordinates and Visualization -----
# takes a lot of time to run 

# the original version
tsne.map2 <- SubCellBarCode::tsneVisualization(protein.data = dfr2,
                                            markerProteins = r.markers2,
                                            dims = 3,
                                            theta = c(0.01),
                                            perplexity = c(60))

# the triplicate version 
tsne.map3 <- tsneVisualization3(protein.data = dfr3,
                                markerProteins = r.markers3,
                                dims = 3,
                                theta = c(0.01),
                                perplexity = c(60),
                                figure_path = figs_path)




## ---- Step 4 SVM classifier and probability of cluster membership ----



###  -------- Step 4.1 Build model and estimate thresholds  ----------

# In this step the SVM is first tuned using the markerproteins and e071 package
# Then the best paramateters are used to train the model and predict the classification
# probability for all the 15 compartments of eacn non-marker protein. 
# Taner used a multiclassification scheme. 

# there might be minor prob differences for same replicates, same proteins, across seeds due
# different shuffling and initialization of the kernel. Same applies when we use the duplicate
# and triplicate analysis since the data is shuffled differently under the same seed


# duplicates 
set.seed(4)
cls <- SubCellBarCode::svmClassification(markerProteins = r.markers2,
                                         protein.data = dfr2,
                                         markerprot.df = markerProteins)

# extract the compartment probabilities of the test set
test2.A <- cls[[1]]$svm.test.prob.out
test2.B <- cls[[2]]$svm.test.prob.out
head(test2.A)

# extract the compartment probabilities of the whole data
all2.A <- cls[[1]]$all.prot.pred
all2.B <- cls[[2]]$all.prot.pred



# estimate compartment classification thresholds 
t.c.df2 <- SubCellBarCode::computeThresholdCompartment(test.repA = test2.A, test.repB = test2.B)

# for some reasonn C4 has NA as a threshold, so I might add the value myself and see differences
t.c.df2_man <- t.c.df2
t.c.df2_man[t.c.df2_man["Compartment"] == "C4","OptedThreshold"] <- 0.515

# triplicates 
set.seed(4)
cls3 <- svmClassification3(markerProteins = r.markers3,
                           protein.data = dfr3,
                           markerprot.df = SubCellBarCode::markerProteins)

# extract the compartment probabilities of the test set
test3.A <- cls3[[1]]$svm.test.prob.out
test3.B <- cls3[[2]]$svm.test.prob.out
test3.C <- cls3[[3]]$svm.test.prob.out
head(test3.A)

# extract whole data probabilities 
all3.A <- cls3[[1]]$all.prot.pred
all3.B <- cls3[[2]]$all.prot.pred
all3.C <- cls3[[3]]$all.prot.pred


# estimate compartment classification thresholds
t.c.df3 <- computeThresholdCompartment3(test.repA = test3.A,
                                        test.repB = test3.B,
                                        test.repC = test3.C)




### ------ Step 4.2 Apply thresholds to compartment level classificatons  -------


# duplicates
c.cls.df2 <- SubCellBarCode::applyThresholdCompartment(all.repA = all2.A, all.repB = all2.B,
                                      threshold.df = t.c.df2)
c.cls.df2_man <- SubCellBarCode::applyThresholdCompartment(all.repA = all2.A, all.repB = all2.B,
                                                       threshold.df = t.c.df2_man)
table(c.cls.df2$svm.pred)
table(c.cls.df2_man$svm.pred) # nothing changed though... possibly no category there 



# triplicates 
c.cls.df3 <- applyThresholdCompartment3(all.repA = all3.A, all.repB = all3.B,
                                       all.repC = all3.C,
                                       threshold.df = t.c.df3)

table(c.cls.df3$svm.pred) # since different proteins were picked, we observe c4 memberships 


 
### ------ Step 4.3 Compute Neighborhood Classification Thresholds ---------

# duplicates 
t.n.df2 <- SubCellBarCode::computeThresholdNeighborhood(test.repA = test2.A, test.repB = test2.B)

# triplicates 
t.n.df3 <- computeThresholdNeighborhood3(test.repA = test3.A,
                                         test.repB = test3.B,
                                         test.repC = test3.C)


### ------ Step 4.4 Apply threshold to neighborhood level classifications --------

# duplicates 
n.cls.df2 <- SubCellBarCode::applyThresholdNeighborhood(all.repA = all2.A, all.repB = all2.B, 
                                       threshold.df = t.n.df2)

table(n.cls.df2$svm.pred.all)



n.cls.df3 <- applyThresholdNeighborhood3(all.repA = all3.A,
                                        all.repB = all3.B,
                                        all.repC = all3.C,
                                        threshold.df = t.n.df3)
table(n.cls.df3$svm.pred.all)




### ----- Step 4.5 Merge compartment and Neighborhood classification -------
# we will use the in-pipeline functon, no changes are needed here 


# duplicates 
cls.df2 <- SubCellBarCode::mergeCls(compartmentCls = c.cls.df2, neighborhoodCls = n.cls.df2)


# triplicates 
cls.df3 <- SubCellBarCode::mergeCls(compartmentCls = c.cls.df3, neighborhoodCls = n.cls.df3)

