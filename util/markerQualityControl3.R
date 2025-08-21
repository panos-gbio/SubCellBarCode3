library(SubCellBarCode)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(stats)



# coveredProteins <- c.prots_t
# protein.data=dfr_t
# figure_path <- fig_path

# The function subsets from the protein markers of the original SCBC pipelne the proteins that exist in the dataset
# The calculations are based on proteins without missing values. 
# A replicate wise correlation is performed, where pearson should be > 0.8 for each protein across replicates.
# Then proteins are correlated to the origal 5CL mass-spec measurements and the minimum spearman and pearson are kept
# The threshold here is pearson >0.8 and spearman >0.6 
# The proteins within the 'Markerproteins' list that pass the QC, are the marker-proteins for this specific experiment 
# and will be used to train the SVM.

markerQualityControl3 <- function(coveredProteins, protein.data, figure_path){
  
  # the function returns a vector with the gene names 
  
  # replicate-wise correlation marker QC
  ## Subset proteomic dataset using 5CL protein marker list from SCBC.
  m.prot.df <- protein.data[coveredProteins, ]
  m.prot.df <- m.prot.df[complete.cases(m.prot.df),]
  
  if( nrow(m.prot.df) < 1)
    stop('Make sure your inputs are correct. Type ?markerQualityControl')
  
  # The subsetting of columns is based on string match of the A/B/C replicate name 
  fracs.df.A <- m.prot.df[, grepl("\\.A\\.", colnames(m.prot.df))]
  fracs.df.B <- m.prot.df[, grepl("\\.B\\.", colnames(m.prot.df))]
  fracs.df.C <- m.prot.df[, grepl("\\.C\\.", colnames(m.prot.df))]
  
  # pairwise replicate wise correlations for all protein markers (same protein different replicates)
  replicate.df <- sapply(rownames(fracs.df.A),
                         function (x) {
                           RepAB <- cor(unlist(fracs.df.A[x, ]), unlist(fracs.df.B[x, ]),
                                        method="pearson") %>% as.numeric()
                           
                           RepBC <- cor(unlist(fracs.df.B[x, ]), unlist(fracs.df.C[x, ]),
                                        method="pearson") %>% as.numeric()
                           
                           RepAC <- cor(unlist(fracs.df.A[x, ]), unlist(fracs.df.C[x, ]),
                                        method="pearson") %>% as.numeric()
                           
                           data.frame(RepAB, RepBC, RepAC)
                         }
  ) %>%
    t() %>%
    as.data.frame()
  
  
  # Convert to dataframe with doubles - all the replicate-wise pearson coefficients
  replicate.df <- data.frame(Protein = rownames(replicate.df),
                             repAB = as.vector(unlist(replicate.df$RepAB)),
                             repBC = as.vector(unlist(replicate.df$RepBC)),
                             repAC = as.vector(unlist(replicate.df$RepAC)))
  
  
  
  # Density plots for the protein markers with patchwork 
  p1 <- ggplot(replicate.df, aes(x = repAB)) +
    geom_density(alpha = .7, fill = "deepskyblue") +
    theme_minimal() +
    theme(text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.text.x = element_text(face = "bold", color="black"),
          axis.text.y = element_text(face = "bold", color="black")) +
    geom_vline(xintercept = 0.8, linetype="dashed", color = "red") +
    labs(title = "AB-Wise Marker QC",
         tag = "A",
         y = "Density",
         x = "Pearson Corr.")
  
  p2 <- ggplot(replicate.df, aes(x = repBC)) +
    geom_density(alpha = .7, fill = "deepskyblue") +
    theme_minimal() +
    theme(text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.text.x = element_text(face = "bold", color="black"),
          axis.text.y = element_text(face = "bold", color="black")) +
    geom_vline(xintercept = 0.8, linetype="dashed", color = "red") +
    labs(title = "BC-Wise Marker QC",
         tag = "",
         y = "",
         x = "Pearson Corr.")
  
  p3 <- ggplot(replicate.df, aes(x = repAC)) +
    geom_density(alpha = .7, fill = "deepskyblue") +
    theme_minimal() +
    theme(text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.text.x = element_text(face = "bold", color="black"),
          axis.text.y = element_text(face = "bold", color="black")) +
    geom_vline(xintercept = 0.8, linetype="dashed", color = "red") +
    labs(title = "AC-Wise Marker QC",
         tag = "",
         y = "",
         x = "Pearson Corr.")
  
  # Keep figures together 
  p4 <- p1+p2+p3+plot_annotation(title = "Replicate-wise correlation of protein markers")
  
  ggsave(paste0(figure_path,"Replicate_wiseQC.pdf"), plot = p4, dpi = 600, width = 12, height = 6)
  
  # remove replicate-wise correlations below 0.8 for all three replicates 
  rep.prots <- replicate.df %>%
    dplyr::filter(repAB > 0.8 & repBC > 0.8 & repAC > 0.8)
  
  message("Number of removed replicate-wise proteins: ", length(setdiff(replicate.df$Protein, rep.prots$Protein)))
  
  # sample-wise correlation marker QC: 
  # proteins that passed the replicate wise QC
  prot.names <- rep.prots$Protein
  
  #subset the 5CL protein markers with those that passed the QC
  markerProteins <- SubCellBarCode::markerProteins[prot.names,][,3:7]
  
  #subset the proteomic dataset with those that passed the QC
  m.prot.df <- m.prot.df[prot.names,]
  
  # sanity check
  # which(!rownames(markerProteins) == prot.names)
  
  # function for to calculate sample-wise correlation for each protein,
  
  # Correlation of each replicate vs 5CL abundance profile (which was the average of the values in each fraction)
  # Checks how well the protein-markers of the sample are conrcondant with the 5CL mass-pec signals (benchmarked markers)
  prot.cor <- function(df, marker.df, cor.method = c("spearman", "pearson")){
    unlist(lapply(rownames(m.prot.df), function(x){
      p.cor <- cor(t(df[x,]), t(marker.df[x,]), method = cor.method)
      names(p.cor) <- x
      p.cor
    }))}
  
  pearson.corA <- prot.cor(df = m.prot.df[grepl("\\.A\\.",
                                                colnames(m.prot.df))],
                           marker.df = markerProteins,
                           cor.method = "pearson")
  
  pearson.corB <- prot.cor(df = m.prot.df[grepl("\\.B\\.",
                                                colnames(m.prot.df))],
                           marker.df = markerProteins,
                           cor.method = "pearson")
  
  pearson.corC <- prot.cor(df = m.prot.df[grepl("\\.C\\.",
                                                colnames(m.prot.df))],
                           marker.df = markerProteins,
                           cor.method = "pearson")
  
  
  # Data frame with all the replicate vs 5CL pearsons and minimum of the three
  pearson.cor <- data.frame(RepA = pearson.corA, RepB = pearson.corB, RepC = pearson.corC)
  pearson.cor$MinP.Cor <- apply(pearson.cor, 1, min)
  
  
  # Repeat the same with Spearman correlation 
  spear.corA <- prot.cor(df = m.prot.df[grepl("\\.A\\.",
                                              colnames(m.prot.df))],
                         marker.df = markerProteins,
                         cor.method = "spearman")
  
  spear.corB <- prot.cor(df = m.prot.df[grepl("\\.B\\.",
                                              colnames(m.prot.df))],
                         marker.df = markerProteins,
                         cor.method = "spearman")
  
  spear.corC <- prot.cor(df = m.prot.df[grepl("\\.C\\.",
                                              colnames(m.prot.df))],
                         marker.df = markerProteins,
                         cor.method = "spearman")
  
  spear.cor <- data.frame(RepA = spear.corA, RepB = spear.corA, RepC = spear.corC)
  spear.cor$MinS.cor <- apply(spear.cor, 1, min)
  
  df_cor <- data.frame(Protein = rownames(spear.cor),
                       Pearson = pearson.cor$MinP.Cor,
                       Spearman = spear.cor$MinS.cor)
  
  # Take the colors from the 5CL based on the markers that passed the QC in this run (prot.names)
  cols <- SubCellBarCode::markerProteins[prot.names,][8]
  Color <- cols$Colour
  
  p5 <- ggplot(df_cor, aes(x = Pearson, y = Spearman)) +
    geom_point(colour = Color, size = 2) +
    geom_hline(yintercept = 0.6, linetype="dashed", color = "red") +
    geom_vline(xintercept = 0.8, linetype="dashed", color = "red") +
    labs(title = "Sample-Wise Correlation Marker QC",
         tag = "B",
         y = "Spearman Corr.",
         x = "Pearson Corr.") +
    theme_minimal() +
    theme(text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, color = "black", size = 14),
          axis.text.x = element_text(face = "bold", color="black", size = 10),
          axis.text.y = element_text(face = "bold", color="black", size = 10))
  
  sample.removed.prot <- df_cor[df_cor$Pearson < 0.8 | df_cor$Spearman < 0.599,]
  sample.removed.prot <- as.character(sample.removed.prot$Protein)
  
  ggsave(paste0(figure_path,"Sample_wiseQC_5CL.pdf"), plot = p5, dpi = 600, width = 10, height = 6)
  
  
  message("Number of removed sample-wise proteins: ",
          length(sample.removed.prot))
  
  robustMarkerProteins <- setdiff(prot.names, sample.removed.prot)
  
  message("Number of total removed marker proteins: ",
          length(sample.removed.prot) + length(setdiff(replicate.df$Protein, rep.prots$Protein)))
  p4/p5
  
  
  # check if there is a depletion after marker qc
  compartment.size <- c(358, 351, 252, 174, 192, 121, 231, 198, 242,
                        132, 220, 215, 341, 69, 269)
  
  compartments <- c("S1", "S2", "S3", "S4", "N1", "N2", "N3", "N4",
                    "C1", "C2", "C3", "C4", "C5", "M1", "M2")
  r.marker.df <- SubCellBarCode::markerProteins[robustMarkerProteins, ]
  
  coverageCompWise <- lapply(seq_len(length(compartments)), function(x){
    temp.df <- r.marker.df[r.marker.df$Compartments == compartments[x], ]
    values <- list(Compartments = compartments[x],
                   ProteinCoverage = 100 * ((dim(temp.df)[1]) /compartment.size[x]))
  })
  
  r.cov.df <- as.data.frame(do.call("rbind", coverageCompWise))
  
  non.enriched.loc <- r.cov.df[r.cov.df$ProteinCoverage < 20, ]
  if(nrow(non.enriched.loc) == 1){
    warning("There is not enough enrichment at: ",
            as.character(non.enriched.loc$Compartments),
            "\nWe recommend you to perform the fractionation, again.")
  }else if(nrow(non.enriched.loc) > 1){
    comp <- paste(as.character(non.enriched.loc$Compartments),
                  collapse = ",")
    warning("There are not enough enrichments at: ",
            comp, "\nWe recommend you to perform the fractionation.")
  }
  
  return(robustMarkerProteins)
}
