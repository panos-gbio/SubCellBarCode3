# load Libraries
library(SubCellBarCode)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(stats)

# # SCBC Parameters from Vignette
# protein.data <- dfr3
# markerProteins <- r.markers3
# figure_path <- figs_path
# dims = 3
# theta = c(0.1, 0.2)
# perplexity = c(50, 60)




tsneVisualization3 <- function(protein.data,
                              markerProteins, dims,
                              theta,
                              perplexity,
                              figure_path){
  
  # subset dataset by protein markers that passed the QC
  tsne.df <- protein.data[markerProteins, ]
  annotation.df <- SubCellBarCode::markerProteins[rownames(tsne.df), ] # here we extract the data from the 5CL
  
  if( ! identical(rownames(tsne.df), rownames(annotation.df)) )
    stop('Make sure your markerProteins are covered in 3365 Proteins and
            there is no extra entry.')
  
  if(dims == 3){
    
    message("Optimization process started. This may take some time!")
    
    lst1 <- unlist(lapply(theta, function(x)
      lapply(perplexity, function(y) sort(x + y))))
    lst2 <- unlist(lapply(theta, function(x)
      lapply(perplexity, function(y){
        tsne.mpa <- Rtsne::Rtsne(tsne.df, dims=3, theta=x, perplexity=y)
        tsne.mpa$itercosts[1]})))
    
    lst.df <- data.frame(Parameters = lst1, Cost = lst2)
    lst.df <- lst.df[order(lst.df$Cost, decreasing = FALSE), ]
    
    min.theta.perp <- unlist(strsplit(as.character(lst.df[1, 1]),
                                      split = ".", fixed = TRUE))
    
    message("Optimization was performed.")
    
    #get the optimization parameters
    theta.val <- as.numeric(paste(0,
                                  as.character(min.theta.perp[2]), sep = "."))
    
    perplexity.val <- as.numeric(min.theta.perp[1])
    
    cat("Theta value: ", theta.val)
    cat("\nPerplexity value: ", perplexity.val)
    
    
    rtsne.map <- Rtsne::Rtsne(tsne.df,
                              dims = 3,
                              theta= theta.val,
                              perplexity = perplexity.val)
    
    
    #plot 3D tsne-map and save using base R mechanisms 
    # Open a PDF device
    pdf(paste0(figure_path,"tsne_3D.pdf"), width = 12, height = 5, pointsize = 10)  # adjust size as needed
    op <- par(mfrow = c(1, 3), mar = c(4, 4, 1, 1))  # 1 row, 3 columns
    on.exit(par(op))
    
    d <- data.frame(x=rtsne.map$Y[, 1],
                    y=rtsne.map$Y[, 2],
                    z=rtsne.map$Y[, 3])
    
    scatterplot3d::scatterplot3d(d$y, d$x, d$z,
                                 pch=20,
                                 color=annotation.df$Colour,
                                 grid = TRUE,
                                 box = FALSE)
    # legend("topright",
    #        legend = c("S1", "S2", "S3", "S4","N1", "N2", "N3", "N4",
    #                   "C1", "C2", "C3", "C4", "C5", "M1", "M2"),
    #        pch = 16,
    #        cex = 0.75,
    #        col = c("gold", "orange", "salmon", "tomato2",
    #                "grey90","grey70", "grey50", "grey30",
    #                "lightblue", "aquamarine", "cyan", "deepskyblue2",
    #                "turquoise3", "burlywood4", "tan4"))

    scatterplot3d::scatterplot3d(d$x, d$y, d$z,
                                 pch=20,
                                 color=annotation.df$Colour,
                                 grid = TRUE,
                                 box = FALSE)
    # legend("topright",
    #        legend = c("S1", "S2", "S3", "S4","N1", "N2", "N3", "N4",
    #                   "C1", "C2", "C3", "C4", "C5", "M1", "M2"),
    #        pch = 16,
    #        cex = 0.75,
    #        col = c("gold", "orange", "salmon", "tomato2",
    #                "grey90", "grey70", "grey50", "grey30",
    #                "lightblue","aquamarine", "cyan","deepskyblue2",
    #                "turquoise3","burlywood4","tan4"))
    # 
    scatterplot3d::scatterplot3d(d$z, d$y, d$x,
                                 pch=20,
                                 color = annotation.df$Colour,
                                 grid = TRUE,
                                 box = FALSE)
    legend("topright",
           legend = c("S1", "S2", "S3", "S4","N1", "N2", "N3", "N4",
                      "C1", "C2", "C3", "C4", "C5", "M1", "M2"),
           pch = 16,
           cex = 1,
           pt.cex = 2,
           bty = "n",
           title = "Compartments",
           col = c("gold", "orange", "salmon", "tomato2",
                   "grey90", "grey70", "grey50", "grey30",
                   "lightblue", "aquamarine", "cyan", "deepskyblue2",
                   "turquoise3", "burlywood4", "tan4"))
    
    
    par(op)   # restore settings
    dev.off() # close PDF device and save file
    
    tsneMap.df <- data.frame(Proteins = rownames(tsne.df), d[, seq_len(3)])
    
  }else if(dims == 2){
    
    message("Optimization process started. This may take some time!")
    
    lst1 <- unlist(lapply(theta,
                          function(x) lapply(perplexity, function(y) sort(x + y))))
    lst2 <- unlist(lapply(theta,
                          function(x) lapply(perplexity, function(y){
                            tsne.mpa <- Rtsne::Rtsne(tsne.df, dims = 2, theta=x, perplexity=y)
                            tsne.mpa$itercosts[1]})))
    
    lst.df <- data.frame(Parameters = lst1, Cost = lst2)
    lst.df <- lst.df[order(lst.df$Cost, decreasing = FALSE), ]
    
    
    min.theta.perp <- unlist(strsplit(as.character(lst.df[1, 1]),
                                      split = ".", fixed = TRUE))
    
    message("Optimization was performed.")
    
    #get the optimization parameters
    theta.val <- as.numeric(paste(0,
                                  as.character(min.theta.perp[2]), sep = "."))
    
    perplexity.val <- as.numeric(min.theta.perp[1])
    
    cat("Theta value: ", theta.val)
    cat("\nPerplexity value: ", perplexity.val)
    rtsne.map <- Rtsne::Rtsne(tsne.df,
                              dims = 2,
                              theta = theta.val,
                              perplexity = perplexity.val)
    
    #plot 2D tsne-map
    d <- data.frame(x=rtsne.map$Y[, 1],
                    y=rtsne.map$Y[, 2])
    
    plot(d$x, d$y, col = annotation.df$Colour, type = "p" , pch = 20)
    legend("topright",
           inset=c(-0.04,0),
           legend = c("S1", "S2", "S3", "S4", "N1", "N2", "N3", "N4",
                      "C1", "C2", "C3", "C4", "C5", "M1", "M2"),
           pch = 16,
           cex = 0.75,
           xpd = TRUE,
           col = c("gold", "orange", "salmon", "tomato2",
                   "grey90", "grey70", "grey50", "grey30",
                   "lightblue", "aquamarine", "cyan", "deepskyblue2",
                   "turquoise3", "burlywood4", "tan4"))
    
    tsneMap.df <- data.frame(Proteins = rownames(tsne.df), d[, seq_len(2)])
    
  } else{
    
    stop("Plesae select the correct dimension. It is either 2 or 3.")
    
  }
}