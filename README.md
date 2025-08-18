# SubCellBarCode 3 
A repository with source code adapted from Taner Arslan's SCBC bioinformatic pipeline. The functions are based on the SubCellBarCode package of Bioconductor. <br>

## Project Structure 
```
SubCellBarCode3/             # Project main Root
│
├── data/                   # TMT-15plex MS-based proteomics
│   ├── raw
│   ├── processed 
|
├── figs                    # project figures saved here 
├── src                     # source code for the analysis 
├── util                    # helper functions for the SCBC3 pipeline 
|
├── renv                    # Package manager for Rstudio 
│   
```

## Changes 
- The pipeline has been adapted for triplicates TMT-sets (15-plex) and requires the original package. The new functions are saved in the util directory and are imported in the R scripts located at the src. Installations, simulated data, plots and an analysis walkthrough are in the src directory. 
- The changes involve markerprotein correlation analysis, tSNE visualization, and the SVM probabilties calculation across the three subcellular replicates for each of the 15 original clusters of the SCBC2 pipeline.
- A de novo markerprotein analysis for the leukemia cell lines has not been done yet but it is highly recommended. 
- Lets hope it does not crash a thousand times 🤘🤘

## Usage 
- First clone the repository.
- Make an Rstudio project and use the Renv to solve the enviroment.
- Check in the '0.installations' script in src that you can load the libraries. If not (unlikely), install what is missing there. 
- Have fun! 