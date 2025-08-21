library(SubCellBarCode)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(stats)
library(tidyr)

# paths

data_path <- paste0(getwd(), "/data/raw/")
save_path <- paste0(getwd(), "/data/processed/")
fig_path <- paste0(getwd(), "/figs/")

# we will use the in-package dataset and generate triplicate MS-data
df <- loadData(protein.data = hcc827Ctrl)
cat(dim(df));str(df);head(df)



df_t <- df %>%
  mutate(
    FS1.C.HCC827 = rowMeans(dplyr::select(df,FS1.A.HCC827,FS1.B.HCC827)), 
    .after = FS1.B.HCC827
  ) %>%
  mutate(
    FS2.C.HCC827 = rowMeans(dplyr::select(df,FS2.A.HCC827,FS2.B.HCC827)),
    .after = FS2.B.HCC827
  )  %>%
  mutate(
    FP1.C.HCC827 = rowMeans(dplyr::select(df,FP1.A.HCC827,FP1.B.HCC827)),
    .after = FP1.B.HCC827
  ) %>%
  mutate(
    FP2.C.HCC827 = rowMeans(dplyr::select(df,FP2.A.HCC827,FP2.B.HCC827)),
    .after = FP2.B.HCC827
  ) %>%
  mutate(
    FP3.C.HCC827 = rowMeans(dplyr::select(df,FP3.A.HCC827,FP3.B.HCC827)),
    .after = FP3.B.HCC827
  )

# save the simulated triplicate 
write.csv(df_t, file = paste0(save_path, "triplicates.csv"))

# load the triplicate and convert first column  into rownames (happens with readr and read.table package)
df3 <- read_csv(file = paste0(save_path, "triplicates.csv"))
df3 <- df3 %>% column_to_rownames("...1")


##  ------ Normalization and NAs ---------------

# check if data is mean-subtracted across the TMT-channels for all proteins 
if (all(df3 > 0)) {
  print("The mean across all channels of each protein is not subtracted")
} else {
  print("The log intensity is mean subtracted across all the TMT channels")
}


# calculate NAs across columns 
lapply(df3, function(x) {
  c("Missing_Values" = sum(is.na(x)))
}) %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(),
               names_to = "samples",
               values_to = "missing_values") %>% 
ggplot(aes(y=samples,x=missing_values)) + geom_col()


# distribution of NAs across proteins row-by-row
apply(df3, 1, function(x){
  c(count = sum(is.na(x)))
}) %>% 
  as.data.frame() %>% 
  rename(count = ".") %>% 
  ggplot(aes(x=count)) + 
  geom_histogram(binwidth = 1, fill="gold", color="black")


## ---------------
