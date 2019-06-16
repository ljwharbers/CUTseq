#Calculating spearman coefficient between segmented files

require(data.table)
require(viridis)
require(forcats)
source("E:/SciLife Lab/save_and_plot.R")

#setwd
setwd("E:/SciLife Lab/CUTseq/Revision2/")

#load in data
files <- list.files("Data/SKBR3_serialDilution/50kb/tsvfiles/", full.names = T, recursive = F)
files1 <- list.files("../Data/SKBR3 deepseq/forCor/50kb/", full.names = T, recursive = F)

files <- c(files, files1)
########### 
# #ROBOT ONLY
# counts <- fread('../Data/BICRO155/readcounts_final.txt', select = c(1:2, 5))
# 
# #formatting
# counts$`Reads after UMI filtering` <- as.numeric(gsub(",", "", counts$`Reads after UMI filtering`))
# counts <- counts[order(counts$Barcode)]
# 
# # #subsetting for 300K+ reads
# files <- files[counts$`Reads after UMI filtering` >= 300000]
# 
# files1 <- list.files("../Data/BICRO167/tsvfiles/100Kbp/", full.names = T, recursive = F)
# counts <- fread('../Data/BICRO167/readcounts_final.txt', select = c(1:2, 5))
# 
# #formatting
# counts$`Reads after UMI filtering` <- as.numeric(gsub(",", "", counts$`Reads after UMI filtering`))
# counts <- counts[order(counts$Barcode)]
# 
# # #subsetting for 300K+ reads
# files1 <- files1[counts$`Reads after UMI filtering` >= 300000]
# 
# files <- c(files, files1)
###########
data <- fread(files[1], select = 1:4)

for(i in files){
  x <- fread(i, select = 5)
  setnames(x, gsub(".*\\/|.tsv|\\.dedup.*", "", i))
  #setnames(x, gsub(".*\\/", "2_", i))
  data <- cbind(data, x)
}

#correlate columns
toCor <- data[, 5:ncol(data)]
res <- cor(toCor, method = "spearman")
res <- melt(res)

#plot matrix
plt <- ggplot(res, aes(x = Var1 , y = fct_rev(Var2), fill = value)) +
  geom_tile() +
  #geom_text(aes(label = format(value, digits = 3))) +
  scale_y_discrete("") + 
  scale_x_discrete("", position = "top") +
  scale_fill_viridis("Spearman's\ncorrelation", option="inferno", begin = 0.75, direction = -1) +
  theme_minimal() +
  theme(axis.text.x.top = element_text(angle = 90))
  
save_and_plot(plt, "./Plots/SKBR3_serialdilution_correlations/spearman_50kb_noLabel", height = 8, width = 8)
