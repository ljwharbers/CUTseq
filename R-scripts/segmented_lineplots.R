#libraries
source("E:/SciLife Lab/save_and_plot.R")
require(data.table)
require(rtracklayer)
require(naturalsort)

#set sample name
samplename <- "Figure_4"

#load in readcount data (only robot exp)
#counts <- fread("E:/SciLife Lab/CUTseq/Data/BICRO167/readcounts_final.txt")

#set+create output directory
output <- paste0("E:/SciLife Lab/CUTseq/Revision2/Plots/profileplots/", samplename, "/")
if(!dir.exists(output)){
  dir.create(output, recursive = TRUE)}

#setwd
setwd("E:/SciLife Lab/CUTseq/Data/JohanBRCA_tsv+bed/tsvfiles/")

#load and sort segmented files
files <- list.files("./", full.names = TRUE, recursive = FALSE)
files <- naturalsort(files)
#counts <- counts[order(counts[, 2])]

########
#subsetting
annotation <- fread("E:/SciLife Lab/CUTseq/Revision2/Data/Figure4_selectionNewlabels.txt")

#############
# #extra selection
x <- numeric(0)

for(i in 1:nrow(annotation)){
  x[i] <- which(grepl(annotation$`Library ID`[i], files) & grepl(annotation$Barcode[i], files))
}
x <- sort(x)

files <- files[x]
annotation$runOrder <- as.numeric(gsub("XZ", "", annotation$`Library ID`))
setorder(annotation, runOrder, Barcode)

# #subsetting for 300K+ reads
# files <- files[counts$V3 >= 300000]
# counts <- counts[counts$V3 >= 300000]

#for robot experiment
#load in one file to preload breaks/switch data
tmp <- fread(files[1])

#get x-axis tick locations, switch locations
ticks <- numeric(0)
for(k in 1:22){
  ticks[k] <- sum(tmp$chromosome == k) / 2
}

#get switch point locations
chrSwitch <- c(0, which(tmp$chromosome != dplyr::lag(tmp$chromosome)))

#set labels
labels <- c("1", rep("", 3), "5", rep("", 5), "11", rep("", 8), "20", "", "")
rm(tmp)


plotlist <- list()
counter <- 1
for(i in files){
  #load in data
  data <- fread(i)
  
  #change name and factor levels for forced order
  names(data)[5] <- "value"
  data$feature <- factor(data$feature, levels = data$feature)
  
  #data lower or higher than (-)2, set to (-)2
  data$value[data$value >= 2] <- 2
  data$value[data$value <= -2] <- -2
  
  
  i <- gsub("\\./", "", i)
  i <- gsub("\\.bam\\.tsv", "", i)
  
  # #only for the downsampled SKBR3
  # i <- gsub("XZ174BC158", "", i)
  # i <- gsub("\\.perc-0.50", "22M", i)
  # i <- gsub("\\.perc-0.25", "11M", i)
  # i <- gsub("\\.perc-0.10", "4.3M", i)
  # i <- gsub("\\.perc-0.05", "2.2M", i)
  # i <- gsub("\\.perc-0.01", "400K", i)
  # i <- gsub("\\.perc-0.005", "200K", i)

  
  # #optional extra gsub for 96 well libraries
  i <- gsub("\\.dedup.*", "", i)
  #i <- paste0(counts[counter, 1], " ", counts[counter, 5])
  
  i <- paste0("KI_", annotation[counter, 6])

  #plotting
  plotlist[[counter]] <- ggplot(data, aes(x = feature, y = value, group = 1)) +
                                geom_path() +
                                geom_vline(xintercept = chrSwitch, linetype = "dashed", size = 0.0000001) +
                                ggtitle(i) +
                                theme_minimal() +
                                scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 2)) +
                                scale_x_discrete(breaks = data$feature[chrSwitch + ticks], labels = labels) +
                                theme(text = element_text(size=11),
                                      #axis.text.x = element_blank(),
                                      axis.title.x = element_blank(),
                                      axis.title.y = element_blank(),
                                      line = element_blank(),
                                      axis.line = element_line(),
                                      plot.background = element_rect(fill = "transparent", color = NA))
  
  print(paste0("sample ", counter, ": Done"))
  counter <- counter + 1
}

page <- 1
for(k in seq(1, length(plotlist), 33)){
  #save as eps
  #postscript(file=paste0(output, samplename, "_page", page, "_postscript.eps"),
   #          onefile = TRUE, paper="special", height=11, width=8, pointsize=8, horizontal=FALSE)
  png(file=paste0(output, samplename, "_page", page, ".png"), height=11, width=8, units = "in", res = 400)
  #check if length is not more than plotlist and generate grid.arrange
  if(k + 32 <= length(plotlist)){
    plots <- grid.arrange(grobs = plotlist[k: (k+32)], ncol = 3, nrow = 11)
    dev.off()
  }else{
    plots <- grid.arrange(grobs = plotlist[k: length(plotlist)], ncol = 3, nrow = 11)
    dev.off()
  }
  print(paste0("page ", page, ": Done"))
  page <- page + 1
}

#visually check
# page <- 1
# for(k in seq(1, length(plotlist), 33)){
#   #check if length is not more than plotlist and generate grid.arrange
#   if(k + 32 <= length(plotlist)){
#     grid.arrange(grobs = plotlist[k: (k+32)], ncol = 3, nrow = 11)
# 
#   }else{
#     grid.arrange(grobs = plotlist[k: length(plotlist)], ncol = 3, nrow = 11)
# 
#   }
#   print(paste0("page ", page, ": Done"))
#   page <- page + 1
# }


