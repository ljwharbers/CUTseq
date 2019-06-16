#libraries
source("E:/SciLife Lab/save_and_plot.R")
require(data.table)
require(rtracklayer)
require(naturalsort)

#set sample name
samplename <- "Figure_4_corrected"

#set+create output directory
output <- paste0("E:/SciLife Lab/CUTseq/Revision2/Plots/profileplots_website/", samplename, "/")
if(!dir.exists(output)){
  dir.create(output, recursive = T)
  dir.create(paste0(output, "/png/"), recursive = T)
  dir.create(paste0(output, "/eps/"), recursive = T)
}

#setwd into directory of bed+tsv files
setwd("E:/SciLife Lab/CUTseq/Data/JohanBRCA_tsv+bed/")

#load and sort segmented files
files.tsv <- list.files("./tsvfiles/", pattern = ".tsv$", full.names = TRUE, recursive = FALSE)
files.tsv <- naturalsort(files.tsv)
files.bed <- list.files("./bedfiles/", pattern = ".bed$", full.names = TRUE, recursive = FALSE)
files.bed <- naturalsort(files.bed)

if(!all.equal(gsub(".tsv", "", files.tsv), gsub(".bed", "", files.bed))) 
  stop("input files are not equal")

#load in possible extra annotation files/selection etc
#annotation <- fread("E:/SciLife Lab/CUTseq/Data/BICRO167/readcounts_final.txt")
#annotation <- annotation[order(annotation$Barcode)]
annotation <- fread("E:/SciLife Lab/CUTseq/Revision2/Data/Figure4_selectionNewlabels_corrected.txt")

#############
# #extra selection
x <- numeric(0)

for(i in 1:nrow(annotation)){
  #print(paste0(i,"_", which(grepl(annotation$`Library ID`[i], files.tsv) & grepl(annotation$Barcode[i], files.tsv))))
  x[i] <- which(grepl(annotation$`Library ID`[i], files.tsv) & grepl(annotation$Barcode[i], files.tsv))
}
x <- sort(x)

files.tsv <- files.tsv[x]
files.bed <- files.bed[x]
annotation$runOrder <- as.numeric(gsub("XZ", "", annotation$`Library ID`))
setorder(annotation, runOrder, Barcode)

#############
#load in full bins 
bins <- fread("E:/SciLife Lab/References/hg19_100kbBins_nochr.bed")
setnames(bins, 1:3, c("chromosome", "start", "end"))

#add 1 for consistency between datasets
bins$start <- bins$start + 1

#remove X+Y
bins <- subset(bins, grepl("[0-9]", bins$chromosome))
bins$chromosome <- as.numeric(bins$chromosome)

#setorder
setorder(bins, chromosome, start)

#get x-axis tick locations, switch locations
ticks <- numeric(0)
for(k in 1:22){
  ticks[k] <- nrow(bins[bins$chromosome == k]) / 2
}

#get switch point locations
chrSwitch <- c(0, which(bins$chromosome != dplyr::lag(bins$chromosome)))
labels <- as.character(1:22)

setnames(bins, 1:3, c("chromosome", "start", "end"))
bins$feature <- paste0(bins$chromosome, ":", bins$start, "-", bins$end)

#check if file lengths are equal
if(!length(files.tsv) == length(files.bed)) stop("list of tsv and bed files are not equal in length")

#loop through files
for(i in 1:length(files.tsv)){
  
  #load in data
  data.tsv <- fread(files.tsv[i])
  data.bed <- fread(files.bed[i], header = F, select = 1:5)
  
  setnames(data.tsv, 5, "value")
  setnames(data.bed, c("chromosome", "start", "end", "feature", "value"))
  
  to.add <- subset(bins, !bins$feature %in% data.tsv$feature)
  to.add[, value := NA]
  
  data.tsv <- rbind(data.tsv, to.add)
  data.bed <- rbind(data.bed, to.add)
  
  #make chr numeric
  data.tsv$chromosome <- as.numeric(data.tsv$chromosome)
  data.bed$chromosome <- as.numeric(data.bed$chromosome)
  
  #order
  setorder(data.tsv, chromosome, start)
  setorder(data.bed, chromosome, start)
  
  #change name and factor levels for forced order
  data.tsv$feature <- factor(data.tsv$feature, levels = data.tsv$feature)
  
  data.bed <- data.bed[, c(4, 1:3, 5)]
  data.bed$feature <- factor(data.bed$feature, levels = data.bed$feature)
  
  #data lower or higher than (-)4, set to (-)4
  data.tsv$value[data.tsv$value >= 4] <- 4
  data.tsv$value[data.tsv$value <= -4] <- -4
  
  #set called column for fill color
  data.tsv[data.tsv$value >= log2(2.5/2), color := "red"]
  data.tsv[data.tsv$value <= log2(1.5/2), color := "blue"]
  data.tsv[(data.tsv$value > log2(1.5/2)) & (data.tsv$value < log2(2.5/2)) , color := "black"]
  
  #get sample name
  name <- files.tsv[i]
  name <- gsub(".*/", "", name)
  name <- gsub("\\.bam|\\.tsv", "", name)

  #for robots
  #name <- paste0(annotation[i, 1], "_", annotation[i, 5])
  
  #for johanBRCA
  name <- paste0("KI_", annotation[i, 6])
  
  #plot and save segmented line plot
  plt <-
    ggplot(data.tsv, aes(x = feature, y = value, group = 1, color = color)) +
    geom_point(data = data.bed, aes(x = feature, y = value, group = 1), color = "black", size = 0.1) +
    geom_point(size = 0.8) +
    geom_vline(xintercept = chrSwitch) +
    ggtitle(name) +
    theme_minimal() +
    scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 2)) +
    scale_x_discrete(drop = F, breaks = bins$feature[chrSwitch + ticks], labels = labels) +
    scale_color_identity() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          line = element_blank(),
          axis.line.y = element_line(),
          plot.background = element_rect(fill = "transparent", color = NA))
  
  png(file = paste0(output, "png/", name, ".png"), height = 6, width = 12, units = "in", res = 300)
  print(plt)
  dev.off()
  
  postscript(file = paste0(output, "eps/", name, "_postscript.eps"), height = 6, width = 12, onefile = TRUE, paper="special", pointsize=8, horizontal=TRUE)
  print(plt)
  dev.off()
  
  print(paste0("Sample ", i, " out of ", length(files.tsv), " done."))
}


