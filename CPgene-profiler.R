old <- Sys.time()

library(dplyr)
library(tidyverse)
library(pheatmap)
library(UpSetR)
library(scales)
library(ape)
library(Biostrings, quietly=T, warn.conflicts=F)
library(reshape2)
library(gridExtra)

############--------CONSTANT VARIABLES------##########

blast_db = "/home/datta/Desktop/CPE/R_package/ARG-annot-20200406T084430Z-001/ARG-annot/ARG-annot_CPGene_DB/ARG-annot_CPGene.DB"
#blast_db = args[1]
blastn = "/usr/bin/blastn" 
evalue = 1e-6
format = "\'6 qseqid sseqid pident nident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\'"
num_threads = 4
fastalocation = "/home/datta/Desktop/CPE/R_package/ARG-annot-20200406T084430Z-001/ARG-annot/test/Assemblies"
#fastalocation = args[2]

############--------CLEANING IF PREVIOUS OUTPUT FILE EXISTS------##########
files.to.delete <- list.files(pattern = "txt|pdf|csv|png|jpg")
file.remove(files.to.delete)

############--------GO TO ASSEMBLY FILES LOCATION------##########

setwd(fastalocation)
getwd()

colnames <- c("assemblyName",
              "qseqid",
              "sseqid",
              "pident",
              "nident",
              "length",
              "mismatch",
              "gapopen",
              "qstart",
              "qend",
              "sstart",
              "send",
              "evalue",
              "bitscore",
              "qlen", 
              "slen")



############--------CLEANING IF PREVIOUS OUTPUT FILE EXISTS------##########
files.to.delete <- list.files(pattern = "txt|pdf|csv|png|jpg|tiff")
file.remove(files.to.delete)

############--------1) BLASTN ASSEMBLY ON ARG-ANNOT DATA------##########

#List of all the fasta files
#https://github.com/garridoo/atsphere_wgs/blob/master/assembly.stats.R

files <- list.files(pattern = "*.fasta$", recursive = F)
files
## Read in all files using a for loop and perform BLASTN and save the results in a output

datalist = list()
for(i in 1:length(files)){
 # files[i]
  min.length <- 0
  seqs <- readDNAStringSet(files[i])
  #seqs <- seqs[width(seqs) >= min.length]
  #names(seqs) <- gsub("^", "contig_", 1:length(seqs))
  names(seqs)
  
  # calculate assembly statistics from contig lengths
  
  lengths.table <- sort(width(seqs))
  lengths <- data.frame(ctg_number=1:length(seqs), acc=cumsum(lengths.table))
  tot.length <- lengths[dim(lengths)[1], 2]
 # tot.length
  
  n50.idx <- which(lengths[, 2] >= tot.length * .50)[1]
  n50 <- lengths.table[n50.idx]
 # n50
  
  n90.idx <- which(lengths[, 2] >= tot.length * .10)[1]
  n90 <- lengths.table[n90.idx]
#  n90
  
  datalist[[i]] <- cbind(files[i],n50,n90,tot.length) 

input = files[i]

  blast_out <- paste(input,system2(command = blastn,
          args = c("-db", blast_db,
                   "-query", input,
                   "-outfmt", format,
                   "-evalue", evalue,
                   "-ungapped",
                   "-num_threads", num_threads),
          wait = TRUE,
          stdout = TRUE ), sep = "\t")
write.table(blast_out,"blastResults.txt",row.names = FALSE,col.names = FALSE, append = TRUE, quote = FALSE) ##OUTPUT1

}


big_data  <- do.call(rbind.data.frame, datalist)
#big_data

big_data$V1 <- as.character(big_data$V1)
big_data$n50 <- as.numeric(as.character(big_data$n50))
big_data$n90 <- as.numeric(as.character(big_data$n90))
big_data$tot.length <- as.numeric(as.character(big_data$tot.length))

names(big_data) <- c("assemblyName","N50","N90","assemblySize")
write.table(big_data,"assemblyStats.txt",row.names = FALSE,col.names = FALSE, append = TRUE, quote = FALSE) ##OUTPUT2

# Color Palette
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

N50_plot <- ggplot(big_data,mapping = aes(x = assemblySize, y = N50)) + geom_point(size=3, colour="#0072B2") + 
  scale_x_continuous(labels = unit_format(unit = "MB", scale = 1e-6), breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(labels = unit_format(unit = "MB", scale = 1e-6), breaks = scales::pretty_breaks(n = 10)) +
  ggtitle("Assembly Size vs N50 plot") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) + 
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

N90_plot <- ggplot(big_data,mapping = aes(x = assemblySize, y = N90)) + geom_point(size=3, colour="#D55E00") + 
  scale_x_continuous(labels = unit_format(unit = "MB", scale = 1e-6), breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(labels = unit_format(unit = "MB", scale = 1e-6), breaks = scales::pretty_breaks(n = 10)) +
  ggtitle("Assembly Size vs N90 plot") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) + 
  theme(plot.margin = unit(c(1,1,1,1), "cm"))


blastResults_df <- read.table("blastResults.txt",sep = "\t", fill = TRUE)
#blastResults_df 
blastResults_df <- blastResults_df[complete.cases(blastResults_df), ]
############--------ADD COLNAMES TO BLASTN OUTPUT TO A TABLE ------##########
names(blastResults_df) <- colnames

############--------FIND CP GENE CONTAINING CONTIGS WITH CPGENE MATCHING 100% ID AND COVERAGE ------##########
blastResults_df$cov <- blastResults_df$nident / blastResults_df$slen
#head(blastResults_df)

#----TRIMMING STRINGS IN THE TABLE
#Samp_CP_df$AssemblyFile <- gsub("_assembly.fasta", "", Samp_CP_df$AssemblyFile)
blastResults_df$sseqid <- gsub("\\(Bla\\)", "", blastResults_df$sseqid)
blastResults_df$sseqid <- gsub("[:_].*", "", blastResults_df$sseqid)
head(blastResults_df)

## 1) Can create a table with AssemblyName | CP contig | Length | CPgene- coordinates from BLAST | Circular/Linear Contig - DONE!

blastResults.filt <- blastResults_df[blastResults_df$cov==1 & blastResults_df$pident==100, c("assemblyName" , "qseqid", "sseqid", "qlen","qstart","qend")] # subsets only the gene name with if CP gene has 100% cov
write.table(blastResults.filt,"blastResults.filt.txt",row.names = FALSE,col.names = FALSE, append = TRUE, quote = FALSE) # ##OUTPUT3 contains AssemblyName | CP contig | Length | CPgene- coordinates from BLAST | Circular/Linear Contig
#blastResults.filt

## 2)  Plot length in chart for each CPgene - Save all the plot in single PDF - DONE!

NDMqlen <- subset(blastResults_df, blastResults_df$cov==1 & blastResults_df$pident==100 & grepl("NDM", sseqid ))[,"qlen", drop=FALSE] 
names(NDMqlen) <- c("NDMqlen")
NDMqlen
table(NDMqlen)
class(NDMqlen)

KPCqlen <- subset(blastResults_df, blastResults_df$cov==1 & blastResults_df$pident==100 & grepl("KPC", sseqid ))[,"qlen", drop=FALSE] 
names(KPCqlen) <- c("KPCqlen")
#KPCqlen

OXAqlen <- subset(blastResults_df, blastResults_df$cov==1 & blastResults_df$pident==100 & grepl("OXA", sseqid ))[,"qlen", drop=FALSE] 
names(OXAqlen) <- c("OXAqlen")
#OXAqlen

##Generate the tables histogram tables for the CPgene lengths

CPcontlen <- function(cl)
{  
#https://stackoverflow.com/a/27889176

min1 <- min(cl)
#min1
max1 <- max(cl) + 10000
#max1
br <- seq(min1,max1, by=10000)
#br
ranges = paste(head(br,-1), br[-1], sep="-")
#ranges

#range(NDMqlen$NDMqlen)

freq   = hist(cl, breaks=br, include.lowest=TRUE, plot=FALSE)
#hist(cl, breaks=br, labels = TRUE, include.lowest=TRUE, plot=TRUE)
data.frame(range = ranges, frequency = freq$counts)

subset(data.frame(range = ranges, frequency = freq$counts), data.frame(range = ranges, frequency = freq$counts)$frequency > 0)
}

N_Hist <- CPcontlen(NDMqlen$NDMqlen)
names(N_Hist) <- c("ContigSize Range","CPContig_Number")
#N_Hist
K_Hist <- CPcontlen(KPCqlen$KPCqlen)
names(K_Hist) <- c("ContigSize Range","CPContig_Number")
#K_Hist
O_Hist <- CPcontlen(OXAqlen$OXAqlen)
names(O_Hist) <- c("ContigSize Range","CPContig_Number")
#O_Hist


line="\n\n========++++++++NDM Contig Size Distribution+++++++++++============\n"
write(line,file="CPContigSizeDist.txt",append=TRUE)
write.table(N_Hist,file="CPContigSizeDist.txt", row.names = FALSE, quote = FALSE, append=TRUE)

line="\n\n========++++++++KPC Contig Size Distribution+++++++++++============\n"
write(line,file="CPContigSizeDist.txt",append=TRUE)
write.table(K_Hist,file="CPContigSizeDist.txt",row.names = FALSE, quote = FALSE, append=TRUE)

line="\n\n========++++++++OXA Contig Size Distribution+++++++++++============\n"
write(line,file="CPContigSizeDist.txt",append=TRUE)
write.table(O_Hist,file="CPContigSizeDist.txt", row.names = FALSE, quote = FALSE, append=TRUE)


##OUTPUT4

N <- ggplot() + geom_histogram(aes(NDMqlen$NDMqlen),fill='orange', color='white') +  xlab("Contig Length") + ylab("Number of Contigs") + ggtitle("NDM Carbapenamase Contig Length Distribution") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) +  theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous(labels = unit_format(unit = "KB", scale = 1e-3), breaks = scales::pretty_breaks(n = 15))
K <- ggplot() + geom_histogram(aes(KPCqlen$KPCqlen),fill='black', color='white') + xlab("Contig Length") + ylab("Number of Contigs") + ggtitle("KPC Carbapenamase Contig Length Distribution") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) +  theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous(labels = unit_format(unit = "KB", scale = 1e-3), breaks = scales::pretty_breaks(n = 10))
O <- ggplot() + geom_histogram(aes(OXAqlen$OXAqlen),fill='black', color='white') + xlab("Contig Length") + ylab("Number of Contigs") + ggtitle("OXA Carbapenamase Contig Length Distribution") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) +  theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous(labels = unit_format(unit = "KB", scale = 1e-3), breaks = scales::pretty_breaks(n = 10))

#----CREATE A MATRIX 
blastResults_Matrix <- table(blastResults_df[blastResults_df$cov==1 & blastResults_df$pident==100, c("assemblyName" , "sseqid")]) # subsets only the gene name with if CP gene has 100% cov
#class(blastResults_Matrix)
#blastResults_Matrix

# # 3) Report how many Assemblies have:
# #     a) Number of samples with co-cocarriage of different genes in different contigs eg: NDM-1#Contig1, KPC-2Contig5  - - Done!
# blastResults.filt 
dupAssemblies_logical <- duplicated(blastResults.filt[,c("assemblyName")]) | duplicated(blastResults.filt[,c("assemblyName")],fromLast=TRUE)
dupAssemblies <- blastResults.filt[dupAssemblies_logical,]
length(unique(dupAssemblies$assemblyName)) ## This many assemblies have multiple genes or same genes in various places
dupAssemblies

##SAMECP_SAMECONTIG - Done!

SameCP_SameContig <- dupAssemblies %>% 
  group_by(assemblyName,qseqid, sseqid) %>% 
  filter(n() > 1) 
# print(n = Inf)

#SameCP_SameContig

##SAMECP_DIFFCONTIG - groupby AssemblyName and sseqid - Done!

# For example, batch0_01032019_ENT1675_assembly.fasta and KPC-2 are grouped as one entity (something like a key) and filtered if they have distinct qseqid more than 1. 
#In this case of course, it have 2 qseqid for the record so it is displayed

SameCP_DiffContig <- dupAssemblies %>%
  group_by(assemblyName,sseqid) %>%
  mutate(key = n_distinct(qseqid)) %>%
  filter(key > 1)  %>%
  select(-key) 
# %>% print(n = Inf)

#SameCP_DiffContig


##DIFFCP_SAMECONTIG

DiffCP_SameContig <- dupAssemblies %>% 
  mutate(key = paste0(sseqid)) %>%
  group_by(assemblyName,qseqid) %>%
  filter(n_distinct(key) > 1) %>%
  select(-key) 
# %>% print(n = Inf)

#DiffCP_SameContig

##DIFFCP_DIFFCONTIG

DiffCP_DiffContig <- dupAssemblies %>% 
  mutate(key = paste0(qseqid)) %>%
  group_by(assemblyName) %>%
  filter(n_distinct(key) > 1 | n() == 1) %>%
  select(-key) %>%
  mutate(key2 = paste0(sseqid)) %>%
  group_by(assemblyName) %>%
  filter(n_distinct(key2) > 1 | n() == 1) %>%
  select(-key2) %>%
  print(n = Inf)

#DiffCP_DiffContig %>%
  print(n = Inf)

write.table(DiffCP_DiffContig,"DiffCP_DiffContig.txt",row.names = FALSE,col.names = TRUE, quote = FALSE) ## Save DiffCP_DiffContig output to txt file ##OUTPUT5
write.table(DiffCP_SameContig,"DiffCP_SameContig.txt",row.names = FALSE,col.names = TRUE, quote = FALSE) ## Save DiffCP_SameContig output to txt file ##OUTPUT6
write.table(SameCP_DiffContig,"SameCP_DiffContig.txt",row.names = FALSE,col.names = TRUE, quote = FALSE) # Output contains SameCP_DiffContig ##OUTPUT7
write.table( SameCP_SameContig,"SameCP_SameContig.txt",row.names = FALSE,col.names = TRUE, quote = FALSE) # Output contains SameCP_DiffContig ##OUTPUT8

SameCP_SameContig_Count <- length(unique(SameCP_SameContig$assemblyName))
SameCP_DiffContig_Count <- length(unique(SameCP_DiffContig$assemblyName))
DiffCP_DiffContig_Count <- length(unique(DiffCP_DiffContig$assemblyName))
DiffCP_SameContig_Count <- length(unique(DiffCP_SameContig$assemblyName))

SameCP_SameContig_Count
SameCP_DiffContig_Count
DiffCP_DiffContig_Count
DiffCP_SameContig_Count

cat("Assemblies with SameCP_SameContig_Count are : ", SameCP_SameContig_Count,file="Co-carriage_Report.txt",append = TRUE)
cat("\nAssemblies with SameCP_DiffContig_Count are : ", SameCP_DiffContig_Count,file="Co-carriage_Report.txt",append = TRUE)
cat("\nAssemblies with DiffCP_DiffContig_Count are : ", DiffCP_DiffContig_Count,file="Co-carriage_Report.txt",append = TRUE)
cat("\nAssemblies with DiffCP_SameContig_Count are : ", DiffCP_SameContig_Count,file="Co-carriage_Report.txt",append = TRUE)



#--> For Future 
# #     e) co-cocarriage of different variant genes in different contigs eg: NDM-1#Contig1, NDM-7#Contig5
# #     f) co-cocarriage of different variant genes in same contigs eg: NDM-1#Contig1, NDM-7#Contig1

#----CREATE A PRESENCE/ABSENCE MATRIX 
blastResults_BinaryMatrix <- as.matrix((blastResults_Matrix > 0) + 0) ## Convert to binary matrix
#class(blastResults_BinaryMatrix)
#blastResults_BinaryMatrix
#class(blastResults_Matrix )
#blastResults_BinaryMatrix


blastResults_MatrixLong <- melt(blastResults_Matrix) ##long format
heatMap_ggplot <- ggplot(blastResults_MatrixLong, aes(sseqid, assemblyName)) +
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient(low = "#F4B41A", high = "#143D59") +  theme(axis.text.y = element_blank()) + xlab("Carbapenamase Genes") + ylab("Assembly") + ggtitle("Carbapenamase Gene Profile Heatmap ")  ##https://www.tailorbrands.com/blog/logo-color-combinations
heatMap_ggplot 

############--------3) FIND INTERSECTIONS OF CPGENES IN VARIOUS ASSEMBLIES USING UPSETR------##########
write.csv(blastResults_BinaryMatrix, file = "cp_presence-abence_matrix.csv", row.names = TRUE,quote=FALSE) ##OUTPUT9

upsetdf <- read.csv(file = "cp_presence-abence_matrix.csv", check.names = FALSE)
upset_plot <- upset(upsetdf, order.by = "degree", nsets = 40, number.angles = 0, point.size = 1.5, line.size = 1, 
                    mainbar.y.label = "Sample Count", sets.x.label = "Carbapenamase Gene Set Size", sets.bar.color = "red",
                    text.scale = c(1.2, 1.3, 1, 1, 1.2, 1.6)) ##OUTPUT10

##SAVE UPSETR PLOT INTO A PDF

#----DRAW PHEATMAP
#PH <- pheatmap(blastResults_BinaryMatrix, color = colorRampPalette(c("white", "skyblue"))(50),fontsize= 7,  main = "Carbapenamase Gene HeatMap", show_rownames = FALSE) # labels_row = NULL)

#https://stackoverflow.com/questions/20500706/saving-multiple-ggplots-from-ls-into-one-and-separate-files-in-r
#MyPlots = list(N50_plot, N90_plot, , )

#png(file="filename.png", width = 1200, height = 800, units = "px",, res = 60)
tiff("filename.tiff", width = 1500, height = 2000, units = 'px', res = 150)
heatMap_ggplot
dev.off()


pdf("CPgene-profiler.pdf",width = 8, height = 10)

N50_plot
N90_plot
N
K 
O
upset_plot
heatMap_ggplot
#grid::grid.newpage()
#sprintf("Assemblies with SameCP_SameContig_Count are %i\nSameCP_DiffContig_Count%i-DiffCP_DiffContig_Count%i-DiffCP_SameContig_Count%i.pdf", SameCP_SameContig_Count, SameCP_DiffContig_Count, DiffCP_DiffContig_Count, DiffCP_SameContig_Count)
#sprintf("Assemblies with SameCP_SameContig_Count are %i\nSameCP_DiffContig_Count%i-DiffCP_DiffContig_Count%i-DiffCP_SameContig_Count%i.pdf", SameCP_SameContig_Count, SameCP_DiffContig_Count, DiffCP_DiffContig_Count, DiffCP_SameContig_Count)
#grid.table(DiffCP_DiffContig)

#grid.newpage()
dev.off()


# Tables in PDF # # https://stackoverflow.com/a/45353185
# library(gridExtra)
# pdf("data_output.pdf", height=11, width=11)
# grid.table(DiffCP_DiffContig)
# grid::grid.newpage()
# grid.table(DiffCP_SameContig)
# grid::grid.newpage()
# grid.table(SameCP_DiffContig)
# grid::grid.newpage()
# grid.table(SameCP_SameContig)
# grid::grid.newpage()
# dev.off()

cat("####--- PLOTTING DONE! Check \"pheatmap.pdf\" -------","\n")


cat("####---INTERSECTIONS OF CPGENES IS CHECKED! PLOTTING DONE! Check \"UpSetRPlot.pdf\" -------","\n")
#cat("\n\n\n ------ Check the following output files: \n 1) \"cp_presence-abence_matrix.csv\" \n 2) \"heatmap\" \n 3) \"CPGeneAssignment.txt\" \n 4) \"UpSetRPlot.pdf\" -------","\n")

# 4) Give more options to users to redefine the charts and blast tool - at the end simple
# 5) Remove passing assembly location from user - just the user have to the location where assemblies exist - at the end simple
# 7) Basic Stats of all assembly files - Total fasta, Genome Assembly Size plot, N50 plot


# ELAPSED TIME

new <- Sys.time() - old # calculate difference
print(new) # print in nice format

# ELAPSED TIME
