# We want to find the CPgenes from a bacterial assembly and plot a heatmap of CP genes from all the assemblies based on the categories of carbapenamases
# First step is to blast the Assembly against CP gene and find the list of the contigs with CP genes
# Second, the CP genes are assigned to the Assembly file and a table is created
# From the table a heatmap is drawn
old <- Sys.time()

###########--------LIBRARIES------##########

library(tidyverse)
library(pheatmap)
library(UpSetR)

############--------CONSTANT VARIABLES------##########

blast_db = "/home/datta/Desktop/CPE/R_package/ARG-annot-20200406T084430Z-001/ARG-annot/ARG-annot_CPGene_DB/ARG-annot_CPGene.DB"
blastn = "/usr/bin/blastn"
evalue = 1e-6
format = "\'6 qseqid sseqid pident nident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\'"
num_threads = 4


colnames <- c("qseqid",
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

############--------ASSEMBLY FILES LOCATION------##########

# TODO - Request path from user
setwd("/home/datta/Desktop/CPE/R_package/ARG-annot-20200406T084430Z-001/ARG-annot/test/Assemblies")
#setwd("/home/datta/Desktop/CPE/NP_assemblies")

############--------CLEANING IF ALREADY FILE EXISTS------##########
unlink("CPGeneAssignment.txt")

############--------1) BLASTN ASSEMBLY ON ARG-ANNOT DATA------##########

#List of all the fasta files
files <- list.files()

lapply(files, function(assemblyFile) { 

#input = "/home/datta/Desktop/CPE/NP_assemblies/batch14_02082019_ENT1646_assembly.fasta"
input = assemblyFile

blast_out <- system2(command = blastn, 
                     args = c("-db", blast_db, 
                              "-query", input, 
                              "-outfmt", format, 
                              "-evalue", evalue,
                              "-ungapped",
                              "-num_threads", num_threads),
                     wait = TRUE,
                     stdout = TRUE)

#blast_out

############--------BLASTN OUTPUT TO A TABLE ------##########

tidy_blast <- blast_out %>%
  as_tibble() %>% 
  separate(col = value, 
           into = colnames,
           sep = "\t",
           convert = TRUE)

############--------FIND CP GENE CONTAINING CONTIGS WITH CPGENE MATCHING 100% ID AND COVERAGE ------##########
tidy_blast$cov <- tidy_blast$nident / tidy_blast$slen
#tidy_blast
cp <- subset(tidy_blast$sseqid, tidy_blast$cov==1 & tidy_blast$pident==100) # subsets only the gene name with if CP gene has 100% cov
#sprintf("Assembly:%s CPgene: %s ", input, cp)

############--------OUTPUT ASSEMBLY NAME AND CORRESPONDING CP GENE IN ASSEMBLY TO A FILE ------##########
CPassignment = paste(input,cp,sep=",")
write(CPassignment,file="CPGeneAssignment.txt",append=TRUE)
})

############--------2) MAKE HEATMAP OF CPGENE PRESENCE/ABSENCE IN THE ASSEMBLIES------##########
##Part 2. collect all the list of the CP genes for all assemblies, make unique list of the CPgenes and create a presence-absence matrix

#----LOAD THE BLASTN RESULTS OUTPUT FROM ABOVE
Samp_CP_df <- read.table("CPGeneAssignment.txt", header = FALSE, sep = ",")
colnames(Samp_CP_df)<-c("AssemblyFile","CPGene")

Samp_CP_df[Samp_CP_df==""]<-NA #remove columns with NA

#----TRIMMING STRINGS IN THE TABLE
Samp_CP_df$AssemblyFile <- gsub("_assembly.fasta", "", Samp_CP_df$AssemblyFile)
Samp_CP_df$CPGene <- gsub("\\(Bla\\)", "", Samp_CP_df$CPGene)
Samp_CP_df$CPGene <- gsub("[:_].*", "", Samp_CP_df$CPGene)
#Samp_CP_df$CPGene <- gsub("^$", "NoCP", Samp_CP_df$CPGene)

#----CREATE A PRESENCE/ABSENCE MATRIX 
binarymat <- data.frame(unclass(table(Samp_CP_df[1:2])),check.names = FALSE) ##Generate a df for convenience
binarymat <- as.matrix((binarymat > 0) + 0) ## Convert to binary matrix

#----DRAW BASE HEATMAP
png(filename="basehm.png",height=8,width=12,res=400,units="in")
heatmap(t(binarymat),Rowv=NA,Colv=NA,na.rm=T,scale="none",col=terrain.colors(100),
        xlab="",ylab="",main="Carbapenamase Gene HeatMap", labCol = FALSE)
dev.off()

#----DRAW PHEATMAP
xx <- pheatmap(binarymat, color = colorRampPalette(c("yellow", "red"))(50),fontsize= 7,  main = "Carbapenamase Gene HeatMap", show_rownames = FALSE) # labels_row = NULL)

#Save heatmap
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(xx, "heatmap.pdf")

############--------3) FIND INTERSECTIONS OF CPGENES IN VARIOUS ASSEMBLIES USING UPSETR------##########
write.csv(binarymat, file = "upset_input.csv", row.names = TRUE,quote=FALSE)

upsetdf <- read.csv(file = "upset_input.csv", check.names = FALSE)
upset_plot <- upset(upsetdf, order.by = "degree", nsets = 40, number.angles = 0, point.size = 1.5, line.size = 1, 
      mainbar.y.label = "Sample Count", sets.x.label = "Carbapenamase Gene Set Size", sets.bar.color = "red",
      text.scale = c(1.2, 1.3, 1, 1, 1.2, 1.6))

##SAVE UPSETR PLOT INTO A PDF

pdf(file="upset_plot.pdf") # or other device
upset_plot
dev.off()

# ELAPSED TIME
new <- Sys.time() - old # calculate difference
print(new) # print in nice format
