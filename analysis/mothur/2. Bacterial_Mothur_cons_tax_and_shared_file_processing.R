#Based on a script originally created by Nick Sard on September 18th, 2017 
#JK

#ABOUT: This script was written to take the cons.taxonomy file produced by mothur in bacterial metabarcoding
#analyses and make 2 tables containing 3 columns (OTU#, Count, Phylum classification name or Lowest classification
#name) so it can be used to give names to the columns in the .shared file containing the read counts per sample.  
#If you have a shared file with more OTUs than PAST can handle (10000) that you want to analyze in PAST, you can also open the 
#subsample.shared file, remove singleton OTUs, truncate it to the top 10000 OTUs if necessary and then save it as a file.
#It also combines OTUs into unique Phyla and sums their read counts by Phylum.

#Processing cons.tax file

#loading necessary libraries
library(tidyverse)
library(tidyr)
library(Hmisc)

#setting working directory and reading in data
setwd("C:/")

#listing files
list.files("Input/")

#loading in the bacterial cons.taxonomy file which contains 3 columns named OTU, Size & Taxonomy (it only goes to Genus at the lowest; 6 categories: Domain; Phylum; Class; Order; Family; Genus)
tax <- read.table("Input/stability.trim.contigs.good.unique.good.filter.unique.pick.precluster.pick.opti_mcc.0.03.cons.taxonomy", header = T, sep = "\t", stringsAsFactors = F)
head(tax[,1:2])
head(tax)

#getting rid of the numbers representing classification certainty in the cons.tax file and putting this into a new column called tax2
tax$tax2 <- gsub(pattern = "\\([0-9]*\\)", replacement = "", x = tax$Taxonomy)
head(tax)

#separating the classification in column tax2 out into their own 6 columns
tax <- separate(data = tax, col = tax2, into = c("domain","phylum","class","order","family","lowest","extra"), sep = ";", fill = "right")
head(tax)

#removing that "extra" column
table(tax$extra)
tax$extra <- NULL

#Make a new table with the first column as the OTU#, the second as the count #, and the third as the lowest classification (what it's 
#lowest classification was to; could be any level as the lowest classification obtained is repeated for each of the subsequent categories)
LowestTax <- tax[,c("OTU", "Size", "lowest")]
head(LowestTax)
table(LowestTax$lowest)
tail(LowestTax)

#writing the abbreviated lowest taxonomy table to file
write.table(x = LowestTax, file = "Output/Bacteria.all.OTU97.lowest.taxonomy.abbreviated.txt", append = F, quote = F, sep = "\t", row.names = F, col.names = T)

#Make a table with phylum as the taxonomic level instead of lowest
PhylumTax <- tax[,c("OTU", "Size", "phylum")]
head(PhylumTax)
tail(PhylumTax)

#writing the abbreviated Phylum taxonomy table to file
write.table(x = PhylumTax, file = "Output/Bacteria.all.OTU97.Phylum.taxonomy.abbreviated.txt", append = F, quote = F, sep = "\t", row.names = F, col.names = T)


#########################

#Processing shared file

#Before starting, make sure you have replaced the Mothur sample IDs in the column called Group with your descriptive sample names
#in the subsample.shared file using the names in "RiverNames_01112021".

#loading in the bacterial shared file which has the OTUs at the top
df <- read.table("Input/stability.trim.contigs.good.unique.good.filter.unique.pick.precluster.pick.opti_mcc.0.03.subsample.shared_NAME.txt", header = T, sep = "\t", stringsAsFactors = F)

#get rid of the label and numOtus columns
df1 <- df %>% select(-label, -numOtus)
head(df1[,1:2])

#make a copy of df1 to remove singleton Otus from
df1_NS <- df1

#Set variables for for loop
Otu_num <- "NONE"
reads <- NULL
col_num <- 2

#ID and remove OTUs with only 1 read from df1_NS (singletons)
for (column in 2:ncol(df1)) {
  Otu_num <- colnames(df1)[col_num]
  #print(Otu_num)
  reads <- sum(as.numeric(as.character(df1[, col_num])))
  #print(reads)
  if (reads == 1){
    df1_NS <- df1_NS %>% select(-Otu_num)
  }
  
  #Reset variables
  Otu_num <- "NONE"
  reads <- NULL
  col_num <- col_num + 1
}

head(df1_NS)[,1:5]
tail(df1_NS)[,1:5]

#Calculate # of columns in new df
l <- ncol(df1_NS)
l

#truncate the file to the top 10,000 OTUs if necessary (or a desired #, would have to change comparison value for l) 
if (l > 10001) {
  df1_NS_trunc <- df1_NS %>% select(-(10002:ncol(df1_NS)))
}else {
  df1_NS_trunc <- df1_NS
}

head(df1_NS_trunc)[,1:5]
tail(df1_NS_trunc)[,1:5]

#Writing the truncated shared table to file
write.table(x = df1_NS_trunc, file = "Output/stability.trim.contigs.good.unique.good.filter.unique.pick.precluster.pick.opti_mcc.shared_NS_truncated_10000.txt", append = F, quote = F, sep = "\t", row.names = F, col.names = T)

##################################

#link OTU number and lowest classifcation level name in df1_NS_trunc and output the file

#Make a copy of df1_NS_trunc to modify
df1_NS_trunc_Lowest <- df1_NS_trunc

#Set variables for for loop
OTU1 <- "NONE"
col_num1 <- 2
lowest_name <-"NONE"

#loop through the truncated shared file and find matches to the appropriate abbreviated taxonomy file; replace the Otu# with the lowest tax 
#name in the output file

for (column1 in 2:ncol(df1_NS_trunc_Lowest)) {
  OTU1 <- colnames(df1_NS_trunc_Lowest)[col_num1]
  print(OTU1)
  lowest_name <- filter(LowestTax, OTU == OTU1)[[3]]
  colnames(df1_NS_trunc_Lowest)[colnames(df1_NS_trunc_Lowest)== OTU1] <- lowest_name
  print(lowest_name)
  
  #Reset variables
  OTU1 <- "NONE"
  col_num1 <- col_num1 + 1
  lowest_name <- "NONE"
}

head(df1_NS_trunc_Lowest)[,1:5]
tail(df1_NS_trunc_Lowest)[,1:5]

write.table(x = df1_NS_trunc_Lowest, file = "Output/stability.trim.contigs.good.unique.good.filter.unique.pick.precluster.pick.opti_mcc.shared_NS_truncated_10000_Lowest.txt", append = F, quote = F, sep = "\t", row.names = F, col.names = T)


#link OTU number and Phylum name in df1_NS_trunc and output the file; also creates a Unique_Phlya file that only contains the Phyla found in your subsampled shared file for use below

#Make a copy of df1_NS_trunc to modify
df1_NS_trunc_Phyla <- df1_NS_trunc

#Set variables for for loop
OTU2 <- "NONE"
col_num2 <- 2
phylum_name <-"NONE"
Unique_Phyla <- character()

#loop through the truncated shared file and find matches to the appropriate abbreviated taxonomy file; replace the Otu# with the 
#phylum name in the output file

for (column2 in 2:ncol(df1_NS_trunc_Phyla)) {
  OTU2 <- colnames(df1_NS_trunc_Phyla)[col_num2]
  print(OTU2)
  phylum_name <- filter(PhylumTax, OTU == OTU2)[[3]]
  colnames(df1_NS_trunc_Phyla)[colnames(df1_NS_trunc_Phyla)== OTU2] <- phylum_name
  print(phylum_name)
  if(phylum_name %nin% Unique_Phyla){
    Unique_Phyla <- c(Unique_Phyla, phylum_name)
  } 
  #Reset variables
  OTU2 <- "NONE"
  col_num2 <- col_num2 + 1
  phylum_name <- "NONE"
}

head(df1_NS_trunc_Phyla)[,1:5]
tail(df1_NS_trunc_Phyla)[,1:5]
Unique_Phyla

write.table(x = df1_NS_trunc_Phyla, file = "Output/stability.trim.contigs.good.unique.good.filter.unique.pick.precluster.pick.opti_mcc.shared_NS_truncated_10000_Phyla.txt", append = F, quote = F, sep = "\t", row.names = F, col.names = T)

###################################

#Make a dataframe with the phyla collapsed and read count summed for each sample from df1_NS_trunc_Phyla

#make a new dataframe to receive the condensed phyla
df1_Condensed_Phyla <- matrix(data = NA, nrow = nrow(df1_NS_trunc_Phyla), ncol = 1)
df1_Condensed_Phyla[,1] <- df1_NS_trunc_Phyla$Group
colnames(df1_Condensed_Phyla)[1] <- "Group"
df1_Condensed_Phyla <- as.data.frame(df1_Condensed_Phyla)

#set the empty variables for the loop
Phylum <- "NONE"
Total <- NULL
dftemp <- NULL

#Loop through the list and for each phylum collect matching columns from the df1_NS_trunc_Phyla dataframe and sum them for each sample and add
#the summed column to df1_Condensed_Phyla
for (taxon in Unique_Phyla) {
  Phylum <- taxon
  print(Phylum)
  dftemp <- df1_NS_trunc_Phyla[ , grepl(Phylum, names(df1_NS_trunc_Phyla))]
  z <- ncol(dftemp)
  print(z)
  if (is.null(z)) {
    Total <- dftemp
    Total <- data.frame(Total)
    names(Total) <- "X"
  }
  else {
    Total <- rowSums(dftemp[, 1:ncol(dftemp)]) 
    Total <- data.frame(Total)                                   
    names(Total) <- "X" 
  }
  
  df1_Condensed_Phyla$X <- Total$X 
  colnames(df1_Condensed_Phyla)[colnames(df1_Condensed_Phyla)=="X"] <- Phylum
  
  #reset variables
  Phylum <- "NONE"
  Total <- NULL
  dftemp <- NULL
  z <- NULL
}

#number of columns of table should be equal to the number of items in Unique_Phlya + 1
ncol(df1_Condensed_Phyla)
length(Unique_Phyla)

#write out the table with the summed read counts for each unique Phylum
write.table(x = df1_Condensed_Phyla, file = "Output/stability.trim.contigs.good.unique.good.filter.unique.pick.precluster.pick.opti_mcc.shared_NS_truncated_10000_Phyla_combined.txt", append = F, quote = F, sep = "\t", row.names = F, col.names = T)

#fin!