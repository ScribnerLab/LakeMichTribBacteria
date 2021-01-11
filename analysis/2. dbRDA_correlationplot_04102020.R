#written by jared homola
#visual correlation matrix edited by gabrielle sanfilippo - April 10, 2020
library(readxl)
library(vegan)
library(corrplot)
library(Hmisc)
library(ggrepel)
library(tidyverse)

setwd("/Users/gabriellesanfilippo/Desktop/Scribner Lab/Bacterial Water Filter Seq_Fall2019/dbRDA")

### Load data ###
# Lat/Longs
latLong <- read_excel("2019Bacteria_RiverCoordinates_11720.xlsx") %>% 
  mutate(River = recode(River, `Pere Marquette` = "PereMarquette")) %>% 
  select(River, Lat, Long)

# Transformed environmental data
envTrans <- read_excel("bacteria_river_sites_with_reduced_landscape_variables_raw_and_transformed.xlsx",
                       sheet = "transformed") %>% 
  filter(River != "Ontonagon",
         River != "Sturgeon", 
         River != "Black") %>% 
  select(-COMID) %>% 
  mutate(River = recode(River, `Pere Marquette` = "PereMarquette")) %>% 
  left_join(latLong, by = "River") %>% 
  mutate(River = as.factor(River),
         Lat = as.numeric(Lat)) %>% 
  mutate(Lat_T = log(Lat + 0.01),
         Long_T = log(-1 * Long + 0.01)) %>% 
  select(-Lat, -Long)

# Time 1 OTUs - must have separate files for April samples and June samples. Here, April samples are called "Time 1".
otu1 <- read_excel("OTUmatrix_time1.xlsx") %>% 
  column_to_rownames(var = "Group") %>% 
  as_tibble(rownames = NA)

# Time 2 OTUs - need a separate file for June ("time 2") samples
otu2 <- read_excel("OTUmatrix_time2.xlsx") %>% 
  column_to_rownames(var = "Group") %>% 
  as_tibble(rownames = NA)

#examine explanatory variables
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
pdf("corrplot_04102020.pdf")
corrplot(cor(as.matrix(envTrans[2:11])), 
         method = "color",
         col=col(200),
         addCoef.col = "black",
         tl.col = "black", tl.srt = 40, tl.cex = 0.7,
         type = "upper",
         outline = TRUE,
         outline.color = "white",
         number.cex=0.8)
dev.off()