library(tidyverse)
library(pheatmap)
dat <- read_delim("G:/My Drive/Side projects/SturgeonRiverBacteria/data/OTUmatrices_04242020_rel.csv", delim = ",")
datfact <- read_delim("G:/My Drive/Side projects/SturgeonRiverBacteria/data/OTUmatrices_04242020Factors.csv", delim = ",", col_types = "fffff") %>%
  separate(Group, into = c(NA, "River", "Time")) %>%
  mutate(Time = recode(Time, "2" = "June"),
         Time = recode(Time, "1" = "April")) %>%
  unite("Group", River:Time) %>%
  column_to_rownames("Group")

selectList <- dat %>%
  pivot_longer(cols = -Group,
               names_to = "taxa",
               values_to = "count") %>%
  group_by(taxa) %>%
  summarize(n_taxa = sum(count)) %>%
  mutate(prop = n_taxa / sum(n_taxa)) %>%
  filter(prop > 0.002)

dat.filtered <- dat %>%
  select(selectList$taxa) %>%
  mutate(site = dat$Group) %>%
  column_to_rownames("site") %>%
  mutate_if(is.numeric, ~log10(.+0.001)) %>%
  t()

siteNames <- dat %>%
  select(Group) %>%
  separate(Group, into = c(NA, "River", "Time")) %>%
  mutate(Time = recode(Time, "2" = "June"),
         Time = recode(Time, "1" = "April")) %>%
  unite("Group", River:Time) %>%
  .$Group

colnames(dat.filtered) <- siteNames

### Pheatmap
library(RColorBrewer)
mat_colors <- list(UrbanAgriNLCD = brewer.pal(3, "Set1")[1:2],
                   ForestNLCD = brewer.pal(3, "Set2")[1:2],
                   Subbasin = brewer.pal(3, "Set3"),
                   SamplingTime = brewer.pal(3, "Dark2")[1:2])
names(mat_colors$UrbanAgriNLCD) <- unique(datfact$UrbanAgriNLCD)
names(mat_colors$ForestNLCD) <- unique(datfact$ForestNLCD)
names(mat_colors$Subbasin) <- unique(datfact$Subbasin)
names(mat_colors$SamplingTime) <- unique(datfact$SamplingTime)

pheatmap(dat.filtered,
         annotation_col = datfact,
         annotation_colors = mat_colors,
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean"
         )

## Modeling
library(pvclust)
result.site <- pvclust(dat.filtered,
                       method.dist="euclidean",
                       method.hclust="complete",
                       nboot=1000)

result.taxa <- pvclust(t(dat.filtered),
                       method.dist="euclidean",
                       method.hclust="complete",
                       nboot=1000)

par(mar=c(1,1,1,1))
plot(result.site,
     labels = rownames(datfact),
     print.pv = 2)

plot(result.taxa,
     labels = rownames(dat.filtered),
     print.pv = 2)
