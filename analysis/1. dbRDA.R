####### Lake Michigan bacterial community analysis ########
## Jared Homola
## March 2020

library(readxl)
library(vegan)
library(corrplot)
library(Hmisc)
library(ggrepel)
library(tidyverse)

setwd("C:/Users/HP/Desktop/SturgeonRiverBacteria")


### Load data ###
# Lat/Longs
latLong <- read_excel("2019Bacteria_RiverCoordinates_11720.xlsx") %>% 
  mutate(River = recode(River, `Pere Marquette` = "PereMarquette")) %>% 
  select(River, Lat, Long)

# Raw environmental data
envRaw <- read_excel("bacteria_rivers_reduced_landscape_vars_raw_&_transformed.xlsx",
                     sheet = "raw_data") %>% 
  filter(River != "Ontonagon",
         River != "Sturgeon", 
         River != "Black") %>% 
  select(-COMID) %>% 
  mutate(River = recode(River, `Pere Marquette` = "PereMarquette")) %>% 
  left_join(latLong, by = "River") %>% 
  mutate(River = as.factor(River),
         Lat = as.numeric(Lat))

# Transformed environmental data
envTrans <- read_excel("bacteria_rivers_reduced_landscape_vars_raw_&_transformed.xlsx",
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

# Time 1 OTUs
otu1 <- read_excel("OTUmatrix_time1.xlsx") %>% 
  column_to_rownames(var = "Group") %>% 
  as_tibble(rownames = NA)

# Time 2 OTUs
otu2 <- read_excel("OTUmatrix_time2.xlsx") %>% 
  column_to_rownames(var = "Group") %>% 
  as_tibble(rownames = NA)

################################################################################

### Analysis 1 - Examine explanatory variables ###
## Is there collinearity?
corrplot(cor(as.matrix(envRaw[2:11])), 
         method = "number",
         type = "upper")

## Remove N, P, Forest

## Are data normally distributed?
ggplot(envTrans, aes(x = UMD_T)) +
  geom_histogram(alpha = 0.2) +
  geom_density() +
  theme_bw(base_size = 18)

## Are the data normal
shapiro.test(envRaw$N_areasqkm)
shapiro.test(envRaw$UMD)
shapiro.test(envRaw$UDOR)
shapiro.test(envRaw$N_total_n_yield)
shapiro.test(envRaw$N_total_p_yield)
shapiro.test(envRaw$N_urban)
shapiro.test(envRaw$N_forest)
shapiro.test(envRaw$N_agriculture)

shapiro.test(envTrans$N_areasqkm_T)
shapiro.test(envTrans$UMD_T)
shapiro.test(envTrans$UDOR_T)
shapiro.test(envTrans$N_urban_T)
shapiro.test(envTrans$N_agriculture_T)
shapiro.test(envTrans$Lat_T)
shapiro.test(envTrans$Long_T)


################################################################################
### Analysis 3 - Distance-based redundancy analysis
## Time 1
# Parse tibbles
rda.otu1.tib <- otu1 %>% 
  rownames_to_column("River") %>% 
  left_join(envTrans, by = "River") %>% 
  column_to_rownames("River")

otu1.dat <- rda.otu1.tib[1:5688]
otu1.env <- rda.otu1.tib[5689:ncol(rda.otu1.tib)] %>% 
  scale() %>%  ## Data of different units need to be scaled
  as_tibble()

# Run analysis
rda.otu1 <- dbrda(otu1.dat ~ N_areasqkm_T + UMD_T + UDOR_T +
                  N_urban_T + N_agriculture_T + Lat_T + Long_T,
                  data = otu1.env,
                  distance = "bray")

# Check for multiple collinearity issues
# Values >10 are problematic
vif.cca(rda.otu1)

# Evaluate model significance
anova(rda.otu1, permutations = 10000)
anova(rda.otu1, by = 'axis', permutations = 10000)
anova(rda.otu1, by = 'terms', permutations = 10000)


## Time 2
# Parse tibbles
rda.otu2.tib <- otu2 %>% 
  rownames_to_column("River") %>% 
  left_join(envTrans, by = "River") %>% 
  column_to_rownames("River")

otu2.dat <- rda.otu2.tib[1:5688]
otu2.env <- rda.otu2.tib[5689:ncol(rda.otu2.tib)] %>% 
  scale() %>%  ## Data of different units need to be scaled
  as_tibble()

# Run analysis
rda.otu2 <- dbrda(otu2.dat ~ N_areasqkm_T + UMD_T + UDOR_T +
                    N_urban_T + N_agriculture_T + Lat_T + Long_T,
                  data = otu2.env,
                  distance = "bray")

# Check for multiple collinearity issues
# Values >10 are problematic
vif.cca(rda.otu2)

# Evaluate model significance
anova(rda.otu2, permutations = 10000)
anova(rda.otu2, by = 'axis', permutations = 10000)
anova(rda.otu2, by = 'terms', permutations = 10000)


