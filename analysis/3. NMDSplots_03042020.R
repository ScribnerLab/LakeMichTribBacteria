#creating NMDS plots for bacterial community data in vegan
#written by Gabrielle Sanfilippo - 04 March 2020

#This script was written to create three different NMDS plots of bacterial data
#ordinated using time (April vs. June collection) and geographic location in
#the Lake Michigan basin (North, East, West). In order to run the script for state ordination,
#the file "bacteria_state" (indicates where each river is located) must be used.
#to ordinate using time of collection, "bacteria_time" must be loaded in. 
#Ensure that 

#load vegan
library(vegan)
library(ggplot2)
#setting working directory and reading in data
setwd("")
list.files("")

#1-Bray-Curtis Dissimilarity Matrix 

#read in the Bacterial data  NOTE: make sure to have row.names = 1 so you can calculate the distances & have the sample names on your tree
#ensure that all 
bacteria.spp <- read.table(file = "Bacterial_OTU_Matrix", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
bacteria.spp_groups <- read.csv("bacteria_state.csv", row.names = 1)

#relative abundance
bacteria.spp.rel <- decostand(bacteria.spp, method = "total")

#calculate distance
bacteria.spp_distmat <- vegdist(bacteria.spp.rel, method = "bray")
bacteria.spp_distmat <- as.matrix(bacteria.spp_distmat, labels = T)
write.csv(bacteria.spp_distmat, "bacteria_spp_distmat.csv")

#Run NMDS using metaMDS
bacteria.spp_NMS <- metaMDS(bacteria.spp_distmat, distance = "bray", k = 3, maxit = 999, trymax = 250, wascores = TRUE)

#Shepards test/goodness of fit
#pdf(file = "shepardsplot_may_state.pdf")
#goodness(bacteria.spp_NMS)
#myplot <- stressplot(bacteria.spp_NMS, p.col="gray28", l.col = "black")
#dev.off()

#plot points in ordination
pdf(file = "NMDS_plot_stateMayl.pdf")
plot(bacteria.spp_NMS, "sites") #produces distance
orditorp(bacteria.spp_NMS, "sites")
with(bacteria.spp_groups, levels(State))
colvec <- c("black", "grey28")
pchvec <- c(24, 19, 15)
plot(bacteria.spp_NMS)
with(bacteria.spp_groups, points(bacteria.spp_NMS, 
                                 display = "sites",
                                 col = "black",
                                 pch = pchvec[State],
                                 bg = colvec[State]))
ordihull(bacteria.spp_NMS,
         bacteria.spp_groups$State,
         display = "sites",
         draw = c("polygon"),
         col = NULL,
         border = c("black", "black", "black"),
         lty = c(1, 2),
         lwd = 2.5)

#treat=c(rep("Time1",5),rep("Time2",5))
#ordiplot(bacteria.spp_NMS, type="n")
#ordihull(bacteria.spp_NMS, groups=treat,draw="polygon", col="grey90", label=F)
#orditorp(bacteria.spp_NMS,display="species",col="red",air=0.01)
#orditorp(bacteria.spp_NMS, display = "sites", col=c(rep("blue",5),rep("red",5)),air=0.01,cex=1.25) #gives points labels
dev.off()

bacteria.spp_NMS

#extract NMDS data for use in ggplot
#data.scores = as.data.frame(scores(bacteria.spp_NMS))
#add columns to data frame
#data.scores$Sample = bacteria.spp$Group
#data.scores$Time = bacteria.spp$Time
#data.scores$State = bacteria.spp$State

#head(data.scores)

#plot with ggplot
#pdf(file = "nmds_plot_gg.pdf")
#nmds_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
#geom_point(size = 9, aes(shape = State, colour = Time))+ 
#theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
#axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
#legend.text = element_text(size = 12, face ="bold", colour ="black"), 
#legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
#axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
#legend.title = element_text(size = 14, colour = "black", face = "bold"), 
#panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
#legend.key=element_blank()) + 
#labs(x = "NMDS1", colour = "Time", y = "NMDS2", shape = "State")  + 
#scale_colour_manual(values = c("purple", "green")) 

#nmds_plot
#ggsave("nmds.svg")

#dev.off()