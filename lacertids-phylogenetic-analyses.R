### 
library(SlicerMorphR)
library(geomorph)
library(tidyr)
library(paleomorph)
library(Morpho)
library(SamplingStrata)
library(ggplot2)
library(RColorBrewer)
library(geiger)
library(hot.dots)


### Load landmark data from Slicer and metadata

origin <- read.table("~/lizards/Analyses/metadata-lacertid.txt",h=T)

SM.log.file = "G:/slicer/GPA/2024-12-06_16_45_56/analysis.log"
SM.log <- parser(SM.log.file, forceLPS = T)

#SM.log$LM2 <- SM.log$LM[-c(45,47),,]

# Get mean shape, Csize, GPA aligned coords
SM.output <- read.csv(file = paste(SM.log$output.path, 
                                   SM.log$OutputData, 
                                   sep = "/"))
# Load PC scores
SlicerMorph.PCs <- read.table(file = paste(SM.log$output.path, 
                                           SM.log$pcScores, 
                                           sep="/"), 
                              sep = ",", header = TRUE, row.names = 1)



PD <- SM.output[,2] # Get Procrustes distances

# Number of landmarks, skipped LMS stored in SM.log$skipped
if (!SM.log$skipped) {
  no.LM <- SM.log$no.LM 
} else {
  no.LM <- SM.log$no.LM - length(SM.log$skipped.LM)
}
PD
no.LM

# Convert to 3D array
Coords <- arrayspecs(SM.output[, -c(1:3)], 
                     p = no.LM, 
                     k = 3)

dimnames(Coords) <- list(1:no.LM, 
                         c("x","y","z"),
                         SM.log$ID)


### Mirror half

sym <- matrix(nrow = 65, ncol = 2)
sym[,1] <- c(4:15,18,20:36, 38:44,46:54,56:60,62:69,71:73,75:77)
sym[,2] <- c(78:142)


#array(unlist(SM.log$LM), dim=c(dim(SM.log$LM[[1]]), length(SM.log$LM)))


SMextended <- array(data = NA, dim = c(nrow(sym)*2 + 12, 3, length(dimnames(SM.log$LM)[[3]])))

SMextended[c(1:77),,] <- SM.log$LM


LMsym <- mirrorfill(SMextended, l1 = as.integer(c(1,2,3,16,17,19,45,37,55,61,70,74)),l2 =sym)

######## GPA

gpa <- gpagen(LMsym)
gpa$coords <- gpa$coords[1:77,,]
dimnames(gpa$coords)[[3]] <- dimnames(SM.log$LM)[[3]]


df1 <- origin[ origin$Specimen_id %in% SM.log$ID,]

df2 <- df1[match(SM.log$ID, df1$Specimen_id),]
rownames(df2) <- df2$Specimen_id


####### Importing and preparing the Phyologeny

fulltree<-read.tree("~/data-metadata/squamates_Title_Science2024_ultrametric_constrained.tre") #
La_tree <- extract.clade(fulltree,node = 10623)

# Add taxa

La_tree<-bind.tip(La_tree,"Podarcis_filfolensis_lija",where=which(La_tree$tip.label=="Podarcis_filfolensis"),
                   position=0.5) # Salvi et al. 2014  https://doi.org/10.1111/mec.12668

La_tree<-bind.tip(La_tree,"Podarcis_thais",where=which(La_tree$tip.label=="Podarcis_peloponnesiacus"),
               position=1.72) # Kiourtsoglou et al. 2021  https://onlinelibrary.wiley.com/doi/full/10.1111/jzs.12540
La_tree<-bind.tip(La_tree,"Podarcis_tiliguerta_C",where=which(La_tree$tip.label=="Podarcis_tiliguerta"),
                   position=4) # Add for exact estimates of divergence time!!!

La_tree<-bind.tip(La_tree,"Podarcis_hispanicus_AM",where=which(La_tree$tip.label=="Podarcis_hispanicus"),
                  position=2.5) # Add for exact estimates of divergence time!!!

La_tree<-bind.tip(La_tree,"Podarcis_ionicus",where=which(La_tree$tip.label=="Podarcis_tauricus"),
                  position=4.91) # Psonis et al. 2021 https://www.sciencedirect.com/science/article/pii/S1055790321000543 

La_tree<-bind.tip(La_tree,"Podarcis_muralis_SA",where=which(La_tree$tip.label=="Podarcis_muralis"),
                  position=5.76) 

La_tree<-bind.tip(La_tree,"Podarcis_muralis_Greece",where=which(La_tree$tip.label=="Podarcis_muralis"),
                  position=4.05) 

La_tree<-bind.tip(La_tree,"Zootoca_carniolica",where=which(La_tree$tip.label=="Zootoca_vivipara"),
                  position=4.5) 

La_tree<-bind.tip(La_tree,"Gastropholis_echinata",where=which(La_tree$tip.label=="Gastropholis_prasina"),
                  position=6.127028) # Tonin et al. 2016

La_tree<-bind.tip(La_tree,"Acanthodactylus_lineomaculatus",where=which(La_tree$tip.label=="Acanthodactylus_erythrurus"),
                  position=5.99) # from Ranchillac et al. 2023 https://doi.org/10.24072/pcjournal.301

La_tree$tip.label[La_tree$tip.label == "Podarcis_muralis"] <- "Podarcis_muralis_Brongniardii"
La_tree$tip.label[La_tree$tip.label == "Podarcis_hispanicus"] <- "Podarcis_hispanicus_GAL"
La_tree$tip.label[La_tree$tip.label == "Podarcis_tiliguerta"] <- "Podarcis_tiliguerta_S"
La_tree$tip.label[La_tree$tip.label == "Eremias_kavirensis"] <- "Eremias_fasciata" # From Orlova et al. 2023 https://doi.org/10.11646/zootaxa.5369.3.2 


#### Note: I couldn't find a phyologeny for Eremias nigrocellata, Philochortus hardeggeri


#df2$Species[df2$Species == "Acanthodactylus_lineomaculatus"] <- "Acanthodactylus_erythrurus" # Check phylogeny when adding related sp A. erythrurus !!!

df2$Species[df2$Species == "Pseuderemias_striatus"] <- "Pseuderemias_smithii" # Check phylogeny if ever adding another Pseuderemias species !!!

##Philochortus has phylogeny for P. spinalis only !!! 


#tree2<-read.tree("G:/macroevo/squam_shl_new_Consensus_9755.tre") ## Tree from Tonin et al. 2016

## Extract Podarcis
# Podarcis_tree <- extract.clade(La_tree,node = 534)
# plotTree(Podarcis_tree,fsize =0.5, lwd =1.5, ftype="i", node.numbers = T)
# 
# Philo.tree <- keep.tip(tree2, tree2$tip.label[grep(x = tree2$tip.label, pattern = "Philochortus")]) 
# plotTree(Philo.tree,fsize =0.5, lwd =1.5, ftype="i", node.numbers = T)

dimnames(gpa$coords)[[3]] <- df2$Species

LacertShape<-two.d.array(gpa$coords)

producer <- df2$producer
names(producer) <- df2$Species

head(LacertShape)
red.data<-treedata(La_tree,LacertShape)
red.data2 <- treedata(La_tree,producer)


## Exploratory analyses - Polymorphospace #####
# PCA.phylo <- gm.prcomp(phy = red.data$phy, gpa$coords[,,df2$Species != "Darevskia_obscura"
#                                                       & df2$Species != "Eremias_nigrocellata"
#                                                       & df2$Species != "Philochortus_hardeggeri"])
# par(font.lab = 2)
# plot(PCA.phylo, axis1 = 1 ,axis2 = 2 ,phylo = TRUE, main = "",pch = 21, bg = "darkgreen",
#      phylo.par = list(tip.labels = T, 
#                      tip.txt.cex = 0.5,
#                       edge.color = "grey", edge.width = 1, node.txt.cex = 0.1, node.cex = 0.5)#,time.plot = T
# )
# abline(h = 0, lty = 2, col = "lightgrey")
# abline(v = 0, lty = 2, col = "lightgrey")
# 
# #Phylogenetic PCA - GLS centering and project as in Revell (2009)
# PCA.w.phylo <- gm.prcomp(phy = red.data$phy,gpa$coords[,,df2$Species != "Darevskia_obscura"
#                                                        & df2$Species != "Eremias_nigrocellata"
#                                                        & df2$Species != "Philochortus_hardeggeri"],GLS = TRUE, transform = TRUE)
# 
# # Can be used to assess adaptive radiation. 
# 
# 
# ## PCA independent of phylogeny
# par(font.lab = 2)
# plot(PCA.w.phylo, axis1 = 1 ,axis2 = 4 ,phylo = TRUE, main = "",pch = 21, bg = "blue",
#      phylo.par = list(tip.labels = T, 
#                       tip.txt.cex = 0.8, edge.color = "grey", edge.width = 1, node.txt.cex = 0.1, node.cex = 0.5)# ,time.plot = T
#      )
# 
# 
# ## PaCA
# PaCA.w.phylo <- gm.prcomp(phy = red.data$phy,gpa$coords[,,df2$Species != "Darevskia_obscura"
#                                                         & df2$Species != "Eremias_nigrocellata"
#                                                         & df2$Species != "Philochortus_hardeggeri"],align.to.phy = TRUE)
# plot(PaCA.w.phylo, phylo = TRUE, main = "",pch = 21, bg = "blue",
#      phylo.par = list(tip.labels = TRUE, 
#                       tip.txt.cex = 0.8, edge.color = "grey", edge.width = 1, node.txt.cex = 0.1, node.cex = 0.1)# ,time.plot = T
# )
# 
# 
# 
# msh <- mshape(gpa$coords[,,df2$Species != "Darevskia_obscura"
#                          & df2$Species != "Eremias_nigrocellata"
#                          & df2$Species != "Philochortus_hardeggeri"])
# 


### Generate wrapped mesh for mean shape
# findMeanSpec(gpa$coords)
refmesh <- read.ply("G:/macroevo/99814-mean-shape.ply")
ref3d <- warpRefMesh(refmesh, SM.log$LM[,,"99814"], msh, color = "lightgray", centered = FALSE)


# reconstruct ancester states
#ancester <- warpRefMesh(refmesh, SM.log$LM[,,"99814"], arrayspecs(PCA.w.phylo$ancestors,60,3)[,,1], color = "lightgray", centered = FALSE)


# deformGrid3d(PCA.phylo$shapes$shapes.comp1$min,PCA.phylo$shapes$shapes.comp1$max,col1 = "lightgrey", col2 ="orange", size = 0.003, lcol = "orange",lwd = 4)
# #wire3d(defPC1min,col=4)
# shade3d(ref3d,col = skin1 , alpha = 0.4
# )
# rgl.bg( sphere = FALSE, fogtype = "none", color = "black" )

#### Evolutionary rates
partition <- read.table("~/data-metadata/modules.txt", h =T)

#partition <- partition[partition$name != "corner-postfront-postrbi",]
EMR<-compare.multi.evol.rates(A=gpa$coords[complete.cases(partition$module),,df2$Species != "Darevskia_obscura"
                                           & df2$Species != "Eremias_nigrocellata"
                                           & df2$Species != "Philochortus_hardeggeri"],gp=partition$module[complete.cases(partition$module)], 
                              Subset=TRUE, phy= red.data$phy, iter=9999,print.progress = T)
summary(EMR)
plot(EMR, xlim =c(0.5,1.7))

hist(EMR$random.sigma, breaks = 100, xlim = c(0.9, 1.8),
     xlab = "Rate Ratios",
     main = "")
abline(v = 1.55545, lwd= 0.5, col = "blue")



mod.la <- phylo.modularity(gpa$coords[complete.cases(partition$module),,df2$Species != "Darevskia_obscura"
                            & df2$Species != "Eremias_nigrocellata"
                            & df2$Species != "Philochortus_hardeggeri"],partition.gp = partition$module[complete.cases(partition$module)],
                 phy= red.data$phy, iter=9999,print.progress = T)


#### Phylogenetic signal
PS.shape <- physignal(A=gpa$coords[,,df2$Species != "Darevskia_obscura"
                                   & df2$Species != "Eremias_nigrocellata"
                                   & df2$Species != "Philochortus_hardeggeri"],phy=red.data$phy,iter=9999,print.progress = T)
summary(PS.shape)


PS.shape2 <- physignal.z(A=gpa$coords[,,df2$Species != "Darevskia_obscura"
                                      & df2$Species != "Eremias_nigrocellata"
                                      & df2$Species != "Philochortus_hardeggeri"],
                         lambda = "front", PAC.no = 50 ,
                         phy=red.data$phy,iter=9999,print.progress = T)
summary(PS.shape2)

plot(PS.shape2)

### Compare phylogenetic signals between modules

meso <- gpa$coords[complete.cases(partition$module) & partition$module == "mes",,df2$Species != "Darevskia_obscura"
                   & df2$Species != "Eremias_nigrocellata"
                   & df2$Species != "Philochortus_hardeggeri"]
NC <- gpa$coords[complete.cases(partition$module) & partition$module == "NC",, df2$Species != "Darevskia_obscura"
                 & df2$Species != "Eremias_nigrocellata"
                 & df2$Species != "Philochortus_hardeggeri"]

PS.meso <- physignal.z(A = meso, phy = red.data$phy, 
                       lambda = "front", PAC.no = 7,
                       iter=9999, verbose = T)

PS.NC <- physignal.z(A = NC, phy = red.data$phy, 
                       lambda = "front", PAC.no = 7,
                       iter=9999, verbose = T)

PS.list <-list(PS.meso, PS.NC)
names(PS.list) <- c("mesoderm", "NC")

PS.Z <- compare.physignal.z(PS.list)

summary(PS.Z)

plot(PS.Z)


PS.Z2 <- compare.physignal.z(PS.meso,PS.NC)
PS.Z2

##### Disparity - not a great method
# disp.meso <- morphol.disparity(meso ~ 1)
# disp.meso
# disp.meso/nrow(meso[,,1])
# 
# 
# ###### 
# 
# disp.nc <- morphol.disparity(NC ~ 1)
# disp.nc
# disp.nc/nrow(NC[,,1])




##### Per-landmark variance

var.la <- per_lm_rates(shape.data = gpa$coords[,,df2$Species != "Darevskia_obscura"
                                                  & df2$Species != "Eremias_nigrocellata"
                                                  & df2$Species != "Philochortus_hardeggeri"], phy = red.data$phy)

spheres3d(msh, col = var.la$Rate_Colors, radius =  0.005)
shade3d(ref3d,col = skin1, alpha = 0.5
)
rgl.bg( sphere = FALSE, fogtype = "none", color = "black" )
#snapshot3d("~/lizards/Analyses/Lacertid-rates-vent-hotdots.png")

# Reorder colors according to sorted Log_Variance
sorted_indices <- order(var.la$Log_Rates)
sorted_log_variance <- var.la$Log_Rates[sorted_indices]
sorted_colors <- var.la$Rate_Colors[sorted_indices]

# Create a color ramp function
color_ramp <- colorRampPalette(sorted_colors)

# Generate a sequence of colors
color_sequence <- color_ramp(100)


# Plot the legend
par(mar = c(5, 4, 4, 2) + 1)

image(1, seq(min(sorted_log_variance), max(sorted_log_variance), length.out = 100), t(matrix(seq(min(sorted_log_variance), max(sorted_log_variance), length.out = 100), ncol = 1)), col = color_sequence, axes = FALSE, xlab = "", ylab = "")
axis(4, at = seq(min(sorted_log_variance), max(sorted_log_variance), length.out = 10), labels = round(seq(min(sorted_log_variance), max(sorted_log_variance), length.out = 10), 2), las = 1)
mtext("Log Variance", side = 4, line = 3)

# Only showing extreme ticks
image(1, seq(min(sorted_log_variance), max(sorted_log_variance), length.out = 100), t(matrix(seq(min(sorted_log_variance), max(sorted_log_variance), length.out = 100), ncol = 1)), col = color_sequence, axes = FALSE, xlab = "", ylab = "")
axis(4, at = c(min(sorted_log_variance), median(sorted_log_variance), max(sorted_log_variance)), labels = round(c(min(sorted_log_variance), median(sorted_log_variance), max(sorted_log_variance)), 2))
mtext("Log Variance", side = 4, line = 3)

variance.la <- per_lm_variance(shape.data = gpa$coords[,,df2$Species != "Darevskia_obscura"
                                                       & df2$Species != "Eremias_nigrocellata"
                                                       & df2$Species != "Philochortus_hardeggeri"])
fac.mod <- partition$module
fac.mod[is.na(fac.mod)] <- "none"
fac.mod <- as.factor(fac.mod)
plot(var.la$Per_Lm_Rates, variance.la$Per_Lm_Variance, pch = 21, bg = c("orange","blue","grey")[fac.mod],
     xlab = "Landmark evolutionary rate", ylab = "Landmark variance")
abline()

lm.var.rates <- lm(var.la$Per_Lm_Rates ~  partition$module*variance.la$Per_Lm_Variance)
summary(lm.var.rates)

lm(var.la$Per_Lm_Rates[partition$module =="mes"] ~  variance.la$Per_Lm_Variance[partition$module =="mes"])


Pvar <- readRDS("G:/results/figures/Fig2orSM/lmOrder.rds") # Export from testing correlations with P. muralis

d.rates <- data.frame(lm = partition$name,   
                       variance_la_Per_Lm_Variance = variance.la$Log_Variance,
                       var_la_Per_Lm_Rates = var.la$Log_Rates,
                       var = variance.la$Per_Lm_Variance,
                       rate = var.la$Per_Lm_Rates,
                       partition_module = partition$module,
                       Podarcis_var <- Pvar$var
                       )


ggplot(d.rates[complete.cases(d.rates$partition_module),], aes( x = var_la_Per_Lm_Rates, y = variance_la_Per_Lm_Variance, color = partition_module)) +
  geom_point(cex = 2) +
  geom_smooth(method = "lm", se = T, alpha = 0.05) +
  scale_color_manual(values = c("mes" = "orange2", "NC" = "cornflowerblue")) +
  labs(title = "", x = "log landmark evolutionary rate", y = "log landmark variance") + theme_classic() +
  labs(colour = "Module") +
  theme(
    axis.title = element_text(size = 16),  # Increase axis label size
    axis.text = element_text(size = 14)   # Increase plot title size
    )

#### Bar plot

d.rates$partition_module[is.na(d.rates$partition_module)] <- "unasigned"
d.rates$partition_module<- factor(d.rates$partition_module, levels = c("NC","mes","unasigned"))
d.rates <- d.rates[order(d.rates$partition_module,d.rates$Podarcis_var),]

colors <- ifelse(d.rates$partition_module == 'mes', 'orange1', ifelse(d.rates$partition_module == 'NC', 'cornflowerblue', 'white'))

par(mar=c(10,5,4,4))
# variance barplot
barplot(d.rates$var, names.arg = d.rates$lm, las = 2, col = colors, 
        xlab = '', ylab = '', 
        main = '')
# rates barplot
barplot(d.rates$rate, names.arg = d.rates$lm, las = 2, col = colors, 
        xlab = '', ylab = '', 
        main = '')
dev.off()
plot.new()
legend('top', legend = c('Neural crest', 'Mesoderm','Unasigned'), fill = c( 'lightblue', 'orange1','white'), bty = "n")
