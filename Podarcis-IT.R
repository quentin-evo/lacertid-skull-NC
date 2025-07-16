library(SlicerMorphR)
library(rgl)
library(geomorph)
library(tidyr)
library(paleomorph)
library(Morpho)
library(SamplingStrata)
library(ggplot2)
library(RColorBrewer)
library(hot.dots)

#load("G:/Bonn-lizards.RData")

origin <- read.table("~/data-metadata/metadata-Podarcis.txt",h=T)
#origin$type <- as.factor(origin$darkness)
#origin$colour <- as.factor(origin$colour)
#origin<-separate(origin, "darkness",into = c("colour","sex"), sep = "-")
origin$phenot_ordi_SpecCheck <- NA
origin$phenot_ordi_SpecCheck[complete.cases(origin$spec)] <- var.bin(
  log(origin$spec[complete.cases(origin$spec)]),bins = 5, iter.max = 100)



SM.log.file = "~/data-metadata/2024-12-03_15_49_50/analysis.log"
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

#### Check outliers
plotOutliers(gpa$coords)



##### Exploring the data - Preliminary analyses ######
# 
# 
# gdf2 <- geomorph.data.frame(coords = gpa$coords, size = gpa$Csize, phenotype = df2$phenotype, sex = df2$sex, SVL = df2$SVL,
#                             spec = df2$spec, phenotO = df2$phenotype_ordinal, lin = df2$lineage, 
#                             loc = df2$Location)
# 
# PCA <- gm.prcomp(gdf2$coords)
# 
# 
# m1 <- procD.lm(gdf2$coords ~ log(size) , iter= 9999, data = gdf2,RRPP = T)
# 
# # pairwise
# m0 <- procD.lm(coords ~ log(size) , iter= 9999, data = gdf2,RRPP = T)
# PW <- pairwise(m1, m0, groups = gdf2$colour)
# summary(PW, test.type = "dist", confidence = 0.95)
# 
# summary(m1)
# 
# ## Plot PCA
# 
# 
# p1 <- plot(PCA,axis1 =1,axis2 = 2,  pch =  c(21,22)[as.factor(gdf2$sex)], 
#      bg = c(
#        "orange","brown", "orange3","lightblue","lightgrey",
#        "green","green3","darkolivegreen1","darkolivegreen","darkolivegreen3"
#        
#             )[as.factor(gdf2$type)], cex = 1.5,
#      main = "", xlab = "PC1 (26.4%)", ylab = "PC2 (11.4%)")
# 
# 
# picknplot.shape(p1)
# 
# 
# # gla<-levels(gdf.af$group)
# # for(i in 1:length(gla)) {
# #   x<-plotaf$pc.scores[,c(1,2)][gdf.af$group == gla[i], ]
# #   CHP2<-chull(x)
# #   CHP2<-c(CHP2,CHP2[1])
# #   points(x[CHP2,],type='l',col= i)
# # }
# legend("topleft",legend =  levels(as.factor(gdf2$colour)), cex = 1)
# legend("topleft",pch = c(21), pt.bg =c(
#   "orange","brown", "orange3","lightblue","lightgrey",
#   "green","green3","darkolivegreen1","darkolivegreen","darkolivegreen3"
#   
# ),
#        legend =  levels(as.factor(gdf2$type)), cex = 1)
# 
# ## Plot ref to target
# msh <- mshape(gpa$coords)
# plotRefToTarget(M1 = msh, M2 = PCA$shapes$shapes.comp1$min, method = "points", #magn =3
#                 )
# 
# 
# # Method with surface
# findMeanSpec(gdf2$coords)
# refmesh <- read.ply("~/lizards/Analyses/100429-half.ply")
# ref3d <- warpRefMesh(refmesh, SM.log$LM[,,"100429"], mshIT, color = "lightgray", centered = FALSE)
# 
# shapePC1min <- plotRefToTarget(M1 = msh, M2 = PCA$shapes$shapes.comp1$min, mesh = ref3d, method = "surface")
# 
# l1 <- define.links(gdf2$coords[,,3 ])
# 
# plotRefToTarget(M1 = msh, M2 = PCA$shapes$shapes.comp1$max, method = "vector")
# 
# ############ With morpho
# 
# 
# defPC1min <- tps3d(ref3d, refmat = msh, tarmat = PCA$shapes$shapes.comp1$min)
# 
# deformGrid3d(msh,PCA$shapes$shapes.comp1$min,col1 = "lightgrey", col2 ="darkgreen", size = 0.002)
# require(rgl)
# #wire3d(defPC1min,col=4)
# shade3d(ref3d,col = skin1#, alpha = 0.5
#         )
# #hmapPC1min <- meshDist(ref3d, defPC1min, rampcolors = c("blue","white","yellow"),add = T)
# 
# 
# defPC1max <- tps3d(ref3d, refmat = msh, tarmat = PCA$shapes$shapes.comp1$max)
# deformGrid3d(msh,PCA$shapes$shapes.comp1$max,col1 = "lightgrey", col2 ="darkgreen", size = 0.002, lwd = 2)
# #hmapPC1max<- meshDist(ref3d, defPC1max, rampcolors = c("blue","white","yellow"),add = T)
# shade3d(ref3d,col = skin1 #, alpha = 0.5
#         )
# 
# ### PC2
# deformGrid3d(msh,PCA$shapes$shapes.comp2$min,col1 = "lightgrey", col2 ="darkgreen", size = 0.002)
# require(rgl)
# #wire3d(defPC1min,col=4)
# shade3d(ref3d,col = skin1#, alpha = 0.5
#         )
# deformGrid3d(msh,PCA$shapes$shapes.comp2$max,col1 = "lightgrey", col2 ="darkgreen", size = 0.002)
# require(rgl)
# #wire3d(defPC1min,col=4)
# shade3d(ref3d,col = skin1#, alpha = 0.5
# )
# 
# #####################################
# 
# 
# ######### Common alloletry
# 
# 
# names(gdf2$size) <- dimnames(gdf2$coords)[[3]]
# twobsize <- two.b.pls(A1 = gdf2$size, A2 = gdf2$coords,iter = 9999)
# 
# summary(twobsize)
# 
# p2bsize <- plot(twobsize)
# 
# p1 <- picknplot.shape(p2bsize)
# 
# plot(p2bsize$plot.args, pch =  c(21,22)[as.factor(gdf2$sex)], 
#      bg = c(         "orange","brown", "orange3","lightblue",
#          "green","green3","darkolivegreen1","darkolivegreen","darkolivegreen3"
#        )[as.factor(gdf2$type)]
# , cex = 1.3, ylab = "PC1",xlab = "Csize")
# 
# legend("bottomright",pch = c(21), pt.bg =c(
#   "orange","brown", "orange3","lightblue","lightgrey",
#   "green","green3","darkolivegreen1","darkolivegreen","darkolivegreen3"
#   
# ),
# legend =  levels(as.factor(gdf2$type)), cex = 1)
# 
#    
# # 
# # preds <- shape.predictor(gdf2$coords, gdf2$size,
# #                           Intercept = FALSE,
# #                           method = "PLS",
# #                        pred1 = min(p2bsize$plot.args$x), pred2 = max(p2bsize$plot.args$x))
# # 
# # plotRefToTarget(M1 = msh, M2 =p1$shapes[[1]], mesh = ref3d, method = "point",links = l1)
# # plotRefToTarget(M1 = msh, M2 =p1$shapes[[2]], mesh = ref3d, method = "point",links = l1)
# # 
# # plotRefToTarget(M1 = msh, M2 = preds$pred1, mesh = ref3d, method = "TPS",links = l1)
# 
# gdf2$type_sex <- paste(gdf2$type,gdf2$sex, sep = "_")
#       
# traj1 <- procD.lm(gdf2$coords ~ log(size)*type_sex, data = gdf2, iter = 1999)
# 
# reveal.model.designs(traj1)
# TA1 <- trajectory.analysis(traj1, groups = gdf2$type_sex, traj.pts = gdf2$size)



###########################################################
################## Males & IT #############################
###########################################################

######## GPA

gdfmIT <- geomorph.data.frame(coords = gpa$coords[,, df2$sex == "M" & df2$lineage == "IT"],
                              size = gpa$Csize[df2$sex == "M" & df2$lineage == "IT"], 
                              phenotype = df2$phenotype[df2$sex == "M" & df2$lineage == "IT"], 
                              SVL = df2$SVL[df2$sex == "M" & df2$lineage == "IT"], 
                              origin = df2$origin[df2$sex == "M" & df2$lineage == "IT"],
                              locality = df2$Location[df2$sex == "M" & df2$lineage == "IT"],
                              ID = df2$Specimen_id[df2$sex == "M" & df2$lineage == "IT"],
                              spec = df2$spec[df2$sex == "M" & df2$lineage == "IT"],
                              P_ordi = df2$phenotype_ordinal[df2$sex == "M" & df2$lineage == "IT"],
                              P_ordi_spec =df2$phenot_ordi_SpecCheck[df2$sex == "M" & df2$lineage == "IT"])

plotAllSpecimens(gdfmIT$coords)


PCAmIT.all <- gm.prcomp(gdfmIT$coords)

# Set color palette
pal1 <- brewer.pal(5,"YlGn")

# Generate mesh for consensus shape
findMeanSpec(gdf2$coords)
refmesh <- read.ply("~/lizards/Analyses/reference-mesh-100429.ply")
ref3d <- warpRefMesh(refmesh, SM.log$LM[,,"100429"], mshIT, color = "lightgray", centered = FALSE)


means.raw <- aggregate(two.d.array(gdfmIT$coords) ~ as.factor(gdfmIT$P_ordi), FUN=mean)
means.mn <- matrix(as.numeric(means.raw[1,-1]), ncol=3, byrow=T)

mshIT <- mshape(gdfmIT$coords)


### Exploratory analyses of IT males only #######

# p1.all <- plot(PCAmIT.all,axis1 =1,axis2 = 2,  pch =  c(21,22)[as.factor(gdfmIT$origin)], 
#                bg = pal1[gdfmIT$P_ordi], cex = 1.5,
#                main = "", xlab = "PC1 (14.8%)", ylab = "PC2 (12.5%)")
# 
# gl.w<-levels(as.factor(gdfmIT$P_ordi))
# for(i in 1:length(gl.w)) {
#   x<- PCAmIT.all$x[gdfmIT$P_ordi == gl.w[i], ]
#   CHP<-chull(x)
#   CHP<-c(CHP,CHP[1])
#   points(x[CHP,],type='l',col= pal1[i], lwd = 2)
# }

# 
# # Shape changes in specimens of score x (change raw line in means.raws[])
# deformGrid3d(mshIT,
#              matrix(as.numeric(means.raw[5,-1]), ncol=3, byrow=T),
#              col1 = "lightgrey", col2 ="darkgreen", size = 0.002)
# #wire3d(defPC1min,col=4)
# shade3d(ref3d,col = skin1#, alpha = 0.5
# )


### Ancestral vs. Nigri only

# gdfmIT$type <- NA
# gdfmIT$type[gdfmIT$phenotype == "ancestral" | gdfmIT$phenotype == "nigriventris"] <- "typical"
# gdfmIT$type[!complete.cases(gdfmIT$type)]  <- "other"
# 
# typical.males <- coords.subset(m0$GM$residuals, gdfmIT$type)
# 
# gdfmIT$typeNcap <- NA
# gdfmIT$typeNcap[gdfmIT$type == "typical" | gdfmIT$origin == "captivity"] <- "typeNcap"
# gdfmIT$typeNcap[!complete.cases(gdfmIT$typeNcap)]  <- "other"
# 
# typicNcap <- coords.subset(m0$GM$residuals, gdfmIT$typeNcap)
# capNint <- coords.subset(typicNcap$typeNcap, gdfmIT$phenotype[gdfmIT$typeNcap == "typeNcap"])
# 
# typeNcap.raw <-  coords.subset(gdfmIT$coords, gdfmIT$typeNcap)
# m.mIT.typeNcap <- procD.lm(typeNcap.raw$typeNcap ~ log(gdfmIT$size[gdfmIT$typeNcap == "typeNcap"]) + 
#                     gdfmIT$origin[gdfmIT$typeNcap == "typeNcap"] + 
#                     gdfmIT$P_ordi[gdfmIT$typeNcap == "typeNcap"] , iter= 9999,RRPP = T)
# 
# 
# m.mIT.typeNcap0 <- procD.lm(typeNcap.raw$typeNcap ~ log(gdfmIT$size[gdfmIT$typeNcap == "typeNcap"]) + 
#                               gdfmIT$origin[gdfmIT$typeNcap == "typeNcap"] , iter= 9999,RRPP = T)
#   
# PWtypeNcap <- pairwise(m.mIT.typeNcap, m.mIT.typeNcap0, groups = as.factor(gdfmIT$P_ordi[gdfmIT$typeNcap == "typeNcap"])#,covariate = log(gdfmIT$size)
# )
# summary(PWtypeNcap, test.type = "dist", confidence = 0.95) # 3 and 5 don't realy differ
# 
# # Plot PCA
# dots <- list()
# retx <- TRUE
# scale. <- dots$scale. <- FALSE
# center <- TRUE
# tol <- dots$tol
# 
# x.w.typ <- two.d.array(typical.males$typical)
# 
# PCAmIT <- prcomp(x.w.typ, center = center, scale. = scale., retx = retx, 
#                     tol = NULL)
# 
# PCAmIT2 <- gm.prcomp(typicNcap$typeNcap)
# 
# par(mfrow = c(1,2))
# p1 <- plot(PCAmIT,axis1 =1,axis2 = 2,  pch =  c(21,22)[as.factor(gdfmIT$origin[gdfmIT$type == "typical"])], 
#            bg = c(
#              "#FFFFCC","#006837" )[as.factor(gdfmIT$P_ordi[gdfmIT$type == "typical"])], cex = 1.5,
#            main = "", xlab = "PC1 (19.6%)", ylab = "PC2 (11.1%)")
# #text(p1$PC.points[,1], p1$PC.points[,2], labels = gdfmIT$locality, cex = 0.5, col = "grey")
# legend("bottomright",pch = c(21,21,21,22), pt.bg =c(
#   "#FFFFCC","#006837","black","black" ),
#   legend =  c("Ances","nigri","captive","wild"), cex = 0.8, pt.cex = 1)
# 
# plot(PCAmIT2,axis1 =1,axis2 = 2,  pch =  c(21,22)[as.factor(gdfmIT$origin[gdfmIT$typeNcap == "typeNcap"])], 
#      bg = c( "#FFFFCC","green","#006837")[as.factor(gdfmIT$P_ordi[gdfmIT$typeNcap == "typeNcap"])], cex = 1.5,
#      main = "", xlab = "PC1 (19.6%)", ylab = "PC2 (11.1%)")
# #text(p1$PC.points[,1], p1$PC.points[,2], labels = gdfmIT$locality, cex = 0.5, col = "grey")
# legend("bottomright",pch = c(21,21,21,22), pt.bg =c(
#   "#FFFFCC","#006837","black","black" ),
#   legend =  c("Ances","nigri","captive","wild"), cex = 0.8, pt.cex = 1)
# 
# ### Plot "intermediate" into tangent space of typical
# 
# 
# pcdata1 <- PCAmIT$x
# 
# plot(pcdata1[,1],pcdata1[,2],  pch =  c(21,22)[as.factor(gdfmIT$origin[gdfmIT$type == "typical"])], 
#      bg = c(
#        "#FFFFCC","#006837" )[as.factor(gdfmIT$P_ordi[gdfmIT$type == "typical"])], cex = 1.5,
#      main = "", xlab = "PC1 (19.6%)", ylab = "PC2 (11.1%)",xlim = c(-0.03,0.03))
# 
# gl.w<-levels(as.factor(gdfmIT$P_ordi[gdfmIT$type == "typical"]))
# for(i in 1:length(gl.w)) {
#   x<- pcdata1[gdfmIT$P_ordi[gdfmIT$type == "typical"] == gl.w[i], ]
#   CHP<-chull(x)
#   CHP<-c(CHP,CHP[1])
#   points(x[CHP,],type='l',col= c(
#     "#C2E699","#006837" )[i], lwd = 2)
# }
# 
# td.int <- two.d.array(capNint$intermediate)
# 
# int.res <- predict(PCAmIT, td.int)
# points(int.res[,1],int.res[,2],pch = 21, cex = 1.5,bg = "#78C679")
# 
# CHP<-c(chull(int.res ),chull(int.res )[1])
# points(int.res[CHP,],type='l',col= "#78C679", lwd = 1.2)
# 


##### Testing for differences among groups in mean shape and allometric trajectories ########
##### MANCOVA ######

## All IT males

### Full model
m.mIT <- procD.lm(coords ~ origin + log(size)*as.factor(P_ordi) , data = gdfmIT, iter= 9999,RRPP = T)

sumIT <- summary(m.mIT)



# pairwise comparisons
m.red <- procD.lm(coords ~ origin + log(size) + as.factor(P_ordi) , data = gdfmIT, iter= 9999,RRPP = T) # reduced model 1
m.redSizeCap <- procD.lm(coords ~ origin + log(size), data = gdfmIT, iter= 9999,RRPP = T) # reduced model 2

t1 <- anova.lm.rrpp(m.red,m.mIT)
summary(t1)

m0 <- procD.lm(coords ~ origin , iter= 9999, data = gdfmIT,RRPP = T) #Null model

PW <- pairwise(m.red, m0, groups = as.factor(gdfmIT$P_ordi)) # Paiwise comparisons for Procrustes distances among Ordinal scores

PW2 <- pairwise(m.mIT, m.red, groups = as.factor(gdfmIT$P_ordi) ,covariate = log(gdfmIT$size)) # # Paiwise comparisons for slopes    
PW3 <- pairwise(m.mIT, m.redSizeCap, groups = as.factor(gdfmIT$P_ordi)) # To compare slopes    
sum.angles <- summary(PW2, test.type = "VC", confidence = 0.95)
sum.mean <- summary(PW3, test.type = "dist", confidence = 0.95)

# Export tables for later analyses with Lacertidae
write.table(sum.angles$summary.table,file = "G:/results/pairwise-podarcis-angles.txt")
write.table(sum.mean$summary.table,file = "G:/results/pairwise-podarcis-dist.txt", sep = ",")
write.csv(sumIT$table, file = "G:/results/Podarcis-IT-RRPP-output.txt")

##### Exploratory - Plot allometry Males IT #####

# plotAllomPred <- plotAllometry(m.mIT, size = gdfmIT$size, logsz = T, method = "RegScore",
#                                             pch = c(21,22)[as.factor(gdfmIT$origin)], 
#                                             bg =pal1[as.factor(gdfmIT$P_ordi)],
#                                             main = "")
# legend("topright",pch = c(21,21,21,21,21,22,22), pt.bg =c(
#   pal1,"black","black"
# ),
# legend =  c(levels(as.factor(gdfmIT$P_ordi)),levels(as.factor(gdfmIT$origin))), cex = 0.8)
# 
# 
# preds <- shape.predictor(m.mIT$GM$fitted, 
#                          x= plotAllomPred$RegScore, Intercept = F,
#                          predmin = min(plotAllomPred$RegScore),
#                          predmax = max(plotAllomPred$RegScore))
# 
# deformGrid3d(msh,preds$predmax,col1 = "lightgrey", col2 ="darkgreen", size = 0.002)
# require(rgl)
# #wire3d(defPC1min,col=4)
# shade3d(ref3d,col = skin1#, alpha = 0.5
# )
# 

##### Same analyses but taking with residuals from regression of shape on Size & Origin instead of raw data

coord.res <- m.redSizeCap$GM$residuals
for( i in 1:length(gdfmIT$ID)){
  coord.res[,, i] <- m.redSizeCap$GM$residuals[,, i] + mshIT
}

PCAm.resSize <- gm.prcomp(coord.res)
pPCAm.resSize<- plot(PCAm.resSize,axis1 =1,axis2 = 2,  pch =  21, 
     bg = pal1[as.factor(gdfmIT$P_ordi)], cex = 1.5,
     main = "", xlab =  "PC1  (16.1%)", ylab = "PC2  (9.7%)")

pal2 <- pal1
pal2[1] <- 'lightgrey'

gl.w<-c("1",NULL,NULL,NULL,"5")
for(i in 1:length(gl.w)) {
  x<- PCAm.resSize$x[gdfmIT$P_ordi == gl.w[i], ]
  CHP<-chull(x)
  CHP<-c(CHP,CHP[1])
  points(x[CHP,],type='l',col= pal2[i], lwd = 2)
}

#text(pPCAm.resSize$PC.points[,1], pPCAm.resSize$PC.points[,2], labels = gdfmIT$ID, cex = 1, col = "grey")
# there seems to be nothing wrong with PS31

# 
# deformGrid3d(msh,PCAm.resSize$shapes$shapes.comp2$min,col1 = "lightgrey", col2 ="darkgreen", size = 0.002)
# shade3d(ref3d,col = skin1#, alpha = 0.5
# )
# 
# deformGrid3d(msh,PCAm.resSize$shapes$shapes.comp2$max,col1 = "lightgrey", col2 ="darkgreen", size = 0.002)
# shade3d(ref3d,col = skin1#, alpha = 0.5
# )


## PCA residuals from regession of shape on Origin only

coord.res.m0 <- m0$GM$residuals
for( i in 1:length(gdfmIT$ID)){
  coord.res.m0[,, i] <- coord.res.m0[,, i] + mshIT
}

PCAm0.resSize <- gm.prcomp(coord.res.m0)
pPCAm0.resSize<- plot(PCAm0.resSize,axis1 =1,axis2 = 2,  pch =  21, 
                     bg = pal1[as.factor(gdfmIT$P_ordi)], cex = 1.5,
                     main = "", xlab =  "PC1  (15.4%)", ylab = "PC2  (11.3%)")


gl.w<-c("1",NULL,NULL,NULL,"5")
for(i in 1:length(gl.w)) {
  x<- PCAm0.resSize$x[gdfmIT$P_ordi == gl.w[i], ]
  CHP<-chull(x)
  CHP<-c(CHP,CHP[1])
  points(x[CHP,],type='l',col= pal2[i], lwd = 2)
}

#text(pPCAm.resSize$PC.points[,1], pPCAm.resSize$PC.points[,2], labels = gdfmIT$ID, cex = 1, col = "grey")
# there seems to be nothing wrong with PS31
# 
# 
# deformGrid3d(msh,PCAm0.resSize$shapes$shapes.comp2$min,col1 = "lightgrey", col2 ="darkgreen", size = 0.002)
# shade3d(ref3d,col = skin1#, alpha = 0.5
# )
# 
# deformGrid3d(msh,PCAm0.resSize$shapes$shapes.comp2$max,col1 = "lightgrey", col2 ="darkgreen", size = 0.002)
# shade3d(ref3d,col = skin1#, alpha = 0.5
# )

# 
# #### 
# names(gdfmIT$size) <- gdfmIT$ID
# tbO <- two.b.pls(A1 = gdfmIT$size, A2 = gdfmIT$coords,iter = 9999)
# 
# summary(twobsize)
# 
# tb0plot <- plot(twobsize)
# 
# plot(tb0plot$plot.args, pch =  c(21,22)[as.factor(gdfmIT$origin)], 
#      bg = c(
#        "orange3","lightgrey", "darkolivegreen1" )[as.factor(gdfmIT$phenotype)]
#      , cex = 1.3, ylab = "PC1",xlab = "Csize")
# 
# legend("bottomright",pch = c(21), pt.bg =c(
#   "orange","brown", "orange3","lightblue","lightgrey",
#   "green","green3","darkolivegreen1","darkolivegreen","darkolivegreen3"
#   
# ),
# legend =  levels(as.factor(gdf2$type)), cex = 1)


# # Shape changes along PC1 of fitted values
# 
# deformGrid3d(msh,PredLmax$shapes[[2]],
#              col1 = "lightgrey", col2 ="darkgreen", size = 0.002, lwd = 2)
# shade3d(ref3d,col = skin1 #, alpha = 0.5
# )
# 
# rgl.snapshot('~/lizards/Analyses/PreLmax-lat.png', fmt = 'png')


###################### With Ordinal scores  ########################################


# 
# pd <-position_jitter(h=0.15,w=0.15)

#### Exploratory - Compare Ordinal scores vs. spec values 
# ggplot(data = origin, aes(x = phenotype_ordinal,y = log(spec))) +
#   geom_point(position = pd) + xlab("Ordinal score") + ylab("log(Spec greenness)")
# 
# ### Assigned ordinal scores vs. transformed spec data into ordinal
# ggplot(data = origin, aes(x = phenotype_ordinal,y = phenot_ordi_SpecCheck)) +
#   geom_point(position = pd) + xlab("Ordinal score") + ylab("Ordinal transformed greenness")



### Plot captive only

## Exploratory - Ordinal scores from spec
# p.PCAmITord <- plot(PCAm.resSize$x[gdfmIT$origin == "captivity",1],PCAm.resSize$x[gdfmIT$origin == "captivity",2],
#      pch =  21, 
#      bg = pal1[as.factor(gdfmIT$P_ordi_spec[gdfmIT$origin == "captivity"])],
#  cex = 1.5,
#      main = "", xlab = "PC1 (17.5%)", ylab = "PC2 (9.34%)")
# 
# 
# legend("bottomright",pch = 21, pt.bg =pal1,
# legend =  levels(as.factor(gdfmIT$P_ordi_spec)), cex = 0.8, pt.cex = 2)

# ### Spec
# library(classInt)
# int=classIntervals(gdfmIT$spec[gdfmIT$origin == "captivity"],5)
# #int=classIntervals(dat$DOC,n,style="fixed",fixedBreaks=seq(0,6,l=2));#if you want specific breaks
# col=pal1[findInterval(gdfmIT$spec[gdfmIT$origin == "captivity"],int$brks,all.inside = T)]
# plot(PCAm.resSize$x[gdfmIT$origin == "captivity",1],PCAm.resSize$x[gdfmIT$origin == "captivity",2],
#      pch =  21, 
#      bg = col,
#      cex = 1.5,
#      main = "Continuous spec data", xlab = "PC1 (17.5%)", ylab = "PC2 (9.34%)") 
# legend("bottomright",pch = 21, pt.bg =c(col[gdfmIT$spec[gdfmIT$origin == "captivity"] == min(gdfmIT$spec[gdfmIT$origin == "captivity"])],
#                                         col[gdfmIT$spec[gdfmIT$origin == "captivity"] == max(gdfmIT$spec[gdfmIT$origin == "captivity"])]),
#        legend =  c(round(min(gdfmIT$spec[gdfmIT$origin == "captivity"]),2),
#                    round(max(gdfmIT$spec[gdfmIT$origin == "captivity"]),2)), cex = 0.8, pt.cex = 2)



### CVA analysis 

## With residuals from regression on Size & Origin

cva.1 <- CVA(dataarray = coord.res, groups = gdfmIT$P_ordi, cv = TRUE)

cva.1$Var
par(font.lab = 2)
plot(cva.1$CVscores[,1:2], bg=pal1[gdfmIT$P_ordi], pch=21, typ="p",asp=1, cex = 1.2, #xlim = c(-6,10),
     xlab = paste("CV1 (",round(cva.1$Var[1,2],1), "%)", sep = ""),
     ylab = paste(paste("CV2 (",round(cva.1$Var[2,2],1), "%)", sep = "")))

# add 95% ellipses
require(car)
for(ii in 1:length(levels(as.factor(gdfmIT$P_ordi)))){
  dataEllipse(cva.1$CVscores[as.factor(gdfmIT$P_ordi)==levels(as.factor(gdfmIT$P_ordi))[ii],1],
              cva.1$CVscores[as.factor(gdfmIT$P_ordi)==levels(as.factor(gdfmIT$P_ordi))[ii],2], 
              add=TRUE,levels=.95, col=pal2[ii], center.pch = NULL, plot.points = F)
  }


if (require(lattice)) {
  histogram(~cva.1$CVscores[,1]|as.factor(gdfmIT$P_ordi),
            layout=c(1,length(levels(as.factor(gdfmIT$P_ordi)))),
            xlab=paste("1st canonical axis", paste(round(cva.1$Var[1,2],1),"%")))
  histogram(~cva.1$CVscores[,2]|gdfmIT$P_ordi, layout=c(1,length(levels(as.factor(gdfmIT$P_ordi)))),
            xlab=paste("2nd canonical axis", paste(round(cva.1$Var[2,2],1),"%")))
} 

cvvis5 <- mean(cva.1$CVscores[gdfmIT$P_ordi == "5", 1])*matrix(cva.1$CVvis[,1],nrow(cva.1$Grandm),ncol(cva.1$Grandm))+cva.1$Grandm
cvvis1 <- mean(cva.1$CVscores[gdfmIT$P_ordi == "1", 1])*matrix(cva.1$CVvis[,1],nrow(cva.1$Grandm),ncol(cva.1$Grandm))+cva.1$Grandm


deformGrid3d(cvvis1, cvvis5,col1 = pal1[1], col2 = pal1[4], size = 0.003 , )
shade3d(ref3d,col = skin1, alpha = 0.5
)
rgl.bg( sphere = FALSE, fogtype = "none", color = "black" )

snapshot3d("~/lizards/Analyses/CV1.png")



### CVA incorporating size information (not residuals) 

cvaM0 <- CVA(dataarray = coord.res.m0, groups = gdfmIT$P_ordi, cv = TRUE, round = 1000)

plot(cvaM0$CVscores[,1:2], bg=pal1[gdfmIT$P_ordi], pch=21, typ="p",asp=1, cex = 1.5, 
     xlim = c(-5,6),
     ylim = c(-4,5),
     xlab = paste("CV1 (",round(cvaM0$Var[1,2],1), "%)", sep = ""),
     ylab = paste(paste("CV2 (",round(cvaM0$Var[2,2],1), "%)", sep = "")))

# add 95% ellipses
require(car)
for(ii in 1:length(levels(as.factor(gdfmIT$P_ordi)))){
  dataEllipse(cvaM0$CVscores[as.factor(gdfmIT$P_ordi)==levels(as.factor(gdfmIT$P_ordi))[ii],1],
              cvaM0$CVscores[as.factor(gdfmIT$P_ordi)==levels(as.factor(gdfmIT$P_ordi))[ii],2], 
              add=TRUE,levels=.95, col=pal2[ii], center.pch = NULL, plot.points = F)
}


if (require(lattice)) {
  histogram(~cvaM0$CVscores[,1]|as.factor(gdfmIT$P_ordi),
            layout=c(1,length(levels(as.factor(gdfmIT$P_ordi)))),
            xlab=paste("1st canonical axis", paste(round(cva.1$Var[1,2],1),"%")))
  histogram(~cvaM0$CVscores[,2]|gdfmIT$P_ordi, layout=c(1,length(levels(as.factor(gdfmIT$P_ordi)))),
            xlab=paste("2nd canonical axis", paste(round(cva.1$Var[2,2],1),"%")))
} 

cvvisM05 <- mean(cvaM0$CVscores[gdfmIT$P_ordi == "5", 1])*matrix(cvaM0$CVvis[,1],nrow(cvaM0$Grandm),ncol(cvaM0$Grandm))+cvaM0$Grandm
cvvisM01 <- mean(cvaM0$CVscores[gdfmIT$P_ordi == "1", 1])*matrix(cvaM0$CVvis[,1],nrow(cva.1$Grandm),ncol(cvaM0$Grandm))+cvaM0$Grandm

deformGrid3d(cvvisM01, cvvisM05,col1 = "lightgrey", col2 ="darkgreen", size = 0.002)
shade3d(ref3d,col = skin1#, alpha = 0.5
)



## Exploratory - Try pooling scores 1 & 3 - With size & captive vs. wild residuals

# gdfmIT$P_ordi1.3 <- gdfmIT$P_ordi
# gdfmIT$P_ordi1.3[gdfmIT$P_ordi1.3 == "3"] <- "1"
# cva.2 <- CVA(dataarray = coord.res, groups = gdfmIT$P_ordi1.3, cv = TRUE)
# 
# cva.2$Var
# 
# plot(cva.2$CVscores[,1:2], bg=pal1[as.factor(gdfmIT$P_ordi)], pch=21, typ="p",asp=1, cex = 1.5)
# 
# for(ii in 1:length(levels(as.factor(gdfmIT$P_ordi1.3)))){
#   dataEllipse(cva.2$CVscores[as.factor(gdfmIT$P_ordi1.3)==levels(as.factor(gdfmIT$P_ordi1.3))[ii],1],
#               cva.2$CVscores[as.factor(gdfmIT$P_ordi1.3)==levels(as.factor(gdfmIT$P_ordi1.3))[ii],2], 
#               add=TRUE,levels=.95, col=pal2[ii], center.pch = NULL, plot.points = F)
# }
# text(cva.2$CVscores[,1], cva.2$CVscores[,2], labels = gdfmIT$locality, cex = 0.8, col = "grey")


##### Exploratory - CVA with males IT + La Spezia
# gdfmIT.S <- geomorph.data.frame(coords = gpa$coords[,, df2$sex == "M" & df2$lineage == "IT" | df2$phenotype == "green-Spezia"],
#                               size = gpa$Csize[df2$sex == "M" & df2$lineage == "IT" | df2$phenotype == "green-Spezia"], 
#                               origin = df2$origin[df2$sex == "M" & df2$lineage == "IT" | df2$phenotype == "green-Spezia"],
#                               ID = df2$Specimen_id[df2$sex == "M" & df2$lineage == "IT" | df2$phenotype == "green-Spezia"],
#                               lin = df2$lineage[df2$sex == "M" & df2$lineage == "IT" | df2$phenotype == "green-Spezia"],
#                               P_ordi = df2$phenotype_ordinal[df2$sex == "M" & df2$lineage == "IT" | df2$phenotype == "green-Spezia"])
# gdfmIT.S$lin.score <- paste(gdfmIT.S$P_ordi, gdfmIT.S$lin,sep="_")
# 
# m.redSizeCap2 <- procD.lm(coords ~ origin + log(size), data = gdfmIT.S, iter= 9999,RRPP = T)
# 
# coord.res2 <- m.redSizeCap2$GM$residuals
# for( i in 1:length(gdfmIT.S$ID)){
#   coord.res2[,, i] <- m.redSizeCap2$GM$residuals[,, i] + msh
# }
# 
# cva.IT.S <- CVA(dataarray = coord.res2, groups = gdfmIT.S$lin.score, cv = TRUE)
# 
# cva.IT.S$Var
# par(font.lab = 2)
# plot(cva.IT.S$CVscores[,1:2], bg=pal1[gdfmIT.S$P_ordi], pch=c(22,21)[as.factor(gdfmIT.S$lin)], typ="p",asp=1, cex = 1.2, #xlim = c(-6,10),
#      xlab = paste("CV1 (",round(cva.IT.S$Var[1,2],1), "%)", sep = ""),
#      ylab = paste(paste("CV2 (",round(cva.IT.S$Var[2,2],1), "%)", sep = "")))
# 
# for(ii in 1:length(levels(as.factor(gdfmIT.S$lin.score)))){
#   dataEllipse(cva.IT.S$CVscores[as.factor(gdfmIT.S$lin.score)==levels(as.factor(gdfmIT.S$lin.score))[ii],1],
#               cva.IT.S$CVscores[as.factor(gdfmIT.S$lin.score)==levels(as.factor(gdfmIT.S$lin.score))[ii],2], 
#               add=TRUE,levels=.95, col=c("lightgrey", "#C2E699",   "#78C679",   "#31A354","#31A354",   "#006837")[ii]
#               , center.pch = NULL, plot.points = F)
# }
# 
# 
# cvvisMS05 <- max(cva.IT.S$CVscores[, 1])*matrix(cva.IT.S$CVvis[,1],nrow(cva.IT.S$Grandm),ncol(cva.IT.S$Grandm))+cva.IT.S$Grandm
# cvvisMS01 <- min(cva.IT.S$CVscores[, 1])*matrix(cva.IT.S$CVvis[,1],nrow(cva.IT.S$Grandm),ncol(cva.IT.S$Grandm))+cva.IT.S$Grandm
# 
# deformGrid3d(cvvisMS01, cvvisMS05,col1 = pal1[1], col2 =pal1[4], size = 0.003 , )
# shade3d(ref3d,col = skin1#, alpha = 0.5
# )
# rgl.bg( sphere = FALSE, fogtype = "none", color = "black" )

############ Modularity test

#### IT
partition <- read.table("~/data-metadata/modules.txt", h =T)
modul1 <- modularity.test(gdfmIT$coords[complete.cases(partition$module),,],partition.gp = partition$module[complete.cases(partition$module)], iter = 9999)
plot(modul1, xlim(0.5,1.2))

### Same but considering landmark on the Parietal bone that are along the midline  as NC

partion2 <- partition
partion2$module[partion2$name == "par-ant-process" | partion2$name == "partiet-ant-lat-corner" |
                partion2$name == "par-ant-lat" | partion2$name == "par-and-mid"] <- NA
modul2 <- modularity.test(gdfmIT$coords[complete.cases(partion2$module),,],partition.gp = partition$module[complete.cases(partion2$module)], iter = 9999)
summary(modul2)

par(mfrow = c(2,1), font.lab = 2)
hist(modul2$random.CR,xlim = c(0.5, 1.5), main ="", xlab = "CR Coefficents")
abline(v = modul2$CR)
box()
hist(modul1$random.CR,xlim = c(0.8, 1.1), main ="", xlab = "CR Coefficents")
abline(v = modul1$CR)
box()
dev.off()


##### Pairwise distances NC vs. Mesoderm - IT

m.mIT.NC <- procD.lm(gdfmIT$coords[partition$module == "NC" & complete.cases(partition$module),,]  ~ origin + log(size)*as.factor(P_ordi) , 
                       data = gdfmIT, iter= 9999,RRPP = T)
m.mIT.NCred <- procD.lm(gdfmIT$coords[partition$module == "NC" & complete.cases(partition$module),,]  ~ origin + log(size) ,
                          data = gdfmIT, iter= 9999,RRPP = T)
PW.NC <- pairwise(m.mIT.NC, m.mIT.NCred, groups = as.factor(gdfmIT$P_ordi))

summary(PW.NC, confidence = 0.95, test.type = "dist") 

m.mIT.mes <- procD.lm(gdfmIT$coords[partition$module == "mes" & complete.cases(partition$module),,] ~ origin + log(size)*as.factor(P_ordi) , data = gdfmIT, iter= 9999,RRPP = T)
m.mIT.mesRed <- procD.lm(gdfmIT$coords[partition$module == "mes" & complete.cases(partition$module),,] ~ origin + log(size) , data = gdfmIT, iter= 9999,RRPP = T)
PW.mes <- pairwise(m.mIT.mes, m.mIT.mesRed, groups = as.factor(gdfmIT$P_ordi))
summary(PW.mes, confidence = 0.95, test.type = "dist")

# ##### Exploratory - Pairwise distances NC vs. Mesoderm - IT + Spezia
# 
# m.mIT.S.NC <- procD.lm(gdfmIT.S$coords[partition$module == "NC" & complete.cases(partition$module),,]  ~ origin + log(size)*as.factor(lin.score) , 
#                        data = gdfmIT.S, iter= 9999,RRPP = T)
# m.mIT.S.NCred <- procD.lm(gdfmIT.S$coords[partition$module == "NC" & complete.cases(partition$module),,]  ~ origin + log(size) ,
#                           data = gdfmIT.S, iter= 9999,RRPP = T)
# PW.NC <- pairwise(m.mIT.S.NC, m.mIT.S.NCred, groups = as.factor(gdfmIT.S$lin.score))
# 
# summary(PW.NC, confidence = 0.95, test.type = "dist") 
# 
# m.mIT.S.mes <- procD.lm(gdfmIT.S$coords[partition$module == "mes" & complete.cases(partition$module),,] ~ origin + log(size)*as.factor(lin.score) , data = gdfmIT.S, iter= 9999,RRPP = T)
# m.mIT.S.mesRed <- procD.lm(gdfmIT.S$coords[partition$module == "mes" & complete.cases(partition$module),,] ~ origin + log(size) , data = gdfmIT.S, iter= 9999,RRPP = T)
# PW.mes <- pairwise(m.mIT.S.mes, m.mIT.S.mesRed, groups = as.factor(gdfmIT.S$lin.score))
# summary(PW.mes, confidence = 0.95, test.type = "dist") 



####### Phenotypic distance plot ####

dist <- c(0.015
, 0.011
 )
Origin <- c("NC","Mesoderm")
df <- data.frame(Origin,dist)

par(font.lab = 2, cex.lab = 1.4, oma = c(0,1,0,0))
barplot(df$dist, names.arg = c("Neural crest","Mesoderm"), width = 0.8, ylab = "Phenotypic distance",
        xlab = "Module",
        col = c("lightblue","orange1"), ylim = c(0,0.02), xlim = c(0,2),axis.lty = 1,
        cex.axis = 1,
        cex.names = 1)
abline(h =0)

##### Per-landmark variance

var.mIT <- per_lm_variance(shape.data = gdfmIT$coords)

spheres3d(mshIT, col = var.mIT$Variance_Colors, radius =  0.005)
shade3d(ref3d,col = skin1, alpha = 0.5
)
rgl.bg( sphere = FALSE, fogtype = "none", color = "black" )
#snapshot3d("~/lizards/Analyses/Podarcis-hotdots-vent.png")

str(var.mIT)
var.mIT$Log_Variance

# Reorder colors according to sorted Log_Variance
sorted_indices <- order(var.mIT$Log_Variance)
sorted_log_variance <- var.mIT$Log_Variance[sorted_indices]
sorted_colors <- var.mIT$Variance_Colors[sorted_indices]

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


var.table <- data.frame(lm = as.factor(partition$name), module = partition$module, var = var.mIT$Per_Lm_Variance)
summary(var.table)

var.table$module[is.na(var.table$module)] <- "unasigned"
var.table$module <- factor(var.table$module, levels = c("NC","mes","unasigned"))
var.table <- var.table[order(var.table$module, var.table$var), ]

lmOrder <- var.table
# saveRDS(var.table, "G:/results/figures/Fig2orSM/lmOrder.rds") # Export from downstream analyses with Lacertids

colors <- ifelse(var.table$module == 'mes', 'orange1', ifelse(var.table$module == 'NC', 'cornflowerblue', 'white'))

par(mar=c(10,5,4,4))
barplot(var.table$var, names.arg = var.table$lm, las = 2, col = colors, cex.names = 0.5,
        xlab = '', ylab = '', 
        main = '')


legend('topright', legend = c('mes', 'NC', 'unasigned'), fill = c('orange1', 'lightblue', 'white'))
