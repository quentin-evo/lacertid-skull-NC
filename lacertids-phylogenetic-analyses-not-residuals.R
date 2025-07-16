
###### CVa Phylo not corrected for anything

cva.raw <- CVA(dataarray = La.gdfmIT$coords, groups = La.gdfmIT$P_ordi, cv = TRUE)

La.pred.raw <- predict(cva.raw,newdata = coord.Lacertid)

plot(cva.raw$CVscores[,1:2], bg=pal1[La.gdfmIT$P_ordi], pch=21, typ="p",asp=1, cex = 1.2, 
     xlim = c(-10,10),
     ylim = c(-8,8),
     xlab = paste("CV1 (",round(cva.raw$Var[1,2],1), "%)", sep = ""),
     ylab = paste(paste("CV2 (",round(cva.raw$Var[2,2],1), "%)", sep = "")))
points(La.pred.raw[,1:2], pch=21, bg = alpha("grey", 0.2), typ="p",asp=1, cex = 0.8)


red.tree.raw<-treedata(La_tree,La.pred.raw)
red.tree.raw$data <- red.tree.raw$data[!duplicated(rownames(red.tree.raw$data)), ]

### BM model
Lac_BM.raw <-fitContinuous(phy = red.tree.raw$phy, dat = red.tree.raw$data[, "CV 1"], model = "BM")

### Early burst model
fitEB_gs.raw <-fitContinuous(phy = red.tree.raw$phy, dat = red.tree.raw$data[, "CV 1"], model ="EB")

#### OU model
fitOU_gs.raw <-fitContinuous(phy = red.tree.raw$phy, dat = red.tree.raw$data[, "CV 1"], model ="OU", ,bounds=list(alpha=c(0,50)))
fitOU_gs.raw

aic_gs<-setNames(c(AIC(Lac_BM.raw),
                   AIC(fitEB_gs.raw),AIC(fitOU_gs.raw)),
                 c("BM","EB","OU"))
par(mfrow = c(2,2))
plotTree.barplot(tree = red.tree$phy,red.tree$data[, "CV1"],
                 args.plotTree=list(fsize=0.5),
                 args.barplot=list(#col=tip.cols,
                   xlab="CVA 1",
                   cex.lab=0.8))
plotTree.barplot(tree = red.tree.raw$phy,red.tree.raw$data[, "CV 1"],
                 args.plotTree=list(fsize=0.5),
                 args.barplot=list(#col=tip.cols,
                   xlab="CVA 1",
                   cex.lab=0.8))

cv1.raw <- red.tree.raw$data[,"CV 1"]

par(bty="n") # set plot box
# Set priors
priorOU.raw<-make.prior(red.tree.raw$phy,
                    dists=list(dalpha="dhalfcauchy",
                               dsig2="dhalfcauchy",
                               dk="cdpois",dtheta="dnorm"),
                    param=list(dalpha=list(scale=1),
                               dsig2=list(scale=1),
                               dk=list(lambda=10, kmax=50),
                               dsb=list(bmax=1, prob=1),
                               dtheta=list(mean=mean(cv1.raw),
                                           sd=1.5*sd(cv1.raw))),
                    plot.prior=TRUE)


# set up to MCMC run
mcmcOU.raw<-bayou.makeMCMC(red.tree.raw$phy,cv1.raw,
                       prior=priorOU,plot.freq=NULL,file.dir=NULL,
                       ticker.freq=100)

# run MCMC
raw.rjMCMC<-mcmcOU.raw$run(200000)

plot(raw.rjMCMC)

raw.rjMCMC.burn<-set.burnin(lac.rjMCMC,0.20)

plot(raw.rjMCMC.burn)

lac.rjMCMC.res <- summary(lac.rjMCMC.burn)
lac.rjMCMC.res$statistics



# Plot output

par(mar=c(1.1,1.1,3.1,0.1),mfrow=c(1,2)) # ^set panels
plotSimmap.mcmc(raw.rjMCMC.burn,edge.type="regimes",
                lwd=2,pp.cutoff=0.25,cex=0.6) # regimes shifts with PP>0.5

mtext("(a)",adj=0,line=1)

plotSimmap.mcmc(raw.rjMCMC.burn,edge.type="theta",
                pp.labels = T,
                pp.col = "lightgrey",
                pp.cex= 0.5,
                lwd=1,pp.cutoff=0.25,cex=0.6
                # ,legend_settings=
                #   list(x=0.2*max(nodeHeights(red.tree$phy)),
                #        y=0.7*Ntip(red.tree$phy))
) # mean theta
mtext("(b)",adj=0,line=1)



##### CVa Phylo scores with species data size corrected beside Podarcis 
#### It seems like one cannot do that as the variance is shifted to the upward around the Lacertid average size

LaSizeRes <- procD.lm(gpa.all$coords[,, df.all2$Species %in% la.species] ~ log(gpa.all$Csize[df.all2$Species %in% la.species]))

plot(LaSizeRes, type = "diagnostic")
msh.la <- mshape(gpa.all$coords[,, df.all2$Species %in% la.species])

all.res <- LaSizeRes$GM$residuals


for( i in 1:length(length(gpa.all$Csize[df.all2$Species %in% la.species]))){
  all.res[,, i] <-LaSizeRes$GM$residuals[,, i] + msh.la
}


La.pred.SizRes <- predict(cva.la,newdata = all.res)


plot(cva.la$CVscores[,1:2], bg=pal1[La.gdfmIT$P_ordi], pch=21, typ="p",asp=1, cex = 1.2, 
     xlim = c(-25,25),
     ylim = c(-12,12),
     xlab = paste("CV1 (",round(cva.la$Var[1,2],1), "%)", sep = ""),
     ylab = paste(paste("CV2 (",round(cva.la$Var[2,2],1), "%)", sep = "")))
points(La.pred.SizRes[,1:2], pch=21, bg = alpha("grey", 0.2), typ="p",asp=1, cex = 0.8)

##### CVa Phylo scores with species data size corrected altogether with Podarcis 

linS <- paste( df.all2$lineage, df.all2$sex, sep = "_")
allSizeRes <- procD.lm(gpa.all$coords[,, df.all2$Species %in% la.species | linS == "IT_M"] ~ log(gpa.all$Csize[df.all2$Species %in% la.species | linS == "IT_M"]))


msh.all <- mshape(gpa.all$coords[,, df.all2$Species %in% la.species | linS == "IT_M"])

all.res2 <- allSizeRes$GM$residuals


for( i in 1:length(gpa.all$Csize[df.all2$Species %in% la.species | linS == "IT_M"])){
  all.res2[,, i] <-allSizeRes$GM$residuals[,, i] + msh.all
}

sel <- df.all2[df.all2$Species %in% la.species |linS == "IT_M" ,]
all.res2P <- all.res2[,, sel$lineage == "IT"]

cva.la2 <- CVA(dataarray = all.res2P, groups = sel$phenotype_ordinal[sel$lineage == "IT"], cv = TRUE)
La.pred.SizRes2 <- predict(cva.la2,newdata = all.res2)


plot(cva.la2$CVscores[,1:2], bg=pal1[as.factor(sel$phenotype_ordinal)], pch=21, typ="p",asp=1, cex = 1.2, 
     xlim = c(-25,25),
     ylim = c(-12,12),
     xlab = paste("CV1 (",round(cva.la$Var[1,2],1), "%)", sep = ""),
     ylab = paste(paste("CV2 (",round(cva.la$Var[2,2],1), "%)", sep = "")))
points(La.pred.SizRes2[,1:2], pch=21, bg = alpha("grey", 0.2), typ="p",asp=1, cex = 0.8)