library(SlicerMorphR)
library(geomorph)
library(tidyr)
library(paleomorph)
library(Morpho)
library(SamplingStrata)
library(ggplot2)
library(RColorBrewer)
library(hot.dots)
library(phytools)
library(geiger)
library(bayou)
library(viridisLite)
library(car)

# Previous objects were generated with lacertids-phylogenetic-analyses.R

rate.la <- per_lm_rates(shape.data = coord.Lacertid[,,sel.spec != "Darevskia_obscura"
                                                    & sel.spec != "Eremias_nigrocellata"
                                                    & sel.spec != "Philochortus_hardeggeri"], phy = red.tree$phy)
variance.la <- per_lm_variance(shape.data = coord.Lacertid[,,sel.spec != "Darevskia_obscura"
                                                           & sel.spec != "Eremias_nigrocellata"
                                                           & sel.spec != "Philochortus_hardeggeri"])

ances <- mshape(gdfmIT$coords[,,gdfmIT$P_ordi =="1"])
nigri <- mshape(gdfmIT$coords[,,gdfmIT$P_ordi =="5"])

ANDist <- (sqrt(rowSums((ances - nigri)^2))) # Euclidian distance between landmarks of Ancestral and Nigriventris 5


d.rates <- data.frame(lm = partition$name,   
                      var_la = variance.la$Per_Lm_Variance,
                      rates_la = rate.la$Per_Lm_Rates,
                      partition_module = partition$module,
                      lmdistP = ANDist
)


lmDistVar <- lm(var_la ~ partition_module*lmdistP, data = d.rates[complete.cases(d.rates$partition_module),])
summary(lmDistVar)
plot(lmDistVar, which = 4)
anova(lmDistVar)

lmDistR <- lm(rates_la ~ partition_module*lmdistP, data = d.rates[complete.cases(d.rates$partition_module),])
summary(lmDistR)
anova(lmDistR)


fac.mod <- as.factor(d.rates$partition_module)

par(mfrow = c(1,2))

plot(d.rates$lmdistP, d.rates$var_la, pch = 21, bg = c("orange","blue","grey")[fac.mod], cex = 1,
     xlab = "Euclidean distance in nigriventris", ylab = "Variance in Lacertidae", xlim = c(0, 0.005))

for (i in 1:2) {
  subset_data <- d.rates[d.rates$partition_module == levels(as.factor(d.rates$partition_module))[i] & complete.cases(d.rates$partition_module),]
  
  # Fit a linear model
  lm_fit <- lm(var_la ~ lmdistP, data = subset_data)
  
  # Add the regression line to the plot
  abline(lm_fit, col = c("orange2","cornflowerblue")[i], lwd = 2)
}

###Evolutionary rates
plot(d.rates$lmdistP, d.rates$rates_la, pch = 21, bg = c("orange","blue","grey")[fac.mod],cex = 1,
     xlab = "Euclidean distance in nigriventris", ylab = "Evolutionary rate in Lacertidae", xlim = c(0, 0.005))

for (i in 1:2) {
  subset_data <- d.rates[d.rates$partition_module == levels(as.factor(d.rates$partition_module))[i] & complete.cases(d.rates$partition_module),]
  
  # Fit a linear model
  lm_fit <- lm(rates_la ~ lmdistP, data = subset_data)
  
  # Add the regression line to the plot
  abline(lm_fit, col = c("orange2","cornflowerblue")[i], lwd = 2)
}
legend('topleft', legend = c('Neural crest', 'Mesoderm','Unasigned'), pch = 21, pt.bg =  c( 'blue', 'orange1','white'), bty = "n", pt.cex = 1, cex = 0.8)



###Cook's D 

cooks.distance(lmDistVar)

plot(cooks.distance(lmDistVar),type="b",pch=18,col="red")

N = length(complete.cases(d.rates$partition_module))
k = 1
cutoff = 4/ (N-k-1)
abline(h=cutoff,lty=2)

cooks.distance(lmDistR)
plot(lmDistR, which = 4)
plot(cooks.distance(lmDistR),type="b",pch=18,col="red")

N = length(complete.cases(d.rates$partition_module))
k = 1
cutoff = 4/ (N-k-1)
abline(h=cutoff,lty=2)

### Without outlier
lmDistVar2 <- lm(var_la ~ partition_module*lmdistP, data = d.rates[complete.cases(d.rates$partition_module) & 
                                                                     d.rates$partition_module != 	"par-post-lat",])
summary(lmDistVar2)

lmDistR2 <- lm(rates_la ~ partition_module*lmdistP, data = d.rates[complete.cases(d.rates$partition_module) & 
                                                                     d.rates$partition_module != 	"postorbi-post",])
summary(lmDistR2)
