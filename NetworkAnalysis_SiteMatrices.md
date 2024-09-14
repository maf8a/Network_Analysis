# NetworkAnalysis_SiteMatrices
***

**Purpose:** calculate nestedness and modularity of attachment site matrices

**Use package bipartite to:**
1. visualize the networks \
2. compute network indices to summarize structure \
3. statisically test for differences between observed and random networks \

**use this vignette:** https://cran.r-project.org/web/packages/bipartite/vignettes/Intro2bipartite.pdf

**analysis performed:** October 2020


***


### Load required packages

```{r load packages}
library(igraph)
library(bipartite)
```
***

### Load the data
***

read in full data matrix

```{r loaddata}
quantScores <- read.csv2( file ="quantScores_forNetworkAnalysis.csv", header = TRUE, row.names = 1)

quantScores_matrix <- data.matrix(quantScores_reduced)#turn it into a matrix format

visweb(quantScores_matrix)#visualize the matrix
```


### Make new matrix for each attachment site
***

```{r siteMatrices}
quantScores_F_reduced <- dplyr::select(quantScores_reduced, contains("_F"))
quantMatrix_F_reduced <- data.matrix(quantScores_F_reduced)

quantScores_D_reduced <- dplyr::select(quantScores_reduced, contains("_D"))
quantMatrix_D_reduced <- data.matrix(quantScores_D_reduced)

quantScores_L4_reduced <- dplyr::select(quantScores_reduced, contains("_L4"))
quantMatrix_L4_reduced <- data.matrix(quantScores_L4_reduced)

quantScores_L5_reduced <- dplyr::select(quantScores_reduced, contains("_L5"))
quantMatrix_L5_reduced <- data.matrix(quantScores_L5_reduced)

quantScores_E_reduced <- dplyr::select(quantScores_reduced, contains("_E"))
quantMatrix_E_reduced <- data.matrix(quantScores_E_reduced)
```

change the colnames of the matrices to exclude the site names
```{r colnames}
colnames(quantMatrix_F_reduced) <- c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12")
colnames(quantMatrix_D_reduced) <- c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12")
colnames(quantMatrix_L4_reduced) <- c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12")
colnames(quantMatrix_L5_reduced) <- c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12")
colnames(quantMatrix_E_reduced) <- c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12")
```

visualize the matrices (by default sorted into maximally nested configuration)
```{r visualizeMatrices}
visweb(quantMatrix_F_reduced, labsize = 2.8)
visweb(quantMatrix_D_reduced, labsize = 2.5)
visweb(quantMatrix_L4_reduced, labsize = 1.8)
visweb(quantMatrix_L5_reduced, labsize = 1)
visweb(quantMatrix_E_reduced, labsize = 2)
```


### NULL MODEL
***


specify null model to compare against (I want to use swsh_samp from vegan package). The output WNODAmatrix from nest.smdm() is the same as the NODF statistic from nestednodf()

```{r nullmodel}
nmF <- vegan::nullmodel(quantMatrix_F_reduced, "swsh_samp")#set up the null model approach
nullsF <- simulate(nmF, nsim = 1000)#do the actual simulations of null models
str(nullsF)#array of 1000 networks "behind" each other

#look at a few of the simulated matrices
nullsF[, ,3]
nullsF[, ,11]


#repeat for other sites
nmD <- vegan::nullmodel(quantMatrix_D_reduced, "swsh_samp")#set up the null model approach
nullsD <- simulate(nmD, nsim = 1000)#do the actual simulations of null models


nmL4 <- vegan::nullmodel(quantMatrix_L4_reduced, "swsh_samp")#set up the null model approach
nullsL4 <- simulate(nmL4, nsim = 1000)#do the actual simulations of null models


nmL5 <- vegan::nullmodel(quantMatrix_L5_reduced, "swsh_samp")#set up the null model approach
nullsL5 <- simulate(nmL5, nsim = 1000)#do the actual simulations of null models


nmE <- vegan::nullmodel(quantMatrix_E_reduced, "swsh_samp")#set up the null model approach
nullsE <- simulate(nmE, nsim = 1000)#do the actual simulations of null models
```


### MODULARITY
***

Compute modularity and plot result for each attachment site. This uses function DIRTLPAwb+ (Beckett 2016) as default

```{r modularity}
QF_obs <- DIRT_LPA_wb_plus(quantMatrix_F_reduced)$modularity#computes only modularity (for later comparison with null)
modF <- computeModules(quantMatrix_F_reduced)#creates a richer moduleWeb class object for plotting
plotModuleWeb(modF, labsize = 0.8)

QD_obs <- DIRT_LPA_wb_plus(quantMatrix_D_reduced)$modularity
modD <- computeModules(quantMatrix_D_reduced)
plotModuleWeb(modD)

QL4_obs <- DIRT_LPA_wb_plus(quantMatrix_L4_reduced)$modularity
modL4 <- computeModules(quantMatrix_L4_reduced)
plotModuleWeb(modL4)

QL5_obs <- DIRT_LPA_wb_plus(quantMatrix_L5_reduced)$modularity
modL5 <- computeModules(quantMatrix_L5_reduced)
plotModuleWeb(modL5)

QE_obs <- DIRT_LPA_wb_plus(quantMatrix_E_reduced)$modularity
modE <- computeModules(quantMatrix_E_reduced)
plotModuleWeb(modE)
```


Statistically evaluate modularity values by comparing matrices to random (null) matrices. Output specifically the modularity index for each matrix, to compare against observed data (and plot). The modularity index (Q) can be found in "likelihood" attribute in output of computeModules().

```{r modularityStats}
par(mar = c(5.1, 4.1, 4.1, 2.1))

#DIRT_LPA_wb_plus() and accessing $modularity. 
QF_nulls <- apply(nullsF, 3, function(x) DIRT_LPA_wb_plus(x)$modularity)#THIS TAKES A LONG TIME

#compute z scores
z_QF <- (QF_obs - mean(QF_nulls))/sd(QF_nulls)

#calculate p value
p_QF <- 2*pnorm(-abs(z_QF))

#plot observed and null-modeled nestedness
plot(density(QF_nulls), xlim = c(0, 0.6), lwd = 2, main = "Modularity F", xlab = paste("Q =", round(QF_obs, digits = 3), "    z =", round(z_QF, digits = 3), "    p =", round(p_QF, digits = 3)) , ylab = "Density")
abline(v = QF_obs, col = "red", lwd = 2)#plot observed modularity
rel_QF <- (QF_obs - mean(QF_nulls))/mean(QF_nulls)#calculate relative modularity (Fortuna et al 2008)...record as M*



#repeat for other sites
QD_nulls <- apply(nullsD, 3, function(x) DIRT_LPA_wb_plus(x)$modularity)#THIS TAKES A LONG TIME (~25 MIN)
z_QD <- (QD_obs - mean(QD_nulls))/sd(QD_nulls)
p_QD <- 2*pnorm(-abs(z_QD))
plot(density(QD_nulls), xlim = c(0, 0.6), lwd = 2, main = "Modularity D", xlab = paste("Q =", round(QD_obs, digits = 3), "    z =", round(z_QD, digits = 3), "    p =", round(p_QD, digits = 5)) , ylab = "Density")
abline(v = QD_obs, col = "red", lwd = 2)
rel_QD <- (QD_obs - mean(QD_nulls))/mean(QD_nulls)


QL4_nulls <- apply(nullsL4, 3, function(x) DIRT_LPA_wb_plus(x)$modularity)#THIS TAKES A LONG TIME 
z_QL4 <- (QL4_obs - mean(QL4_nulls))/sd(QL4_nulls)
p_QL4 <- 2*pnorm(-abs(z_QL4))
plot(density(QL4_nulls), xlim = c(0, 0.6), lwd = 2, main = "Modularity L4", xlab = paste("Q =", round(QL4_obs, digits = 3), "    z =", round(z_QL4, digits = 3), "    p =", round(p_QL4, digits = 10)) , ylab = "Density")
abline(v = QL4_obs, col = "red", lwd = 2)
rel_QL4 <- (QL4_obs - mean(QL4_nulls))/mean(QL4_nulls)


QL5_nulls <- apply(nullsL5, 3, function(x) DIRT_LPA_wb_plus(x)$modularity)#THIS TAKES A LONG TIME 
z_QL5 <- (QL5_obs - mean(QL5_nulls))/sd(QL5_nulls)
p_QL5 <- 2*pnorm(-abs(z_QL5))
plot(density(QL5_nulls), xlim = c(0, 0.6), lwd = 2, main = "Modularity L5", xlab = paste("Q =", round(QL5_obs, digits = 3), "    z =", round(z_QL5, digits = 3), "    p =", round(p_QL5, digits = 13)) , ylab = "Density")
abline(v = QL5_obs, col = "red", lwd = 2)
rel_QL5 <- (QL5_obs - mean(QL5_nulls))/mean(QL5_nulls)


QE_nulls <- apply(nullsE, 3, function(x) DIRT_LPA_wb_plus(x)$modularity)#THIS TAKES A LONG TIME 
z_QE <- (QE_obs - mean(QE_nulls))/sd(QE_nulls)
p_QE <- 2*pnorm(-abs(z_QE))
plot(density(QE_nulls), xlim = c(0, 0.6), lwd = 2, main = "Modularity E", xlab = paste("Q =", round(QE_obs, digits = 3), "    z =", round(z_QE, digits = 3), "    p =", round(p_QE, digits = 15)) , ylab = "Density")
abline(v = QE_obs, col = "red", lwd = 2)
rel_QE <- (QE_obs - mean(QE_nulls))/mean(QE_nulls)
```



### NESTEDNESS
***

Calculate nestedness of each matrix using WNODA index (Pinheiro et al 2019)

```{r nestedness}
wnodaF <- nest.smdm(quantMatrix_F_reduced, weighted = TRUE, decreasing = "abund")#use decreasing = "abund" to invoke WNODA. Default is to sort before calculating
wnodaF_obs <- wnodaF$WNODAmatrix#nestedness for entire matrix (WNODA statistic)


#repeat for other matrices
wnodaD <- nest.smdm(quantMatrix_D_reduced, weighted = TRUE, decreasing = "abund")#use decreasing = "abund" to invoke WNODA. Default is to sort before calculating
wnodaD_obs <- wnodaD$WNODAmatrix


wnodaL4 <- nest.smdm(quantMatrix_L4_reduced, weighted = TRUE, decreasing = "abund")#use decreasing = "abund" to invoke WNODA. Default is to sort before calculating
wnodaL4_obs <- wnodaL4$WNODAmatrix


wnodaL5 <- nest.smdm(quantMatrix_L5_reduced, weighted = TRUE, decreasing = "abund")#use decreasing = "abund" to invoke WNODA. Default is to sort before calculating
wnodaL5_obs <- wnodaL5$WNODAmatrix


wnodaE <- nest.smdm(quantMatrix_E_reduced, weighted = TRUE, decreasing = "abund")#use decreasing = "abund" to invoke WNODA. Default is to sort before calculating
wnodaE_obs <- wnodaE$WNODAmatrix
```



Statistically evaluate nestedness by comparing to null. Output specifically the WNODA index for each matrix, to compare against observed data.

```{r nestednessStats}
wnodaF_nulls <- apply(nullsF, 3, function(x)nest.smdm(x, weighted = TRUE, decreasing = "abund")$WNODAmatrix)#to get all the WNODA statistics, just leave out the "$WNODAmatrix" from this line

#compute z scores
z_wnodaF <- (wnodaF_obs - mean(wnodaF_nulls))/sd(wnodaF_nulls)

#calculate p value
p_wnodaF <- 2*pnorm(-abs(z_wnodaF))

#plot observed and null-modelled nestedness
plot(density(wnodaF_nulls), xlim = c(0, 60), lwd = 2, main = "Nestedness F", xlab = paste("WNODA =", round(wnodaF_obs, digits = 3), "    z =", round(z_wnodaF, digits = 3), "    p =", round(p_wnodaF, digits = 5)) , ylab = "Density")
abline(v = wnodaF_obs, col = "red", lwd = 2)#plot observed nestedness
rel_wnodaF <- (wnodaF_obs - mean(wnodaF_nulls))/mean(wnodaF_nulls)


#repeat with other sites
wnodaD_nulls <- apply(nullsD, 3, function(x)nest.smdm(x, weighted = TRUE, decreasing = "abund")$WNODAmatrix)
z_wnodaD <- (wnodaD_obs - mean(wnodaD_nulls))/sd(wnodaD_nulls)
p_wnodaD <- 2*pnorm(-abs(z_wnodaD))
plot(density(wnodaD_nulls), xlim = c(0, 60), lwd = 2, main = "Nestedness D", xlab = paste("WNODA =", round(wnodaD_obs, digits = 3), "    z =", round(z_wnodaD, digits = 3), "    p =", round(p_wnodaD, digits = 11)) , ylab = "Density")
abline(v = wnodaD_obs, col = "red", lwd = 2)
rel_wnodaD <- (wnodaD_obs - mean(wnodaD_nulls))/mean(wnodaD_nulls)

wnodaL4_nulls <- apply(nullsL4, 3, function(x)nest.smdm(x, weighted = TRUE, decreasing = "abund")$WNODAmatrix)
z_wnodaL4 <- (wnodaL4_obs - mean(wnodaL4_nulls))/sd(wnodaL4_nulls)
p_wnodaL4 <- 2*pnorm(-abs(z_wnodaL4))
plot(density(wnodaL4_nulls), xlim = c(0, 60), lwd = 2, main = "Nestedness L4", xlab = paste("WNODA =", round(wnodaL4_obs, digits = 3), "    z =", round(z_wnodaL4, digits = 3), "    p =", round(p_wnodaL4, digits = 3)) , ylab = "Density")
abline(v = wnodaL4_obs, col = "red", lwd = 2)
rel_wnodaL4 <- (wnodaL4_obs - mean(wnodaL4_nulls))/mean(wnodaL4_nulls)

wnodaL5_nulls <- apply(nullsL5, 3, function(x)nest.smdm(x, weighted = TRUE, decreasing = "abund")$WNODAmatrix)
z_wnodaL5 <- (wnodaL5_obs - mean(wnodaL5_nulls))/sd(wnodaL5_nulls)
p_wnodaL5 <- 2*pnorm(-abs(z_wnodaL5))
plot(density(wnodaL5_nulls), xlim = c(0, 60), lwd = 2, main = "Nestedness L5", xlab = paste("WNODA =", round(wnodaL5_obs, digits = 3), "    z =", round(z_wnodaL5, digits = 3), "    p =", round(p_wnodaL5, digits = 7)) , ylab = "Density")
abline(v = wnodaL5_obs, col = "red", lwd = 2)
rel_wnodaL5 <- (wnodaL5_obs - mean(wnodaL5_nulls))/mean(wnodaL5_nulls)

wnodaE_nulls <- apply(nullsE, 3, function(x)nest.smdm(x, weighted = TRUE, decreasing = "abund")$WNODAmatrix)
z_wnodaE <- (wnodaE_obs - mean(wnodaE_nulls))/sd(wnodaE_nulls)
p_wnodaE <- 2*pnorm(-abs(z_wnodaE))
plot(density(wnodaE_nulls), xlim = c(0, 60), lwd = 2, main = "Nestedness E", xlab = paste("WNODA =", round(wnodaE_obs, digits = 3), "    z =", round(z_wnodaE, digits = 3), "    p =", round(p_wnodaE, digits = 13)) , ylab = "Density")
abline(v = wnodaE_obs, col = "red", lwd = 2)
rel_wnodaE <- (wnodaE_obs - mean(wnodaE_nulls))/mean(wnodaE_nulls)
```
 

