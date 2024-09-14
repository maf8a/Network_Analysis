# NetworkAnalysis_ConsensusMatrices
***

**Purpose:** calculate nestedness and modularity of consensus matrices 

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

read in and visualize full data matrix 

```{r loaddata}
quantScores <- read.csv2( file ="quantScores_forNetworkAnalysis.csv", header = TRUE, row.names = 1)

#turn it into a matrix format
quantScores_matrix <- data.matrix(quantScores)#turn it into a matrix format

visweb(quantScores_matrix)#visualize the matrix
```

Make the matrix binary
```{r binary}
quantScores_bin<- as.data.frame(ifelse(quantScores > 0.5, 1, 0))
```
***

## MATRIX 1: PARASITE/SITE
***


Prepare the matrix

```{r prepParSite}
#make separate DF for each site
quantScores_bin_F <- dplyr::select(quantScores_bin, contains("_F"))
quantScores_bin_D <- dplyr::select(quantScores_bin, contains("_D"))
quantScores_bin_E <- dplyr::select(quantScores_bin, contains("_E"))
quantScores_bin_L4 <- dplyr::select(quantScores_bin, contains("_L4"))
quantScores_bin_L5 <- dplyr::select(quantScores_bin, contains("_L5"))


#copy df, then add column which is sum of each row
quantScores_bin_F_sum <- quantScores_bin_F
quantScores_bin_D_sum <- quantScores_bin_D
quantScores_bin_E_sum <- quantScores_bin_E
quantScores_bin_L4_sum <- quantScores_bin_L4
quantScores_bin_L5_sum <- quantScores_bin_L5


quantScores_bin_F_sum$Fsum <- rowSums(quantScores_bin_F_sum)
quantScores_bin_D_sum$Dsum <- rowSums(quantScores_bin_D_sum)
quantScores_bin_E_sum$Esum <- rowSums(quantScores_bin_E_sum)
quantScores_bin_L4_sum$L4sum <- rowSums(quantScores_bin_L4_sum)
quantScores_bin_L5_sum$L5sum <- rowSums(quantScores_bin_L5_sum)


#add binary column indicating whether each row has any non-zero cells
quantScores_bin_F_sum$Fbin <- ifelse(quantScores_bin_F_sum$Fsum > 0, 1, 0)
quantScores_bin_D_sum$Dbin <- ifelse(quantScores_bin_D_sum$Dsum > 0, 1, 0)
quantScores_bin_E_sum$Ebin <- ifelse(quantScores_bin_E_sum$Esum > 0, 1, 0)
quantScores_bin_L4_sum$L4bin <- ifelse(quantScores_bin_L4_sum$L4sum > 0, 1, 0)
quantScores_bin_L5_sum$L5bin <- ifelse(quantScores_bin_L5_sum$L5sum > 0, 1, 0)


#turn df into matrix so cbind is easier
quantScores_bin_F_sum_matrix <- data.matrix(quantScores_bin_F_sum)
quantScores_bin_D_sum_matrix <- data.matrix(quantScores_bin_D_sum)
quantScores_bin_E_sum_matrix <- data.matrix(quantScores_bin_E_sum)
quantScores_bin_L4_sum_matrix <- data.matrix(quantScores_bin_L4_sum)
quantScores_bin_L5_sum_matrix <- data.matrix(quantScores_bin_L5_sum)


#make a dataframe which combines all the bin columns
quantScores_siteBin <- as.data.frame(cbind("F"= quantScores_bin_F_sum_matrix[,14],
                                           "D"= quantScores_bin_D_sum_matrix[,14],
                                           "E"= quantScores_bin_E_sum_matrix[,14],
                                           "L4"= quantScores_bin_L4_sum_matrix[,14],
                                           "L5"= quantScores_bin_L5_sum_matrix[,14]))


#visualize matrix (by default, maximally nested configuration)
visweb(quantScores_siteBin, labsize = 2.8)
```



Null model

```{r nullsParSite}
nmsiteBin <- vegan::nullmodel(quantScores_siteBin, "quasiswap")#set up the null model approach
nullssiteBin <- simulate(nmsiteBin, nsim = 1000)#do the actual simulations of null models
```

Nestedness

```{r nestednessParSite}
par(mar = c(5.1, 4.1, 4.1, 2.1))

nodfsiteBin <- nest.smdm(quantScores_siteBin, weighted = FALSE, decreasing = "fill")#use decreasing = "fill" to invoke nodf. Default is to sort before calculating
nodfsiteBin_obs <- nodfsiteBin$NODFmatrix

nodfsiteBin_nulls <- apply(nullssiteBin, 3, function(x)nest.smdm(x, weighted = FALSE, decreasing = "fill")$NODFmatrix)#takes a long time
z_nodfsiteBin <- (nodfsiteBin_obs - mean(nodfsiteBin_nulls))/sd(nodfsiteBin_nulls)
p_nodfsiteBin <- 2*pnorm(-abs(z_nodfsiteBin))
plot(density(nodfsiteBin_nulls), xlim = c(0, 60), lwd = 2, main = "Nestedness siteBin", xlab = paste("nodf =", round(nodfsiteBin_obs, digits = 3), "    z =", round(z_nodfsiteBin, digits = 3), "    p =", round(p_nodfsiteBin, digits = 24)) , ylab = "Density")
abline(v = nodfsiteBin_obs, col = "red", lwd = 2)
rel_nodfsiteBin <- (nodfsiteBin_obs - mean(nodfsiteBin_nulls))/mean(nodfsiteBin_nulls)
```

Modularity

```{r modParSite}
#plot modular network
quantScores_siteBin_Mod <- computeModules(quantScores_siteBin)
plotModuleWeb(quantScores_siteBin_Mod, labsize = 0.8)

Q_siteBin_obs <- DIRT_LPA_wb_plus(quantScores_siteBin)$modularity

Q_siteBin_nulls <- apply(nullssiteBin, 3, function(x) DIRT_LPA_wb_plus(x)$modularity)#takes a long time
z_Q_siteBin <- (Q_siteBin_obs - mean(Q_siteBin_nulls))/sd(Q_siteBin_nulls)
p_Q_siteBin <- 2*pnorm(-abs(z_Q_siteBin))
plot(density(Q_siteBin_nulls), xlim = c(0, 0.6), lwd = 2, main = "Modularity siteBin", xlab = paste("Q =", round(Q_siteBin_obs, digits = 3), "    z =", round(z_Q_siteBin, digits = 3), "    p =", round(p_Q_siteBin, digits = 25)) , ylab = "Density")
abline(v = Q_siteBin_obs, col = "red", lwd = 2)#plot observed modularity
rel_Q_siteBin <- (Q_siteBin_obs - mean(Q_siteBin_nulls))/mean(Q_siteBin_nulls)
```
***

## MATRIX 2: PARASITE/HOST
***

Prepare the matrix

```{r parHost}
#make separate DF for each host
quantScores_bin_H1 <- dplyr::select(quantScores_bin, contains("H1_"))
quantScores_bin_H2 <- dplyr::select(quantScores_bin, contains("H2_"))
quantScores_bin_H3 <- dplyr::select(quantScores_bin, contains("H3_"))
quantScores_bin_H4 <- dplyr::select(quantScores_bin, contains("H4_"))
quantScores_bin_H5 <- dplyr::select(quantScores_bin, contains("H5_"))
quantScores_bin_H6 <- dplyr::select(quantScores_bin, contains("H6_"))
quantScores_bin_H7 <- dplyr::select(quantScores_bin, contains("H7_"))
quantScores_bin_H8 <- dplyr::select(quantScores_bin, contains("H8_"))
quantScores_bin_H9 <- dplyr::select(quantScores_bin, contains("H9_"))
quantScores_bin_H10 <- dplyr::select(quantScores_bin, contains("H10_"))
quantScores_bin_H11 <- dplyr::select(quantScores_bin, contains("H11_"))
quantScores_bin_H12 <- dplyr::select(quantScores_bin, contains("H12_"))


#copy df, then add column which is sum of each row
quantScores_bin_H1_sum <- quantScores_bin_H1
quantScores_bin_H2_sum <- quantScores_bin_H2
quantScores_bin_H3_sum <- quantScores_bin_H3
quantScores_bin_H4_sum <- quantScores_bin_H4
quantScores_bin_H5_sum <- quantScores_bin_H5
quantScores_bin_H6_sum <- quantScores_bin_H6
quantScores_bin_H7_sum <- quantScores_bin_H7
quantScores_bin_H8_sum <- quantScores_bin_H8
quantScores_bin_H9_sum <- quantScores_bin_H9
quantScores_bin_H10_sum <- quantScores_bin_H10
quantScores_bin_H11_sum <- quantScores_bin_H11
quantScores_bin_H12_sum <- quantScores_bin_H12

quantScores_bin_H1_sum$H1sum <- rowSums(quantScores_bin_H1_sum)
quantScores_bin_H2_sum$H2sum <- rowSums(quantScores_bin_H2_sum)
quantScores_bin_H3_sum$H3sum <- rowSums(quantScores_bin_H3_sum)
quantScores_bin_H4_sum$H4sum <- rowSums(quantScores_bin_H4_sum)
quantScores_bin_H5_sum$H5sum <- rowSums(quantScores_bin_H5_sum)
quantScores_bin_H6_sum$H6sum <- rowSums(quantScores_bin_H6_sum)
quantScores_bin_H7_sum$H7sum <- rowSums(quantScores_bin_H7_sum)
quantScores_bin_H8_sum$H8sum <- rowSums(quantScores_bin_H8_sum)
quantScores_bin_H9_sum$H9sum <- rowSums(quantScores_bin_H9_sum)
quantScores_bin_H10_sum$H10sum <- rowSums(quantScores_bin_H10_sum)
quantScores_bin_H11_sum$H11sum <- rowSums(quantScores_bin_H11_sum)
quantScores_bin_H12_sum$H12sum <- rowSums(quantScores_bin_H12_sum)

#add binary column indicating whether each row has any non-zero cells
quantScores_bin_H1_sum$H1bin <- ifelse(quantScores_bin_H1_sum$H1sum > 0, 1, 0)
quantScores_bin_H2_sum$H2bin <- ifelse(quantScores_bin_H2_sum$H2sum > 0, 1, 0)
quantScores_bin_H3_sum$H3bin <- ifelse(quantScores_bin_H3_sum$H3sum > 0, 1, 0)
quantScores_bin_H4_sum$H4bin <- ifelse(quantScores_bin_H4_sum$H4sum > 0, 1, 0)
quantScores_bin_H5_sum$H5bin <- ifelse(quantScores_bin_H5_sum$H5sum > 0, 1, 0)
quantScores_bin_H6_sum$H6bin <- ifelse(quantScores_bin_H6_sum$H6sum > 0, 1, 0)
quantScores_bin_H7_sum$H7bin <- ifelse(quantScores_bin_H7_sum$H7sum > 0, 1, 0)
quantScores_bin_H8_sum$H8bin <- ifelse(quantScores_bin_H8_sum$H8sum > 0, 1, 0)
quantScores_bin_H9_sum$H9bin <- ifelse(quantScores_bin_H9_sum$H9sum > 0, 1, 0)
quantScores_bin_H10_sum$H10bin <- ifelse(quantScores_bin_H10_sum$H10sum > 0, 1, 0)
quantScores_bin_H11_sum$H11bin <- ifelse(quantScores_bin_H11_sum$H11sum > 0, 1, 0)
quantScores_bin_H12_sum$H12bin <- ifelse(quantScores_bin_H12_sum$H12sum > 0, 1, 0)

#turn df into matrix so cbind is easier
quantScores_bin_H1_sum_matrix <- data.matrix(quantScores_bin_H1_sum)
quantScores_bin_H2_sum_matrix <- data.matrix(quantScores_bin_H2_sum)
quantScores_bin_H3_sum_matrix <- data.matrix(quantScores_bin_H3_sum)
quantScores_bin_H4_sum_matrix <- data.matrix(quantScores_bin_H4_sum)
quantScores_bin_H5_sum_matrix <- data.matrix(quantScores_bin_H5_sum)
quantScores_bin_H6_sum_matrix <- data.matrix(quantScores_bin_H6_sum)
quantScores_bin_H7_sum_matrix <- data.matrix(quantScores_bin_H7_sum)
quantScores_bin_H8_sum_matrix <- data.matrix(quantScores_bin_H8_sum)
quantScores_bin_H9_sum_matrix <- data.matrix(quantScores_bin_H9_sum)
quantScores_bin_H10_sum_matrix <- data.matrix(quantScores_bin_H10_sum)
quantScores_bin_H11_sum_matrix <- data.matrix(quantScores_bin_H11_sum)
quantScores_bin_H12_sum_matrix <- data.matrix(quantScores_bin_H12_sum)


#make a dataframe which combines all the bin columns
quantScores_hostBin <- data.frame(cbind("H1"= quantScores_bin_H1_sum_matrix[,7], 
                                        "H2"= quantScores_bin_H2_sum_matrix[,7],
                                        "H3"= quantScores_bin_H3_sum_matrix[,7],
                                        "H4"= quantScores_bin_H4_sum_matrix[,7],
                                        "H5"= quantScores_bin_H5_sum_matrix[,7],
                                        "H6"= quantScores_bin_H6_sum_matrix[,7],
                                        "H7"= quantScores_bin_H7_sum_matrix[,7],
                                        "H8"= quantScores_bin_H8_sum_matrix[,7],
                                        "H9"= quantScores_bin_H9_sum_matrix[,7],
                                        "H10"= quantScores_bin_H10_sum_matrix[,7],
                                        "H11"= quantScores_bin_H11_sum_matrix[,7],
                                        "H12"= quantScores_bin_H12_sum_matrix[,7]))

#visualize matrix (by default, maximally nested configuration)
visweb(quantScores_hostBin, labsize = 2)
```



Null model

```{r nullsParHost}

nmhostBin <- vegan::nullmodel(quantScores_hostBin, "quasiswap")#set up the null model approach
nullshostBin <- simulate(nmhostBin, nsim = 1000)#do the actual simulations of null models
```

Nestedness

```{r nestParHost}
par(mar = c(5.1, 4.1, 4.1, 2.1))

#nestedness:
nodfhostBin <- nest.smdm(quantScores_hostBin, weighted = FALSE, decreasing = "fill")#use decreasing = "fill" to invoke nodf. Default is to sort before calculating
nodfhostBin_obs <- nodfhostBin$NODFmatrix

nodfhostBin_nulls <- apply(nullshostBin, 3, function(x)nest.smdm(x, weighted = FALSE, decreasing = "fill")$NODFmatrix)#takes a long time
z_nodfhostBin <- (nodfhostBin_obs - mean(nodfhostBin_nulls))/sd(nodfhostBin_nulls)
p_nodfhostBin <- 2*pnorm(-abs(z_nodfhostBin))
plot(density(nodfhostBin_nulls), xlim = c(0, 80), lwd = 2, main = "Nestedness hostBin", xlab = paste("nodf =", round(nodfhostBin_obs, digits = 3), "    z =", round(z_nodfhostBin, digits = 3), "    p =", round(p_nodfhostBin, digits = 4)) , ylab = "Density")
abline(v = nodfhostBin_obs, col = "red", lwd = 2)
rel_nodfhostBin <- (nodfhostBin_obs - mean(nodfhostBin_nulls))/mean(nodfhostBin_nulls)
```

Modularity

```{r modParHost}
#plot modular network
quantScores_hostBin_Mod <- computeModules(quantScores_hostBin)
plotModuleWeb(quantScores_hostBin_Mod, labsize = 0.6)

Q_hostBin_obs <- DIRT_LPA_wb_plus(quantScores_hostBin)$modularity

Q_hostBin_nulls <- apply(nullshostBin, 3, function(x) DIRT_LPA_wb_plus(x)$modularity)#takes a long time
z_Q_hostBin <- (Q_hostBin_obs - mean(Q_hostBin_nulls))/sd(Q_hostBin_nulls)
p_Q_hostBin <- 2*pnorm(-abs(z_Q_hostBin))

par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(density(Q_hostBin_nulls), xlim = c(0, 0.6), lwd = 2, main = "Modularity hostBin", xlab = paste("Q =", round(Q_hostBin_obs, digits = 3), "    z =", round(z_Q_hostBin, digits = 3), "    p =", round(p_Q_hostBin, digits = 3)) , ylab = "Density")
abline(v = Q_hostBin_obs, col = "red", lwd = 2)#plot observed modularity
rel_Q_hostBin <- (Q_hostBin_obs - mean(Q_hostBin_nulls))/mean(Q_hostBin_nulls)
```
