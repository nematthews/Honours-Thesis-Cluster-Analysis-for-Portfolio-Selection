# Attempt HRP algorithm and checks 
# 18/06/2021
# Author: Nina Matthews
# Project: Honours Thesis: Cluster Analysis for Portfolio Construction
# Partner: Siphesihle Cele
# Ref Source: 
# Kipnis: https://quantstrattrader.com/2017/05/22/the-marcos-lopez-de-prado-hierarchical-risk-parity-algorithm/


### Note:
# We intend to work with Marcos Lopaz de Prado's 2015 work on HRP, hence
# we will stay as close as possible to his coding conventions such as variable and function names.
rm(list=ls())
library(StatMatch)
library(reshape2)
library(ggplot2)

###################################################################
################### HRP Algorithm  ###############################
###################################################################


######## STAGE 1: Tree Clustering #############

# Functions:
# 1.1 plotCorrMatrix
# 1.2 correlDist
# 1.2 Euc_dist
# 1.3 cluster_fn

# 1.1  Plot heatmap of corr matrix
plotCorrMatrix <- function(path, corr, labels = None){
        # Heatmap
        heatmap()
}

# 1.2 distance matrix based on correlation Fn
correlDist <- function(corr){
        # Distance matrix based on the returns correlation matrix
        # where distances are: 0 <= d[i,j] <= 1
        dist <- sqrt((1-corr)/2) # based on equation on pg.5\
        return(dist)
        
        #### CHECKED: works ####
}

# 1.3 Euclidean Distance Fn
Euc_dist <- function(distance){
        # where euclidean dist needs 0 <= ~d <= sqrt(N)
        Euclidean_d <- as.matrix(dist(distance,diag = TRUE, upper = TRUE))
        return(Euclidean_d)
        
        #### CHECKED: works ####
}

# 1.4 Cluster cols: (i*,j*) = argmin{di, dj}
cluster_fn <- function(dist){
        # Takes Euclidean distance to create cluster structure
        object_name <- as.dist(dist)
        clus <- hclust(object_name, method = "single")
        
        #### CHECKED: works ####
        
} 

######## STAGE 2: Quasi-Diagonalization #############
# Functions:
# 2.1 getQuasiDiag

# 2.1 Ordering obtained from resulting cluster hierarchy
getQuasiDiag <- function(link){
        # Sort cluster items by distance
        clust_order <- link$order
        return(clust_order)
}


######## STAGE 3: Recursive Bisection #############
# Functions:
# 3.1 getIVP (3)
# 3.2 getClusterVar  (2)
# 3.3 getRecBipart (1)
# Extra to put into getRecBipart

# Process Notes
# 

### Inverse Variance Portfolio (IVP)

getIVP <- function(cov) {
        # inv-var weights
        invDiag <- 1/diag(as.matrix(cov))
        w <- invDiag/sum(invDiag)
        return(w)
}

getClusterVar <- function(cov, cItems) {
        # compute cluster var per clsuter
        covSlice <- cov[cItems, cItems]
        # get weights from IVP for calculation: wVw
        w <- getIVP(covSlice)
        cVar <- t(w) %*% as.matrix(covSlice) %*% w
        return(cVar)
}

getRecBipart <- function(cov, sortIx) {
        assign("w", value = rep(1, ncol(cov)), envir = .GlobalEnv)
        
        # recursive call on recCall, only returns when all len(clusters) = 1
        recCall(cov, sortIx)
        return(w)
}

recCall <- function(cov, sortIx) {
        ###### Cluster Split ############################
        # get index values for the slip, truncate to insure ints
        subIdx <- 1:trunc(length(sortIx)/2)
        # use index vals to split ordered list
        cItems0 <- sortIx[subIdx]
        cItems1 <- sortIx[-subIdx]
        
        ######## ClusterVar call ##########################
        # Call ClusterVar for subcluster varaince on each half
        cVar0 <- getClusterVar(cov, cItems0)
        cVar1 <- getClusterVar(cov, cItems1)
        
        ######## Weights update ############################
        #calculate weighting factor 
        alpha <- 1 - cVar0/(cVar0 + cVar1)
        # updating of weights 
        w[cItems0] <<- w[cItems0] * alpha # weight 1
        w[cItems1] <<- w[cItems1] * (1-alpha) # weight 2
        
        ######## recursive calls #########################
        # recursive calls: conditioned on cluster length > 1
        # First half
        if(length(cItems0) > 1) {
                recCall(cov, cItems0)
        }
        # Second half
        if(length(cItems1) > 1) {
                recCall(cov, cItems1)
        }
}

###################################################################
##### CHECK 1: ####### 3 x 3 cor matrix used by MLdP ############
###################################################################

cor_d <- matrix(nrow  =  3, ncol =3, c(1,0.7,0.2,0.7,1,-0.2,0.2,-0.2,1), 
                byrow = TRUE)


## 1.2 Initial D matrix
dist1 <- correlDist(cor_d)

# 1.3 Euclidean distance matrix
dist2 <- Euc_dist(dist1)

# Clusters
link <- cluster_fn(dist = dist2)
order <- link$order

plot(as.dendrogram(link))

# Cluster order
sortIx <- getQuasiDiag(link)

# Cor heat map
#meltcor <- melt(cor_d)

###################################################################
##### CHECK 2: Full data from Lopez ############
###################################################################

Lopezdata <- read.csv("x_output.csv")

corLop <- round(cor(Lopezdata, method = "pearson"),9)

#### HEAT PLOT Check:


corLop_melted <- melt(corLop)


ggplot(data = corLop_melted, aes(x=Var1, y=Var2, fill=value)) + 
        geom_tile()+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())

###################################################################
##### CHECK 3: Full data from Lopez (retrieved by Kipnis) ############
###################################################################

cor <- as.matrix(read.csv("corMat.csv"))
rownames(cor) <- c("V1","V2", "V3", "V4","V5","V6","V7","V8","V9","V10" )
colnames(cor) <-  c("V1","V2", "V3", "V4","V5","V6","V7","V8","V9","V10" )

# HEAT PLOT Correlation

cor_melted <- melt(cor)


ggplot(data = cor_melted, aes(x=Var1, y=Var2, fill=value)) + 
        geom_tile()



cov <- as.matrix(read.csv("covMat.csv"))
rownames(cov) <- c("V1","V2", "V3", "V4","V5","V6","V7","V8","V9","V10" )
colnames(cov) <-  c("V1","V2", "V3", "V4","V5","V6","V7","V8","V9","V10" )

# HEAT PLOT Covariance
cov_melted <- melt(cov)


ggplot(data = cov_melted, aes(x=Var1, y=Var2, fill=value)) + 
        geom_tile()+
        theme(axis.title.x = element_blank(),
                axis.title.y = element_blank())

### Checking clustering with Lopez data: #####

## Initial D matrix
dist1Lop <- correlDist(cor)

# Euclidean distance matrix
dist2Lop <- Euc_dist(dist1Lop)

# Clusters
linkLop <- cluster_fn(dist = dist2Lop)
sortIx <- linkLop$order

plot(as.dendrogram(linkLop))

# Cluster order
sortIx <- getQuasiDiag(linkLop)

#### Step 2: Quasi-Diagonalization ########

# Reorganize the covar matrix based on cluster order:
#  9  2 10  1  7  3  6  4  5  8

# Diag Matrix
library(seriation)

library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "Spectral"))(10000)
diag <- hmap(cor, method = "HC_single", col = coul) 

######## STAGE 3: Recursive Bisection #############

outTest <- getRecBipart(cov, sortIx)

###################################################################
#               ########### END  ############
###################################################################


