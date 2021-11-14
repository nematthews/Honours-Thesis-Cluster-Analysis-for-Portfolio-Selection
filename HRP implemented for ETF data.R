# Attempt HRP algorithm and checks 
# 18/06/2021
# Author: Nina Matthews
# Project: Honours Thesis: Cluster Analysis for Portfolio Construction
# Partner: Siphesihle Cele


### Note:
# We intend to work with Marcos Lopaz de Prado's 2015 work on HRP, hence
# we will stay as close as possible to his coding conventions such as variable and function names.
rm(list=ls())
library(StatMatch)
library(reshape2)
library(ggplot2)
library(seriation)
library(RColorBrewer)
library(graphics)

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
# 3.1 getRecBipart 
# 3.2 recCall
# 3.3 getClusterVar  
# 3.4 getIVP 


getRecBipart <- function(cov, sortIx) {
        assign("w", value = rep(1, ncol(cov)), envir = .GlobalEnv)
        
        # recursive call on recCall, only returns when all len(clusters) = 1
        recCall(cov, sortIx)
        return(w)
}

recCall <- function(cov, sortIx) {
        ###### Cluster Split ############################
        # get index values for the slip, truncate to insure ints
        Idx_split <- (1:trunc(length(sortIx)/2))
        # use index vals to split ordered list
        cItems0 <- sortIx[Idx_split]
        cItems1 <- sortIx[-Idx_split]
        
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

### (Sub) Cluster Variance 

getClusterVar <- function(cov, cItems) {
        # Slices diag covar matrix for assets under each subclust
        diag_slice <- cov[cItems, cItems]
        # get weights from IVP for calculation: wVw
        w <- getIVP(diag_slice)
        cVar <- t(w) %*% as.matrix(diag_slice) %*% w
        return(cVar)
}

### Inverse Variance Portfolio (IVP)
getIVP <- function(slice) {
        # inv-var weights
        inv_diag <- 1/diag(as.matrix(slice))
        # For division of trace: sum diag
        w <- inv_diag/sum(inv_diag)
        return(w)
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

plot(as.dendrogram(link),main = "Example Hierarchical Dendrogram", ylab = "Cluster Leaves Distances", xlab = "Assets")

# Cluster order
sortIx <- getQuasiDiag(link)

# Cor heat map
#meltcor <- melt(cor_d)

###################################################################
##### CHECK 2: Full data from Lopez ############
###################################################################

### Correlation from Lopez gen
Lopezdata.cor <- as.matrix(read.csv("x_output.csv"))
Lopezdata.cor <- Lopezdata.cor[,-1]
rownames(Lopezdata.cor) <- c("Asset 1","Asset 2", "Asset 3", "Asset 4","Asset 5","Asset 6","Asset 7","Asset 8","Asset 9","Asset 10" )
colnames(Lopezdata.cor) <-  c("Asset 1","Asset 2", "Asset 3", "Asset 4","Asset 5","Asset 6","Asset 7","Asset 8","Asset 9","Asset 10" )

#### HEAT PLOT Check:
corLop_melted <- melt(Lopezdata.cor)


ggplot(data = corLop_melted, aes(x=Var1, y=Var2, fill=value)) + 
        geom_tile()+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank()) +
        theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10))

###################################################################
##### CHECK 3: Full data from Lopez (retrieved by Kipnis) ############
###################################################################

# Correct
cor <- as.matrix(read.csv("corMat.csv"))
rownames(cor) <- c("Asset 1","Asset 2", "Asset 3", "Asset 4","Asset 5","Asset 6","Asset 7","Asset 8","Asset 9","Asset 10" )
colnames(cor) <-  c("Asset 1","Asset 2", "Asset 3", "Asset 4","Asset 5","Asset 6","Asset 7","Asset 8","Asset 9","Asset 10" )

# HEAT PLOT Correlation

cor_melted <- melt(cor)


ggplot(data = cor_melted, aes(x=Var1, y=Var2, fill=value)) + 
        geom_tile() + theme(axis.title.x = element_blank(),
                          axis.title.y = element_blank()) +
        theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10))



cov <- as.matrix(read.csv("covMat.csv"))
rownames(cov) <- c("Asset 1","Asset 2", "Asset 3", "Asset 4","Asset 5","Asset 6","Asset 7","Asset 8","Asset 9","Asset 10" )
colnames(cov) <-  c("Asset 1","Asset 2", "Asset 3", "Asset 4","Asset 5","Asset 6","Asset 7","Asset 8","Asset 9","Asset 10" )
# HEAT PLOT Covariance
cov_melted <- melt(cov)


ggplot(data = cov_melted, aes(x=Var1, y=Var2, fill=value)) +
        geom_tile()+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank()) +
theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10))

### Checking clustering with Lopez data: #####

## Initial D matrix
dist1Lop <- correlDist(cor)

# Euclidean distance matrix
dist2Lop <- Euc_dist(dist1Lop)

# Clusters
linkLop <- cluster_fn(dist = dist2Lop)
sortIx <- linkLop$order

plot(as.dendrogram(linkLop), main = "Hierarchical Clustering Dendrogram", ylab = "Cluster Leaves Distances", xlab = "Assets")


# Cluster order
sortIx <- getQuasiDiag(linkLop)

#### Step 2: Quasi-Diagonalization ########

# Reorganize the covar matrix based on cluster order:
#  9  2 10  1  7  3  6  4  5  8

# Diag Matrix

coul <- colorRampPalette(brewer.pal(8, "Spectral"))(100000)
# need to supply the cor matrix as HRP uses cor dist for clustering
diag <- hmap(cor, method = "HC_single", col = coul) 

######## STAGE 3: Recursive Bisection #############

outTest <- getRecBipart(cov, sortIx)


###################################################################
#               ########### EFT Data  ############
###################################################################

library(readxl)
ETF_cov <- as.matrix(read_excel("DataPostCov2.xlsx"))

ETF_cor <- cov2cor(ETF_cov)

rownames(ETF_cor) <- c("ETFSAP", "ETFSWX", "ETFT40", "ETFGLD", "ETFPLD", "ETFPLT", "ASHINF", "ASHMID", "ASHT40", "CSPROP", "SMART","NFEMOM", "NFGOVI", "NFILBI", "MAPPSG", "GIVISA", "NFSH40", "NFTRCI", "NGPLD", "STX40", "STXDIV",  "STXRAF" ,"SYGEU"  )
colnames(ETF_cor) <-  c("ETFSAP", "ETFSWX", "ETFT40", "ETFGLD", "ETFPLD", "ETFPLT", "ASHINF", "ASHMID", "ASHT40", "CSPROP", "SMART","NFEMOM", "NFGOVI", "NFILBI", "MAPPSG", "GIVISA", "NFSH40", "NFTRCI", "NGPLD", "STX40", "STXDIV",  "STXRAF" ,"SYGEU"  )


# HEAT PLOT Covariance
cor_melted <- melt(ETF_cor, id.vars = NULL)


ggplot(data = cor_melted, aes(x=Var1, y=Var2, fill=value), title = "Unclustered Correlations") + 
        geom_tile()+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())+
        theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10))


## Initial D matrix
dist1Lop <- correlDist(ETF_cor)

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

coul <- colorRampPalette(brewer.pal(8, "Spectral"))(10000)
diag <- hmap(ETF_cor, method = "HC_single", col = coul) 

######## STAGE 3: Recursive Bisection #############

outTest <- getRecBipart(ETF_cov, sortIx)

###################################################################
#               ########### END  ############
###################################################################


