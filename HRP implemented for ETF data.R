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


