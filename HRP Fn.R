##### HRP Function ######


#library(seriation)
#library(RColorBrewer)


HRP_Fn <- function(corr,cov){
######## STAGE 1: Tree Clustering #############

# Functions:
# 1.2 correlDist
# 1.2 Euc_dist
# 1.3 cluster_fn

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

## Initial D matrix
dist1Lop <- correlDist(corr)

# Euclidean distance matrix
dist2Lop <- Euc_dist(dist1Lop)

# Clusters
linkLop <- cluster_fn(dist = dist2Lop)
sortIx <- linkLop$order

# Cluster order
sortIx <- getQuasiDiag(linkLop)

# Diag Matrix
# coul <- colorRampPalette(brewer.pal(8, "Spectral"))(10000)

# Recursive Bisec
outTest <- getRecBipart(cov, sortIx)

return(outTest)

}

 cor1 <- as.matrix(read.csv("corMat.csv"))
 cov1 <- as.matrix(read.csv("covMat.csv"))

ClusterTest <- HRP_Fn(corr = cor1, cov = cov1)
