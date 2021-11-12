# 13/10/2021
# Author: Nina Matthews
# Project: Honours Thesis: Cluster Analysis for Portfolio Construction
# Partner: Siphesihle Cele
# Supervisor: Tim Gebbie

############################################################################# PACKAGES
#############################################################################

rm(list=ls())

# load required libraries
library(zoo)
library(xts)
library(timeSeries) 
library(rbenchmark)
library(nloptr) # for SQP
library(quadprog) # for QP
library(ggplot2)


############################################################################# DATA 
#############################################################################

## Load previously prepared data
load(file = "PT-TAA.RData")

## 1. Checking for missing data
# check for missing data
head(is.na(tsGRet))
# define tickers of interest
Entities = colnames(tsGRet)
# remove the money market asset (we will compute excess returns!)
Entities <- Entities[-c(grep('STEFI',Entities))]
Entities <- Entities[-c(grep('ALSI',Entities))]
# reference out the risk-free asset returns
tsRRF <- tsGRet[,'STEFI']
# reference out the "market portfolio"
tsMKT <- tsGRet[,'ALSI']
# reference out the tickers of interest
tsGRet <- tsGRet[,Entities] # BOND + EQTY. INDEX. PORTFOLIO


############################################################################# IS / OOS split  
#############################################################################

#### DOC: Insample SR vs OS ####
# 60/40 split

### Training Data
tsGRet <- tsGRet[-1,]
tsRRF <- tsRRF[-1,]
tsMKT <- tsMKT[-1,]
# establish 60 split:
# dim(tsGRet)[1]*0.6 = 88.2
ISGret <- head(tsGRet,88)
OOSGRet <- tail(tsGRet,59)

ISRFR <- head(tsRRF,88)
OOSRFR <- tail(tsRRF,59)

ISMKT <- head(tsMKT,88)
OOSMKT <- tail(tsMKT,59)

#############
# Increased window size: 20:80 split

bigISGret <- head(tsGRet,29)
bigOOSGRet <- tail(tsGRet,118)

bigISRFR <- head(tsRRF,29)
bigOOSRFR <- tail(tsRRF,118)

bigISMKT <- head(tsMKT,29)
bigOOSMKT <- tail(tsMKT,118)

############################################################################# 
# Iteration and Min/Max BT Length
#############################################################################
### iv) Estimate of The Maximum of the Sample
getMax<- function(N) {
        if(N < 2)
                stop("N must be greater than 1")
        #Euler-Mascheroni constant
        euler.const<-(-digamma(1))
        a<-(1-1/N)
        b<-(1-(1/N)*( exp(-1)))
        exp.max<-(1-euler.const)*qnorm(a) + (euler.const)*qnorm(b)
        # Estimate of Maximum
        return(round(exp.max, 6))
}

### v) Get Minimum Backtest Length

getMinBackTestLen<- function(N,ex.maxSharpe=1) {
        if(is.na(ex.maxSharpe))
                stop("NA not valid")
        if(ex.maxSharpe==0)
                stop("Division by zero not valid")
        len<-(getMax(N)/ex.maxSharpe)^2
        #Return Integer
        return (round(len, 0))
}
#Upper bound of Minimum Backtest Length
getMinBTL_UB<- function(N,ex.maxSharpe=1) {
        if(is.na(ex.maxSharpe))
                stop("NA not valid")
        if(ex.maxSharpe==0)
                stop("Division by zero not valid")
        len<-(2*log(N)/ex.maxSharpe^2)
        #Return Integer
        return (round(len, 0))
}


############################################################################# IS analysis
# 60:40
#############################################################################


## 2. Compute the Geometric mean
# without correcting for missing data (NA)
mIS <- colMeans(ISGret)
rfrIS <- colMeans(ISRFR)
sIS <- colStdevs(ISGret)
# include the missing data (NA)
mIS <- colMeans(ISGret, na.rm=TRUE)
rfrIS <- colMeans(ISRFR, na.rm=TRUE)
EMktIS <- colMeans(ISMKT, na.rm=TRUE)
VMktIS <- colStdevs(ISMKT, na.rm=TRUE)
sIS <- colStdevs(ISGret, na.rm=TRUE)
cIS <- var(ISGret, na.rm=TRUE)
# omit missing data rows (clean data)
IScltsGRet = na.omit(ISGret)
# can use !is.na
mIS <- colMeans(IScltsGRet)
# visualise the covariance matrix
heatmap(cIS)
# visualise the correlation matrix
rho1 <- cov2cor(cIS)

par(mar = c(5,4,4,4))
plot(sIS, mIS,
     ylab="Expected Return [%]",
     xlab ="Volatility [%]",
     main="IS Monthly Hist. Risk & Return",
     # plot.type="s",
     ylim = c(0, 0.025), xlim = c(0, 0.12))

# turn on the grid
grid()
# label points
text(sIS, mIS,labels=names(mIS), cex= 1, pos = 4)


############################################################################# Optimal SR Max IS 60:40
#############################################################################

## Plot the efficient frontier.
# create the range of risk aversion parameters
lambda <- seq(from=0,to=1,length.out=80)
# Fully Invested
A <- matrix (1, nrow=length(mIS))
b <- 1
meq <- 1
# No short-selling
A <- cbind(1, diag(length(mIS)))
b <- c(b, rep(0, length(mIS)))
# initialise the weights
Wts <- matrix(NA,length(lambda),length(mIS))
# 8. Find the weight vector for each return level
for (i in 1:length(lambda)) {
        f <- mIS * lambda[i] # This moves the solution up along the efficient frontier
        H <- cIS # the covariance matrix
        sol <- solve.QP(H, f, A, b, meq=1)
        Wts[i,] <- sol$solution 
}

# use replicate and element-wise multiplication to find the Expected Returns
ISERet2 <- rowSums(Wts * t(replicate(nrow(Wts),mIS)))
# use matrix multiplication to find the risk
# loop over each asset pair and compute risk each weigth vector
# preallocate zero (could use NA)
ISERisk2 <- numeric(nrow(Wts))
ISERisk2[] <- NA
# pre-compute covariance matrix
IS.Sigma2 <- matrix(cIS,nrow(cIS),nrow(cIS))
# compute the portfolio volatility
for (i in 1:nrow(Wts)) {ISERisk2[i] <- as.numeric(Wts[i,]) %*% IS.Sigma2 %*% as.numeric(Wts[i,])}
# add the return and risk columns to data frame
ISdfA2 <- cbind(sqrt(ISERisk2),ISERet2)
# compute the sharpe ration
IS.SR <- (ISERet2 - rfrIS) / sqrt(ISERisk2)
# add risk return curves to plot
points(ISdfA2,pch = 19, type = "l", lwd=1, col ="blue")
# annotate
text(0.025,0.013,labels='IS Efficient Frontier',pos = 4, col ="blue")

## 9. Find The Sharpe Ratio maximising portfolio
# Initial values as fully invested equally weighted portfolio
Ones0 <- seq(1,1,length.out = length(mIS))
# equally weighted portfolio
Wts0 <- Ones0 / length(Ones0)
# unit vector
e <- rep(1,length(Wts0)) # useful matrix (ones)
# initialise the weights
Wts <- matrix(NA,1,length(mIS))
# 8. Maximise the Sharpe Ratio  
fn0 <- function(x) {return(-(x%*% mIS - rfrIS)/ sqrt(x %*% cIS %*% x) )}
# Fully Invested + Return Target
heq0 <- function(x) {return(x %*% e - 1)} # fully invested
# Use SQP to solve for the tangency portfolio  
soln <- slsqp(Wts0, fn = fn0, gr = NULL, # target returns
              lower = rep(0,length(Wts0)), # no short-selling
              upper = rep(1,length(Wts0)), # no leverage
              heq = heq0, # fully invested constraint function
              control = list(xtol_rel = 1e-8),
              nl.info = TRUE) # SQP
Wts <- soln$par
# print the weight matrix
print(Wts)
### 103 ITERATIONS
SR.wts <- as.matrix(Wts)



###############################################################################
############## IS Equally weighted portfolio 60:40  ###################### 
###############################################################################


# Initialise weights as an equally weighted portfolio
Wts0 <- as.vector(seq(1,1,length.out = length(mIS)) / length(mIS))

# Returns for any portfolio with equal weights
IS.EquiRet <-sum(Wts0*mIS)

## Risk for equal port
# remove row/col names
IS.Sigma2 <- matrix(cIS,nrow(cIS),nrow(cIS))
#IS.EquiRisk2 <- matrix(NA,length(Wts0))

IS.EquiRisk2 <- t(Wts0) %*% IS.Sigma2 %*% Wts0

###############################################################################
############## IS HRP portfolio 60:40 ###################### 
###############################################################################
source("HRP Fn.R")


# get correlation matrix
IS.VarMat <- var(IScltsGRet)
# annualized
IS.corMat <- cov2cor(IS.VarMat)


#### CLUSTERING ###

IS.HRP.wts <- HRP_Fn(corr = IS.corMat, cov = IS.VarMat)

ISm <- matrix(mIS,nrow = length(mIS), ncol = 1)

# Returns for HRP portfolio
IS.HRP.Ret <-sum(IS.HRP.wts%*%mIS)
# IS.Sigma2 <- matrix(covar,nrow(covar),nrow(covar))
IS.HRP.risk2 <- t(IS.HRP.wts) %*% IS.Sigma2 %*% IS.HRP.wts

points(sqrt(IS.HRP.risk2),IS.HRP.Ret, col = "darkgreen",cex = 2,lwd = 3)
points(0.02,0.298)

legend('bottomright', legend = c("IS Sharpe Ratio Maximizing", 'IS HRP', "Equally weighted", "IS Efficiency Frontier", "IS Sharpe Ratio"), 
       col = c('magenta', 'darkgreen',"orange","blue", "red"), lwd = c(3,3,3,1,1), lty = c(NA, NA, NA,"solid", "dashed"), cex = 0.9, pch = c(1,1,1, NA,NA))

###############################################################################

## 10. Compute return and risk for each weight vector
IS.ERetPSR <- Wts %*% mIS
# use matrix multiplication to find the risk
IS.ERiskPSR <- Wts %*% IS.Sigma2 %*% Wts 
# add risk return curves to plot
points(sqrt(IS.ERiskPSR),IS.ERetPSR,pch = 1, type = "p", col ="magenta", cex = 2, lwd = 3)
points(sqrt(IS.EquiRisk2),IS.EquiRet, col = "orange",cex = 2, lwd = 3)
# include text
#text(sqrt(IS.ERiskPSR), IS.ERetPSR,labels='IS Maximal Sharpe Ratio', col = 'magenta', cex= 1, pos = 2)
text(0.06,0.025, "SR", col = "red")
# 
# # 11. Plot SML (Security Market Line)
# IS.SR0 <- ((IS.ERetPSR-rfrIS) / sqrt(IS.ERiskPSR))
# IS.eq = function(IS.x){return(as.numeric(rfrIS + c(IS.SR0) * IS.x))}
# IS.x <- seq(from=0,to=0.30,length.out=20)
# points(IS.x, IS.eq(IS.x), type = "l", col="red")
# text(0.09,0.20,labels='Market Line', srt = 45, col = 'red', pos = 4)

# plot the Sharpe Ratio against risk levels
#par(mfrow=c(6,4),mar=c(1,1,1,1), oma=c(3,1,0,0))
par(new = T)
plot(sqrt(ISERisk2),IS.SR, axes=F, type ="l", lty = 2, col="red",xlab=NA, ylab=NA,xlim = c(0, 0.12))
axis(side = 4)

mtext(expression(Sharpe ~ Ratio: ~ ~ frac(mu-R[f],sigma )), side = 4, col ="red", line = 3)


#**************************************************************

############################################################################# IS analysis 80:20
#############################################################################

#**************************************************************


## 2. Compute the Geometric mean
# without correcting for missing data (NA)
big.mIS <- colMeans(bigISGret)
big.rfrIS <- colMeans(bigISRFR)
big.sIS <- colStdevs(bigISGret)
# include the missing data (NA)
big.mIS <- colMeans(bigISGret, na.rm=TRUE)
big.rfrIS <- colMeans(bigISRFR, na.rm=TRUE)
big.EMktIS <- colMeans(bigISMKT, na.rm=TRUE)
big.VMktIS <- colStdevs(bigISMKT, na.rm=TRUE)
big.sIS <- colStdevs(bigISGret, na.rm=TRUE)
big.cIS <- var(bigISGret, na.rm=TRUE)
# omit missing data rows (clean data)
big.IScltsGRet = na.omit(bigISGret)
# can use !is.na
big.mIS <- colMeans(big.IScltsGRet)
# visualise the covariance matrix
heatmap(big.cIS)
# visualise the correlation matrix
rho1 <- cov2cor(big.cIS)

# par(mar = c(5,4,4,4))
# plot(big.sIS,big.mIS,
#      ylab="Expected Return [%]",
#      xlab ="Volatility [%]",
#      main="IS Monthly Hist. Risk & Return",
#      # plot.type="s",
#      ylim = c(0, 0.025), xlim = c(0, 0.12))
# 
# # turn on the grid
# grid()
# # label points
# text(sIS, mIS,labels=names(mIS), cex= 1, pos = 4)

#**************************************************************

############################################################################# Optimal SR Max IS 60:40
#############################################################################

#**************************************************************

## Plot the efficient frontier.
# create the range of risk aversion parameters
biglambda <- seq(from=0,to=1,length.out=80)
# Fully Invested
big.A <- matrix (1, nrow=length(big.mIS))
big.b <- 1
big.meq <- 1
# No short-selling
big.A <- cbind(1, diag(length(big.mIS)))
big.b <- c(big.b, rep(0, length(big.mIS)))
# initialise the weights
big.Wts <- matrix(NA,length(biglambda),length(big.mIS))
# 8. Find the weight vector for each return level
for (i in 1:length(biglambda)) {
        f <- big.mIS * biglambda[i] # This moves the solution up along the efficient frontier
        H <- big.cIS # the covariance matrix
        sol <- solve.QP(H, f, big.A, big.b, meq=1)
        big.Wts[i,] <- sol$solution
}

# use replicate and element-wise multiplication to find the Expected Returns
big.ISERet2 <- rowSums(big.Wts * t(replicate(nrow(big.Wts),big.mIS)))
# use matrix multiplication to find the risk
# loop over each asset pair and compute risk each weigth vector
# preallocate zero (could use NA)
big.ISERisk2 <- numeric(nrow( big.Wts))
big.ISERisk2[] <- NA
# pre-compute covariance matrix
big.IS.Sigma2 <- matrix(big.cIS,nrow(big.cIS),nrow(big.cIS))
# compute the portfolio volatility
for (i in 1:nrow(big.Wts)) {big.ISERisk2[i] <- as.numeric(big.Wts[i,]) %*% big.IS.Sigma2 %*% as.numeric(big.Wts[i,])}
# add the return and risk columns to data frame
big.ISdfA2 <- cbind(sqrt(big.ISERisk2),big.ISERet2)
# compute the sharpe ration
big.IS.SR <- (big.ISERet2 - big.rfrIS) / sqrt(big.ISERisk2)
# # add risk return curves to plot
# points(ISdfA2,pch = 19, type = "l", lwd=1, col ="blue")
# # annotate
# text(0.025,0.013,labels='IS Efficient Frontier',pos = 4, col ="blue")

## 9. Find The Sharpe Ratio maximising portfolio
# Initial values as fully invested equally weighted portfolio
big.Ones0 <- seq(1,1,length.out = length(big.mIS))
# equally weighted portfolio
big.Wts0 <- big.Ones0 / length(big.Ones0)
# unit vector
big.e <- rep(1,length(big.Wts0)) # useful matrix (ones)
# initialise the weights
big.Wts <- matrix(NA,1,length(big.mIS))
# 8. Maximise the Sharpe Ratio  
big.fn0 <- function(x) {return(-(x%*% big.mIS - big.rfrIS)/ sqrt(x %*% big.cIS %*% x) )}
# Fully Invested + Return Target
big.heq0 <- function(x) {return(x %*% big.e - 1)} # fully invested
# Use SQP to solve for the tangency portfolio  
big.soln <- slsqp(big.Wts0, fn = big.fn0, gr = NULL, # target returns
              lower = rep(0,length(big.Wts0)), # no short-selling
              upper = rep(1,length(big.Wts0)), # no leverage
              heq = big.heq0, # fully invested constraint function
              control = list(xtol_rel = 1e-8),
              nl.info = TRUE) # SQP
big.Wts <- big.soln$par
# print the weight matrix
#### 14 iterations
print(big.Wts)

# ## 10. Compute return and risk for each weight vector
# big.IS.ERetPSR <- big.Wts %*% big.mIS
# # use matrix multiplication to find the risk
# big.IS.ERiskPSR <- big.Wts %*% big.IS.Sigma2 %*% big.Wts 
# # add risk return curves to plot
# points(sqrt(IS.ERiskPSR),IS.ERetPSR,pch = 1, type = "p", col ="magenta", cex = 2, lwd = 3)
# points(sqrt(IS.EquiRisk2),IS.EquiRet, col = "orange",cex = 2, lwd = 3)
# # include text
# #text(sqrt(IS.ERiskPSR), IS.ERetPSR,labels='IS Maximal Sharpe Ratio', col = 'magenta', cex= 1, pos = 2)
# text(0.06,0.025, "SR", col = "red")
# # 

############################################################################# OOS analysis  
#############################################################################

## 2. Compute the Geometric mean
# without correcting for missing data (NA)
mOOS <- colMeans(OOSGRet)
rfrOOS<- colMeans(OOSRFR)
sOOS <- colStdevs(OOSGRet)
# include the missing data (NA)
mOOS <- colMeans(OOSGRet, na.rm=TRUE)
rfrOOS <- colMeans(OOSRFR, na.rm=TRUE)
EMktOOS <- colMeans(OOSMKT, na.rm=TRUE)
VMktOOS <- colStdevs(OOSMKT, na.rm=TRUE)
sOOS <- colStdevs(OOSGRet, na.rm=TRUE)
cOOS <- var(OOSGRet, na.rm=TRUE)
# omit missing data rows (clean data)
OOScltsGRet = na.omit(OOSGRet)
# can use !is.na
mOOS <- colMeans(OOScltsGRet)
# visualise the covariance matrix
heatmap(cOOS)
# visualise the correlation matrix
rho1 <- cov2cor(cOOS)

plot(sOOS, mOOS,
     ylab="Expected Return",
     xlab ="Volatility",
     main="OOS Monthly Hist. Risk & Return",
     # plot.type="s",
     ylim = c(0, 0.026), xlim = c(0, 0.14))

# turn on the grid
grid()
# label points
text(sOOS, mOOS,labels=names(mOOS), cex= 1, pos = 4)


############################################################################# Optimal SR Max OOS 
#############################################################################

## Plot the efficient frontier.
# create the range of risk aversion parameters
lambda <- seq(from=0,to=1,length.out=80)
# Fully Invested
A <- matrix (1, nrow=length(mOOS))
b <- 1
meq <- 1
# No short-selling
A <- cbind(1, diag(length(mOOS)))
b <- c(b, rep(0, length(mOOS)))
# initialise the weights
OOS.Wts <- matrix(NA,length(lambda),length(mOOS))
# 8. Find the weight vector for each return level
for (i in 1:length(lambda)) {
        f <- mOOS * lambda[i] # This moves the solution up along the efficient frontier
        H <- cOOS # the covariance matrix
        OOS.sol <- solve.QP(H, f, A, b, meq=1)
        OOS.Wts[i,] <- OOS.sol$solution
}

# use replicate and element-wise multiplication to find the Expected Returns
OOSERet2 <- rowSums(OOS.Wts * t(replicate(nrow(OOS.Wts),mOOS)))
# use matrix multiplication to find the risk
# loop over each asset pair and compute risk each weigth vector
# preallocate zero (could use NA)
OOSERisk2 <- numeric(nrow(OOS.Wts))
OOSERisk2[] <- NA
# pre-compute covariance matrix
OOS.Sigma2 <- matrix(cOOS,nrow(cOOS),nrow(cOOS))
# compute the portfolio volatility
for (i in 1:nrow(OOS.Wts)) {OOSERisk2[i] <- as.numeric(OOS.Wts[i,]) %*% OOS.Sigma2 %*% as.numeric(OOS.Wts[i,])}
# add the return and risk columns to data frame
OOSdfA2 <- cbind(sqrt(OOSERisk2),OOSERet2)
# compute the sharpe ration
OOS.SR <- (OOSERet2 - rfrOOS) / sqrt(OOSERisk2)
# add risk return curves to plot
points(OOSdfA2,pch = 19, type = "l", lwd=1, col ="blue")
# annotate
text(0.02,0.012,labels='IS Efficient Frontier',pos = 4, col ="blue")

## 9. Find The Sharpe Ratio maximising portfolio
# Initial values as fully invested equally weighted portfolio
Ones0 <- seq(1,1,length.out = length(mOOS))
# equally weighted portfolio
Correct.Wts0 <- Ones0 / length(Ones0)
# unit vector
e <- rep(1,length(Correct.Wts0)) # useful matrix (ones)
# initialise the weights
Correct.Wts <- matrix(NA,1,length(mOOS))
# 8. Maximise the Sharpe Ratio
fn0 <- function(x) {return(-(x%*% mOOS - rfrOOS)/ sqrt(x %*% cOOS %*% x) )}
# Fully Invested + Return Target
heq0 <- function(x) {return(x %*% e - 1)} # fully invested
# Use SQP to solve for the tangency portfolio
cor.soln <- slsqp(Correct.Wts0, fn = fn0, gr = NULL, # target returns
              lower = rep(0,length(Correct.Wts0)), # no short-selling
              upper = rep(1,length(Correct.Wts0)), # no leverage
              heq = heq0, # fully invested constraint function
              control = list(xtol_rel = 1e-8),
              nl.info = TRUE) # SQP
Correct.Wts <- cor.soln$par
# print the weight matrix
### 110 ITERATIONS
print(Correct.Wts)


###############################################################################
############## OOS Equally weighted portfolio  ###################### 
###############################################################################


# Initialise weights as an equally weighted portfolio
Wts0 <- as.vector(seq(1,1,length.out = length(mOOS)) / length(mOOS))


# Returns for any portfolio with equal weights
OOS.EquiRet <-sum(Wts0*mOOS)

## Risk for equal port
# remove row/col names
OOS.Sigma2 <- matrix(cOOS,nrow(cOOS),nrow(cOOS))
#IS.EquiRisk2 <- matrix(NA,length(Wts0))

OOS.EquiRisk2 <- t(Wts0) %*% OOS.Sigma2 %*% Wts0

points(sqrt(OOS.EquiRisk2),OOS.EquiRet, col = "orange",cex = 2, lwd = 3, pch = 4)

####################
############## OOS HRP portfolio  ###################### 
###############################################################################

###############################################################################
############## IS HRP portfolio  ###################### 
###############################################################################
source("HRP Fn.R")


# get correlation matrix
OOS.VarMat <- var(OOScltsGRet)
# annualized
OOS.corMat <- cov2cor(OOS.VarMat)


#### CLUSTERING ###

OOS.HRP.wts <- HRP_Fn(corr = OOS.corMat, cov = OOS.VarMat)

OOSm <- matrix(mOOS,nrow = length(mOOS), ncol = 1)

# Returns for HRP portfolio
OOS.HRP.Ret <-sum(OOS.HRP.wts%*%mOOS)
# IS.Sigma2 <- matrix(covar,nrow(covar),nrow(covar))
OOS.HRP.risk2 <- t(OOS.HRP.wts) %*% OOS.Sigma2 %*% OOS.HRP.wts

points(sqrt(OOS.HRP.risk2),OOS.HRP.Ret, col = "darkgreen",cex = 2,lwd = 3, pch = 1)
#points(0.02,0.298)

#####################################################################################################

### Use weights from IS:

IS.HRP.wts <- as.matrix(IS.HRP.wts, nrow= length(IS.HRP.wts), ncol =1)

#### CLUSTERING ###

mOOS <- matrix(mOOS,nrow = length(mOOS), ncol = 1)

# Returns for HRP portfolio
OOS.HRP.Ret <-t(IS.HRP.wts)%*%mOOS
# IS.Sigma2 <- matrix(covar,nrow(covar),nrow(covar))
OOS.HRP.risk2 <- t(IS.HRP.wts) %*% OOS.Sigma2 %*% IS.HRP.wts

points(sqrt(OOS.HRP.risk2),OOS.HRP.Ret, col = "darkgreen",cex = 2, lwd = 3, pch = 4)

legend('bottomright', legend = c("IS Sharpe Ratio","OOS Sharpe Ratio", "IS HRP",'OOS HRP', "Equally weighted", "IS Efficient Frontier", "IS Sharpe Ratio"), 
       col = c('magenta','magenta', 'darkgreen',"darkgreen","orange","blue", "red"), lwd = c(3,3,3,3,3,1,1), lty = c(NA,NA,NA, NA, NA,"solid", "dashed"), cex = 0.9, pch = c(1,4,1,4,4, NA,NA))
####################

## 10. Compute return and risk for each weight vector
Correct.OOS.ERetPSR <- Correct.Wts %*% mOOS
# use matrix multiplication to find the risk
Correct.OOS.ERiskPSR <- Correct.Wts %*% OOS.Sigma2 %*% Correct.Wts 
# add risk return curves to plot
points(sqrt(Correct.OOS.ERiskPSR),Correct.OOS.ERetPSR,pch = 1, type = "p", col ="magenta", cex = 2, lwd = 3)


## 10. Compute return and risk for each weight vector
OOS.ERetPSR <- Wts %*% mOOS
# use matrix multiplication to find the risk
OOS.ERiskPSR <- Wts %*% OOS.Sigma2 %*% Wts 
# add risk return curves to plot
points(sqrt(OOS.ERiskPSR),OOS.ERetPSR, type = "p", col ="magenta", cex = 2, lwd =3, pch = 4)
points(sqrt(OOS.EquiRisk2),OOS.EquiRet, col = "orange",cex = 2, lwd = 3, pch = 4)

text(0.04,0.025, "SR", col = "red")
# include text
#text(sqrt(OOS.ERiskPSR), OOS.ERetPSR,labels='SOO Maximal Sharpe Ratio', col = 'magenta', cex= 1, pos = 4)

#***************************************
#* 80:20 Stats:

## 10. Compute return and risk for each weight vector
big.OOS.ERetPSR <- big.Wts %*% mOOS
# use matrix multiplication to find the risk
big.OOS.ERiskPSR <- big.Wts %*% OOS.Sigma2 %*% big.Wts 
# add risk return curves to plot
points(sqrt(big.OOS.ERiskPSR),big.OOS.ERetPSR, type = "p", col ="magenta", cex = 2, lwd =3, pch = 5)


#***************************************


# plot the Sharpe Ratio against risk levels
par(new = T)
plot(sqrt(OOSERisk2),OOS.SR, axes=F, type ="l", lty = 2, col="red",xlab=NA, ylab=NA,xlim=c(0.0,0.15))
axis(side = 4)

mtext(expression(Sharpe ~ Ratio: ~ ~ frac(mu-R[f],sigma )), side = 4, col ="red", line = 3)

