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

# establish 60 split:
# dim(tsGRet)[1]*0.6 = 88.8
ISGret <- head(tsGRet,72)
OOSGRet <- tail(tsGRet,76)

ISRFR <- head(tsRRF,72)
OOSRFR <- tail(tsRRF,76)

ISMKT <- head(tsMKT,72)
OOSMKT <- tail(tsMKT,76)


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

plot(sIS, mIS,
     ylab="Expected Return [%]",
     xlab ="Volatility [%]",
     main="IS Monthly Hist. Risk & Return",
     # plot.type="s",
     ylim = c(0, 0.04), xlim = c(0, 0.12))

# turn on the grid
grid()
# label points
text(sIS, mIS,labels=names(mIS), cex= 1, pos = 4)


############################################################################# Optimal SR Max IS 
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
text(0.025,0.013,labels='Efficient Frontier',pos = 4, col ="blue")

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
              control = list(xtol_rel = 1e-8)) # SQP
Wts <- soln$par
# print the weight matrix
print(Wts)

###############################################################################
############## IS Equally weighted portfolio  ###################### 
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
############## IS HRP portfolio  ###################### 
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

legend('bottomright', legend = c("Sharpe Ratio Maximizing", 'HRP', "Equally weighted", "Efficiency Frontier", "SR vs Risk"), 
       col = c('magenta', 'darkgreen',"orange","blue", "red"), lwd = c(3,3,3,1,1), lty = c('solid', 'solid', "solid","solid", "dashed"))

###############################################################################

## 10. Compute return and risk for each weight vector
IS.ERetPSR <- Wts %*% mIS
# use matrix multiplication to find the risk
IS.ERiskPSR <- Wts %*% IS.Sigma2 %*% Wts 
# add risk return curves to plot
points(sqrt(IS.ERiskPSR),IS.ERetPSR,pch = 1, type = "p", col ="magenta", cex = 2, lwd = 3)
points(sqrt(IS.EquiRisk2),IS.EquiRet, col = "orange",cex = 2, lwd = 3)
# include text
text(sqrt(IS.ERiskPSR), IS.ERetPSR,labels='IS Maximal Sharpe Ratio', col = 'magenta', cex= 1, pos = 2)
text(0.06,0.025, "SR vs Risk", col = "red")

## 11. Plot SML (Security Market Line)
# IS.SR0 <- ((IS.ERetPSR-rfrIS) / sqrt(IS.ERiskPSR))
# IS.eq = function(IS.x){return(as.numeric(rfrIS + c(IS.SR0) * IS.x))}
# IS.x <- seq(from=0,to=0.30,length.out=20)
# points(IS.x, IS.eq(IS.x), type = "l", col="red")
#text(0.09,0.20,labels='Market Line', srt = 45, col = 'red', pos = 4)

# plot the Sharpe Ratio against risk levels
par(new = T)
plot(sqrt(ISERisk2),IS.SR, axes=F, type ="l", lty = 2, col="red",xlab=NA, ylab=NA,xlim = c(0, 0.12))
axis(side = 4)

mtext(expression(Sharpe ~ Ratio: ~ ~ frac(mu-R[f],sigma )), side = 4, col ="red", line = 3)
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
     ylab="Expected Return [%]",
     xlab ="Volatility [%]",
     main="OOS Monthly Hist. Risk & Return",
     # plot.type="s",
     ylim = c(0, 0.026), xlim = c(0, 0.14))

# turn on the grid
grid()
# label points
text(sOOS, mOOS,labels=names(mOOS), cex= 1, pos = 4)


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

points(sqrt(OOS.EquiRisk2),OOS.EquiRet, col = "orange",cex = 2, lwd = 3)

####################
############## OOS HRP portfolio  ###################### 
###############################################################################

### Use weights from IS:

IS.HRP.wts <- as.matrix(IS.HRP.wts, nrow= length(IS.HRP.wts), ncol =1)

#### CLUSTERING ###

mOOS <- matrix(mOOS,nrow = length(mOOS), ncol = 1)

# Returns for HRP portfolio
OOS.HRP.Ret <-t(IS.HRP.wts)%*%mOOS
# IS.Sigma2 <- matrix(covar,nrow(covar),nrow(covar))
OOS.HRP.risk2 <- t(IS.HRP.wts) %*% OOS.Sigma2 %*% IS.HRP.wts

points(sqrt(OOS.HRP.risk2),OOS.HRP.Ret, col = "darkgreen",cex = 2, lwd = 3)

legend('bottomright', legend = c("Sharpe Ratio Maximizing", 'HRP', "Equally weighted", "Efficiency Frontier", "SR vs Risk"), 
       col = c('magenta', 'darkgreen',"orange","blue", "red"), lwd = c(3,3,3,1,1), lty = c('solid', 'solid', "solid","solid", "dashed"))
####################

## 10. Compute return and risk for each weight vector
OOS.ERetPSR <- Wts %*% mOOS
# use matrix multiplication to find the risk
OOS.ERiskPSR <- Wts %*% OOS.Sigma2 %*% Wts 
# add risk return curves to plot
points(sqrt(OOS.ERiskPSR),OOS.ERetPSR,pch = 13, type = "p", col ="magenta", cex = 2)
points(sqrt(OOS.EquiRisk2),OOS.EquiRet, col = "orange",cex = 2)
# include text
text(sqrt(OOS.ERiskPSR), OOS.ERetPSR,labels='SOO Maximal Sharpe Ratio', col = 'magenta', cex= 1, pos = 4)

## 11. Plot SML (Security Market Line)
OOS.SR0 <- ((OOS.ERetPSR-rfrIS) / sqrt(OOS.ERiskPSR))
OOS.eq = function(OOS.x){return(as.numeric(rfrOOS + c(OOS.SR0) * OOS.x))}
OOS.x <- seq(from=0,to=0.30,length.out=20)
points(OOS.x, OOS.eq(OOS.x), type = "l", col="red")
text(0.09,0.20,labels='Market Line', srt = 45, col = 'red', pos = 4)

# plot the Sharpe Ratio against risk levels
par(new = T)
plot(sqrt(OOS.ERisk2),OOS.SR, axes=F, type ="l", lty = 2, col="blue",xlab=NA, ylab=NA,xlim=c(0.05,0.28))
axis(side = 4)

mtext(expression(Sharpe ~ Ratio: ~ ~ frac(mu-R[f],sigma )), side = 4, col ="red", line = 3)
