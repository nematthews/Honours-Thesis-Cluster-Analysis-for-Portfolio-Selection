# 09/10/2021
# Author: Nina Matthews
# Project: Honours Thesis: Cluster Analysis for Portfolio Construction
# Partner: Siphesihle Cele
# Supervisor: Tim Gebbie

#############################################################################
#############################################################################


### Doc Summary:
# Using PT assignment data to test portfolio construction and rolling windows

### Functions for: ####
# 1. Equally Weighted Port
# 2. SR Maximizing Port
# 3. Buy-Hold
# 4. HRP Port
# 5. Constant Mix Port

### Simulation Windows ###
# Window of length of half of data (rounded if odd number)
# Overlapping Rolling Window
# Growing Window

#############################################################################
############## Needed packages ###################### 
#############################################################################

rm(list=ls())
#library("writexl") # uncomment to write certain results to excel docs
library(zoo)
library(xts)
library(timeSeries) 
library(rbenchmark)
library(nloptr) # for SQP
library(quadprog) # for QP
library(ggplot2)
library(dplyr)
library(lubridate)


###############################################################################
############## DATA format requirements ###################### 
###############################################################################
load(file = "PT-TAA.RData")

## 1. Checking for missing data
######## tsGRet = monthly geometric returns ########
head(is.na(tsGRet))

# define tickers of interest
Entities = colnames(tsGRet)
# remove the money market asset (we will compute excess returns!)
Entities <- Entities[-c(grep('STEFI',Entities))]
Entities <- Entities[-c(grep('ALSI',Entities))]


tsGRet = na.omit(tsGRet)
# reference out the risk-free asset returns
tsRRF <- tsGRet[,'STEFI']
# reference out the "market portfolio"
tsMKT <- tsGRet[,'ALSI']
# reference out the "bonds"
tsBonds <- tsGRet[,'ALBI'] 
## Constant mix data:
tsCM <- cbind(tsMKT,tsBonds)
# reference out the tickers of interest
tsGRet <- tsGRet[,Entities] # BOND + EQTY. INDEX. PORTFOLIO
## 2. Compute the Geometric mean
# Expected return for each asset
m <- colMeans(tsGRet) # n x 1 (where n = # assets)
covar <- var(tsGRet) # n x n
s <- colStdevs(tsGRet) # standard deviation of assets for plot


# Plot 2003-2015 monthly
plot(tsGRet,plot.type = c("single"), 
     format = "auto", 
     at=pretty(tsTAA), 
     ylab = "Returns",
     main = "TRI for sectors")


### Covar condition: 52.03
Cond.Covar <- kappa(covar, exact = TRUE)
Cond.Covar
###############################################################################
## Plots Assets by volatility (sd) and return (mean return)#################
###############################################################################


plot(s, m,
     ylab="Expected Ann. Return [%]",
     xlab ="Ann. Volatility [%]",
     main="Monthly Hist. Risk & Return",
     # plot.type="s",
     ylim = c(0, 0.025), xlim = c(0, 0.12))
# turn on the grid
grid()
# label points
text(s, m,labels=names(m), cex= 1, pos = 4)


###############################################################################
## Optimization Fns to be called in rolling windows #################
###############################################################################

#######################################################
# ************** 1. Equally Weighted **************

#######################################################
# ************** 2. SR Maximizing MV **************

RFR  <- colMeans(tsRRF) # Risk Free Rate (STEFI)

SR.fn <- function(m, covar, RFR = RFR){
        
        ## Check inputs:
        # check dimensions and for square matrix 
        if (length(m) != nrow(covar))
                stop("Mean & covar matrices do not conform")
        
        if (nrow(covar) != ncol(covar)) {
                stop("Covar matrix is not square")
        }
        #Check covar semi-definiteness
        else if(any(eigen(covar)$values < 0))
                stop("Covar matrix is not positive semidefinite")
        
        
        # Initialize wts as equal weight
        Ones.vec <- seq(1,1,length.out = length(m))  # 1 x n
        Wts0  <- Ones.vec / length(Ones.vec) # equally weighted portfolio as starting
        
        Wts  <- matrix(NA,1,length(m)) # initialise weights storage
        
        # 8. Maximise the Sharpe Ratio 
        
        # Objective Fn: SR function: 
        # (ReturnP - Rf)/sd [take the neg as it needs to max]
        fn0  <- function(x) {return(-(x%*% m - RFR)/ sqrt(x %*% covar %*% x) )}
        # Fully Invested + Return Target
        heq0 <- function(x) {return(x %*% Ones.vec - 1)} # fully invested
        
        # Use SQP to solve for the tangency portfolio  
        opt <- slsqp(Wts0, fn = fn0, gr = NULL, # target returns
                     lower = rep(0,length(Wts0)), # no short-selling
                     upper = rep(1,length(Wts0)), # no leverage
                     heq = heq0, # fully invested constraint function
                     control = list(xtol_rel = 1e-8)) # SQP
        Wts <- opt$par
        
        return(Wts)
        
}
#### TEST####
SR.wts <- SR.fn(m = m ,covar = covar, RFR = RFR)
SR.wts
#### WORKS ###

#######################################################
# ************** 3. Buy-Hold **************
# Optimize for first month's weights and hold portfolio over full period. No active rebalancing

#######################################################
# ************** 4. HRP **************

source("HRP Fn.R")

covar_HRP <- var(tsGRet, na.rm=TRUE) # n x n
corr_HRP <- cov2cor(covar_HRP)
test_HRP <- HRP_Fn(corr = corr_HRP, cov = covar_HRP)

#######################################################
# ************** 5. Constant Mix **************

## Rebalance to same initial month 1 weights each month moving forward

###############################################################################
################## Simulation Windows  #################
###############################################################################

###### initialize storage and inputs:
Window  <- 36 # (Can change between 3 to 4 years: 36 -> 48)
N <- dim(tsGRet)
#######################################################
# ************** Overlapping Window **************

#1. Equally Weighted
Overlap_tsERet        <- tsGRet[,1]*0
names(Overlap_tsERet) <- "Equally Weighted CM"

# 2. SR storage
Overlap_tsWts         <- tsGRet* 0 # weight per asset at each time step
Overlap_tsPRet        <-  tsGRet[,1]*0 # store Return that we would get if we used said w's
names(Overlap_tsPRet) <- "SR Maximizing"

# 3. Buy-Hold with SR (uses the first window 1:74)
m.1     <- colMeans(tsGRet[1 :(Window-1),], na.rm = T)
covar.1 <- var(tsGRet[1:(Window-1),], na.rm = T)
# Mnth 1 initial weights
BHWts <- SR.fn(m = m.1, covar = covar.1,RFR = RFR)
# Overlap-window weights beginning of months
Overlap_tsBH_Wts0    <- tsGRet* 0  # for i-th month
# Overlap-window weights end of months
Overlap_tsBH_WtsEnd   <- tsGRet* 0 # for i-th month
# Insert initial optimized weights
Overlap_tsBH_Wts0[Window,] <- BHWts
# BH Portfolio total realised Returns per rolled month
Overlap_tsBHRet <- Overlap_tsPRet 
names(Overlap_tsBHRet) <- "SR Maximizing Buy-Hold"

#4. HRP storage
Overlap_HRP_Wts   <- tsGRet* 0 # weight per asset at each time step
Overlap_HRP_PRet  <-  tsGRet[,1]*0 # store Return that we would get if we used said w's
names(Overlap_HRP_PRet) <- "HRP"

#5. Constant Mix Port
# Initialise weights
CMWts <- matrix(c(0.6,0.4),1,2)
Overlap_tsCMRet <- Overlap_tsPRet # CM Returns per rolled month
names(Overlap_tsCMRet) <- "Balanced Portfolio 60:40 CM"

# Step forwards by a month using loop
tot <- dim(tsGRet)

# loop deals with a single month at a time starting with month (Window)
for (i in Window:(tot[1]-1)){
        
        #### need new stats of new window each time
        # Overlaping window
        m.i <- colMeans(tsGRet[(1+i-Window) :(i-1),], na.rm = T)
        covHRP.i <- cov(tsGRet[(1+i-Window) :(i-1),])
        covar.i <- var(tsGRet[(1+i-Window) :(i-1),], na.rm = T)
        corr.i <- cov2cor(covar.i)
        
        #### Call opt using above inputs, store returned wts
        #1.  Equal
        EWts <- rep(1/N[2], length.out = N[2])
        
        #2. SR
        Overlap_tsWts[i,] <- SR.fn(m = m.i, covar = covar.i, RFR = RFR)
        
        # 3. BH 
        # Insert initial optimized weights using SR max
        Overlap_tsBH_Wts0[Window,] <- BHWts
        
        #4. HRP
        Overlap_HRP_Wts[i,] <- HRP_Fn(corr = corr.i, cov = covHRP.i)
        
        #5. Constant Mix
        # Weights calculated outside of for loop
        
        #### Calc + store realised returns, wts * actual market returns
        
        #1. Equally weighted realised returns
        Overlap_tsERet[i] <- EWts %*% t(tsGRet[i,])
        
        #2 SR realised returns
        Overlap_tsPRet[i] <- Overlap_tsWts[i,] %*% t(tsGRet[i,])
        
        #3. BH realised returns
        ## Realised Port returns (Sum across all assets) 
        Overlap_tsBHRet[i] <- Overlap_tsBH_Wts0[i,] %*% t(tsGRet[i,])
        
        #4. HRP realised returns
        Overlap_HRP_PRet[i] <- Overlap_HRP_Wts[i,] %*% t(tsGRet[i,])
        
        #5. CM realised returns
        Overlap_tsCMRet[i] <- CMWts %*% t(tsCM[i,])
        
        
        #######  Update BH Weights ######
        # Calc month end weight
        Overlap_tsBH_WtsEnd[i,] <- (tsGRet[i,]*Overlap_tsBH_Wts0[i,])+Overlap_tsBH_Wts0[i,]
        # Calc month (i+1) weights
        Overlap_tsBH_Wts0[(i+1),] <- Overlap_tsBH_WtsEnd[i,]/sum(Overlap_tsBH_WtsEnd[i,])
}

#######################################################
# ************** Growing Window **************
# 1. Ew under rolling, will be equivalent
#2. SR storage
Grow_tsWts         <- tsGRet* 0 # weight per asset at each time step
Grow_tsPRet        <-  tsGRet[,1]*0 # store Return that we would get if we used said w's
names(Grow_tsPRet) <- "SR Maximizing"

# 3. Buy-Hold with SR (uses the first window 1:74)
# Grow window weights beginning of months
Grow_tsBH_Wts0        <- tsGRet* 0  # for i-th month
# Grow window weights end of months
Grow_tsBH_WtsEnd        <- tsGRet* 0 # for i-th month
# Insert initial optimized weights
Grow_tsBH_Wts0[Window,] <- BHWts
# BH Portfolio total Returns per rolled month
Grow_tsBHRet <- Grow_tsPRet 
names(Grow_tsBHRet) <- "SR Maximizing Buy-Hold"

#4. HRP storage
Grow_HRP_Wts         <- tsGRet* 0 # weight per asset at each time step
Grow_HRP_PRet        <-  tsGRet[,1]*0 # store Return that we would get if we used said w's
names(Grow_HRP_PRet) <- "HRP"

# #5. Constant Mix
# CM under rolling, will be equivalent


for (i in Window:(tot[1]-1)){
        
        #### need new stats of new window each time
        # Growing Window
        m.i <- colMeans(tsGRet[1:(i-1),], na.rm=TRUE)       
        covHRP.i <- cov(tsGRet[1:(i-1),])
        covar.i <- var(tsGRet[1:(i-1),], na.rm=TRUE)
        corr.i <- cov2cor(covHRP.i)
        
        #### Call opt using above inputs, store returned wts
        #1. Equal
        # Calculated under rolling window, they will be equal
        
        #2. SR
        Grow_tsWts[i,] <- SR.fn(m = m.i, covar = covar.i, RFR = RFR)
        
        # 3. BH: Initial weights calc outside
        # Insert initial optimized weights
        Grow_tsBH_Wts0[Window,] <- BHWts
        
        #4. HRP
        Grow_HRP_Wts[i,] <- HRP_Fn(corr = corr.i, cov = covHRP.i)
        
        #5. Constant-Mix
        # Calculated under rolling window, they will be equal
        
        #### Calc + store realised returns, wts * actual next mnths returns
        #1 Equally weighted realised returns
        # Calculated under rolling window, they will be equa
        #2. SR realised returns
        Grow_tsPRet[i] <- Grow_tsWts[i,] %*% t(tsGRet[i,])
        #3. BH realised returns
        ## Realised Port returns (Sum across all assets) 
        Grow_tsBHRet[i] <- Grow_tsBH_Wts0[i,] %*% t(tsGRet[i,])
        #4. HRP realised returns
        Grow_HRP_PRet[i] <- Grow_HRP_Wts[i,] %*% t(tsGRet[i,])
        #5. CM realised returns
        #Calculated under rolling window, they will be equal
        
        # Calc month end weight
        Grow_tsBH_WtsEnd[i,] <- (tsGRet[i,]*Grow_tsBH_Wts0[i,])+Grow_tsBH_Wts0[i,]
        # Calc month i+1 weights
        Grow_tsBH_Wts0[(i+1),] <- Grow_tsBH_WtsEnd[i,]/sum(Grow_tsBH_WtsEnd[i,])
        
        
}

###############################################################################
########## Output check  #################
###############################################################################

### Export GReturns from 4 months post window
GRet_4 <- as.data.frame(tsGRet[Window:(Window+4),])

# Insert desired file path below:
# write_xlsx(GRet_4,"/Users/Ninamatthews/Desktop/THESIS/Excel Data Check/GRet_4.xlsx")


#######################################################
# ************** Overlap Window **************

# SR Maximizing Port:
# Remove 0 padding: 2009-11-30 to end
clean_Overlap_tsWts <- Overlap_tsWts[Window:i,] # row 1 fine
clean_Overlap_HRP_Wts <- Overlap_HRP_Wts[Window:i,] # row 1 fine
clean_Overlap_tsPRet <- Overlap_tsPRet[Window:i,] ## row 1 fine

# Pull first 4 months
Overlap_tsWts_4 <- as.data.frame(clean_Overlap_tsWts[1:4,])
Overlap_tsPRet_4 <- as.data.frame(clean_Overlap_tsPRet[1:4,])

### Export Data for checking: Overlaping Window 
# Insert desired file path below:
#write_xlsx(Overlap_tsWts_4,"/Users/Ninamatthews/Desktop/THESIS/Excel Data Check/Overlap_tsWts_4.xlsx")
#write_xlsx(Overlap_tsPRet_4,"/Users/Ninamatthews/Desktop/THESIS/Excel Data Check/Overlap_tsPRet_4.xlsx")

#######################################################
# ************** Growing Window **************

# SR Maximizing Port:
# Remove 0 padding: 2009-11-30 to end
clean_Grow_tsWts <- Grow_tsWts[Window:i,] # row 1 fine
clean_Grow_tsPRet <- Grow_tsPRet[Window:i,] # row 1 fine


# Pull first 4 months
# ts Real data 4 mnths
Grow_tsWts_4 <- as.data.frame(clean_Grow_tsWts[1:4,])
Grow_tsPRet_4 <- as.data.frame(clean_Grow_tsPRet[1:4,])

### Export Data for checking: Growing Window
# Insert desired file path below:
#write_xlsx(Grow_tsWts_4,"/Users/Ninamatthews/Desktop/THESIS/Excel Data Check/Grow_tsWts_4.xlsx")
#write_xlsx(Grow_tsPRet_4,"/Users/Ninamatthews/Desktop/THESIS/Excel Data Check/Grow_tsPRet_4.xlsx")

###############################################################################
########## Plot Equity Curves using Wealth Index  #################
###############################################################################

#par(mfrow=c(2,1))

#######################################################
# ************** Overlaping Window PLOT **************

## ALSI
tsALSI <- tsMKT[Window:i,]
# 0 initial return
tsALSI[1,]=0
### ALBI
tsBonds <- tsBonds[Window:i,]
tsBonds[1,]=0

#Compute the wealth index
O_tsIndx <- merge(colCumprods(exp(Overlap_tsERet)),colCumprods(exp(Overlap_tsPRet)))
O_tsIndx <- merge(O_tsIndx,colCumprods(exp(Overlap_tsBHRet)))
O_tsIndx <- merge(O_tsIndx,colCumprods(exp(Overlap_HRP_PRet)))
O_tsIndx <- merge(O_tsIndx,colCumprods(exp(Overlap_tsCMRet)))
# remove zeros from the first half
O_tsIndx <- O_tsIndx[Window:i,]
# Add ALSI
O_tsALL <- merge(O_tsIndx,colCumprods(exp(tsALSI)))
# Add ALBI
O_tsALL <- merge(O_tsALL,colCumprods(exp(tsBonds)))
# print header
head(O_tsALL)


## Visualise the Equity Curves
# plot the merge indices
plot(O_tsALL,plot.type = "s", col = c("orange", "magenta", "blue", "darkgreen", "purple","red", "lightblue") ,lwd = c(1,1,1,3,1,1,1),at = "chic", format = "%Y %b", ylim = c(1,5), xlab = "Time", ylab = "Index")
text(as.POSIXct("2009-02-28"),3.4,"Global Financial Crisis",pos = 2, cex = 1.1, srt = 90)
text(as.POSIXct("2014-01-31"),4.9,"QE tappering",pos = 2, cex = 1.1, srt = 90)
text(as.POSIXct("2012-08-31"),4.9,"Marikana massacre",pos = 2, cex = 1.1, srt = 90)
text(as.POSIXct("2015-07-31"),2.5,"Anitdumping Policy",pos = 2, cex = 1.1, srt = 90)
## Add event lines
abline(v = as.POSIXct("2009-02-28"),lwd = 3,col = "black") #GFC
abline(v = as.POSIXct("2014-01-31"),lwd = 3,col = "black") # QE tappering
abline(v = as.POSIXct("2012-08-31"), lwd = 3,col = "black") # Marikana massacre
abline(v = as.POSIXct("2015-07-31"), lwd = 3,col = "black") # Antidumping
# title
title(main = "Overlapping Rolling Window Equity Curve")
# legend
# EQW, SR, BH, HRP, CM
legend("topleft",c(names(O_tsALL)[1:5], "ALSI (Equity)", "ALBI (Bonds)"),col = c("orange","magenta","blue","darkgreen", "purple","red", "lightblue"), lwd = c(2,2,2,3,2,2,2), lty = c('solid', 'solid', 'solid', 'solid', 'solid','solid','solid'), bty = "o")

#######################################################
# ************** Growing Window PLOT **************


#Compute the wealth index
G_tsIndx <- merge(colCumprods(exp(Overlap_tsERet)),colCumprods(exp(Grow_tsPRet)))
G_tsIndx <- merge(G_tsIndx,colCumprods(exp(Grow_tsBHRet)))
G_tsIndx <- merge(G_tsIndx,colCumprods(exp(Grow_HRP_PRet)))
G_tsIndx <- merge(G_tsIndx,colCumprods(exp(Overlap_tsCMRet)))
# remove the padded zero from the first half
G_tsIndx <- G_tsIndx[Window:i,]
G_tsIndx <- as.timeSeries(G_tsIndx)
# Add ALSI
G_tsALL <- merge(G_tsIndx,colCumprods(exp(tsALSI)))
# Add ALBI
G_tsALL <- merge(G_tsALL,colCumprods(exp(tsBonds)))
# print header
head(G_tsALL)

## Visualise the Equity Curves
# plot the merge indices
plot(G_tsALL,plot.type = "s", col = c("orange", "magenta", "blue", "darkgreen", "purple", "red","lightblue"), lwd = c(1,1,1,3,1,1,1), at = "chic", format = "%Y %b", ylim = c(1,5), xlab = "Time", ylab = "Index")
text(as.POSIXct("2009-02-28"),3.3,"Global Financial Crisis",pos = 2, cex = 1.1, srt = 90)
text(as.POSIXct("2014-01-31"),4.9,"QE tappering",pos = 2, cex = 1.1, srt = 90)
text(as.POSIXct("2012-08-31"),4.9,"Marikana massacre",pos = 2, cex = 1.1, srt = 90)
text(as.POSIXct("2015-07-31"),2.5,"Anitdumping Policy",pos = 2, cex = 1.1, srt = 90)
## Add event lines
abline(v = as.POSIXct("2009-02-28"),lwd = 3,col = "black") #GFC
abline(v = as.POSIXct("2014-01-31"),lwd = 3,col = "black") # QE tappering
abline(v = as.POSIXct("2012-08-31"), lwd = 3,col = "black") # Marikana massacre
abline(v = as.POSIXct("2015-07-31"), lwd = 3,col = "black") # Antidumping

# title
title(main = "Growing Window Equity Curve")
# legend
# EQW, SR, CM, HRP
legend("topleft",c(names(G_tsALL)[1:5], "ALSI (Equity)", "ALBI (Bonds)"),col = c("orange","magenta","blue","darkgreen", "purple","red", "lightblue"), lwd = c(2,2,2,3,2,2,2), lty = c('solid', 'solid', 'solid', 'solid', 'solid','solid','solid'), bty = "o")


###############################################################################
########## Portfolios Turnover  #################
###############################################################################

# Turn over as its proportional to trading costs: wE - w0 = change in wts
# 35 bases points 
# 0.0035 x delta weights (change in weight?) = % cost of rebalancing 

## NEED vec of W0's and vec We's for each port PER simulation type

#######################################################
# ************** Overlaping Window PLOT **************

## Storage for End of Month Weights
# 1. EW
Overlap_tsE_WtsEnd   <- tsGRet* 0 # for i-th month
# 2. SR
Overlap_tsSR_WtsEnd   <- tsGRet* 0 # for i-th month
# 3. BH
## Already have Wts calculated: Overlap_tsBH_WtsEnd
# 4. HRP
Overlap_tsHRP_WtsEnd   <- tsGRet* 0 # for i-th month
# 5. Constant Mix
Overlap_tsCM_WtsEnd   <- tsCM* 0 # for i-th month

## Storage for Monthly % cost p/Portfolio
#1. Equal 
Overlap_tsE_TOver <- tsGRet[,1]*0 
#2. SR 
Overlap_tsSR_TOver <- tsGRet[,1]*0 
#3. BH 
Overlap_tsBH_TOver <- tsGRet[,1]*0 
#4. HRP 
Overlap_tsHRP_TOver <- tsGRet[,1]*0 
#5. CM 
Overlap_tsCM_TOver <- tsCM[,1]*0 

## Storage for adjusted ports realised returns for trading costs
# 1. EW
Overlap_tsE_adj   <- tsGRet[,1]*0 
# 2. 
Overlap_tsSR_adj  <- tsGRet[,1]*0 
# 3. BH
Overlap_tsBH_adj  <- tsGRet[,1]*0 
# 4. HRP
Overlap_tsHRP_adj   <- tsGRet[,1]*0 
# 5. Constant Mix
Overlap_tsCM_adj <- tsCM[,1]*0 

# loop deals with a single month at a time starting with month (Window)
for (i in Window:(tot[1]-1)){
        
        #### End of Month Weights for all ports
        #1.  Equal
        Overlap_tsE_WtsEnd[i,] <- (tsGRet[i,]*EWts)+EWts
        #2. SR
        Overlap_tsSR_WtsEnd[i,] <- (tsGRet[i,]*Overlap_tsWts[i,])+Overlap_tsWts[i,]
        # 3. BH 
        # Overlap_tsBH_WtsEnd
        #4. HRP
        Overlap_tsHRP_WtsEnd[i,] <- (tsGRet[i,]*Overlap_HRP_Wts[i,])+Overlap_HRP_Wts[i,]
        #5. Constant Mix
        Overlap_tsCM_WtsEnd[i,] <- (tsCM[i,]*CMWts)+CMWts
        
        ######## P Turnover per month ## ISSUE LIES HERE
        #1. Equal 
        Overlap_tsE_TOver[i] <- apply(abs(Overlap_tsE_WtsEnd - EWts),1,sum)
        #2. SR
        Overlap_tsSR_TOver[i] <- apply(abs(Overlap_tsSR_WtsEnd - Overlap_tsWts),1,sum)
        #3. BH
        Overlap_tsBH_TOver[i] <- apply(abs(Overlap_tsBH_WtsEnd - Overlap_tsBH_Wts0),1,sum)
        #4. HRP
        Overlap_tsHRP_TOver[i] <- apply(abs(Overlap_tsHRP_WtsEnd - Overlap_HRP_Wts),1,sum)
        #5. CM
        Overlap_tsCM_TOver[i] <- apply(abs(Overlap_tsCM_WtsEnd-CMWts),1,sum)
        
        
        # #### Percentage cost of rebalancing using 35 basis points p/mnth
        #1.  Equal
        Overlap_tsE_cost[i] <- 0.0035* Overlap_tsE_TOver[i]
        #2. SR
        Overlap_tsSR_cost[i] <- 0.0035* Overlap_tsSR_TOver[i]
        # 3. BH
        Overlap_tsBH_cost[i] <- 0.0035* Overlap_tsBH_TOver[i]
        #4. HRP
        Overlap_tsHRP_cost[i] <- 0.0035* Overlap_tsE_TOver[i]
        #5. Constant Mix
        Overlap_tsCM_cost[i] <- 0.0035* Overlap_tsCM_TOver[i]

        ####### Adjusting returns p/mnth
        #1. Equally weighted realised returns
        Overlap_tsE_adj[i] <- Overlap_tsERet[i] - Overlap_tsE_cost[i]
        
        #2 SR realised returns
        Overlap_tsSR_adj[i]  <- Overlap_tsPRet[i] - Overlap_tsSR_cost[i]
        
        #3. BH realised returns
        ## Realised Port returns (Sum across all assets) 
        Overlap_tsBH_adj[i] <- Overlap_tsBHRet[i] - Overlap_tsBH_cost[i]
        
        #4. HRP realised returns
        Overlap_tsHRP_adj[i] <- Overlap_HRP_PRet[i] - Overlap_tsHRP_cost[i]
        
        #5. CM realised returns
        Overlap_tsCM_adj[i] <- Overlap_tsCMRet[i] - Overlap_tsCM_cost[i]
}


##### Plot Overlap Adjusted for Turnover Costs

#Compute the wealth index
O_adj_tsIndx <- cbind(colCumprods(exp(Overlap_tsE_adj)),colCumprods(exp(Overlap_tsSR_adj)))
O_adj_tsIndx <- merge(O_adj_tsIndx,colCumprods(exp(Overlap_tsBH_adj)))
O_adj_tsIndx <- merge(O_adj_tsIndx,colCumprods(exp(Overlap_tsHRP_adj)))
O_adj_tsIndx <- merge(O_adj_tsIndx,colCumprods(exp(Overlap_tsCM_adj)))
# remove zeros from the first half
O_adj_tsIndx <- O_adj_tsIndx[Window:i,]
# Add ALSI
O_adj_tsALL <- merge(O_adj_tsIndx,colCumprods(exp(tsALSI)))
# Add ALBI
O_adj_tsALL <- merge(O_adj_tsALL,colCumprods(exp(tsBonds)))
# print header
head(O_adj_tsALL)

## Visualise the Equity Curves
# plot the merge indices
plot(O_adj_tsALL,plot.type = "s", col = c("orange", "magenta", "blue", "green", "purple","red", "lightblue") ,at = "chic", format = "%Y %b", ylim = c(1,5), xlab = "Time", ylab = "Index")
abline(v = as.POSIXct("2009-02-28"))
abline(v = as.POSIXct("2014-01-31"))
abline(v = as.POSIXct("2015-07-31"))
# title
title(main = "Overlapping Rolling Window Equity Curve")
# legend
# EQW, SR, BH, HRP, CM
legend("topleft",names(O_adj_tsALL),col = c("orange","magenta","blue","green", "purple","red", "lightblue"), lwd = 2, lty = c('solid', 'solid', 'solid', 'solid', 'solid','solid','solid'), bty = "o")


############################################################################### 


# DO TO list:

### NINA
# 1. Compute and plot turnovers 

### SIPHE
# 1. Deflated SR
# 2. PBO
# 3. 


### Relative mix changes over time depending on how the real asset classes in relative perform

### Checking changes