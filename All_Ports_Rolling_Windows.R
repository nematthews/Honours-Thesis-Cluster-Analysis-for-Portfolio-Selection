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

### Rolling Windows ###
# Window of length of half of data (rounded if odd number)
# Shifting Window
# Growing Window

#############################################################################
############## Needed packages ###################### 
#############################################################################

rm(list=ls())
library("writexl")
library(zoo)
library(xts)
library(timeSeries) 
library(rbenchmark)
library(nloptr) # for SQP
library(quadprog) # for QP
library(ggplot2)


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
# reference out the tickers of interest
tsGRet <- tsGRet[,Entities] # BOND + EQTY. INDEX. PORTFOLIO
## 2. Compute the Geometric mean
# Expected return for each asset
m <- colMeans(tsGRet) # n x 1 (where n = # assets)
covar <- var(tsGRet, na.rm=TRUE) # n x n
s <- colStdevs(tsGRet) # standard deviation of assets for plot


# Plot 2003-2015 monthly
plot(tsGRet,plot.type = c("single"), 
     format = "auto", 
     at=pretty(tsTAA), 
     ylab = "Returns",
     main = "TRI for sectors")

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

# Fn that takes in 74 months of data and gets weights for the 75 month.
# Use the actualized return in the data in month 75


#######################################################
# ************** 1. Equally Weighted **************

#######################################################
# ************** 2. SR Maximizing MV **************

RFR  <- colMeans(tsRRF) # Risk Free Rate (STEFI)
#Mkt.R <- colMeans(tsMKT, na.rm=TRUE) # "Market" port E[R]
#Mkt.V <- colStdevs(tsMKT, na.rm=TRUE) # "Market" port Var


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

# Run optimization in month "1" and use these weights every month

#BHwts <- SR.fn()

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
################## Rolling Window (length = n/2) #################
###############################################################################

# need a for loop to step forward by a month, shifting window

###### initialize storage and inputs:
Window  <- round(nrow(tsGRet)/2) # N/2, using half data as window (74)
N <- dim(tsGRet)

#######################################################
# ************** Shifting Window **************

#1. Equally Weighted
Shift_tsERet        <- tsGRet[,1]*0
names(Shift_tsERet) <- " Shift Equally Wt"

# 2. SR storage
Shift_tsWts         <- tsGRet* 0 # weight per asset at each time step
Shift_tsPRet        <-  tsGRet[,1]*0 # store Return that we would get if we used said w's
names(Shift_tsPRet) <- " Shift SR "

# 3. Buy-Hold with SR (uses the first window 1:74)
# m.1     <- colMeans(tsGRet[1 :Window,], na.rm = T)
# covar.1 <- var(tsGRet[1:Window,], na.rm = T)
# 
# BHWts <- SR.fn(m = m.1, covar = covar.1,RFR = RFR)
# Shift_tsBHRet <- Shift_tsPRet # BH Returns per rolled month
# names(Shift_tsBHRet) <- "Shift Buy-Hold"
# # # Step forwards by a month using loop
# tot <- dim(tsGRet)

#4. HRP storage
Shift_HRP_Wts         <- tsGRet* 0 # weight per asset at each time step
Shift_HRP_PRet        <-  tsGRet[,1]*0 # store Return that we would get if we used said w's
names(Shift_HRP_PRet) <- " Grow HRP "

#5. Constant Mix Port
m.1     <- colMeans(tsGRet[1 :Window,], na.rm = T)
covar.1 <- var(tsGRet[1:Window,], na.rm = T)

CMWts <- SR.fn(m = m.1, covar = covar.1,RFR = RFR)
Shift_tsCMRet <- Shift_tsPRet # CM Returns per rolled month
names(Shift_tsCMRet) <- "Shift Constant Mix"

# # Step forwards by a month using loop
tot <- dim(tsGRet)

# loop deals with a single month at a time that is the month after the window
for (i in Window:(tot[1]-1)){
        
        #### need new stats of new window each time
        # Shifting window
        m.i <- colMeans(tsGRet[(1+i-Window) :i,], na.rm = T)
        covar.i <- var(tsGRet[(1+i-Window) :i,], na.rm = T)
        corr.i <- cov2cor(covar.i)
        
        #### Call opt using above inputs, store returned wts
        #1.  Equal
        EWts <- rep(1/N[2], length.out = N[2])
        #2. SR
        Shift_tsWts[i,] <- SR.fn(m = m.i, covar = covar.i, RFR = RFR)
        # 3. BH 
        
        #4. HRP
        Shift_HRP_Wts[i,] <- HRP_Fn(corr = corr.i, cov = covar.i)
        
        #5. Constant Mix
        # Weights calculated outside of for loop
        
        #### Calc + store realised returns, wts * actual next mnths returns
        #1. Equally weighted realised returns
        Shift_tsERet[i] <- EWts %*% t(tsGRet[i,])
        #2 SR realised returns
        Shift_tsPRet[i] <- Shift_tsWts[i,] %*% t(tsGRet[i,])
        #3. BH realised returns
        #Shift_tsBHRet[i] <- BHWts %*% t(tsGRet[i,])
        #4. HRP realised returns
        Shift_HRP_PRet [i] <- Shift_HRP_Wts[i,] %*% t(tsGRet[i,])
        #5. CM realised returns
        Shift_tsCMRet[i] <- CMWts %*% t(tsGRet[i,])
}

#######################################################
# ************** Growing Window **************

#1. Equally Weighted
Grow_tsERet        <-   tsGRet[,1]*0
names(Grow_tsERet) <- " Grow Equally Wt"

#2. SR storage
Grow_tsWts         <- tsGRet* 0 # weight per asset at each time step
Grow_tsPRet        <-  tsGRet[,1]*0 # store Return that we would get if we used said w's
names(Grow_tsPRet) <- " Grow SR "

# #3. Buy-Hold with SR
# m.1     <- colMeans(tsGRet[(1) :Window,], na.rm = T)
# covar.1 <- var(tsGRet[(1) :Window,], na.rm = T)
# 
# BHWts               <- SR.fn(m = m.1, covar = covar.1,RFR = RFR)
# Grow_tsBHRet        <- Grow_tsPRet
# names(Grow_tsBHRet) <- "Grow Buy-Hold"

#4. HRP storage
Grow_HRP_Wts         <- tsGRet* 0 # weight per asset at each time step
Grow_HRP_PRet        <-  tsGRet[,1]*0 # store Return that we would get if we used said w's
names(Grow_HRP_PRet) <- " Grow HRP "

#5. Constant Mix (SR)
m.1     <- colMeans(tsGRet[(1) :Window,], na.rm = T)
covar.1 <- var(tsGRet[(1) :Window,], na.rm = T)

CMWts               <- SR.fn(m = m.1, covar = covar.1,RFR = RFR)
Grow_tsCMRet        <- Grow_tsPRet
names(Grow_tsCMRet) <- "Grow Constant-Mix"


for (i in Window:(nrow(tsGRet)-1)){
        #### need new stats of new window each time
        
        # Growing Window
        m.i <- colMeans(tsGRet[1:i,], na.rm=TRUE)
        covar.i <- var(tsGRet[1:i,], na.rm=TRUE)
        corr.i <- cov2cor(covar.i)
        
        #### Call opt using above inputs, store returned wts
        #1. Equal
        EWts <- rep(1/N[2], length.out = N[2])
        #2. SR
        Grow_tsWts[i,] <- SR.fn(m = m.i, covar = covar.i, RFR = RFR)
        #3. Buy-Hold
        
        #4. HRP
        Grow_HRP_Wts[i,] <- HRP_Fn(corr = corr.i, cov = covar.i)
        #5. Constant-Mix
        # Weights calculated outside of for loop
        
        #### Calc + store realised returns, wts * actual next mnths returns
        
        #1 Equally weighted realised returns
        Grow_tsERet[i] <- EWts %*% t(tsGRet[i,])
        #2. SR realised returns
        Grow_tsPRet[i] <- Grow_tsWts[i,] %*% t(tsGRet[i,])
        #3. BH realised returns
        #Grow_tsBHRet[i] <- BHWts %*% t(tsGRet[i,])
        #4. HRP realised returns
        Grow_HRP_PRet [i] <- Grow_HRP_Wts[i,] %*% t(tsGRet[i,])
        #5. CM realised returns
        Grow_tsCMRet[i] <- CMWts %*% t(tsGRet[i,])
        
}

###############################################################################
########## Output check  #################
###############################################################################

### Export GReturns from 4 months post window
GRet_4 <- as.data.frame(tsGRet[Window:(Window+4),])
write_xlsx(GRet_4,"/Users/Ninamatthews/Desktop/THESIS/Excel Data Check/GRet_4.xlsx")

#######################################################
# ************** Shift Window **************

# Remove 0 padding: 2009-11-30 to end
clean_Shift_tsWts <- Shift_tsWts[Window:i,] # row 1 fine
clean_Shift_tsPRet <- Shift_tsPRet[Window:i,] ## row 1 fine


# Pull first 4 months
Shift_tsWts_4 <- as.data.frame(clean_Shift_tsWts[1:4,])
Shift_tsPRet_4 <- as.data.frame(clean_Shift_tsPRet[1:4,])

### Export Data for checking: Shifting Window 
write_xlsx(Shift_tsWts_4,"/Users/Ninamatthews/Desktop/THESIS/Excel Data Check/Shift_tsWts_4.xlsx")
write_xlsx(Shift_tsPRet_4,"/Users/Ninamatthews/Desktop/THESIS/Excel Data Check/Shift_tsPRet_4.xlsx")

#######################################################
# ************** Growing Window **************

# Remove 0 padding: 2009-11-30 to end
clean_Grow_tsWts <- Grow_tsWts[Window:i,] # row 1 fine
clean_Grow_tsPRet <- Grow_tsPRet[Window:i,] # row 1 fine


# Pull first 4 months
# ts Real data 4 mnths
Grow_tsWts_4 <- as.data.frame(clean_Grow_tsWts[1:4,])
Grow_tsPRet_4 <- as.data.frame(clean_Grow_tsPRet[1:4,])

### Export Data for checking: Growing Window
write_xlsx(Grow_tsWts_4,"/Users/Ninamatthews/Desktop/THESIS/Excel Data Check/Grow_tsWts_4.xlsx")
write_xlsx(Grow_tsPRet_4,"/Users/Ninamatthews/Desktop/THESIS/Excel Data Check/Grow_tsPRet_4.xlsx")

###############################################################################
########## Plot Equity Curves using Wealth Index  #################
###############################################################################

#par(mfrow=c(2,1))

#######################################################
# ************** Shifting Window PLOT **************

#Compute the wealth index
S_tsIndx <- merge(colCumprods(exp(Shift_tsERet)),colCumprods(exp(Shift_tsPRet)))
S_tsIndx <- merge(S_tsIndx,colCumprods(exp(Shift_tsCMRet)))
S_tsIndx <- merge(S_tsIndx,colCumprods(exp(Shift_HRP_PRet)))
# tsIndx <-colCumprods(exp(tsPRet))
# remove the padded zero from the first half
S_tsIndx <- S_tsIndx[Window:i,]
# print header
head(S_tsIndx)

## Visualise the Equity Curves
# plot the merge indices
plot(S_tsIndx,plot.type = "s", col = c("orange", "black", "blue", "green"))
# turn on the grid
grid()
# title
title(main = "Shifting Window Equity Curve")
# legend
# EQW, SR, CM, HRP
legend("topleft",names(S_tsIndx),col = c("orange","black","blue","green"), lwd = 2, lty = c('solid', 'solid', 'solid', 'solid'), bty = "n")

#######################################################
# ************** Growing Window PLOT **************


#Compute the wealth index
G_tsIndx1 <- merge(colCumprods(exp(Grow_tsERet)),colCumprods(exp(Grow_tsPRet)))
G_tsIndx <- merge(G_tsIndx1,colCumprods(exp(Grow_tsCMRet)))
G_tsIndx <- merge(G_tsIndx,colCumprods(exp(Grow_HRP_PRet)))
# tsIndx <-colCumprods(exp(tsPRet))
# remove the padded zero from the first half
G_tsIndx <- G_tsIndx[Window:i,]
# print header
head(G_tsIndx)

## Visualise the Equity Curves
# plot the merge indices
plot(G_tsIndx,plot.type = "s", col = c("orange", "black", "blue", "green"))
# turn on the grid
grid()
# title
title(main = "Growing Window Equity Curve")
# legend
# EQW, SR, CM, HRP
legend("topleft",names(G_tsIndx),col = c("orange","black","blue","green"), lwd = 2, lty = c('solid', 'solid', 'solid', 'solid'),bty = "n")

############################################################################### 

# DO TO list:

### NINA
# 1. Need to change BH to constant mix portfolio (DONE)
# 2. Add proper BH port
# 3. Compute and plot turnovers

### SIPHE
# 1. Geometric Returns
# 2. Winzorise
# 3. Covar