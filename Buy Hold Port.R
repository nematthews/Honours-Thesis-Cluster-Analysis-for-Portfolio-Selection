#### Buy-Hold Port

# 3. Buy-Hold with SR (uses the first window 1:74)
m.1     <- colMeans(tsGRet[1 :Window,], na.rm = T)
covar.1 <- var(tsGRet[1:Window,], na.rm = T)

BHWts <- SR.fn(m = m.1, covar = covar.1,RFR = RFR)
#####################################################################
# Pull first month initial weight
Shift_tsBHWts_first <- as.data.frame(matrix(BHWts, 1, ncol = length(BHWts)))

### Export Data for checking: Shifting Window 
write_xlsx(Shift_tsBHWts_first,"/Users/Ninamatthews/Desktop/THESIS/Excel Data Check/Shift_tsBHWts_first.xlsx")


# Shift window weights beginning of months
Shift_tsBH_Wts0        <- tsGRet* 0  # for i-th month
# Shift window weights end of months
Shift_tsBH_WtsEnd        <- tsGRet* 0 # for i-th month
# Insert initial optimized weights
Shift_tsBH_Wts0[Window,] <- BHWts
# BH Portfolio total Returns per rolled month
Shift_tsBHRet <- Shift_tsPRet 
names(Shift_tsBHRet) <- "Shift Buy-Hold"
# Step forwards by a month using loop
tot <- dim(tsGRet)


for (i in Window:(tot[1]-1)){
        
        #### need new stats of new window each time
        # Shifting window
        m.i <- colMeans(tsGRet[(1+i-Window) :i,], na.rm = T)
        covar.i <- var(tsGRet[(1+i-Window) :i,], na.rm = T)
        corr.i <- cov2cor(covar.i)
        
        #### Call opt using above inputs, store returned wts

        # 3. BH: Initial weights calc outside
        # Insert initial optimized weights
        Shift_tsBH_Wts0[Window,] <- BHWts

        
        #### Calc + store realised returns, wts * actual next mnths returns
       
        #3. BH realised returns
        ## Realised Port returns (Sum across all assets) 
        Shift_tsBHRet[i] <- Shift_tsBH_Wts0[i,] %*% t(tsGRet[i,])
        
        # Calc month end weight
        Shift_tsBH_WtsEnd[i,] <- (tsGRet[i,]*Shift_tsBH_Wts0[i,])+Shift_tsBH_Wts0[i,]
        # Calc month i+1 weights
        Shift_tsBH_Wts0[(i+1),] <- Shift_tsBH_WtsEnd[i,]/sum(Shift_tsBH_WtsEnd[i,])
        
}


