# Siphesihle Cele & Nina Matthews 
# University of Cape Town Statistical Science
# Honours Project
# 2021 November, Cape Town

#clear environment 
rm(list = ls())
#load packages needed for analysis
library(readxl)
library(ggplot2)
library(dplyr)
library(xlsx)
library(writexl)
library(Rcpp)
library(timeSeries)

#This code was used to read data into R environment and perform Exploratory Data Analysis (EDA)
# Exhange traded fund data is implemented in algorithms used.
#Frequency is monthly data spanning 20 years from 2001 to 2021
# For each section, comments are added in order to give context on code
#Code patterns used from other individuals has been referenced accordingly. 

# What we set out to do in this code is the following:
# 1. Read the data into environment 
# 2. remove redundant columns and rows (cleaning)
# 3. Visualize data availability
# 4. Perform truncation were necessary.
# 5. compute necessary variables needed in analysis
# 6. compute winsorizing algorithm
# 7. produce covariance matrix needed for further algorithms in simulation section.


#######################################################################
#####################  Part 1,2: Get Data in Environment & EDA#########

DataExtraction = read_excel("ETF-DATA-2001-2021.xlsx")
#extract ticker names 
tickerNames = (DataExtraction[3,])
#remove the NA col names and retain the ticker names only in the column.
# note, these columns are above the start of observations so no data lost here.

tickerNamesForTRI= tickerNames[!is.na(tickerNames)]

# extract column names from data set 
colnames(DataExtraction)= DataExtraction[5,]
#start data reading from the ticker names and not the empty rows
DataExtraction= DataExtraction[6:nrow(DataExtraction),]

#return dates on data frame
DataExtraction$Dates = openxlsx::convertToDate(DataExtraction$Dates)

#connect data to be plotted and will convert to binary for plotting.
plotData = DataExtraction[,c(1,seq(2,ncol(DataExtraction), by = 3))]
#the data set to be visualized, named with ticker names as you clean.
colnames(plotData)= c("Dates",tickerNames[seq(2,ncol(DataExtraction), by = 3)])

#TRI data that is going to be winsorized, pre-truncation
TRI_Data_Winsorizing <- plotData
#TRI data that will be winsorized post, truncation. Truncation happens in 2015.
TRI_Data_Winsorizing_2015 <- TRI_Data_Winsorizing[c(162:240),]

#Extract the dates from both data sets to check whether truncation was performed correctly.
date = plotData$Dates
date_2015 = TRI_Data_Winsorizing_2015$Dates
 
#Print the dates out for both data sets to see if dates preserved and truncated dates are 
#truncated properly.
date
date_2015

#Extract all data in data frame except dates, these will be included later...
plotData = plotData[,-1]
TRI_Data_Winsorizing_2015 = TRI_Data_Winsorizing_2015[-1]

################################################################################
########################## Part 3,4: Visualization Process #######################

#Create Binary data points for the data frames, assign 1's for places where
#Data exists and 0's where NA's are observed.
plotData = ifelse(plotData == "#N/A N/A", 0,1)
#repeat the same process for truncated data.
TRI_Data_Winsorizing_2015_Preserved <-TRI_Data_Winsorizing_2015
TRI_Data_Winsorizing_2015 = ifelse(TRI_Data_Winsorizing_2015 == "#N/A N/A", 0,1 )

#Visualize the data frames to see whether the binary dataset was created properly.
#Repeat process for both datasets.
plotData = data.frame(date,plotData)
TRI_Data_Winsorizing_2015 = data.frame(date_2015,TRI_Data_Winsorizing_2015)


#view plot data if still correct, print data
plotData
TRI_Data_Winsorizing_2015

# ggplot process.

#Work with preTruncatedData create dataframe to be plotted using ggplot for 
#Visualization.
preTruncData = data.frame(matrix(ncol = 3, nrow = 0))
names(preTruncData) = c("Date", "Check", "Ticker")
for(i in 2:52){
  temp = data.frame(date, plotData[,i], tickerNamesForTRI[i-1])
  names(temp) = c("Date", "Check", "Ticker")
  preTruncData = rbind(preTruncData, temp)
}

#Work with postTruncatedData create dataframe to be plotted using ggplot for 
#Visualization.
postTruncData = data.frame(matrix(ncol = 3, nrow = 0))
names(postTruncData) = c("Date", "Check", "Ticker")
for(j in 2:52){
  temp_2015 = data.frame(date_2015, TRI_Data_Winsorizing_2015[,j], tickerNamesForTRI[j-1])
  names(temp_2015) = c("Date", "Check", "Ticker")
  postTruncData = rbind(postTruncData, temp_2015)
}


preTruncData$Check=as.factor(preTruncData$Check)
#Full data Set
ggplot(preTruncData, aes(x = Ticker, y = Date, col = Check)) + 
  geom_line(size = 2.5) +
  theme(axis.text.x = element_text(colour = 'black', angle = 90, size = 8, hjust = 0.5, vjust = 0.5),axis.title.x=element_blank()) +
  scale_color_gradient(low = "grey", high = "red") +
  ggtitle("ETF Data Visualization") + scale_y_discrete()
  ggsave("plot_long.pdf")

# TRuncated data set
ggplot(postTruncData, aes(x = Ticker, y = Date, col = Check)) + 
  geom_line(size = 2.5) +
  theme(axis.text.x = element_text(colour = 'black', angle = 90, size = 8, hjust = 0.5, vjust = 0.5),axis.title.x=element_blank()) +
  scale_color_gradient(low = "white", high = "navy") +
  ggtitle("Exchange Traded Funds Timeline") 
  ggsave("plot_long.pdf")




################################################################################
################### Part 5,6: Necessary variables/Winsorizing ##################  
  
#As version control, create a temporary data set so, original data is not lost.  
temp = TRI_Data_Winsorizing_2015_Preserved
# Similar as line 142, but remove asset with no Data from data frame.
TRI_Data_Winsorizing_2015_Preserved = temp[,-15] 
# Recode '#N/A N/A' into pure NA observations.
TRI_Data_Winsorizing_2015_Preserved[TRI_Data_Winsorizing_2015_Preserved == "#N/A N/A"] = NA


###############################################################################
############ subset data: Proof of Principle to get Covariance Matrix##########

#Only execute this after attempting full study and failing to get covariance matrix
#Wth all 49 assets.
TRI_TryTest = TRI_Data_Winsorizing_2015_Preserved
#Keeps assets with full data in specified timeline.
TRI_TryTest = TRI_TryTest[,c(7,8,9,10,11,12,17,18,19,24,25,27,28,29,30,32,33,34,35,41,42,46,50)]

##############################################################################


#Create combination variable to create matrix to merge with TRI data
comb = matrix(NA,ncol = ncol(TRI_Data_Winsorizing_2015_Preserved))
#combTs = matrix(NA,ncol = ncol(TRI_TryTest))
for(i in 1:nrow(temp))
{
  a = as.numeric(TRI_Data_Winsorizing_2015_Preserved[i,])
  comb = rbind(comb, a)
}
#Create data frame of the combinatory variable.
comb = as.data.frame(comb[-1,])
colnames(comb) = names(TRI_Data_Winsorizing_2015_Preserved)

#retrieve names
names(comb)
row.names(comb)=date_2015

#Remove virtually empty assets.
comb = comb[,-c(31,40,48)]
#Only run this line of code when you are running the covariance proof of principle code.
#comb = comb[-c(1:53),] 

##################################################################################
##prepare to plot data 

#Convert dataframes into timeSeries data using timeSeries function.
#Alternatively, do not run this code and run data into algorithms as it.
#Problem with second option is that it uses Total return index instead of price index.
TRI_TSD = timeSeries(comb)
plot(TRI_TSD,plot.type = c("single"), 
     format = "auto", 
     at=pretty(TRI_TSD), 
     ylab = "Returns",
     main = "TRI for ETF")


#Converts total return index into Price index
tsIndx = index2wealth(TRI_TSD)
#geometric returns in a continuous time.
tsGRet = diff(log(tsIndx))
#check dimentions of result.
dim(tsIndx)
#plot price index to pick up errors 
plot(tsIndx,plot.type = c("single"), 
     format = "auto", 
     at=pretty(tsIndx), 
     ylab = "Returns",
     main = "TRI for ETF")
#plot geometric returns to pick up outliers.
plot(tsGRet,plot.type = c("single"), 
     format = "auto", 
     at=pretty(tsGRet), 
     ylab = "Returns",
     main = "Continuous Returns")

#leaves me with 3 rows. so no go zone. [This is if we remove assets that have rows with NA's]
#TRI_TSD_NA = na.omit(TRI_TSD)
#dim(TRI_TSD_NA)

##############################################################################
############# Function to calculate Geometric returns ########################

L = matrix(NA, ncol = ncol(comb), nrow = nrow(comb))
 geometric_Returns <- function(data){
#   
   L = matrix(NA, ncol = ncol(data), nrow = nrow(data))
#   
#   
 #Process would be loop through assets.
   for (i in 1:ncol(data)){
#     #loop through each index 
     for (j in 2:(nrow(data))){
#       
#     
       if(is.na(data[j,i])==F & is.na(data[j-1,i])==F){
#         
#         get log differences between assets. should be run, if diff(log) has not been 
         # done in the timeSeries Process
         L[j,i]=log(data[j,i])-log(data[j-1,i])
       }
#     
#       
     }
#     
   }
#   
   return(L)
 }
# 
# returns geometric returns in continuous time 
 g = as.data.frame(geometric_Returns(comb)) 
 
# 
# 




##################################################################################
 # Take the geometric returns in CT as the pre-winsorized data in order to put through winsorizing alg
 #Make data frame.
TRI_Data_Winsorizing_2015_Preserved_pre <- g
 TRI_Data_Winsorizing_2015_preTS <- as.data.frame(tsGRet)
 
# Test case histogram using asset in 5th column
test_1 = hist(TRI_Data_Winsorizing_2015_Preserved_pre[,5])
#Tolerance for convergence
tol <- 0.000001
#function to execute winsorizing on data set, with the following parameters
#data, standard dev, upper bound, lower bound and tolerance
winzorizing_function = function(data, tol){
  
  #loop through all n columns of data set.
  for (i in 1:ncol(data)){
    
    #obtain standard deviation of column we are working on.
    column_Standard_dev = sd(data[,i], na.rm = T)
    diffSD = column_Standard_dev
    #whilst the difference in prev and now SD is greater than tolerance
    while(diffSD > tol){
      # stopping condition irregardless. To ensure we do not run an infinate loop.
      if (diffSD< tol){
        break
      }
      oldSD = sd(data[,i], na.rm = T)
     
      #upper bound of my 99.7% confidence interval for my data
      upper_bound <- mean(data[,i], na.rm = T)+3*(column_Standard_dev)
      #lower bound of my 99.7% confidence interval for my data
      lower_bound <- mean(data[,i], na.rm = T)-3*(column_Standard_dev)
      
      #write a loop that loops through all rows in each column to do winsorizing
      outliers_position = data[,i] > upper_bound | data[,i] < lower_bound
      
      #Print cases, each loop to test process functionality
      #print(i)
      #print(lower_bound)
      #print(upper_bound)
      for (j in 1:nrow(data)){
        
        if (is.na(data[j,i])==F & data[j,i] > upper_bound){
          data[j,i] = upper_bound
          #else replace with lower because the rest would be less than Upper.
        }  else if(is.na(data[j,i])==F & data[j,i] < lower_bound){
          data[j,i] = lower_bound  
        }
        
      }
      newSD= sd(data[,i], na.rm = T)
      
      
      diffSD = newSD-oldSD
      print(oldSD)
      print(newSD)
      print(diffSD)
    } #while loop ends here
    
  } # column loop ends here
  
  return(data)
}

#Put data through winsorizing function
TRI_Data_Winsorizing_2015_Preserved_post = winzorizing_function(TRI_Data_Winsorizing_2015_Preserved_pre, tol)
#Plot post winsorizing histogram to see whether it resembles normal distribution
test_2 = hist(TRI_Data_Winsorizing_2015_Preserved_post[,5], col="lightblue1",
              border="dodgerblue3", xlab = "Returns", main = "ETF5IT SJ Equity Distribution post Winzorising")

#Remove NA row.
TRI_Data_Winsorizing_2015_Preserved_post=TRI_Data_Winsorizing_2015_Preserved_post[-1,]
row.names(TRI_Data_Winsorizing_2015_Preserved_post)=date_2015
colnames(TRI_Data_Winsorizing_2015_Preserved_post)=colnames(comb[-1,])

################################################################################
##################### Part 7: Produce Covariance Matrix ########################

#exponentiate becuase we used logs before.
TRI_Data_Winsorizing_2015_Preserved_post_EX = exp(TRI_Data_Winsorizing_2015_Preserved_post)
#test histogram to make sure it hasnt changed again.
test_6 = hist(TRI_Data_Winsorizing_2015_Preserved_post[,5])

#get covariance wiht NA values
covTRI_withNA = var(TRI_Data_Winsorizing_2015_Preserved_post_EX, na.rm = TRUE)
#replace NA values with 0's 
covTRI_withNA[is.na(covTRI_withNA)]<-0
# Check condition of covariance matrix
cond.cov <- kappa(covTRI_withNA, exact = TRUE)
covTRI_withNA = covTRI_withNA
#create data frame for Cov Matrix
covTRI_withNA = as.data.frame(covTRI_withNA)
#Write into CSV file as new data, will be put into the other algorithms.
write_xlsx(covTRI_withNA,"/Users/Siphesihle/Desktop/Thesis/DataPostCov2.xlsx")


#test, reading the file and see if it works
testDa <- read_excel("DataPostWin.xlsx")

###############################################################################

#The end.





