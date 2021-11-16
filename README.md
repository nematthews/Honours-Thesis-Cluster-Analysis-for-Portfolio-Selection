# REPO: Honours-Thesis-Cluster-Analysis-for-Portfolio-Selection
### R Code 

##  Honours Thesis: Comparing the Hierarchical Risk Parity Algorithm and Mean-Variance Portfolio Selection
### Author: Nina Matthews
### Author: Siphesihle Cele
### Supervisor: A/Prof. Tim Gebbie
### University of Cape Town

#############################################################################
#############################################################################

# Needed Files per Section:

## 3. Data Wrangling Section:
Already prepared data:
Dataset 1:  PT-TAA.RData

Cleaning of Dataset 2:
Main File: Dataset2ETFDataWrangling.r

Requires:
- ETF-DATA-2001-2021.xlsx

## 4. Data Science and function validation:
Main File: HRP - comp and tests -1.R 

Requires:
- x_output.csv
- Lopez_data.csv
- corMat.csv
- covMat.csv

## 5. Compatibility: Algorithms and the Data Section:
Main File: HRP implemented for ETF data.R

Requires: 
- DataPostCov2.xlsx

## 6. Dynamic Backtest Section:
Main File: All_Ports_Rolling_Windows.R 

Requires: 
- HRP Fn.R
- PT-TAA.RData

## 7.  Static Backtest Section:
Main File: IS_SR vs OS_SR.R

Requires:
- PT-TAA.RData

#############################################################################
#############################################################################


## Project Summary:

The Hierarchical Risk Parity Algorithm in Comparison To Mean-Variance for Portfolio Selection.

### Data collection and processing
1. Sourced from Bloomberg
2. Data cleaning
3. Winsorizing
4. Geometrically compounding returns
5. Obtaining Covariance matrix

### Data Science
1. Construction of the HRP algorithm
2. Construction of the Mean-Variance SR Maximising algorithm 

### Functions for: 
 1. Equally Weighted Port
 2. SR Maximising Port
 3. Buy-Hold
 4. HRP Port
 5. Constant Mix Port

### Backtest Simulation
1. Overlapping Rolling window 
2. Growing Window
3. Static IS v OOS 

### Performance measures theory
1. Probabilistic SR
2. Deflated SR
3. Probability of Backtest overfitting

#############################################################################


