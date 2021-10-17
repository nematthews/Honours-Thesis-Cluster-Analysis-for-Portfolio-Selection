# 13/10/2021
# Author: Nina Matthews
# Project: Honours Thesis: Cluster Analysis for Portfolio Construction
# Partner: Siphesihle Cele
# Supervisor: Tim Gebbie

#############################################################################
#############################################################################

#### DOC: Insample SR vs OS ####
# 40/60 split

### Training Data

IS <- head(tsGRet,dim(tsGRet)[1]*0.4)
OOS <- tail(tsGRet,89)

