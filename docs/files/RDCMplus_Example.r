################################################################################################################
################################################################################################################
#this file uses the R scripts from Dr. Andre A. Rupp and Dr. Oliver Wilhelm; Available from:
#http://www.education.umd.edu/EDMS/fac/Rupp/R%20Files%20for%20Mplus%20Input%20File%20Generation.zip
#
#additionally, some scripts have been modified or created by Dr. Jonathan Templin

#this routine: (1) formats the Chapter 9 data in accordance with the R functions expectations
#              (2) writes the Mplus script
#              (3) runs Mplus using the MplusAutomation package
#              (4) collects Mplus output
#              (5) converts Mplus output to LCDM parameters 
#              (6) estimates attribute reliability (Templin and Bradshaw, 2013)
#              (7) creates EAP and MAP estimates for attribute profiles and marginal attributes
################################################################################################################
################################################################################################################

#load MplusAutomation:
if (tryCatch(require("MplusAutomation")) == FALSE){
  install.packages("MplusAutomation")
  library(MplusAutomation)
}

#read in R scripts from Rupp and Wilhelm and Templin:
source("MplusDCM_functions.R")

#ANALYZE EXAMPLE DATA FILE ------------------------------------------------------------------------------

#read in original data and create id variable
ch9data = read.table(file = "ch9data.dat")
colnames(ch9data) = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "TrueClass")
ch9data$id = 1:dim(ch9data)[1]

#read in qmatrix file
qmatrix = read.csv("ch9qmatrix.csv")

#write data set to file for reading by function
write.table(x = ch9data[c("id", "X1", "X2", "X3", "X4", "X5", "X6", "X7")], file = "inputdata.dat", 
            col.names = FALSE, quote = FALSE, row.names = FALSE)

#create output file -- only use columns of qmatrix with 0/1 entries (i.e., no item id column)
myoutput = generate.mplus.input(filename_ds = "inputdata.dat", qmatrix = qmatrix[,2:4], input_file = "ch9mplus.inp", PreAddVars = "id")

#run Mplus
myanalysis = runModels(directory = getwd(), showOutput = TRUE)

#creating function for parsing Mplus output from analysis above
MplusAnalysis = readModels(target = getwd())

#Create LCDM output from Mplus Output
Analysis_Output = Mplus_to_LCDM(MplusAnalysis = MplusAnalysis, qmatrix = qmatrix)

#display LCDM item parameters:
Analysis_Output$LCDM_item_parameters
 
#display attribute reliability estimates
Analysis_Output$attribute_reliability 

#display examinee EAP_estimates for attribute profile
Analysis_Output$EAP_profile_estimates

#display examinee EAP_estimates for marginal attributes
Analysis_Output$EAP_marginal_estimates

#display examinee EAP_estimates for attribute profile
Analysis_Output$MAP_profile_estimates

#display examinee EAP_estimates for marginal attributes
Analysis_Output$MAP_marginal_estimates
