#rm(list = ls())


# ==============================================================================
# QualitativeIndex The function takes as an input a table of values for
# three vulnerability dimensions on a long format. The function calculated the
# percentiles of each dimension using simulating 1000 point from observed
# mean and std. Based on percentiles, each value of each dimension is
# classify as L, I, H, VH, and E. 
# The function returns data frame T that includes original table and qualitative 
# index at the dimension level and data frame of summary of statistics of 
# observed and simulated values
# NOTE that the RCP scenarios must have listed consecutively from RCP26_2050:RCP85_2100
# Smit Vasquez Caballero
# 6/15/2021 
# ==============================================================================

QualitativeIndex <- function(T){
   # change data frame such that each dimension are listed as columns 
   # ================================
  # T$MPA <- factor(Habitat$MPA)
   #library(tidyr)
   # from wide to long along the RCP scenarios
   T <-gather(T,Scenario,Index, RCP26_2050:RCP85_2100)
   # from long to wide along the Dimensions
   T <-spread(T,Dimension,Index)
   
   # create a data frame to collect statistics 
   SStat <- data.frame(Dimension = rep(NA,3), 
                       ObsMean = rep(NA,3),
                       ObsSD = rep(NA,3), 
                       SimMean = rep(NA,3), 
                       SimSD = rep(NA,3), 
                       per20 = rep(NA,3), 
                       per40 = rep(NA,3), 
                       per60 = rep(NA,3), 
                       per80 = rep(NA,3))
   
   
   # create new variable names
   DName <- colnames (T) [ (ncol(T)-2) : ncol(T)]
   NewColName <- paste ("Q", DName, sep="_")
   T[NewColName] <- NA
   
   # Calculate thresholds values and fill out new variables with qualitative index
   # =================================
   for (i in 1:length(DName)) {
     ansmean <- mean(T[[DName[i]]], na.rm = TRUE)
     anssd <- sd(T[[DName[i]]], na.rm = TRUE)
     
     # Simulation
     # =============
     # create simulated data using observed statistics
     set.seed(123) # for reproducibility
     Simdata <- rnorm(n=1000, mean = ansmean, sd = anssd)
     
     # calculate percentiles
     P <- quantile(Simdata, prob = c(.2, .40, .60, .80))
     
     for (j in 1:nrow(T)) {
       if (T[[DName[i]]][j] <= P[1]) {
         T[[NewColName[i]]][j] <- 'Low'
       }
       else if ((T[[DName[i]]][j] > P[1]) & 
                (T[[DName[i]]][j] <= P[2])) {
         T[[NewColName[i]]][j] <- 'Intermediate'  
       }
       else if ((T[[DName[i]]][j] > P[1]) & 
                (T[[DName[i]]][j] <= P[2])) {
         T[[NewColName[i]]][j] <- 'Intermediate'  
       }
       else if ((T[[DName[i]]][j] > P[2]) & 
                (T[[DName[i]]][j] <= P[3])) {
         T[[NewColName[i]]][j] <- 'High'  
       }
       else if ((T[[DName[i]]][j] > P[3]) & 
                (T[[DName[i]]][j] <= P[4])) {
         T[[NewColName[i]]][j] <- 'Very High'  
       }
       else if (T[[DName[i]]][j] > P[4]) {
         T[[NewColName[i]]][j] <- 'Extreme'  
       }
       
     } # end of row loop
     
     # Collect fill in table with statistics
     SStat$Dimension[i] <- DName[i]
     SStat$ObsMean[i] <- ansmean
     SStat$ObsSD[i] <- anssd
     SStat$SimMean[i] <- mean(Simdata)
     SStat$SimSD[i] <- sd(Simdata)
     SStat$per20[i] <- P[1]
     SStat$per40[i] <- P[2]
     SStat$per60[i] <- P[3]
     SStat$per80[i] <- P[4] 
   } # end of column loop 

# create list of return objects
objlist <- list(T,SStat)
return (objlist)   
} # end of function

# ==============================================================================

# set directory and Juan's table 
##setwd("c:/Users/Usuario/Documents/Francesca/IndexCalculation")
##library("readxl")
JuanMatrix <- read_excel(here::here("C:/Users/jbuen/OneDrive/Documents/GitHub/vulnerability-index/postanalysis/xan/tool_inputs/JuanMatrix.xlsx"))

# =====================================
#           Habitat 
# =====================================
Habitat <- read_excel(here::here("inputs/Habitat.xlsx")) # load file

# obtain index per dimension
list1 <- QualitativeIndex(T=Habitat)
THabitat <- list1[[1]]

# Add overall Index based on dimension index
THabitat <- merge(THabitat,JuanMatrix, by.x=c("Q_Ecological adaptive capacity",
                                "Q_Ecological sensitivity", 
                                "Q_Exposure" ), 
           by.y=c("Q_AC",
                  "Q_Sensitivity", 
                  "Q_Exposure"))

# rearrange variables
THabitat <- THabitat[, c("MPA","Habitat","Scenario",
           "Exposure", "Ecological sensitivity", "Ecological adaptive capacity",
           "Q_Exposure","Q_Ecological sensitivity", "Q_Ecological adaptive capacity", 
           "Index")]

# =====================================
#           Users
# =====================================
Users <- read_excel(here::here("inputs/Users.xlsx"))
list2 <- QualitativeIndex(Users)
TUsers <- list2[[1]]
TUsers <- merge(TUsers, JuanMatrix, by.x = c("Q_Ecological adaptive capacity",
                                             "Q_Ecological sensitivity", 
                                             "Q_Exposure"), 
                                    by.y = c("Q_AC",
                                              "Q_Sensitivity", 
                                              "Q_Exposure"))
TUsers <- TUsers[, c("MPA","Users","Scenario",
                     "Exposure", "Ecological sensitivity", "Ecological adaptive capacity",
                     "Q_Exposure","Q_Ecological sensitivity", "Q_Ecological adaptive capacity", 
                     "Index")]


# =====================================
#           Species
# =====================================
Species <- read_excel(here::here("inputs/Species.xlsx"))
list3 <- QualitativeIndex(Species)
TSpecies <- list3[[1]]
TSpecies <- merge(TSpecies, JuanMatrix, by.x = c("Q_Ecological adaptive capacity",
                                             "Q_Ecological sensitivity", 
                                             "Q_Exposure"), 
                by.y = c("Q_AC",
                         "Q_Sensitivity", 
                         "Q_Exposure"))
TSpecies <- TSpecies[, c("MPA","Species","Scenario",
                     "Exposure", "Ecological sensitivity", "Ecological adaptive capacity",
                     "Q_Exposure","Q_Ecological sensitivity", "Q_Ecological adaptive capacity", 
                     "Index")]


# =====================================
#           SST
# =====================================
SST <- read_excel(here::here("inputs/SST.xlsx"))
list4 <- QualitativeIndex(SST)
TSST <- list4[[1]]
   # NOTE: Name of dimension change from ecological to Social
TSST <- merge(TSST, JuanMatrix, by.x = c("Q_Social adaptive capacity",
                                         "Q_Social sensitivity", 
                                         "Q_Ecological Vulnerability"), 
              by.y = c("Q_AC",
                       "Q_Sensitivity", 
                       "Q_Exposure"))
TSST <- TSST[, c("MPA","Scenario",
                 "Ecological Vulnerability", "Social sensitivity", "Social adaptive capacity",
                 "Q_Ecological Vulnerability", "Q_Social sensitivity", "Q_Social adaptive capacity", 
                 "Index")]


# =====================================
#           MWH
# =====================================
MHW <- read_excel(here::here("inputs/MHW.xlsx"))
list5 <- QualitativeIndex(MHW)
TMHW <- list5[[1]]
TMHW <- merge(TMHW, JuanMatrix, by.x = c("Q_Social adaptive capacity",
                                         "Q_Social sensitivity", 
                                         "Q_Ecological Vulnerability"), 
              by.y = c("Q_AC",
                       "Q_Sensitivity", 
                       "Q_Exposure"))
TMHW <- TMHW[, c("MPA","Scenario",
                 "Ecological Vulnerability", "Social sensitivity", "Social adaptive capacity",
                 "Q_Ecological Vulnerability", "Q_Social sensitivity", "Q_Social adaptive capacity", 
                 "Index")]


######################################################################
## save output tables
write.csv (TMHW,     here::here("outputs/TMHW_QI.csv"))
write.csv (TSST,     here::here("outputs/TSST_QI.csv"))
write.csv (TSpecies, here::here("outputs/TSpecies_QI.csv"))
write.csv (TUsers,   here::here("outputs/TUsers_QI.csv"))
write.csv (THabitat, here::here("outputs/THabitat_QI.csv"))

