rm(list=ls())

here::i_am("R/model.R")

options(warn=1)
#' Name: model.R
#' Purpose: Check model computations. This script asks the user to select a file, 
#' then processes all the scenarios and saves the outputs in 
#' several excel files.
#' Author: Juan Luis Herrera Cortijo (juan.luis.herrera.cortijo@gmail.com)
#' Date Created: 08/07/2020
#' Last update: 06/04/2021



source(here::here("shiny/Vulnerability_Index/function_library.R"))

require(magrittr)
require(tidyverse)




# Ask the user to select a MPA file

# data_file <- here::here("data/TemplateCC_New.xlsx")#file.choose()
data_file <- file.choose()

# Inputs

model_definitions_file <- here::here("shiny/Vulnerability_Index/model definitions.xlsx")



# Compute the model for all the scenarios in the file selected

scenarios <- process_data_file(.data_file = data_file,
                               .model_definition_file = model_definitions_file,
                               .normalize = TRUE,
                               .na.rm = TRUE,
                               .drop_correlated_indicators = TRUE)


################################################################################
## MPA RESULTS  ################################################################

## SP vulnerability
MPA_vul <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"data")})
## MPA data coverage 
MPA_cov <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"coverage")})
## MPA indicator contribution
MPA_cont <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution") %>% last})
MPA_cont_IND1 <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"data", "water cond")})
## MPA years time series
MPA_year <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Year") %>% pluck(1,1,"data")})
## MPA quantitative indicators
MPA_quant <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Quantitative") %>% pluck(1,1,"data")})

## Save Excels
writexl::write_xlsx(MPA_vul,here::here("R/MPA_INDEX_result.xlsx"))
writexl::write_xlsx(MPA_cov,here::here("R/MPA_INDEX_coverage.xlsx"))
writexl::write_xlsx(MPA_cont,here::here("R/MPA_INDEX_contribution.xlsx"))
writexl::write_xlsx(MPA_year,here::here("R/MPA_INDEX_years.xlsx"))
writexl::write_xlsx(MPA_quant,here::here("R/MPA_INDEX_quantitative.xlsx"))

################################################################################
#### SP RESULTS ################################################################

## SP vulnerability
SP_vul <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SP INDEX",vars = "Value") %>% pluck(1,1,"data")})
## SP data coverage
SP_cov <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SP INDEX",vars = "Value") %>% pluck(1,1,"coverage")})
## SP indicator contribution
SP_cont <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SP INDEX",vars = "Value") %>% pluck(1,1,"contribution") %>% last})
## SP number of years
SP_year <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SP INDEX",vars = "Year") %>% pluck(1,1,"data")})
## SP quantitative indicators
SP_quant <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SP INDEX",vars = "Quantitative") %>% pluck(1,1,"data")})

## Save Excels
writexl::write_xlsx(SP_vul,here::here("R/SP_INDEX_result.xlsx"))
writexl::write_xlsx(SP_cov,here::here("R/SP_INDEX_coverage.xlsx"))
writexl::write_xlsx(SP_cont,here::here("R/SP_INDEX_contribution.xlsx"))
writexl::write_xlsx(SP_year,here::here("R/SP_INDEX_years.xlsx"))
writexl::write_xlsx(SP_quant,here::here("R/SP_INDEX_quantitative.xlsx"))

################################################################################
#### SST INDEX #################################################################

## SST vulnerability
SST_vul <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SST INDEX",vars = "Value") %>% pluck(1,1,"data")})
## SST data coverage
SST_cov <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SST INDEX",vars = "Value") %>% pluck(1,1,"coverage")})
## SST indicator contribution
SST_cont <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SST INDEX",vars = "Value") %>% pluck(1,1,"contribution") %>% last})
## SST years
SST_year <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SST INDEX",vars = "Year") %>% pluck(1,1,"data")})
## SST quantitative indicators
SST_quant <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SST INDEX",vars = "Quantitative") %>% pluck(1,1,"data")})

## Save Excels
writexl::write_xlsx(SST_vul,here::here("R/SST_INDEX_result.xlsx"))
writexl::write_xlsx(SST_cov,here::here("R/SST_INDEX_coverage.xlsx"))
writexl::write_xlsx(SST_cont,here::here("R/SST_INDEX_contribution.xlsx"))
writexl::write_xlsx(SST_year,here::here("R/SST_INDEX_years.xlsx"))
writexl::write_xlsx(SST_quant,here::here("R/SST_INDEX_quantitative.xlsx"))


################################################################################
#### MHW INDEX #################################################################

## MHW vulnerability
MHW_vul <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MHW INDEX",vars = "Value") %>% pluck(1,1,"data")})
## MHW data coverage
MHW_cov <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MHW INDEX",vars = "Value") %>% pluck(1,1,"coverage")})
## MHW contribution
MHW_cont <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MHW INDEX",vars = "Value") %>% pluck(1,1,"contribution") %>% last})
## MHW years
MHW_year <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MHW INDEX",vars = "Year") %>% pluck(1,1,"data")})
## MHW quantitative indicators
MHW_quant <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MHW INDEX",vars = "Quantitative") %>% pluck(1,1,"data")})

## Save Excels
writexl::write_xlsx(MHW_vul,here::here("R/MHW_INDEX_result.xlsx"))
writexl::write_xlsx(MHW_cov,here::here("R/MHW_INDEX_coverage.xlsx"))
writexl::write_xlsx(MHW_cont,here::here("R/MHW_INDEX_contribution.xlsx"))
writexl::write_xlsx(MHW_year,here::here("R/MHW_INDEX_years.xlsx"))
writexl::write_xlsx(MHW_quant,here::here("R/MHW_INDEX_quantitative.xlsx"))

################################################################################
#### HABITAT RESULTS ###########################################################

## HABITAT vulnerability
HB_vul <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "HB INDEX",vars = "Value") %>% pluck(1,1,"data")})
## HABITAT data coverage
HB_cov <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "HB INDEX",vars = "Value") %>% pluck(1,1,"coverage")})
## HABITAT contribution
HB_cont <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "HB INDEX",vars = "Value") %>% pluck(1,1,"contribution") %>% last})
## HABITAT years
HB_year <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "HB INDEX",vars = "Year") %>% pluck(1,1,"data")})
## HABITAT quantitative indicators
HB_quant <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "HB INDEX",vars = "Quantitative") %>% pluck(1,1,"data")})

## Save Excels
writexl::write_xlsx(HB_vul,here::here("R/HB_INDEX_result.xlsx"))
writexl::write_xlsx(HB_cov,here::here("R/HB_INDEX_coverage.xlsx"))
writexl::write_xlsx(HB_cont,here::here("R/HB_INDEX_contribution.xlsx"))
writexl::write_xlsx(HB_year,here::here("R/HB_INDEX_years.xlsx"))
writexl::write_xlsx(HB_quant,here::here("R/HB_INDEX_quantitative.xlsx"))

################################################################################
#### RECREATIONAL ACTIVITIES INDEX #############################################

## RA vulnerability
RA_vul <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "RA INDEX",vars = "Value") %>% pluck(1,1,"data")})
## RA data coverage
RA_cov <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "RA INDEX",vars = "Value") %>% pluck(1,1,"coverage")})
## RA data contribution
RA_cont <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "RA INDEX",vars = "Value") %>% pluck(1,1,"contribution") %>% last})
## RA years
RA_year <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "RA INDEX",vars = "Year") %>% pluck(1,1,"data")})
## RA quantitative indicators
RA_quant <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "RA INDEX",vars = "Quantitative") %>% pluck(1,1,"data")})

## Save Excels
writexl::write_xlsx(RA_vul,here::here("R/RA_INDEX_result.xlsx"))
writexl::write_xlsx(RA_cov,here::here("R/RA_INDEX_coverage.xlsx"))
writexl::write_xlsx(RA_cont,here::here("R/RA_INDEX_contribution.xlsx"))
writexl::write_xlsx(RA_year,here::here("R/RA_INDEX_years.xlsx"))
writexl::write_xlsx(RA_quant,here::here("R/RA_INDEX_quantitative.xlsx"))


################################################################################
#### PROFESIONAL FISHERS INDEX #################################################

## PF vulnerability
PF_vul <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "PF INDEX",vars = "Value") %>% pluck(1,1,"data")})
## PF data coverage
PF_cov <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "PF INDEX",vars = "Value") %>% pluck(1,1,"coverage")})
## PF contribution
PF_cont <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "PF INDEX",vars = "Value") %>% pluck(1,1,"contribution") %>% last})
## PF years
PF_year <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "PF INDEX",vars = "Year") %>% pluck(1,1,"data")})
## PF quantitative indicators
PF_quant <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "PF INDEX",vars = "Quantitative") %>% pluck(1,1,"data")})

## Save Excels
writexl::write_xlsx(PF_vul,here::here("R/PF_INDEX_result.xlsx"))
writexl::write_xlsx(PF_cov,here::here("R/PF_INDEX_coverage.xlsx"))
writexl::write_xlsx(PF_cont,here::here("R/PF_INDEX_contribution.xlsx"))
writexl::write_xlsx(PF_year,here::here("R/PF_INDEX_years.xlsx"))
writexl::write_xlsx(PF_quant,here::here("R/PF_INDEX_quantitative.xlsx"))

################################################################################
#### LOGS ######################################################################
## CAPTURE THE ERRORS AND SAVE IN TXT FILES

logs <- scenarios %>% compact %>% map2(names(.),~{
  c(paste0("\n\n##### ",.y," #####"), attr(.x,"log") %>% unique())
  
}) %>% unlist

write(logs,file=here::here("R/logs.txt"))

################################################################################
## TRANSFORM INDICATORS TO QUALITATIVE VALUES ##################################
## CONSIDERING DISTRIBUTION OF DATA ACROSS MPAS ################################
## Smit Vasquez Caballero
## Read quantitative values in /inputs and stores qualitative in /outputs
mpa <- substr (names(scenarios) [1], start=1, stop=nchar(names(scenarios) [1])-8)


source(here::here("R/IndexCalculation.R")) ## READS THE SMITS' CODE

# table_thresholds <- rbind (
#                    cbind ("Group"="Habitats", list1[[2]]),  ## Habitats
#                    cbind("Group"="Users", list2[[2]]), ## Users
#                    cbind("Group"="Species", list3[[2]]), ## Species
#                    cbind("Group"="SST", list4[[2]]), ## SST
#                    cbind("Group"="MHW", list5[[2]])) ## MHW

## Tables made by Francesca:
tableEV  <- read_excel (here::here("inputs/EV.xlsx"))
tableMPA <- read_excel (here::here("inputs/MPA.xlsx"))

## Calculate thresholds with Smits' code
ecological_thresholds <- QualitativeIndex (T=tableEV)
social_thresholds     <- QualitativeIndex (T=tableMPA)

## The values of the MPA under analysis:
table_mpa <- data.frame (
  "MPA" = mpa,
  "Scenario" = gsub("[.]", "", paste ("RCP", substr(names (MPA_vul), start=nchar(names (MPA_vul))-7, stop=nchar(names(MPA_vul))), sep="")),
  "EcologicalAdaptiveCapacity" = c(MPA_vul[[1]]$`ECOLOGICAL ADAPTIVE CAPACITY`,
                                   MPA_vul[[2]]$`ECOLOGICAL ADAPTIVE CAPACITY`,
                                   MPA_vul[[3]]$`ECOLOGICAL ADAPTIVE CAPACITY`,
                                   MPA_vul[[4]]$`ECOLOGICAL ADAPTIVE CAPACITY`,
                                   MPA_vul[[5]]$`ECOLOGICAL ADAPTIVE CAPACITY`,
                                   MPA_vul[[6]]$`ECOLOGICAL ADAPTIVE CAPACITY`),
  "EcologicalSensitivity" = c(MPA_vul[[1]]$`ECOLOGICAL SENSITIVITY`,
                              MPA_vul[[2]]$`ECOLOGICAL SENSITIVITY`,
                              MPA_vul[[3]]$`ECOLOGICAL SENSITIVITY`,
                              MPA_vul[[4]]$`ECOLOGICAL SENSITIVITY`,
                              MPA_vul[[5]]$`ECOLOGICAL SENSITIVITY`,
                              MPA_vul[[6]]$`ECOLOGICAL SENSITIVITY`),
  "EcologicalExposure" = c(MPA_vul[[1]]$`EXPOSURE`,
                           MPA_vul[[2]]$`EXPOSURE`,
                           MPA_vul[[3]]$`EXPOSURE`,
                           MPA_vul[[4]]$`EXPOSURE`,
                           MPA_vul[[5]]$`EXPOSURE`,
                           MPA_vul[[6]]$`EXPOSURE`),
  "SocialAdaptiveCapacity" = c(MPA_vul[[1]]$`SOCIAL ADAPTIVE CAPACITY`,
                               MPA_vul[[2]]$`SOCIAL ADAPTIVE CAPACITY`,
                               MPA_vul[[3]]$`SOCIAL ADAPTIVE CAPACITY`,
                               MPA_vul[[4]]$`SOCIAL ADAPTIVE CAPACITY`,
                               MPA_vul[[5]]$`SOCIAL ADAPTIVE CAPACITY`,
                               MPA_vul[[6]]$`SOCIAL ADAPTIVE CAPACITY`),
  "SocialSensitivity" = c(MPA_vul[[1]]$`SOCIAL SENSITIVITY`,
                               MPA_vul[[2]]$`SOCIAL SENSITIVITY`,
                               MPA_vul[[3]]$`SOCIAL SENSITIVITY`,
                               MPA_vul[[4]]$`SOCIAL SENSITIVITY`,
                               MPA_vul[[5]]$`SOCIAL SENSITIVITY`,
                               MPA_vul[[6]]$`SOCIAL SENSITIVITY`),
  "SocialExposure" = c(MPA_vul[[1]]$`ECOLOGICAL VULNERABILITY`,
                       MPA_vul[[2]]$`ECOLOGICAL VULNERABILITY`,
                       MPA_vul[[3]]$`ECOLOGICAL VULNERABILITY`,
                       MPA_vul[[4]]$`ECOLOGICAL VULNERABILITY`,
                       MPA_vul[[5]]$`ECOLOGICAL VULNERABILITY`,
                       MPA_vul[[6]]$`ECOLOGICAL VULNERABILITY`)
  )

names(MPA_vul[[1]])

## calculate QUALITATIVE VALUES FOR THE MPA
breaks.eco.adaptive    <- unlist(c(0, ecological_thresholds [[2]] [1, c("per20", "per40", "per60", "per80")], 1))
breaks.eco.sensitivity <- unlist(c(0, ecological_thresholds [[2]] [2, c("per20", "per40", "per60", "per80")], 1))
breaks.eco.exposure    <- unlist(c(0, ecological_thresholds [[2]] [3, c("per20", "per40", "per60", "per80")], 1))
breaks.social.adaptive    <- unlist(c(0, social_thresholds [[2]] [2, c("per20", "per40", "per60", "per80")], 1))
breaks.social.sensitivity <- unlist(c(0, social_thresholds [[2]] [3, c("per20", "per40", "per60", "per80")], 1))

table_mpa$EcologicalAdaptiveCapacityQ <- cut(table_mpa$EcologicalAdaptiveCapacity, breaks=breaks.eco.adaptive, labels=FALSE, include.lowest = TRUE)
table_mpa$EcologicalSensitivityQ      <- cut(table_mpa$EcologicalSensitivity,      breaks=breaks.eco.sensitivity, labels=FALSE, include.lowest = TRUE)
table_mpa$EcologicalExposureQ         <- cut(table_mpa$EcologicalExposure,         breaks=breaks.eco.exposure, labels=FALSE, include.lowest = TRUE)
table_mpa$SocialAdaptiveCapacityQ     <- cut(table_mpa$SocialAdaptiveCapacity, breaks=breaks.social.adaptive, labels=FALSE, include.lowest = TRUE)
table_mpa$SocialSensitivityQ          <- cut(table_mpa$SocialSensitivity,      breaks=breaks.social.sensitivity, labels=FALSE, include.lowest = TRUE)

number.to.category <- function (x) {
  x [x== 1] <- "Low"
  x [x== 2] <- "Moderate"
  x [x== 3] <- "High"
  x [x== 4] <- "Very high"
  x [x== 5] <- "Extreme"
  return(x)
}

table_mpa$EcologicalAdaptiveCapacityQ <- number.to.category (x=table_mpa$EcologicalAdaptiveCapacityQ)
table_mpa$EcologicalSensitivityQ <- number.to.category (x=table_mpa$EcologicalSensitivityQ)
table_mpa$EcologicalExposureQ <- number.to.category (x=table_mpa$EcologicalExposureQ)
table_mpa$SocialAdaptiveCapacityQ <- number.to.category (x=table_mpa$SocialAdaptiveCapacityQ)
table_mpa$SocialSensitivityQ <- number.to.category (x=table_mpa$SocialSensitivityQ)

## Calculate final index by combination of categories (use JuanMatrix from IndexCalculation.R)
JuanMatrix <- gsub("Intermediate", "Moderate", as.matrix(JuanMatrix))
JuanMatrix <- as.data.frame(JuanMatrix)

Index <- rep (NA, nrow(table_mpa))
for (i in 1:nrow(table_mpa))
{
  one   <- which (toupper(JuanMatrix$Q_AC) == toupper(table_mpa$EcologicalAdaptiveCapacityQ [i]))
  two   <- which (toupper(JuanMatrix$Q_Sensitivity) == toupper(table_mpa$EcologicalSensitivityQ [i]))
  three <- which (toupper(JuanMatrix$Q_Exposure) == toupper(table_mpa$EcologicalExposureQ [i]))
  Index [i] <- JuanMatrix$Index [intersect (one, intersect(two, three))]
}


table_mpa$EcologicalVulnerability <- Index

#########################################################################
table_mpa$SocialExposure <- c(MPA_vul [[1]]$`ECOLOGICAL VULNERABILITY`,
                              MPA_vul [[2]]$`ECOLOGICAL VULNERABILITY`,
                              MPA_vul [[3]]$`ECOLOGICAL VULNERABILITY`,
                              MPA_vul [[4]]$`ECOLOGICAL VULNERABILITY`,
                              MPA_vul [[5]]$`ECOLOGICAL VULNERABILITY`,
                              MPA_vul [[6]]$`ECOLOGICAL VULNERABILITY`)

table_mpa$SocialExposureQ <- table_mpa$EcologicalVulnerability
#########################################################################


## sOCIAL Vulnerability
Index2 <- rep (NA, nrow(table_mpa))
for (i in 1:nrow(table_mpa))
{
  four   <- which (toupper(JuanMatrix$Q_AC) == toupper(table_mpa$SocialAdaptiveCapacityQ [i]))
  five   <- which (toupper(JuanMatrix$Q_Sensitivity) == toupper(table_mpa$SocialSensitivityQ [i]))
  six    <- which (toupper(JuanMatrix$Q_Exposure) == toupper(table_mpa$SocialExposureQ [i]))
  Index2 [i] <- JuanMatrix$Index [intersect (four, intersect(five, six))]
}

table_mpa$SocialVulnerability <- Index2


################################################################################
#### FIGURES ###################################################################
## Use table_mpa for the donuts


## Juan is working on this

## Load donut code and save the figures in figures/
source(here::here("R/Figures.R")) ## Donut



################################################################################
#### FIGURES ##################################################################

## Objects for the histograms
sp_thresholds <- read_excel (here::here("inputs/Species_thresholds.xlsx"), sheet=2)

percentiles.AC.species <- unlist (c(0, sp_thresholds [1, c("Per_20", "Per_40", "Per_60", "Per_80")], 1))
percentiles.S.species <- unlist (c(0, sp_thresholds [2, c("Per_20", "Per_40", "Per_60", "Per_80")], 1))
percentiles.E.species <- unlist (c(0, sp_thresholds [3, c("Per_20", "Per_40", "Per_60", "Per_80")], 1))

make.histogram.SP <- function (Scenario = "2.6_2050")
{
  sheet <- grep(Scenario, names(SP_vul) )
  table_sp <- data.frame ("SP" = SP_vul [[sheet]]$SP, 
              #names(SP_vul[[sheet]]),
              "Scenario" = "RCP26_2050",
              "Exposure" = SP_vul [[sheet]]$EXPOSURE, 
              "Adaptive capacity" = SP_vul [[sheet]]$`ECOLOGICAL ADAPTIVE CAPACITY`, 
              "Sensitivity" = SP_vul [[sheet]]$`ECOLOGICAL SENSITIVITY`,
              "Index"= SP_vul [[sheet]]$Vulnerability)
  
  table_sp$ExposureQ    <- cut (table_sp$Exposure, breaks=percentiles.E.species, labels=FALSE, include.lowest = TRUE)
  table_sp$SensitivityQ <- cut (table_sp$Sensitivity, breaks=percentiles.S.species, labels=FALSE, include.lowest = TRUE)
  table_sp$AdaptiveQ <- cut (table_sp$Adaptive.capacity, breaks=percentiles.AC.species, labels=FALSE, include.lowest = TRUE)
  
  table_sp$ExposureQ    <- number.to.category (table_sp$ExposureQ)
  table_sp$SensitivityQ <- number.to.category (table_sp$SensitivityQ)
  table_sp$AdaptiveQ    <- number.to.category (table_sp$AdaptiveQ)
  
  table_sp$IndexQ <- NA
  for (i in 1:nrow(table_sp))
  {
    one   <- which (toupper(JuanMatrix$Q_AC) == toupper(table_sp$AdaptiveQ [i]))
    two   <- which (toupper(JuanMatrix$Q_Sensitivity) == toupper(table_sp$SensitivityQ [i]))
    three <- which (toupper(JuanMatrix$Q_Exposure) == toupper(table_sp$ExposureQ [i]))
    table_sp$IndexQ [i] <- JuanMatrix$Index [intersect (one, intersect(two, three))]
  }
  return(table_sp)
}

getwd()

sp.scenario1 <- make.histogram.SP(Scenario = "2.6_2050") ; write.csv(sp.scenario1, "sp_scenario1.csv")
sp.scenario2 <- make.histogram.SP(Scenario = "4.5_2050")
sp.scenario3 <- make.histogram.SP(Scenario = "8.5_2050")
sp.scenario4 <- make.histogram.SP(Scenario = "2.6_2100")
sp.scenario5 <- make.histogram.SP(Scenario = "4.5_2100")
sp.scenario6 <- make.histogram.SP(Scenario = "8.5_2100")





#####################################################################################################################
#####################################################################################################################

# Save the results in excel files

#### MPA INDEX ####

MPA_vul <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"data")})
MPA_vul_IND1<- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"data", "FPRES")})
indicator.values.SSTthreat <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "SST threat")})
indicator.values.MHWthreat <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "MHW threat")})
indicator.values.watercond <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "water cond")})
indicator.values.humanpressure <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "human pressure")})
indicator.values.habitatintegritythreats <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "habitat integrity threats")})
indicator.values.speciesintegritythreats <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "species integrity threats")})
indicator.values.hab.redundancy <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "hab.redundancy")})
indicator.values.hab.Recoverypotential <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "hab. Recovery potential")})
indicator.values.sp.Recoverypotential <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "sp. Recovery potential")})
indicator.values.effectiveness <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "effectiveness")})
indicator.values.conservationefforts <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "conservation efforts")})
indicator.values.adaptivemanagement <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "adaptive management")})
indicator.values.Professionalfishingdependency <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "Professional fishing dependency")})
indicator.values.Professionalfishingeffort <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "Professional fishing effort")})
indicator.values.Professionalfishinglocaldependency <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "Professional fishing local dependency")})
indicator.values.Recreationalactivitiesemployment <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "Recreational activities employment")})
indicator.values.Recreationalactivitiesecosystem <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "Recreational activities ecosystem")})
indicator.values.Recreationalactivitiesfacilities <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "Recreational activities facilities")})
indicator.values.Flexibility <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "Flexibility")})
indicator.values.SocialOrganization <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "Social Organization")})
indicator.values.learning <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "learning")})
indicator.values.assets <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "assets")})
indicator.values.agencyandsocioculturalaspects <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution", "agency and socio-cultural aspects")})

## OBTAIN REAL VALUE OF THE INDICATORS (THEY ARE WEIGHTED)
values.SSTthreat <- (as.data.frame (indicator.values.SSTthreat) [-1]) * 1
values.MHWthreat <- as.data.frame (indicator.values.MHWthreat) [-1] * c(1)
values.watercond <- as.data.frame (indicator.values.watercond) [-1] * c(3,3,3)
values.humanpressure <- as.data.frame (indicator.values.humanpressure) [-1] * c(8,8,8,8,8,8,8,8)
values.habitatintegritythreats <- as.data.frame (indicator.values.habitatintegritythreats) [-1] * c(8,8,8,8,8,8,8,8)
values.speciesintegritythreats <- as.data.frame (indicator.values.speciesintegritythreats) [-1] * c(8,8,8,8,8,8,8,8)
values.hab.redundancy <- (1- (as.data.frame (indicator.values.hab.redundancy) [-1] * c(3,3,3))) * 1/3 * 0.148 * 1/3 * 1/3
values.hab.Recoverypotential <- (1- (as.data.frame (indicator.values.hab.Recoverypotential) [-1] * c(5,5,5,5,5))) * 1/5 * 0.19 * 1/3 *1/3
values.sp.Recoverypotential <- (1- (as.data.frame (indicator.values.sp.Recoverypotential) [-1] * c(5,5,5,5,5))) * 1/5 *  0.188  * 1/3 * 1/3
values.effectiveness <- (1- (as.data.frame (indicator.values.effectiveness) [-1] * c(6,6,6,6,6,6))) * 1/6 * 0.146 * 1/3 *1/3
values.conservationefforts <- (1- (as.data.frame (indicator.values.conservationefforts) [-1] * c(6,6,6,6,6,6))) * 1/6 * 0.17 * 1/3 *1/3
values.adaptivemanagement <- (1- (as.data.frame (indicator.values.adaptivemanagement) [-1] * c(2,2))) * 1/2 *0.16 *1/3 *1/3
values.Professionalfishingdependency <- as.data.frame (indicator.values.Professionalfishingdependency) [-1] * c(2,2)
values.Professionalfishingeffort <- as.data.frame (indicator.values.Professionalfishingeffort) [-1] * c(2,2)
values.Professionalfishinglocaldependency <- as.data.frame (indicator.values.Professionalfishinglocaldependency) [-1] * c(3,3,3)
#values.Recreationalactivitiesemployment <- as.data.frame (indicator.values.Recreationalactivitiesemployment) [-1] * c(3,3,3)
values.Recreationalactivitiesecosystem <- as.data.frame (indicator.values.Recreationalactivitiesecosystem) [-1] * c(3,3,3)
values.Recreationalactivitiesfacilities <- as.data.frame (indicator.values.Recreationalactivitiesfacilities) [-1] * c(2,2)
values.Flexibility <- (1 - (as.data.frame (indicator.values.Flexibility) [-1] * c(4,4,4,4,3,3,3) * c(2,2,2,2,2,2,2))) * c( 1/8, 1/8, 1/8, 1/8, 1/6, 1/6, 1/6) * 0.234 * 1/3 * 1/3 
values.SocialOrganization <- (1 - (as.data.frame (indicator.values.SocialOrganization) [-1] * c(8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8) * c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2))) * c( 1/16, 1/16, 1/16, 1/16, 1/16, 1/16, 1/16, 1/16, 1/16, 1/16, 1/16, 1/16, 1/16, 1/16, 1/16, 1/16) * 0.247 * 1/3 * 1/3 
values.learning <- (1- (as.data.frame (indicator.values.learning) [-1] * c(1,1) * c(2,2))) * c( 1/2, 1/2) * 0.214 * 1/3 * 1/3
values.assets <- (1- (as.data.frame (indicator.values.assets) [-1] * c(1,1) * c(2,2))) * c( 1/2, 1/2) * 0.155 * 1/3 * 1/3
values.agencyandsocioculturalaspects <- (1- (as.data.frame (indicator.values.agencyandsocioculturalaspects) [-1] * c(4,4,4,4,3,3,3) * c(2,2,2,2,2,2,2))) * c( 1/8, 1/8, 1/8, 1/8, 1/6, 1/6, 1/6) * 0.15 * 1/3 * 1/3



## EXPORT INDICATORS VALUES IN EXCELL

indicator.values <- cbind(
  as.data.frame(values.hab.redundancy),
  as.data.frame(values.hab.Recoverypotential),
  as.data.frame(values.sp.Recoverypotential),
  as.data.frame(values.effectiveness),
  as.data.frame(values.conservationefforts),
  as.data.frame(values.adaptivemanagement),
  as.data.frame(values.Flexibility),
  as.data.frame(values.SocialOrganization),
  as.data.frame(values.learning),
  as.data.frame(values.assets),
  as.data.frame(values.agencyandsocioculturalaspects))



writexl::write_xlsx(indicator.values,here::here("results/indicator_values.xlsx"))

MPA_cov <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"coverage")})

MPA_cont <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"contribution") %>% last})
MPA_cont_IND1 <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Value") %>% pluck(1,1,"data", "water cond")})

writexl::write_xlsx(MPA_vul,here::here("R/MPA_INDEX_result.xlsx"))
writexl::write_xlsx(MPA_cov,here::here("R/MPA_INDEX_coverage.xlsx"))
writexl::write_xlsx(MPA_cont,here::here("R/MPA_INDEX_contribution.xlsx"))

MPA_year <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Year") %>% pluck(1,1,"data")})
MPA_quant <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MPA INDEX",vars = "Quantitative") %>% pluck(1,1,"data")})

writexl::write_xlsx(MPA_year,here::here("R/MPA_INDEX_years.xlsx"))
writexl::write_xlsx(MPA_quant,here::here("R/MPA_INDEX_quantitative.xlsx"))
#### SP INDEX ####

SP_vul <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SP INDEX",vars = "Value") %>% pluck(1,1,"data")})

SP_cov <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SP INDEX",vars = "Value") %>% pluck(1,1,"coverage")})

SP_cont <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SP INDEX",vars = "Value") %>% pluck(1,1,"contribution") %>% last})

writexl::write_xlsx(SP_vul,here::here("R/SP_INDEX_result.xlsx"))
writexl::write_xlsx(SP_cov,here::here("R/SP_INDEX_coverage.xlsx"))
writexl::write_xlsx(SP_cont,here::here("R/SP_INDEX_contribution.xlsx"))


SP_year <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SP INDEX",vars = "Year") %>% pluck(1,1,"data")})
SP_quant <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SP INDEX",vars = "Quantitative") %>% pluck(1,1,"data")})

writexl::write_xlsx(SP_year,here::here("R/SP_INDEX_years.xlsx"))
writexl::write_xlsx(SP_quant,here::here("R/SP_INDEX_quantitative.xlsx"))


#### SST INDEX ####

SST_vul <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SST INDEX",vars = "Value") %>% pluck(1,1,"data")})

SST_cov <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SST INDEX",vars = "Value") %>% pluck(1,1,"coverage")})

SST_cont <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SST INDEX",vars = "Value") %>% pluck(1,1,"contribution") %>% last})

writexl::write_xlsx(SST_vul,here::here("R/SST_INDEX_result.xlsx"))
writexl::write_xlsx(SST_cov,here::here("R/SST_INDEX_coverage.xlsx"))
writexl::write_xlsx(SST_cont,here::here("R/SST_INDEX_contribution.xlsx"))

SST_year <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SST INDEX",vars = "Year") %>% pluck(1,1,"data")})
SST_quant <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "SST INDEX",vars = "Quantitative") %>% pluck(1,1,"data")})

writexl::write_xlsx(SST_year,here::here("R/SST_INDEX_years.xlsx"))
writexl::write_xlsx(SST_quant,here::here("R/SST_INDEX_quantitative.xlsx"))


#### MHW INDEX ####

MHW_vul <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MHW INDEX",vars = "Value") %>% pluck(1,1,"data")})

MHW_cov <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MHW INDEX",vars = "Value") %>% pluck(1,1,"coverage")})

MHW_cont <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MHW INDEX",vars = "Value") %>% pluck(1,1,"contribution") %>% last})

writexl::write_xlsx(MHW_vul,here::here("R/MHW_INDEX_result.xlsx"))
writexl::write_xlsx(MHW_cov,here::here("R/MHW_INDEX_coverage.xlsx"))
writexl::write_xlsx(MHW_cont,here::here("R/MHW_INDEX_contribution.xlsx"))


MHW_year <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MHW INDEX",vars = "Year") %>% pluck(1,1,"data")})
MHW_quant <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "MHW INDEX",vars = "Quantitative") %>% pluck(1,1,"data")})

writexl::write_xlsx(MHW_year,here::here("R/MHW_INDEX_years.xlsx"))
writexl::write_xlsx(MHW_quant,here::here("R/MHW_INDEX_quantitative.xlsx"))


#### HB INDEX ####

HB_vul <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "HB INDEX",vars = "Value") %>% pluck(1,1,"data")})

HB_cov <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "HB INDEX",vars = "Value") %>% pluck(1,1,"coverage")})

HB_cont <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "HB INDEX",vars = "Value") %>% pluck(1,1,"contribution") %>% last})

writexl::write_xlsx(HB_vul,here::here("R/HB_INDEX_result.xlsx"))
writexl::write_xlsx(HB_cov,here::here("R/HB_INDEX_coverage.xlsx"))
writexl::write_xlsx(HB_cont,here::here("R/HB_INDEX_contribution.xlsx"))

HB_year <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "HB INDEX",vars = "Year") %>% pluck(1,1,"data")})
HB_quant <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "HB INDEX",vars = "Quantitative") %>% pluck(1,1,"data")})

writexl::write_xlsx(HB_year,here::here("R/HB_INDEX_years.xlsx"))
writexl::write_xlsx(HB_quant,here::here("R/HB_INDEX_quantitative.xlsx"))


#### RA INDEX ####

RA_vul <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "RA INDEX",vars = "Value") %>% pluck(1,1,"data")})

RA_cov <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "RA INDEX",vars = "Value") %>% pluck(1,1,"coverage")})

RA_cont <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "RA INDEX",vars = "Value") %>% pluck(1,1,"contribution") %>% last})

writexl::write_xlsx(RA_vul,here::here("R/RA_INDEX_result.xlsx"))
writexl::write_xlsx(RA_cov,here::here("R/RA_INDEX_coverage.xlsx"))
writexl::write_xlsx(RA_cont,here::here("R/RA_INDEX_contribution.xlsx"))

RA_year <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "RA INDEX",vars = "Year") %>% pluck(1,1,"data")})
RA_quant <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "RA INDEX",vars = "Quantitative") %>% pluck(1,1,"data")})

writexl::write_xlsx(RA_year,here::here("R/RA_INDEX_years.xlsx"))
writexl::write_xlsx(RA_quant,here::here("R/RA_INDEX_quantitative.xlsx"))


#### PF INDEX ####

PF_vul <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "PF INDEX",vars = "Value") %>% pluck(1,1,"data")})

PF_cov <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "PF INDEX",vars = "Value") %>% pluck(1,1,"coverage")})

PF_cont <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "PF INDEX",vars = "Value") %>% pluck(1,1,"contribution") %>% last})

writexl::write_xlsx(PF_vul,here::here("R/PF_INDEX_result.xlsx"))
writexl::write_xlsx(PF_cov,here::here("R/PF_INDEX_coverage.xlsx"))
writexl::write_xlsx(PF_cont,here::here("R/PF_INDEX_contribution.xlsx"))

PF_year <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "PF INDEX",vars = "Year") %>% pluck(1,1,"data")})
PF_quant <- scenarios %>% compact %>% map(~{ .x %>% get_results(indexes = "PF INDEX",vars = "Quantitative") %>% pluck(1,1,"data")})

writexl::write_xlsx(PF_year,here::here("R/PF_INDEX_years.xlsx"))
writexl::write_xlsx(PF_quant,here::here("R/PF_INDEX_quantitative.xlsx"))

#### Logs ####

logs <- scenarios %>% compact %>% map2(names(.),~{
  c(paste0("\n\n##### ",.y," #####"), attr(.x,"log") %>% unique())
  
}) %>% unlist

write(logs,file=here::here("R/logs.txt"))

