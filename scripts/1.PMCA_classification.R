# PMCA_algorithm---------------------------------------------------------
#
# PURPOSE: Group the Dx_Code codes at first encounter into Medical Complexity groupings.
# https://github.com/kpwhri/pmca/blob/main/get_PMCA.R 
# 
# LAST UPDATED: July 5, 2024, by Emily Yeo

#clear environment
rm(list = ls())
gc()

#Load Packages
pacman::p_load(knitr, tidyverse, magrittr, lme4, lmerTest, GGally, corrplot, 
               Hmisc, kableExtra, dplyr, plyr, janitor, lubridate, survminer, 
               ggplot2, here, readr, tableone, officer, flextable,finalfit,
               purrr, stringr, lme4, corrplot, pscl, stargazer, MASS, lmerTest,
               pccc, icd, RCurl, reshape2, htmlTable, readxl)

#This section will need to be customized for different environments or data sets 

# SET WORKSPACE AND PATHS #### (edit for personal computers)
basedir <- "/Users/emily/projects/research/ADOR/CHLA_clean/" # home PC
server_in <- "/Volumes/IPHY/ADORLab/Lab\ Projects/CHLAneuro/CHLAclean/"
in_path <- paste(basedir,"clean_inputs/",sep="")
check_path <- paste(basedir, "dx_codes/", sep ="")
out_path <- paste(basedir,"outputs/",sep="")

# First Encounters ####
setwd(server_in)
T_3 <- read_csv("clean_inputs/table3_first_nov24.csv", 
                na = c("NA", "don't know", "","    .",
                       "Missing: Not provided",999.00))

# remove the first column that just has row numbers
T_3 <- T_3[, -1]

# I want to include DOB in T_3. So I am merging that from T_4 below:
# Last encounters
T_4 <- read_csv("clean_inputs/table4_last_nov24.csv", 
                na = c("NA", "don't know", "","    .",
                       "Missing: Not provided",999.00))
T_4 <- T_4[, -1] # remove empty first row
Rid_DOB <- T_4[, c("Research_ID", "DOB")]
unique_Rid_DOB <- unique(Rid_DOB)
merged_T3 <- merge(T_3, unique_Rid_DOB, by = "Research_ID", all.x = TRUE)
rm(Rid_DOB)
rm(unique_Rid_DOB)
rm(T_3)

# Get ICD9 and ICD10 codes with their PMCA systems assignment ####
#setwd(basedir) # local run
setwd(server_in)
icd9 <- read_excel("dx.codes/icd9_pmca_systems.xlsx", sheet = "Sheet1")
icd10 <- read_excel("dx.codes/icd10_pmca_systems.xlsx", sheet = "Sheet1")

# Create list of system codes ####
systems <- data.frame("code" = "cardiac","label" = "cardiac","id"="")
systems <- rbind(systems,data.frame("code" = "cranio","label" = "craniofacial","id"=""))
systems <- rbind(systems,data.frame("label"="dermatological","code"="derm","id"=""))
systems <- rbind(systems,data.frame("label"="endocrinological","code"="endo","id"=""))
systems <- rbind(systems,data.frame("label"="gastrointestinal","code"="gastro","id"=""))
systems <- rbind(systems,data.frame("label"="genetic","code"="genetic","id"=""))
systems <- rbind(systems,data.frame("label"="genitourinary","code"="genito","id"=""))
systems <- rbind(systems,data.frame("label"="hematological","code"="hemato","id"=""))
systems <- rbind(systems,data.frame("label"="immunological","code"="immuno","id"=""))
systems <- rbind(systems,data.frame("label"="malignancy","code"="malign","id"=""))
systems <- rbind(systems,data.frame("label"="mental health","code"="mh","id"=""))
systems <- rbind(systems,data.frame("label"="metabolic","code"="metab","id"=""))
systems <- rbind(systems,data.frame("label"="musculoskeletal","code"="musculo","id"=""))
systems <- rbind(systems,data.frame("label"="neurological","code"="neuro","id"=""))
systems <- rbind(systems,data.frame("label"="pulmonary-respiratory","code"="pulresp","id"=""))
systems <- rbind(systems,data.frame("label"="renal","code"="renal","id"=""))
systems <- rbind(systems,data.frame("label"="ophthalmological","code"="opthal","id"=""))
systems <- rbind(systems,data.frame("label"="otologic","code"="otol","id"=""))
systems <- rbind(systems,data.frame("label"="otolaryngological","code"="otolar","id"=""))
systems <- rbind(systems,data.frame("label"="progressive","code"="progressive","id"=""))

systems <- as.data.frame(systems, stringsAsFactors = FALSE)

# Process the input data ####
df <- merged_T3 

# if DIAG_SEQ_KEY = -1, no diagnosis entered, drop it
df <- subset(df, df$Dx_Display != "-1") # drops 1

# get unique Encounter_ID and Dx_Code set
Enc_Dx <- subset(df,select = c("Encounter_ID","Dx_Code"))
Enc_Dx <- unique(Enc_Dx) # drops from 8919 to 8757 (diff of 162)

# Start processing Dx_Code codes determine if ICD 9 or ICD 10 ####

# remove the dot in the codes
Enc_Dx$Dx_Code <- gsub("\\.","",Enc_Dx$Dx_Code)
# get first chacter
Enc_Dx$first <- substring(Enc_Dx$Dx_Code,1,1)
# if first is A-Z excluding E and V codes -- definitely ICD 10
icd10list <- c('A', 'B', 'C', 'D', 'F', 'G' ,'H', 'I', 'J', 'K', 'L' ,'M', 'N',
               'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'W' ,'X' ,'Y', 'Z')

Enc_Dx$icd10 <- ifelse(Enc_Dx$first %in% icd10list,1,0)

# if first is 0-9 -- definitely ICD 9
icd9list <- c(0:9)
Enc_Dx$icd9 <- ifelse(Enc_Dx$first %in% icd9list,1,0)

# get if E
Enc_Dx$icde <- ifelse(Enc_Dx$first == 'E',1,0)

# unique datset of just Research_ID and Encounter_ID
R_Enc_ID  <- subset(df, select = c("Research_ID",
                                   'Encounter_ID')) %>% unique()

# Count the code types in each encounter
Enc_ID_icd<- Enc_Dx %>% 
  dplyr::group_by(Encounter_ID) %>% 
  dplyr::summarize(icd10ct = sum(icd10),
                   icd9ct = sum(icd9),
                   icdect = sum(icde)) # i548 unique encounters

# if at least one ICD-10 code then mark entire claim;
Enc_ID_icd$icd10codeset <- 0
Enc_ID_icd$icd10codeset <- ifelse(Enc_ID_icd$icd10ct > 0,
                                  1,Enc_ID_icd$icd10codeset)
# if just one code and it is an E-code;
Enc_ID_icd$icd10codeset <- ifelse(Enc_ID_icd$icd10ct == 0 
                                  & Enc_ID_icd$icd9ct == 0 
                                  & Enc_ID_icd$icdect == 1,
                                  1,Enc_ID_icd$icd10codeset)
# all other cases leave as 0

# merge the icd10codeset variable back into the Enc_Dx recordset
Enc_Dx <- subset(Enc_Dx,select = c("Encounter_ID","Dx_Code"))
Enc_Dx <- merge(x=Enc_Dx,y=Enc_ID_icd,by="Encounter_ID",all=TRUE)
Enc_Dx <- subset(Enc_Dx,select= c("Encounter_ID","Dx_Code","icd10codeset"))
Enc_Dx10 <- subset(Enc_Dx,Enc_Dx$icd10codeset == 1)
Enc_Dx9 <- subset(Enc_Dx,Enc_Dx$icd10codeset!= 1)

# Looping the icd codes ####
# as a function
getSystem <- function(icd_base,sys_code,Encounter_ID_base){
  sys_loop <- subset(icd_base,icd_base$system == sys_code)
  foundit <- data.frame("Encounter_ID"=00000,"Dx_Code"="Dx_Code")
  if (nrow(sys_loop) < 1){return(foundit)}
  for (row in 1:nrow(sys_loop)){
    start_code <- tolower(sys_loop$start_code[row])
    no_chars <- sys_loop$no_chars[row]
    tmp <- Encounter_ID_base
    tmp$test <- tolower(substring(tmp$Dx_Code,1,no_chars))
    tmp$found <- ifelse(tmp$test == start_code,1,0)
    tmp2 <- subset(tmp,tmp$found == 1,select=c("Encounter_ID","Dx_Code"))
    foundit <- rbind(foundit,tmp2)
  }
  return(foundit)
}

### ICD 10 claims
matches_10 <- Enc_Dx10 # 8401 icd 10 codes 

### Loop systems
for (row in 1:nrow(systems)){
  mysystem <- systems$code[row]
  print(mysystem)
  matches <- getSystem(icd10,mysystem,Enc_Dx10)
  matches$result <- 1
  colnames(matches)[colnames(matches)== "result"] <- as.character(mysystem)
  colnames((matches))
  matches_10 <- merge(x=matches_10,y=matches,by=c("Encounter_ID","Dx_Code"),
                      all.x=TRUE)
}

# ICD 9 claims
matches_9 <- Enc_Dx9 # 356 

# loop through systems
for (row in 1:nrow(systems)){
  mysystem <- systems$code[row]
  print(mysystem)
  matches <- getSystem(icd9,mysystem,Enc_Dx9)
  matches$result <- 1
  colnames(matches)[colnames(matches)== "result"] <- as.character(mysystem)
  colnames((matches))
  matches_9 <- merge(x=matches_9,y=matches,
                     by=c("Encounter_ID","Dx_Code"),all.x=TRUE)
}

# combine 9 and 10 ICD codes ####
Enc_ID_icd_final <- rbind(matches_9,matches_10) # 9042 

# all NA to 0 
Enc_ID_icd_final[is.na(Enc_ID_icd_final)] <- 0

# compress each claim into one using group_by 
Encounter_ID_sum <- Enc_ID_icd_final %>% 
  dplyr::group_by(Encounter_ID) %>% 
  dplyr::summarize(
    cardiac_t= sum(cardiac),
    cranio_t= sum(cranio),
    derm_t= sum(derm),
    endo_t= sum(endo),
    gastro_t= sum(gastro),
    genetic_t= sum(genetic),
    genito_t= sum(genito),
    hemato_t= sum(hemato),
    immuno_t= sum(immuno),
    malign_t= sum(malign),
    mh_t= sum(mh),
    metab_t= sum(metab),
    musculo_t= sum(musculo),
    neuro_t= sum(neuro),
    pulresp_t= sum(pulresp),
    renal_t= sum(renal),
    opthal_t= sum(opthal),
    otol_t= sum(otol),
    otolar_t= sum(otolar),
    progressive_t= sum(progressive))

Encounter_ID_sum$cardiac_yn <- ifelse(Encounter_ID_sum$cardiac_t >0,1,0)
Encounter_ID_sum$cranio_yn <- ifelse(Encounter_ID_sum$cranio_t >0,1,0)
Encounter_ID_sum$derm_yn <- ifelse(Encounter_ID_sum$derm_t >0,1,0)
Encounter_ID_sum$endo_yn <- ifelse(Encounter_ID_sum$endo_t >0,1,0)
Encounter_ID_sum$gastro_yn <- ifelse(Encounter_ID_sum$gastro_t >0,1,0)
Encounter_ID_sum$genetic_yn <- ifelse(Encounter_ID_sum$genetic_t >0,1,0)
Encounter_ID_sum$genito_yn <- ifelse(Encounter_ID_sum$genito_t >0,1,0)
Encounter_ID_sum$hemato_yn <- ifelse(Encounter_ID_sum$hemato_t >0,1,0)
Encounter_ID_sum$immuno_yn <- ifelse(Encounter_ID_sum$immuno_t >0,1,0)
Encounter_ID_sum$malign_yn <- ifelse(Encounter_ID_sum$malign_t >0,1,0)
Encounter_ID_sum$mh_yn <- ifelse(Encounter_ID_sum$mh_t >0,1,0)
Encounter_ID_sum$metab_yn <- ifelse(Encounter_ID_sum$metab_t >0,1,0)
Encounter_ID_sum$musculo_yn <- ifelse(Encounter_ID_sum$musculo_t >0,1,0)
Encounter_ID_sum$neuro_yn <- ifelse(Encounter_ID_sum$neuro_t >0,1,0)
Encounter_ID_sum$pulresp_yn <- ifelse(Encounter_ID_sum$pulresp_t >0,1,0)
Encounter_ID_sum$renal_yn <- ifelse(Encounter_ID_sum$renal_t >0,1,0)
Encounter_ID_sum$opthal_yn <- ifelse(Encounter_ID_sum$opthal_t >0,1,0)
Encounter_ID_sum$otol_yn <- ifelse(Encounter_ID_sum$otol_t >0,1,0)
Encounter_ID_sum$otolar_yn <- ifelse(Encounter_ID_sum$otolar_t >0,1,0)
Encounter_ID_sum$progressive_yn <- ifelse(Encounter_ID_sum$progressive_t >0,1,0)

# first: join up Encounter_ID with Research_ID
R_Enc_ID_sum <- merge(x=R_Enc_ID,
                      y=Enc_ID_icd_final,
                      by="Encounter_ID",all.x=TRUE) 
# (CHECK THIS) ####
#  Roll up to one record per person, with single flag for each body system, sum 
# across claims per body system, presence of a progressive condition / malignancy. 
#  Calculate final condition determinations.  

# get sum across all claims ####
Research_ID_sum <- R_Enc_ID_sum %>% 
  dplyr::group_by(Research_ID) %>% 
  dplyr::summarize(
    cardiac_claims= sum(cardiac),
    cranio_claims= sum(cranio),
    derm_claims= sum(derm),
    endo_claims= sum(endo),
    gastro_claims= sum(gastro),
    genetic_claims= sum(genetic),
    genito_claims= sum(genito),
    hemato_claims= sum(hemato),
    immuno_claims= sum(immuno),
    malign_claims= sum(malign),
    mh_claims= sum(mh),
    metab_claims= sum(metab),
    musculo_claims= sum(musculo),
    neuro_claims= sum(neuro),
    pulresp_claims= sum(pulresp),
    renal_claims= sum(renal),
    opthal_claims= sum(opthal),
    otol_claims= sum(otol),
    otolar_claims= sum(otolar),
    progressive_claims= sum(progressive))

# indicator yes or no (1 or 0)
Research_ID_sum$cardiac_any <- ifelse(Research_ID_sum$cardiac_claims >0,1,0)
Research_ID_sum$cranio_any <- ifelse(Research_ID_sum$cranio_claims >0,1,0)
Research_ID_sum$derm_any <- ifelse(Research_ID_sum$derm_claims >0,1,0)
Research_ID_sum$endo_any <- ifelse(Research_ID_sum$endo_claims >0,1,0)
Research_ID_sum$gastro_any <- ifelse(Research_ID_sum$gastro_claims >0,1,0)
Research_ID_sum$genetic_any <- ifelse(Research_ID_sum$genetic_claims >0,1,0)
Research_ID_sum$genito_any <- ifelse(Research_ID_sum$genito_claims >0,1,0)
Research_ID_sum$hemato_any <- ifelse(Research_ID_sum$hemato_claims >0,1,0)
Research_ID_sum$immuno_any <- ifelse(Research_ID_sum$immuno_claims >0,1,0)
Research_ID_sum$malign_any <- ifelse(Research_ID_sum$malign_claims >0,1,0)
Research_ID_sum$mh_any <- ifelse(Research_ID_sum$mh_claims >0,1,0)
Research_ID_sum$metab_any <- ifelse(Research_ID_sum$metab_claims >0,1,0)
Research_ID_sum$musculo_any <- ifelse(Research_ID_sum$musculo_claims >0,1,0)
Research_ID_sum$neuro_any <- ifelse(Research_ID_sum$neuro_claims >0,1,0)
Research_ID_sum$pulresp_any <- ifelse(Research_ID_sum$pulresp_claims >0,1,0)
Research_ID_sum$renal_any <- ifelse(Research_ID_sum$renal_claims >0,1,0)
Research_ID_sum$opthal_any <- ifelse(Research_ID_sum$opthal_claims >0,1,0)
Research_ID_sum$otol_any <- ifelse(Research_ID_sum$otol_claims >0,1,0)
Research_ID_sum$otolar_any <- ifelse(Research_ID_sum$otolar_claims >0,1,0)
Research_ID_sum$progressive_any <- ifelse(Research_ID_sum$progressive_claims >0,1,0)

# a row sum gives claims where more than one system was tagged
Research_ID_sum$systems_any =  Research_ID_sum$cardiac_any+
  Research_ID_sum$cranio_any+
  Research_ID_sum$derm_any+
  Research_ID_sum$endo_any+
  Research_ID_sum$gastro_any+
  Research_ID_sum$genetic_any+
  Research_ID_sum$genito_any+
  Research_ID_sum$hemato_any+
  Research_ID_sum$immuno_any+
  ### malign is not a system 
  Research_ID_sum$mh_any+
  Research_ID_sum$metab_any+
  Research_ID_sum$musculo_any+
  Research_ID_sum$neuro_any+
  Research_ID_sum$pulresp_any+
  Research_ID_sum$renal_any+
  Research_ID_sum$opthal_any+
  Research_ID_sum$otol_any+
  Research_ID_sum$otolar_any
### progressive is not a system  

# Assign PMCA ####

#    'Complex Chronic':  1) more than one body system is involved, OR 
#                     2) one or more conditions are progressive, OR 
#                     3) one or more conditions are malignant

Research_ID_sum$cc_3 <- ifelse(Research_ID_sum$systems_any >= 2 | 
                                 Research_ID_sum$progressive_any >= 1 | 
                                 Research_ID_sum$malign_any >= 1,1,0)

# 'Non-complex Chronic': 1) only one body system is involved, AND 
#                        2) the condition is not progressive or malignant
# 

Research_ID_sum$ncc_2 <- ifelse(Research_ID_sum$systems_any == 1 & 
                                  Research_ID_sum$progressive_any < 1 & 
                                  Research_ID_sum$malign_any < 1,1,0)

# 'Non-Chronic': 		1) no body system indicators are present, AND 
#                   2) the condition is not progressive or malignant

Research_ID_sum$nc_1 <- ifelse(Research_ID_sum$systems_any < 1 & 
                                 Research_ID_sum$progressive_any < 1 & 
                                 Research_ID_sum$malign_any < 1,1,0)

# give text designation 
Research_ID_sum$PMCA_category <- ifelse(Research_ID_sum$cc_3 == 1,
                                        "Complex Chronic",
                                        ifelse(Research_ID_sum$ncc_2 == 1,
                                               "Non-complex Chronic",
                                               "Non-chronic"))


# output a file called pmca_end_YYYYMMDD.csv last_date_txt
# variables are Research_ID = Research_ID, PMCA_category text
Research_ID_pmca <- subset(Research_ID_sum, 
                           select=c("Research_ID","PMCA_category"))

#write.csv(Research_ID_pmca,
#          file = "outputs/dataframes/pmca_first_12.05.csv",
#          row.names=FALSE,na="")