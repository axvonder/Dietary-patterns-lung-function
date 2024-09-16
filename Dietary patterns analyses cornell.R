# set working directory and load packages
setwd("//rschfs1x/userrs/R-Z/wav25_RS/Desktop/updated code")
# List of required packages
required_packages <- c("tidyverse", "dplyr", "lme4", "lmerTest", "optimx", "forcats", 
                       "broom", "purrr", "metafor", "stringr", "mgcv", "gamm4", 
                       "survival", "survminer", "coxme", "ggplot2", "openxlsx",
                       "forestplot", "stringer")

# Install packages that are not already installed
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load the required packages
lapply(required_packages, library, character.only = TRUE)

#load data
HABCfinal <- read.csv("HABCfinal.csv")
RASfinal <- read.csv("RASfinal.csv")
#variable transformation & harmonization
nameF <- function(RASdat, RASvar, HABCvar) {
  for (i in seq_along(RASvar)) { #function to rename RAS variables to match HABC
    colname <- RASvar[i]
    new <- HABCvar[i]
    if (colname %in% colnames(RASdat)) {
      colnames(RASdat)[colnames(RASdat) == colname] <- new
    }
  }
  return(RASdat)
}
HABCfinal <- rename(HABCfinal, ID = HABCID)
HABCnumVars <- c("AGE", "HTCM", "FEV1", "FVC", "RATIO", "PACKYR1", "HEI2020score", "DASHscore")
RASnumVars <- c("age", "height", "fev1", "fvc", "ratio", "packyears", "HEI2020score", "DASHscore")
numVars <- HABCnumVars
RASfinal <- nameF(RASfinal, RASnumVars, numVars)
HABCcatVars <- c("RACE", "GENDER", "SMKCAT", "EDUC", "SITE")
RAScatVars <- c("race", "Gender", "SmokingStatus", "education", "InstitutionName")
catVars <- HABCcatVars
RASfinal <- nameF(RASfinal, RAScatVars, catVars)
HABCvars <- c("ID", "QCFEV1", "QCFVC", "TIME", "N", "HT2" ,"AGE2")
RASvars <- c("participantid", "FEV1_score", "FVC_score", "time", "count", "ht2", "age2")
RASfinal <- nameF(RASfinal, RASvars, HABCvars)
HABCfinal <- HABCfinal %>% mutate_at(catVars, as.factor) %>% mutate_at(numVars, as.numeric)
RASfinal <- RASfinal %>% mutate_at(catVars, as.factor) %>% mutate_at(numVars, as.numeric)
HABCfinal$RACE <- fct_recode(HABCfinal$RACE,
                             "White" = "1",
                             "Black" = "2") %>% fct_relevel("White")
HABCfinal$GENDER <- fct_recode(HABCfinal$GENDER,
                              "Male" = "1",
                              "Female" = "2") %>% fct_relevel("Male")
HABCfinal$EDUC <- fct_recode(HABCfinal$EDUC,
                             "< High School" = "1",
                             "High School/GED" = "2",
                             "Postsecondary" = "3") %>% relevel("< High School")
educlevs <- c("< High School", "High School/GED", "Postsecondary")
HABCfinal$EDUC <- factor(HABCfinal$EDUC, educlevs) %>% relevel("< High School")
HABCfinal$SITE <- fct_recode(HABCfinal$SITE,
                             "Memphis" = "1",
                             "Pittsburgh" = "2") %>% fct_relevel("Memphis")
RASfinal$RACE <- fct_recode(RASfinal$RACE,
                            "White" = "1",
                            "Black" = "2",
                            "Hispanic (non-Black)" = "3",
                            "Hispanic (Black)" = "4",
                            "Other" = "5") %>% fct_relevel("White")
racelevs <- c("White", "Black", "Hispanic (non-Black)", "Hispanic (Black)", "Other")
RASfinal$RACE <- factor(RASfinal$RACE, racelevs) %>% fct_relevel("White")
RASfinal$RACE2 <- factor(RASfinal$RACE, levels = c("White", "Black"))
RASfinal$GENDER <- fct_recode(RASfinal$GENDER,
                              "Male" = "0")
smoklevs <- c("Never", "Former", "Intermittent", "Persistent")
HABCfinal$SMKCAT <- factor(HABCfinal$SMKCAT, smoklevs) %>% fct_relevel("Never")
RASfinal$SMKCAT <- factor(RASfinal$SMKCAT, smoklevs) %>% fct_relevel("Never")
HABCfinal <- HABCfinal %>%
  mutate(SMKCAT2 = factor(case_when(
    SMKCAT == "Never" ~ "Never",
    SMKCAT %in% c("Former", "Intermittent", "Persistent") ~ "Ever",
    TRUE ~ NA
  )))
RASfinal <- RASfinal %>%
  mutate(SMKCAT2 = factor(case_when(
    SMKCAT == "Never" ~ "Never",
    SMKCAT %in% c("Former", "Intermittent", "Persistent") ~ "Ever",
    TRUE ~ NA
  )))
dietScores <- c("HEI2020score", "DASHscore")
PFTs <- c("FEV1", "FVC", "RATIO")
HABCHEIquarts <- quantile(HABCfinal$HEI2020score, probs = c(0.25, 0.5, 0.75))
RASHEIquarts <- quantile(RASfinal$HEI2020score, probs = c(0.25, 0.5, 0.75))
assignQuart <- function(score, quartiles) { #function to assign quartiles of HEI score (for healthy/unhealthy analysis)
  cut(score, breaks = c(-Inf, quartiles[1], quartiles[2], quartiles[3], Inf),
      labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)
}
HABCfinal$HEIQuartile <- assignQuart(HABCfinal$HEI2020score, HABCHEIquarts)
RASfinal$HEIQuartile <- assignQuart(RASfinal$HEI2020score, RASHEIquarts)
HABCfinal$HEIhalf <- ifelse(HABCfinal$HEIQuartile == "Q1" | HABCfinal$HEIQuartile == "Q2", "Unhealthy", "Healthy")
RASfinal$HEIhalf <- ifelse(RASfinal$HEIQuartile == "Q1" | RASfinal$HEIQuartile == "Q2", "Unhealthy", "Healthy")
#only participants who will be included in analysis & slice at first observation
HABCbaseline <- HABCfinal %>% group_by(ID) %>%
  filter((any(QCFEV1 %in% c(3, 4) | QCFVC %in% c(3, 4))) & N == 1)
RASbaseline <- RASfinal %>% group_by(ID) %>%
  filter((any(QCFEV1 %in% c(3, 4) | QCFVC %in% c(3, 4))) & N == 1)
rm(educlevs, HABCcatVars, HABCHEIquarts, HABCnumVars, HABCvars, racelevs, RAScatVars, RASHEIquarts, RASnumVars, RASvars, smoklevs, assignQuart, nameF)

###TABLE 1 - baseline characteristics###
catLevs <- list(
  RACE = c("White", "Black", "Hispanic (non-Black)", "Hispanic (Black)", "Other"),
  GENDER = c("Male", "Female"),
  SMKCAT = c("Never", "Former", "Intermittent", "Persistent"),
  EDUC = c("< High School", "High School/GED", "Postsecondary"),
  SITE = c("Memphis", "Pittsburgh", "Altamira Family Med", "Centre de Recherche", "Harbor-UCLA", "London HSC", "London Reg Ca Ctr", "MD Anderson", "Rush Univ Med Ctr", "San Diego, U of CA", "SUNY Stony Brook", "Swedish Medical Ctr", "UCSD - Chula Vista", "Upstate Carolina", "VAMC Jesse Brown", "VAMC Kansas City", "VAMC Minneapolis", "VAMC Puget Sound", "VAMC Washington DC", "Wichita CCOP")
)
meanie <- function(var, df, sub = NULL, cat = NULL) {
  if (!is.null(sub) && !is.null(cat)) { #function to calculate mean, sd, and missingness of numeric variables
    df <- df[df[[sub]] != cat, ]
  }
  mean_val <- mean(df[[var]], na.rm = TRUE)
  sd_val <- sd(df[[var]], na.rm = TRUE)
  miss <- sum(is.na(df[[var]])) / nrow(df) * 100
  return(c(mean = mean_val, sd = sd_val, missing = miss))
}
freqie <- function(var, df) {
  tab <- table(df[[var]], useNA = "ifany")
  total <- sum(tab)
  perc <- 100 * tab / total
  levels <- names(tab)
  freqs <- sprintf("%d (%.2f)", tab, perc)
  if (is.na(levels[length(levels)])) {
    freqs <- freqs[-length(freqs)]
    levels <- levels[-length(levels)]
    missing_value <- sprintf("%.2f", perc[length(perc)])
  } else {
    missing_value <- "0.00"
  }
  return(list(levels = levels, freqs = freqs, missing = missing_value))
}
table1 <- data.frame(characteristic = character(), HABC = character(), missing_HABC = numeric(), 
                     RAS = character(), missing_RAS = numeric(), stringsAsFactors = FALSE) #empty table
for (var in numVars) { #apply numeric variables
  if (var == "PACKYR1") {
    habc <- meanie(var, HABCbaseline, "SMKCAT", "Never")
    ras <- meanie(var, RASbaseline, "SMKCAT", "Never")
  } else {
    habc <- meanie(var, HABCbaseline)
    ras <- meanie(var, RASbaseline)
  }
  table1 <- rbind(table1, data.frame(
    characteristic = var,
    HABC = sprintf("%.2f (±%.2f)", habc["mean"], habc["sd"]),
    missing_HABC = sprintf("%.2f", habc["missing"]),
    RAS = sprintf("%.2f (±%.2f)", ras["mean"], ras["sd"]),
    missing_RAS = sprintf("%.2f", ras["missing"])
  ))
}
for (var in catVars) {
  habcFreqs <- freqie(var, HABCbaseline)
  rasFreqs <- freqie(var, RASbaseline)
  var_levels <- catLevs[[var]]
  table1 <- rbind(table1, data.frame( #"header" row for each categorical variable
    characteristic = var,
    HABC = NA,
    missing_HABC = habcFreqs$missing,
    RAS = NA,
    missing_RAS = rasFreqs$missing
  ))
  for (lev in var_levels) {#add levels of categorical variables (beneath each catLev header)
    habcLevelFreq <- ifelse(lev %in% habcFreqs$levels, habcFreqs$freqs[habcFreqs$levels == lev], "")
    rasLevelFreq <- ifelse(lev %in% rasFreqs$levels, rasFreqs$freqs[rasFreqs$levels == lev], "")
    table1 <- rbind(table1, data.frame(
      characteristic = lev,
      HABC = habcLevelFreq,
      missing_HABC = NA,
      RAS = rasLevelFreq,
      missing_RAS = NA
    ))
  }
}
colnames(table1) <- c("Characteristic", "Health ABC", "missing (%)", "RAS", "missing (%)")
table1$`Health ABC` <- ifelse(table1$`Health ABC` == "", "-", table1$`Health ABC`)
table1$RAS <- ifelse(table1$RAS == "", "-", table1$RAS)
ncount <- data.frame( #add n-count (number of observations in HABCbaseline & RASbaseline)
  Characteristic = NA, `Health ABC` = "n = 2,560", `missing (%)` = NA, RAS = "n = 2,864", `missing (%).1` = NA)
colnames(ncount) <- colnames(table1)
table1 <- rbind(ncount, table1)
rm(catLevs, habcFreqs, rasFreqs, ncount, habc, habcLevelFreq, lev, ras, rasLevelFreq, var, var_levels, freqie, meanie)






###TABLE 2 - dietary pattern by covariates###
table2 <- data.frame(characteristic = character(), HABC_HEI = character(), HABC_DASH = character(),
                     RAS_HEI = character(), RAS_DASH = character(), stringsAsFactors = FALSE) #empty table
t2List <- list() #empty list to store results
meanie2 <- function(var, df) {
  mean_val <- mean(df[[var]], na.rm = TRUE)
  sd_val <- sd(df[[var]], na.rm = TRUE)
  return(sprintf("%.2f (±%.2f)", mean_val, sd_val))
}
pval <- function(dietVar, var, df) {
  if (length(unique(df[[var]])) > 2) {#ANOVA for more than 2 unique values
    test <- aov(df[[dietVar]] ~ df[[var]])
    p_value <- summary(test)[[1]]$'Pr(>F)'[1]
    if (p_value < 0.001) {p_value <- "<0.001"} else {p_value <- round(p_value, 3)}
  } else if (length(unique(df[[var]])) == 2) {#T-test for exactly 2 unique values
    test <- t.test(df[[dietVar]] ~ df[[var]])
    p_value <- test$p.value
    if (p_value < 0.001) {p_value <- "<0.001"} else {p_value <- round(p_value, 3)}
  } else { #If there are less than 2 levels, set p_value to "-" (NA)
    p_value <- "-"
  }
  return(p_value)
}
for (cohort in c("HABC", "RAS")) { #loop to calculate mean (SD) & pval for each diet score by each categorical covariate level
  if (cohort == "HABC") { 
    df <- HABCbaseline
  } else {
    df <- RASbaseline
  }
  for (score in dietScores) {
    xListName <- paste(cohort, score, "List", sep = "")
    xList <- list() 
    for (var in catVars) {
      levs <- levels(as.factor(df[[var]]))
      colName <- paste(cohort, sub("score", "", score), sep = "_")
      for (lev in levs) {
        sub <- df[df[[var]] == lev, ]
        result <- meanie2(score, sub)
        row <- setNames(data.frame(characteristic = paste(var, lev, sep = ": "), result), c("characteristic", colName))
        xList[[length(xList) + 1]] <- row
      }
      p_val <- pval(score, var, df)
      pRow <- setNames(data.frame(characteristic = paste(var, "P-value", sep = ": "), p_val), c("characteristic", colName))
      xList[[length(xList) + 1]] <- pRow
    }
    t2List[[xListName]] <- xList
  }
}
all <- c()
for (name in names(t2List)) { #create a list of all unique characteristics (to shell table)
  for (df in t2List[[name]]) {
    all <- c(all, df$characteristic)
  }
}
all <- unique(all)
orderedLevels <- c("RACE: White", "RACE: Black", "RACE: Hispanic (non-Black)", "RACE: Hispanic (Black)", "RACE: Other", "RACE: P-value",
                    "GENDER: Male", "GENDER: Female", "GENDER: P-value",
                    "SMKCAT: Never", "SMKCAT: Former", "SMKCAT: Intermittent", "SMKCAT: Persistent", "SMKCAT: P-value",
                    "EDUC: < High School", "EDUC: High School/GED", "EDUC: Postsecondary", "EDUC: P-value",
                    "SITE: Memphis", "SITE: Pittsburgh", "SITE: Altamira Family Med", "SITE: Centre de Recherche", "SITE: Harbor-UCLA", "SITE: London HSC", "SITE: London Reg Ca Ctr", "SITE: MD Anderson", "SITE: Rush Univ Med Ctr", "SITE: San Diego, U of CA", "SITE: SUNY Stony Brook", "SITE: Swedish Medical Ctr", "SITE: UCSD - Chula Vista", "SITE: Upstate Carolina", "SITE: VAMC Jesse Brown", "SITE: VAMC Kansas City", "SITE: VAMC Minneapolis", "SITE: VAMC Puget Sound", "SITE: VAMC Washington DC", "SITE: Wichita CCOP", "SITE: P-value")
allOrdered <- factor(all, levels = orderedLevels)
all <- all[order(allOrdered)] #sort
table2 <- data.frame(characteristic = all, 
                     HABC_HEI = rep(NA, length(all)), 
                     HABC_DASH = rep(NA, length(all)),
                     RAS_HEI = rep(NA, length(all)), 
                     RAS_DASH = rep(NA, length(all)), 
                     stringsAsFactors = FALSE)
mapVal <- function(list, charName, value) {
  for (df in list) { #create function to map values from t2List to table2
    if (df$characteristic == charName && !is.null(df[[value]])) {
      return(as.character(df[[value]]))
    }
  }
  return("-") 
}
mapped <- list(HABC_HEI = "HABCHEI2020scoreList", HABC_DASH = "HABCDASHscoreList",
                       RAS_HEI = "RASHEI2020scoreList", RAS_DASH = "RASDASHscoreList")
for (col in names(mapped)) {
  Listy <- mapped[[col]]
  for (i in 1:nrow(table2)) {
    charName <- table2$characteristic[i]
    if (grepl("HEI", col)) { #match columns by dataset & diet score
      val <- ifelse(grepl("HABC", col), "HABC_HEI2020", "RAS_HEI2020")
    } else if (grepl("DASH", col)) {
      val <- ifelse(grepl("HABC", col), "HABC_DASH", "RAS_DASH")
    }
    table2[i, col] <- mapVal(t2List[[Listy]], charName, val)
  }
}
#formatting table
catVarHeads <- function(df) {
  catVarsunique <- unique(sub(":.*", "", df$characteristic))
  for (cat in catVarsunique) {
    idx <- which(startsWith(df$characteristic, paste0(cat, ":"))) #f indices of the category
    if (length(idx) > 0) { #add row before first occurrence of variable
      newRow <- data.frame(characteristic = cat, HABC_HEI = NA, HABC_DASH = NA,
                            RAS_HEI = NA, RAS_DASH = NA, stringsAsFactors = FALSE)
      df <- rbind(df[1:(idx[1]-1), ], newRow, df[idx[1]:nrow(df), ])
    }
  }
  return(df)
}
reName <- function(df) {
  df$characteristic <- gsub(".*: ", "", df$characteristic)
  return(df)
}
table2 <- catVarHeads(table2)
table2 <- reName(table2)
table2 <- table2[-1, ] #first row duplicating, can't figure out why - delete
ncount <- data.frame( #add n-count (number of observations in HABCbaseline & RASbaseline)
  characteristic = NA, HABC_HEI = "n = 2,560", HABC_DASH = NA, RAS_HEI = "n = 2,864", RAS_DASH = NA)
colnames(ncount) <- colnames(table2)
table2 <- rbind(ncount, table2)
rm(HABCbaseline, RASbaseline, df, mapped, ncount, pRow, row, sub, t2List, xList, all, allOrdered, catVars, charName, cohort, col, colName, i, lev, levs, Listy, name, numVars, orderedLevels, p_val, result, score, val, var, xListName, catVarHeads, mapVal, meanie2, pval, reName)




###TABLE 3 - regression; dietary pattern & PFT###
#set functions to run models
digi <- function(value, digits = 2, p_value = FALSE, is_i2 = FALSE) {
  numV <- as.numeric(value) #convert to numeric
  if (is_i2) { #i2 values rounded to the nearest whole number
    return(as.character(round(numV)))
  }
  if (p_value) { 
    if (numV < 0.001) {
      if (numV >= 0.0005) { #if between 0.0005 and 0.001, round to 3 dps
        return("0.001")
      } else {
        return("<0.001")
      }
    } else {
      return(as.character(sprintf("%.3f", numV))) #round pvals to 3 dps
    }
  } else {
    if (abs(numV) < 0.1) { #if <0.1, use scientific notation
      return(formatC(numV, format = "e", digits = digits)) #round to 2 dps (not pvals)
    } else {
      return(as.character(sprintf(paste0("%.", digits, "f"), numV)))
    }
  }
}
rounding <- function(tlist) {
  for (minmax in c("miniList", "maxiList")) {
    for (sub in names(tlist[[minmax]])) {
      for (df in c("cross", "long")) {
        tlist[[minmax]][[sub]][[df]][["B"]] <- digi(tlist[[minmax]][[sub]][[df]][["B"]])
        tlist[[minmax]][[sub]][[df]][["SE"]] <- digi(tlist[[minmax]][[sub]][[df]][["SE"]])
        tlist[[minmax]][[sub]][[df]][["pval"]] <- digi(tlist[[minmax]][[sub]][[df]][["pval"]], p_value = TRUE)
        if (startsWith(sub, "META")) { #only round i2 and Hetpval if the subsection starts with "META" (instead of HABC or RAS)
          tlist[[minmax]][[sub]][[df]][["i2"]] <- digi(tlist[[minmax]][[sub]][[df]][["i2"]], is_i2 = TRUE)
          tlist[[minmax]][[sub]][[df]][["Hetpval"]] <- digi(tlist[[minmax]][[sub]][[df]][["Hetpval"]], p_value = TRUE)
        }
      }
    }
  }
  return(tlist)
}
runModel <- function(model = "mini", PFT, diet, df, race = FALSE) { #minimally-adjusted model
  dat <- get(df, envir = .GlobalEnv)
  cohort <- if (startsWith(df, "RAS")) {
    "RAS"
  } else if (startsWith(df, "HABC")) {
    "HABC"
  }
  if (PFT == "FEV1") { #use QCs associated with PFT
    sub <- quote(QCFEV1 %in% c(3, 4))
  } else if (PFT == "FVC") {
    sub <- quote(QCFVC %in% c(3, 4))
  } else if (PFT == "RATIO") {
    sub <- quote(QCFEV1 %in% c(3, 4) & QCFVC %in% c(3, 4))
  }
  base <- paste(PFT, "~ TIME + AGE + AGE2 + HTCM + HT2",
                if(!race) {"+ RACE"} else {""}, #if stratifying by race, exclude RACE from model (set to TRUE when calling function)
                if(!(cohort == "RAS" || grepl("Male|Female", df))) {"+ GENDER"}, #if using RAS or gender subset data, exclude GENDER from model
                "+", diet, "* TIME + (TIME || ID)")
  full <- if(model == "maxi" && !grepl("Never", df)) {" + SITE + SMKCAT + PACKYR1"} else if(model == "maxi" && grepl("Never", df)) {" + SITE"} else {""} #if dataset contains "Never", don't include smoking variables
  form <- paste(base, full)
  m <- summary(lmer(as.formula(form), 
                    data = dat, REML = TRUE,
                    subset = eval(sub),
                    control = lmerControl(optimizer = "bobyqa", 
                                          optCtrl = list(maxfun = 100000))), ddf = "Kenward-Roger")
  coefs <- m$coefficients
  dietrow <- coefs[rownames(coefs) == diet, ]
  dietrowlong <- coefs[rownames(coefs) == paste0("TIME:", diet), ]
  text <- capture.output(print(m))
  nObs <- as.numeric(str_extract(text[grep("groups:  ID,", text)], "(?<=ID,\\s)\\d+"))
  cross <- c(Name = diet, B = dietrow[["Estimate"]], SE = dietrow[["Std. Error"]], pval = dietrow[["Pr(>|t|)"]],  N = nObs)
  long <- c(Name = paste0("TIME:", diet), B = dietrowlong[["Estimate"]], SE = dietrowlong[["Std. Error"]], pval = dietrowlong[["Pr(>|t|)"]], N = nObs)
  save <- paste0(cohort, "_", PFT, "_", diet)
  listName <- if(model == "mini") {"miniList"} else {"maxiList"}
  results <- list(save = save, cross = cross, long = long)
  return(results)
}
metaFunk <- function(H, R) { #meta-analysis function
  Bcross <- c(as.numeric(H$cross["B"]), as.numeric(R$cross["B"]))
  SEcross <- c(as.numeric(H$cross["SE"]), as.numeric(R$cross["SE"])) #cross-sectional Bs & SEs
  Mcross <- rma(yi = Bcross, sei = SEcross, method = "REML") #meta for cross-sectional
  crossRez <- c(B = as.character(Mcross[["b"]]), SE = Mcross$se, pval = Mcross$pval, i2 = Mcross$I2, Hetpval = Mcross$QEp)
  Blong <- c(as.numeric(H$long["B"]), as.numeric(R$long["B"]))
  SElong <- c(as.numeric(H$long["SE"]), as.numeric(R$long["SE"])) #longitudinal Bs & SEs
  Mlong <- rma(yi = Blong, sei = SElong, method = "REML") #meta for longitudinal
  longRez <- c(B = as.character(Mlong[["b"]]), SE = Mlong$se, pval = Mlong$pval, i2 = Mlong$I2, Hetpval = Mlong$QEp)
  return(list(cross = crossRez, long = longRez))
}
RunMeta <- function(ListX) { #run meta analysis on list of HABC & RAS results
  for(listName in names(ListX)) { 
    currentList <- ListX[[listName]]
    for (name in names(currentList)) { #meta loop to match HABC & RAS
      if (grepl("HABC", name)) { #find HABC/RAS matching counterparts in results list
        counter <- gsub("HABC", "RAS", name) #sub HABC for RAS name
        if (counter %in% names(currentList)) {
          result <- metaFunk(currentList[[name]], currentList[[counter]])
          metaName <- gsub("HABC_", "META_", name) #create a new name for the meta-analyzed result
          currentList[[metaName]] <- result #store meta-analyzed results in pre-defined list (mini or maxi)
        }
      }
    }
    ListX[[listName]] <- currentList #store both mini & maxi meta-analyzed results in list
  }
  return(ListX)
}
getResults <- function(cohorts, race = FALSE) {
  models <- c("mini", "maxi")
  PFTs <- c("FEV1", "FVC", "RATIO")
  dietScores <- c("HEI2020score", "DASHscore")
  miniList <- list() #create mini & maxi lists
  maxiList <- list() 
  for(cohort in cohorts) { #run mini & maxi models for all PFTs, diet scores, both cohorts
    for(PFT in PFTs) {
      for(dietScore in dietScores) {
        for(model in models) { #dynamically call either mini or maxi function
          cat("Running", model, "model:", cohort, ",", PFT, ",", dietScore, "\n")
          result <- runModel(model, PFT, dietScore, cohort, race)
          if (model == "mini") {
            miniList[[result$save]] <- list(cross = result$cross, long = result$long)
          } else {
            maxiList[[result$save]] <- list(cross = result$cross, long = result$long)
          }
        }
      }
    }
  }
  results <- list(miniList = miniList, maxiList = maxiList)
  if (!any(grepl("Female", cohorts))) {
    results <- RunMeta(results) #don't run meta anlaysis if sex = Female
  }
  results <- rounding(results) 
  return(results)
}
filltab <- function(table, list, type) { #function to fill in table with results from list
  for (name in names(list)) {
    crossVals <- list[[name]]$cross #extract cross-sectional
    longVals <- list[[name]]$long #extract longitudinal
    crossi = which(table$PFT == name) #find the row in the table that matches the name
    table[crossi, paste(type, "B", sep = "")] <- crossVals["B"]
    table[crossi, paste(type, "SE", sep = "")] <- crossVals["SE"]
    table[crossi, paste(type, "P", sep = "")] <- crossVals["pval"]
    longi = which(table$PFT == paste(name, "*TIME", sep = "")) #find the row in the table that matches the name
    table[longi, paste(type, "B", sep = "")] <- longVals["B"]
    table[longi, paste(type, "SE", sep = "")] <- longVals["SE"]
    table[longi, paste(type, "P", sep = "")] <- longVals["pval"]
  }
  return(table)
}
cohorts <- c("HABCfinal", "RASfinal")
t3List <- getResults(cohorts)
table3 <- data.frame(PFT = character(), miniB = character(), miniSE = character(), miniP = character(),
                     maxiB = character(), maxiSE = character(), maxiP = character(), stringsAsFactors = FALSE) #empty table
cohortNames <- c("HABC", "RAS", "META")
for (pft in PFTs) { #add row for PFT
  table3 <- rbind(table3, data.frame(PFT = pft, miniB = NA, miniSE = NA, miniP = NA, 
                                     maxiB = NA, maxiSE = NA, maxiP = NA))
  for (cohort in cohortNames) { #add row for cohort
    table3 <- rbind(table3, data.frame(PFT = paste0(cohort, "_", pft), miniB = NA, miniSE = NA, miniP = NA, 
                                       maxiB = NA, maxiSE = NA, maxiP = NA))
    for (score in dietScores) { #add each diet score (+longitudinal) row
      table3 <- rbind(table3, data.frame(PFT = paste0(cohort, "_", pft, "_", score), 
                                         miniB = NA, miniSE = NA, miniP = NA, 
                                         maxiB = NA, maxiSE = NA, maxiP = NA))
      table3 <- rbind(table3, data.frame(PFT = paste0(cohort, "_", pft, "_", paste0(score, "*TIME")), 
                                         miniB = NA, miniSE = NA, miniP = NA, 
                                         maxiB = NA, maxiSE = NA, maxiP = NA))
    }
  }
}
models <- c("mini", "maxi")
for (model in models) { #fill in table with results
  listName <- paste(model, "List", sep = "")
  table3 <- filltab(table3, t3List[[listName]], model)
}


nVals <- function(table, resultsList, prefixes = c("X", "Y")) {
  diet <- "HEI2020score" #doesn't matter which to use, just picking HEI2020
  for (prefix in prefixes) { #loop over each prefix - used for stratified analyses
    for (PFT in PFTs) { #pull n-values from resultsList - repeat for each PFT
      N_mini_HABC <- as.numeric(resultsList[[paste0(prefix, "miniList")]][[paste0("HABC_", PFT, "_", diet)]]$cross["N"])
      N_maxi_HABC <- as.numeric(resultsList[[paste0(prefix, "maxiList")]][[paste0("HABC_", PFT, "_", diet)]]$cross["N"])
      N_mini_RAS <- as.numeric(resultsList[[paste0(prefix, "miniList")]][[paste0("RAS_", PFT, "_", diet)]]$cross["N"])
      N_maxi_RAS <- as.numeric(resultsList[[paste0(prefix, "maxiList")]][[paste0("RAS_", PFT, "_", diet)]]$cross["N"])
      N_mini_META <- N_mini_HABC + N_mini_RAS #sum Nvals for meta-analysis
      N_maxi_META <- N_maxi_HABC + N_maxi_RAS
      rowIndex_HABC <- which(table$PFT == paste0("HABC_", PFT)) #row indices for each cohort
      rowIndex_RAS <- which(table$PFT == paste0("RAS_", PFT))
      rowIndex_META <- which(table$PFT == paste0("META_", PFT))
      #add n-values to central column for each model
      table[rowIndex_HABC, paste0(prefix, "miniSE")] <- paste("n =", N_mini_HABC)
      table[rowIndex_HABC, paste0(prefix, "maxiSE")] <- paste("n =", N_maxi_HABC)
      table[rowIndex_RAS, paste0(prefix, "miniSE")] <- paste("n =", N_mini_RAS)
      table[rowIndex_RAS, paste0(prefix, "maxiSE")] <- paste("n =", N_maxi_RAS)
      table[rowIndex_META, paste0(prefix, "miniSE")] <- paste("n =", N_mini_META)
      table[rowIndex_META, paste0(prefix, "maxiSE")] <- paste("n =", N_maxi_META)
    }
  }
  return(table)
}
table3 <- nVals(table3, t3List, c("")) #add n-values - no prefixes for main table
rm(cohort, cohortNames, cohorts, listName, model, models, pft, score)
#formatting table
rowNames <- function(table, col_name){ #create function to format row names
  table %>%
    mutate("{col_name}" := case_when(
      grepl("FEV1$", .data[[col_name]]) & !grepl("score", .data[[col_name]]) ~ sub("_FEV1$", "", .data[[col_name]]), 
      grepl("FVC$", .data[[col_name]]) & !grepl("score", .data[[col_name]]) ~ sub("_FVC$", "", .data[[col_name]]),
      grepl("RATIO$", .data[[col_name]]) & !grepl("score", .data[[col_name]]) ~ sub("_RATIO$", "", .data[[col_name]]),
      grepl("HEI2020score$", .data[[col_name]]) ~ "HEI2020",
      grepl("DASHscore$", .data[[col_name]]) ~ "DASH",
      grepl("HEI2020score\\*TIME", .data[[col_name]]) ~ "HEI2020 * Time",
      grepl("DASHscore\\*TIME", .data[[col_name]]) ~ "DASH * Time",
      TRUE ~ .data[[col_name]] # Default case to keep original value
    ))
}
table3 <- rowNames(table3, "PFT")












###APPENDIX 3 - 3-way interaction model for dietary pattern & covariates###
#adjust model function for 3-way analysis
way3Model <- function(model = "mini", PFT, diet, df, covariate) { #minimally-adjusted model
  dat <- get(df, envir = .GlobalEnv)
  cohort <- if (startsWith(df, "RAS")) {
    "RAS"
  } else if (startsWith(df, "HABC")) {
    "HABC"
  }
  if (PFT == "FEV1") { #use QCs associated with PFT
    sub <- quote(QCFEV1 %in% c(3, 4))
  } else if (PFT == "FVC") {
    sub <- quote(QCFVC %in% c(3, 4))
  } else if (PFT == "RATIO") {
    sub <- quote(QCFEV1 %in% c(3, 4) & QCFVC %in% c(3, 4))
  }
  race_var <- if(cohort == "RAS" && covariate == "RACE2") { "RACE2" } else { "RACE" }
  smoke_var <- if(covariate == "SMKCAT2") { "SMKCAT2" } else { "SMKCAT" }
  base <- paste(PFT, "~ TIME + AGE + AGE2 + HTCM + HT2",
                "+", race_var,
                if(cohort != "RAS") {"+ GENDER"}, #if using RAS, exclude GENDER from model
                "+", diet, "* TIME *", covariate, "+ (TIME || ID)")
  full <- if(model == "maxi") {paste(" + SITE +", smoke_var, "+ PACKYR1")} else { "" }
  form <- paste(base, full)
  m <- lmer(as.formula(form), 
          data = dat, REML = TRUE,
          subset = eval(sub),
          control = lmerControl(optimizer = "bobyqa", 
                                optCtrl = list(maxfun = 100000)))
  a <- anova(m, type = "III")
  # intTerm <- paste("TIME", covariate, diet, sep = ":")
  # pint <- round(a[rownames(a) == intTerm, "Pr(>F)"], 3)
  pint <- sprintf("%.3f", a[nrow(a), "Pr(>F)"]) #get p-value for 3-way interaction
  save <- paste0(cohort, "_", PFT, "_", diet, "_", covariate)
  result <- list(save = save, pint = pint)
  return(result)
}
way3Results <- function(cohorts, race = FALSE) {
  models <- c("mini", "maxi")
  PFTs <- c("FEV1", "FVC", "RATIO")
  dietScores <- c("HEI2020score", "DASHscore")
  covariates <- c("GENDER", "RACE", "RACE2", "SMKCAT2", "SMKCAT", "HEIhalf")
  miniList <- list() #create mini & maxi lists
  maxiList <- list() 
  for(cohort in cohorts) { #run mini & maxi models for all PFTs, diet scores, both cohorts
    for(PFT in PFTs) {
      for(dietScore in dietScores) {
        for(covariate in covariates) {
          if ((covariate == "GENDER" && startsWith(cohort, "RAS")) ||
              (covariate == "RACE2" && startsWith(cohort, "HABC"))) {
            next # GENDER only applies to HABC cohort; RACE2 only applies to RAS cohort
          }
          for(model in models) { #dynamically call either mini or maxi function
            cat("Running", model, "model:", cohort, ",", PFT, ",", dietScore, ",", covariate, "\n")
            result <- way3Model(model, PFT, dietScore, cohort, covariate)
            if (model == "mini") {
              miniList[[result$save]] <- list(pint = result$pint)
            } else {
              maxiList[[result$save]] <- list(pint = result$pint)
            }
          }
        }
      }
    }
  }
  results <- list(miniList = miniList, maxiList = maxiList)
  return(results)
}
cohorts <- c("HABCfinal", "RASfinal")
ap3List <- way3Results(cohorts)
ap3 <- data.frame(name = character(), miniFEV1 = character(), miniFVC = character(), miniRATIO = character(),
                     maxiFEV1 = character(), maxiFVC = character(), maxiRATIO = character(), stringsAsFactors = FALSE) #empty table
cohortNames <- c("HABC", "RAS")
covariates <- c("GENDER", "SMKCAT2", "SMKCAT", "RACE", "RACE2", "HEIhalf")
for (covariate in covariates) { #add row for covariate
  ap3 <- rbind(ap3, data.frame(name = covariate, miniFEV1 = NA, miniFVC = NA, miniRATIO = NA, 
                               maxiFEV1 = NA, maxiFVC = NA, maxiRATIO = NA))
  for (cohort in cohortNames) { #add row for cohort
    if (!(covariate == "GENDER" && cohort == "RAS") && !(covariate == "RACE2" && cohort == "HABC")) { #skip GENDER for RAS and RACE2 for HABC
      ap3 <- rbind(ap3, data.frame(name = paste0(cohort, "_", covariate), miniFEV1 = NA, miniFVC = NA, miniRATIO = NA, 
                                   maxiFEV1 = NA, maxiFVC = NA, maxiRATIO = NA))
      for (score in dietScores) { #add each diet score (+longitudinal) row
        ap3 <- rbind(ap3, data.frame(name = paste0(cohort, "_", covariate, "_", score), 
                                     miniFEV1 = NA, miniFVC = NA, miniRATIO = NA, 
                                     maxiFEV1 = NA, maxiFVC = NA, maxiRATIO = NA))
      }
    }
  }
}
parse <- function(name) { #function to read/split out components of variable names
  parts <- strsplit(name, "_")[[1]]
  list(cohort = parts[1], covariate = parts[4], score = parts[3])
}
Lists <- list(mini = ap3List$miniList, maxi = ap3List$maxiList)
for (prefix in names(Lists)) {
  current <- Lists[[prefix]]
  for (minmax in names(current)) {
    parsed <- parse(minmax)  
    pattern <- paste0(parsed$cohort, "_", parsed$covariate, "_", parsed$score) #create name pattern to match ap3
    rowi <- grep(pattern, ap3$name) #match row in ap3
    if (length(rowi) > 0) {
      pft <- ifelse(grepl("FEV1", minmax), "FEV1", ifelse(grepl("FVC", minmax), "FVC", "RATIO"))
      ap3[rowi, paste0(prefix, pft)] <- current[[minmax]]$pint #ppdate matching row in ap3 with pint from ap3List (mini/maxi)
    }
  }
}
rm(Lists, parsed, cohort, cohortNames, cohorts, covariate, covariates, current, minmax, pattern, pft, prefix, rowi, score, way3Model, way3Results)
#formatting table
ap3 <- rowNames(ap3, "name")






###APPENDIX 5 - diet—PFT association by race (White, Black)###
#subset datasets for white & black participants
HABCfinal_White <- subset(HABCfinal, RACE == "White")
HABCfinal_Black <- subset(HABCfinal, RACE == "Black")
RASfinal_White <- subset(RASfinal, RACE == "White")
RASfinal_Black <- subset(RASfinal, RACE == "Black")
#run models for each subset
cohortGroups <- list(
  White = c("HABCfinal_White", "RASfinal_White"),
  Black = c("HABCfinal_Black", "RASfinal_Black")
)
ap5List <- list()
for (race in names(cohortGroups)) { #run model function to get results
  prefix <- substr(race, 1, 1) #prefix for races (W, B)
  rez <- getResults(cohortGroups[[race]], race = TRUE)
  ap5List[[paste0(prefix, "miniList")]] <- rez$miniList
  ap5List[[paste0(prefix, "maxiList")]] <- rez$maxiList
}
ap5 <- data.frame(PFT = character(),
                     WminiB = character(), WminiSE = character(), WminiP = character(),
                     WmaxiB = character(), WmaxiSE = character(), WmaxiP = character(),
                     BminiB = character(), BminiSE = character(), BminiP = character(),
                     BmaxiB = character(), BmaxiSE = character(), BmaxiP = character(), stringsAsFactors = FALSE) #empty table
cohortNames <- c("HABC", "RAS", "META")
for (pft in PFTs) { #add row for PFT
  ap5 <- rbind(ap5, data.frame(PFT = pft,
                                     WminiB = NA, WminiSE = NA, WminiP = NA, 
                                     WmaxiB = NA, WmaxiSE = NA, WmaxiP = NA,
                                     BminiB = NA, BminiSE = NA, BminiP = NA, 
                                     BmaxiB = NA, BmaxiSE = NA, BmaxiP = NA))
  for (cohort in cohortNames) { #add row for cohort
    ap5 <- rbind(ap5, data.frame(PFT = paste0(cohort, "_", pft),
                                       WminiB = NA, WminiSE = NA, WminiP = NA, 
                                       WmaxiB = NA, WmaxiSE = NA, WmaxiP = NA,
                                       BminiB = NA, BminiSE = NA, BminiP = NA, 
                                       BmaxiB = NA, BmaxiSE = NA, BmaxiP = NA))
    for (score in dietScores) { #add each diet score (+longitudinal) row
      ap5 <- rbind(ap5, data.frame(PFT = paste0(cohort, "_", pft, "_", score), 
                                         WminiB = NA, WminiSE = NA, WminiP = NA, 
                                         WmaxiB = NA, WmaxiSE = NA, WmaxiP = NA,
                                         BminiB = NA, BminiSE = NA, BminiP = NA, 
                                         BmaxiB = NA, BmaxiSE = NA, BmaxiP = NA))
      ap5 <- rbind(ap5, data.frame(PFT = paste0(cohort, "_", pft, "_", paste0(score, "*TIME")), 
                                         WminiB = NA, WminiSE = NA, WminiP = NA, 
                                         WmaxiB = NA, WmaxiSE = NA, WmaxiP = NA,
                                         BminiB = NA, BminiSE = NA, BminiP = NA, 
                                         BmaxiB = NA, BmaxiSE = NA, BmaxiP = NA))
    }
  }
}
models <- c("Wmini", "Wmaxi", "Bmini", "Bmaxi")
for (model in models) { #fill in table with results
  listName <- paste(model, "List", sep = "")
  ap5 <- filltab(ap5, ap5List[[listName]], model)
}
ap5 <- nVals(ap5, ap5List, prefixes = c("W", "B")) #add n-values
rm(cohortGroups, HABCfinal_Black, HABCfinal_White, RASfinal_Black, RASfinal_White, rez, cohort, cohortNames, listName, model, models, pft, prefix, race, score)
#formatting table
ap5 <- rowNames(ap5, "PFT")





###APPENDIX 6 - dietary pattern by smoking status###
#subset datasets for never & ever smoking participants
HABCfinal_Never <- subset(HABCfinal, SMKCAT == "Never")
HABCfinal_Ever <- subset(HABCfinal, SMKCAT %in% c("Former", "Intermittent", "Persistent"))
RASfinal_Never <- subset(RASfinal, SMKCAT == "Never")
RASfinal_Ever <- subset(RASfinal, SMKCAT %in% c("Former", "Intermittent", "Persistent"))
#run models for each subset
cohortGroups <- list(
  Never = c("HABCfinal_Never", "RASfinal_Never"),
  Ever = c("HABCfinal_Ever", "RASfinal_Ever")
)
ap6List <- list()
for (smk in names(cohortGroups)) { #run model function to get results
  prefix <- substr(smk, 1, 1) #prefix for smoking groups (N, E)
  rez <- getResults(cohortGroups[[smk]])
  ap6List[[paste0(prefix, "miniList")]] <- rez$miniList
  ap6List[[paste0(prefix, "maxiList")]] <- rez$maxiList
}
ap6 <- data.frame(PFT = character(),
                  NminiB = character(), NminiSE = character(), NminiP = character(),
                  NmaxiB = character(), NmaxiSE = character(), NmaxiP = character(),
                  EminiB = character(), EminiSE = character(), EminiP = character(),
                  EmaxiB = character(), EmaxiSE = character(), EmaxiP = character(), stringsAsFactors = FALSE) #empty table
cohortNames <- c("HABC", "RAS", "META")
for (pft in PFTs) { #add row for PFT
  ap6 <- rbind(ap6, data.frame(PFT = pft,
                               NminiB = NA, NminiSE = NA, NminiP = NA, 
                               NmaxiB = NA, NmaxiSE = NA, NmaxiP = NA,
                               EminiB = NA, EminiSE = NA, EminiP = NA, 
                               EmaxiB = NA, EmaxiSE = NA, EmaxiP = NA))
  for (cohort in cohortNames) { #add row for cohort
    ap6 <- rbind(ap6, data.frame(PFT = paste0(cohort, "_", pft),
                                 NminiB = NA, NminiSE = NA, NminiP = NA, 
                                 NmaxiB = NA, NmaxiSE = NA, NmaxiP = NA,
                                 EminiB = NA, EminiSE = NA, EminiP = NA, 
                                 EmaxiB = NA, EmaxiSE = NA, EmaxiP = NA))
    for (score in dietScores) { #add each diet score (+longitudinal) row
      ap6 <- rbind(ap6, data.frame(PFT = paste0(cohort, "_", pft, "_", score), 
                                   NminiB = NA, NminiSE = NA, NminiP = NA, 
                                   NmaxiB = NA, NmaxiSE = NA, NmaxiP = NA,
                                   EminiB = NA, EminiSE = NA, EminiP = NA, 
                                   EmaxiB = NA, EmaxiSE = NA, EmaxiP = NA))
      ap6 <- rbind(ap6, data.frame(PFT = paste0(cohort, "_", pft, "_", paste0(score, "*TIME")), 
                                   NminiB = NA, NminiSE = NA, NminiP = NA, 
                                   NmaxiB = NA, NmaxiSE = NA, NmaxiP = NA,
                                   EminiB = NA, EminiSE = NA, EminiP = NA, 
                                   EmaxiB = NA, EmaxiSE = NA, EmaxiP = NA))
    }
  }
}
models <- c("Nmini", "Nmaxi", "Emini", "Emaxi")
for (model in models) { #fill in table with results
  listName <- paste(model, "List", sep = "")
  ap6 <- filltab(ap6, ap6List[[listName]], model)
}
ap6 <- nVals(ap6, ap6List, prefixes = c("N", "E"))
rm(cohortGroups, HABCfinal_Ever, HABCfinal_Never, RASfinal_Ever, RASfinal_Never, rez, cohort, cohortNames, listName, model, models, pft, prefix, smk, score)
#formatting table
ap6 <- rowNames(ap6, "PFT")









###APPENDIX 7 - dietary pattern by sex###
#subset HABC (RAS is all male)
HABCfinal_Male <- subset(HABCfinal, GENDER == "Male")
HABCfinal_Female <- subset(HABCfinal, GENDER == "Female")
cohortGroups <- list(
  Male = c("HABCfinal_Male", "RASfinal"),
  Female = c("HABCfinal_Female")
)
ap7List <- list()
for (gend in names(cohortGroups)) { #run model function to get results
  prefix <- substr(gend, 1, 1) #prefix for races (M, F)
  rez <- getResults(cohortGroups[[gend]])
  ap7List[[paste0(prefix, "miniList")]] <- rez$miniList
  ap7List[[paste0(prefix, "maxiList")]] <- rez$maxiList
}
ap7 <- data.frame(PFT = character(),
                  MminiB = character(), MminiSE = character(), MminiP = character(),
                  MmaxiB = character(), MmaxiSE = character(), MmaxiP = character(),
                  FminiB = character(), FminiSE = character(), FminiP = character(),
                  FmaxiB = character(), FmaxiSE = character(), FmaxiP = character(), stringsAsFactors = FALSE) #empty table
cohortNames <- c("HABC", "RAS", "META")
for (pft in PFTs) { #add row for PFT
  ap7 <- rbind(ap7, data.frame(PFT = pft,
                               MminiB = NA, MminiSE = NA, MminiP = NA, 
                               MmaxiB = NA, MmaxiSE = NA, MmaxiP = NA,
                               FminiB = NA, FminiSE = NA, FminiP = NA, 
                               FmaxiB = NA, FmaxiSE = NA, FmaxiP = NA))
  for (cohort in cohortNames) { #add row for cohort
    ap7 <- rbind(ap7, data.frame(PFT = paste0(cohort, "_", pft),
                                 MminiB = NA, MminiSE = NA, MminiP = NA, 
                                 MmaxiB = NA, MmaxiSE = NA, MmaxiP = NA,
                                 FminiB = NA, FminiSE = NA, FminiP = NA, 
                                 FmaxiB = NA, FmaxiSE = NA, FmaxiP = NA))
    for (score in dietScores) { #add each diet score (+longitudinal) row
      ap7 <- rbind(ap7, data.frame(PFT = paste0(cohort, "_", pft, "_", score), 
                                   MminiB = NA, MminiSE = NA, MminiP = NA, 
                                   MmaxiB = NA, MmaxiSE = NA, MmaxiP = NA,
                                   FminiB = NA, FminiSE = NA, FminiP = NA, 
                                   FmaxiB = NA, FmaxiSE = NA, FmaxiP = NA))
      ap7 <- rbind(ap7, data.frame(PFT = paste0(cohort, "_", pft, "_", paste0(score, "*TIME")), 
                                   MminiB = NA, MminiSE = NA, MminiP = NA, 
                                   MmaxiB = NA, MmaxiSE = NA, MmaxiP = NA,
                                   FminiB = NA, FminiSE = NA, FminiP = NA, 
                                   FmaxiB = NA, FmaxiSE = NA, FmaxiP = NA))
    }
  }
}
models <- c("Mmini", "Mmaxi", "Fmini", "Fmaxi")
for (model in models) { #fill in table with results
  listName <- paste(model, "List", sep = "")
  ap7 <- filltab(ap7, ap7List[[listName]], model)
}
ap7 <- nVals(ap7, ap7List, prefixes = c("M", "F"))
rm(cohortGroups, HABCfinal_Female, rez, cohort, cohortNames, listName, model, models, pft, prefix, gend, score)
#formatting table
ap7 <- rowNames(ap7, "PFT")






###APPENDIX 8 - dietary pattern among men by smoking status###
HABCfinal_Male_Never <- subset(HABCfinal_Male, SMKCAT == "Never")
HABCfinal_Male_Ever <- subset(HABCfinal_Male, SMKCAT %in% c("Former", "Intermittent", "Persistent"))
RASfinal_Never <- subset(RASfinal, SMKCAT == "Never")
RASfinal_Ever <- subset(RASfinal, SMKCAT %in% c("Former", "Intermittent", "Persistent"))
cohortGroups <- list(
  Never = c("HABCfinal_Male_Never", "RASfinal_Never"),
  Ever = c("HABCfinal_Male_Ever", "RASfinal_Ever")
)
ap8List <- list()
for (smk in names(cohortGroups)) { #run model function to get results
  prefix <- substr(smk, 1, 1) #prefix for smoking groups (N, E)
  rez <- getResults(cohortGroups[[smk]])
  ap8List[[paste0(prefix, "miniList")]] <- rez$miniList
  ap8List[[paste0(prefix, "maxiList")]] <- rez$maxiList
}
ap8 <- data.frame(PFT = character(),
                  NminiB = character(), NminiSE = character(), NminiP = character(),
                  NmaxiB = character(), NmaxiSE = character(), NmaxiP = character(),
                  EminiB = character(), EminiSE = character(), EminiP = character(),
                  EmaxiB = character(), EmaxiSE = character(), EmaxiP = character(), stringsAsFactors = FALSE) #empty table
cohortNames <- c("HABC", "RAS", "META")
for (pft in PFTs) { #add row for PFT
  ap8 <- rbind(ap8, data.frame(PFT = pft,
                               NminiB = NA, NminiSE = NA, NminiP = NA, 
                               NmaxiB = NA, NmaxiSE = NA, NmaxiP = NA,
                               EminiB = NA, EminiSE = NA, EminiP = NA, 
                               EmaxiB = NA, EmaxiSE = NA, EmaxiP = NA))
  for (cohort in cohortNames) { #add row for cohort
    ap8 <- rbind(ap8, data.frame(PFT = paste0(cohort, "_", pft),
                                 NminiB = NA, NminiSE = NA, NminiP = NA, 
                                 NmaxiB = NA, NmaxiSE = NA, NmaxiP = NA,
                                 EminiB = NA, EminiSE = NA, EminiP = NA, 
                                 EmaxiB = NA, EmaxiSE = NA, EmaxiP = NA))
    for (score in dietScores) { #add each diet score (+longitudinal) row
      ap8 <- rbind(ap8, data.frame(PFT = paste0(cohort, "_", pft, "_", score), 
                                   NminiB = NA, NminiSE = NA, NminiP = NA, 
                                   NmaxiB = NA, NmaxiSE = NA, NmaxiP = NA,
                                   EminiB = NA, EminiSE = NA, EminiP = NA, 
                                   EmaxiB = NA, EmaxiSE = NA, EmaxiP = NA))
      ap8 <- rbind(ap8, data.frame(PFT = paste0(cohort, "_", pft, "_", paste0(score, "*TIME")), 
                                   NminiB = NA, NminiSE = NA, NminiP = NA, 
                                   NmaxiB = NA, NmaxiSE = NA, NmaxiP = NA,
                                   EminiB = NA, EminiSE = NA, EminiP = NA, 
                                   EmaxiB = NA, EmaxiSE = NA, EmaxiP = NA))
    }
  }
}
models <- c("Nmini", "Nmaxi", "Emini", "Emaxi")
for (model in models) { #fill in table with results
  listName <- paste(model, "List", sep = "")
  ap8 <- filltab(ap8, ap8List[[listName]], model)
}
ap8 <- nVals(ap8, ap8List, prefixes = c("N", "E"))
rm(cohortGroups, HABCfinal_Male, HABCfinal_Male_Ever, HABCfinal_Male_Never, RASfinal_Ever, RASfinal_Never, rez, cohort, cohortNames, listName, model, models, pft, prefix, smk, score)
#formatting table
ap8 <- rowNames(ap8, "PFT")








###APPENDIX 9 - dietary pattern by diet quality (healthy vs unhealthy)###
HABCfinal_Healthy <- subset(HABCfinal, HEIhalf == "Healthy")
HABCfinal_Unhealthy <- subset(HABCfinal, HEIhalf == "Unhealthy")
RASfinal_Healthy <- subset(RASfinal, HEIhalf == "Healthy")
RASfinal_Unhealthy <- subset(RASfinal, HEIhalf == "Unhealthy")
cohortGroups <- list(
  Healthy = c("HABCfinal_Healthy", "RASfinal_Healthy"),
  Unhealthy = c("HABCfinal_Unhealthy", "RASfinal_Unhealthy")
)
ap9List <- list()
for (hei in names(cohortGroups)) { #run model function to get results
  prefix <- substr(hei, 1, 1) #prefix for diet quality groups (H, U)
  rez <- getResults(cohortGroups[[hei]])
  ap9List[[paste0(prefix, "miniList")]] <- rez$miniList
  ap9List[[paste0(prefix, "maxiList")]] <- rez$maxiList
}
ap9 <- data.frame(PFT = character(),
                  HminiB = character(), HminiSE = character(), HminiP = character(),
                  HmaxiB = character(), HmaxiSE = character(), HmaxiP = character(),
                  UminiB = character(), UminiSE = character(), UminiP = character(),
                  UmaxiB = character(), UmaxiSE = character(), UmaxiP = character(), stringsAsFactors = FALSE) #empty table
cohortNames <- c("HABC", "RAS", "META")
for (pft in PFTs) { #add row for PFT
  ap9 <- rbind(ap9, data.frame(PFT = pft,
                               HminiB = NA, HminiSE = NA, HminiP = NA, 
                               HmaxiB = NA, HmaxiSE = NA, HmaxiP = NA,
                               UminiB = NA, UminiSE = NA, UminiP = NA, 
                               UmaxiB = NA, UmaxiSE = NA, UmaxiP = NA))
  for (cohort in cohortNames) { #add row for cohort
    ap9 <- rbind(ap9, data.frame(PFT = paste0(cohort, "_", pft),
                                 HminiB = NA, HminiSE = NA, HminiP = NA, 
                                 HmaxiB = NA, HmaxiSE = NA, HmaxiP = NA,
                                 UminiB = NA, UminiSE = NA, UminiP = NA, 
                                 UmaxiB = NA, UmaxiSE = NA, UmaxiP = NA))
    for (score in dietScores) { #add each diet score (+longitudinal) row
      ap9 <- rbind(ap9, data.frame(PFT = paste0(cohort, "_", pft, "_", score), 
                                   HminiB = NA, HminiSE = NA, HminiP = NA, 
                                   HmaxiB = NA, HmaxiSE = NA, HmaxiP = NA,
                                   UminiB = NA, UminiSE = NA, UminiP = NA, 
                                   UmaxiB = NA, UmaxiSE = NA, UmaxiP = NA))
      ap9 <- rbind(ap9, data.frame(PFT = paste0(cohort, "_", pft, "_", paste0(score, "*TIME")), 
                                   HminiB = NA, HminiSE = NA, HminiP = NA, 
                                   HmaxiB = NA, HmaxiSE = NA, HmaxiP = NA,
                                   UminiB = NA, UminiSE = NA, UminiP = NA, 
                                   UmaxiB = NA, UmaxiSE = NA, UmaxiP = NA))
    }
  }
}
models <- c("Hmini", "Hmaxi", "Umini", "Umaxi")
for (model in models) { #fill in table with results
  listName <- paste(model, "List", sep = "")
  ap9 <- filltab(ap9, ap9List[[listName]], model)
}
ap9 <- nVals(ap9, ap9List, prefixes = c("H", "U"))
rm(cohortGroups, HABCfinal_Healthy, HABCfinal_Unhealthy, RASfinal_Healthy, RASfinal_Unhealthy, rez, cohort, cohortNames, listName, model, models, pft, prefix, hei, score)
#formatting table
ap9 <- rowNames(ap9, "PFT")





###APPENDIX 10 - dietary pattern component scores & FEV1###
#define component scores
componentScores <- c("totalfruitscore", "wholefruitscore", "totalvegscore", "greensandbeansscore", "wholegrainsscore",
                     "dairyscore", "totalproteinscore", "seafoodandplantscore", "refinedgrainsscore", "fattyacidscore",
                     "sodiumscore", "addsugarscore", "satfatscore", "dashfruitscore", "dashvegscore", "dashnutslegscore", "dashwholegrainsscore",
                     "dashdairyscore", "dashredprocscore", "dashssbscore", "dashsodscore")
models <- c("mini", "maxi")
cohorts <- c("HABCfinal", "RASfinal")
miniList <- list() #create mini & maxi lists
maxiList <- list() 
for(cohort in cohorts) { #run mini & maxi models for all PFTs, diet scores, both cohorts
  for(componentScore in componentScores) {
    for(model in models) { #dynamically call either mini or maxi function
      cat("Running", model, "model:", cohort, ",", "FEV1", ",", componentScore, "\n")
      result <- runModel(model, "FEV1", componentScore, cohort) #doesn't suppress errors or warnings, just messages (was cluttering output)
      if (model == "mini") {
        miniList[[result$save]] <- list(cross = result$cross, long = result$long)
      } else {
        maxiList[[result$save]] <- list(cross = result$cross, long = result$long)
      }
    }
  }
}
results <- list(miniList = miniList, maxiList = maxiList)
results <- RunMeta(results)

Padj <- function(tlist) { #apply FDR adjustment to p-values
  for (minmax in c("miniList", "maxiList")) {
    pvals <- c()  #list to store p-values
    locs <- list()  #list to store locations of p-values (they're in different lists)
    for (sub in names(tlist[[minmax]])) { #first collect p-values from the different results lists
      for (df in c("cross", "long")) {
        pvals <- c(pvals, tlist[[minmax]][[sub]][[df]][["pval"]])
        locs[[length(pvals)]] <- list(minmax, sub, df)
      }
    }
    adjusted_pvals <- p.adjust(pvals, method = "fdr") #adjust p-values using FDR method
    for (i in seq_along(adjusted_pvals)) { #replace original p-values with adjusted p-values
      location <- locs[[i]]
      tlist[[location[[1]]]][[location[[2]]]][[location[[3]]]][["pval"]] <- adjusted_pvals[i]
    }
  }
  return(tlist)
}
results <- Padj(results)
results <- rounding(results)
ap10List <- results
ap10 <- data.frame(PFT = character(), miniB = character(), miniSE = character(), miniP = character(),
                      maxiB = character(), maxiSE = character(), maxiP = character(), stringsAsFactors = FALSE) #empty table
cohortNames <- c("HABC", "RAS", "META")
for (pft in "FEV1") { #add row for PFT
  ap10 <- rbind(ap10, data.frame(PFT = pft,
                                 miniB = NA, miniSE = NA, miniP = NA, 
                                 maxiB = NA, maxiSE = NA, maxiP = NA))
  for (cohort in cohortNames) { #add row for cohort
    ap10 <- rbind(ap10, data.frame(PFT = paste0(cohort, "_", pft),
                                   miniB = NA, miniSE = NA, miniP = NA, 
                                   maxiB = NA, maxiSE = NA, maxiP = NA))
    for (score in componentScores) { #add each diet score (+longitudinal) row
      ap10 <- rbind(ap10, data.frame(PFT = paste0(cohort, "_", pft, "_", score), 
                                    miniB = NA, miniSE = NA, miniP = NA, 
                                    maxiB = NA, maxiSE = NA, maxiP = NA))
      ap10 <- rbind(ap10, data.frame(PFT = paste0(cohort, "_", pft, "_", paste0(score, "*TIME")), 
                                    miniB = NA, miniSE = NA, miniP = NA, 
                                    maxiB = NA, maxiSE = NA, maxiP = NA))
    }
  }
}
models <- c("mini", "maxi")
for (model in models) { #fill in table with results
  listName <- paste(model, "List", sep = "")
  ap10 <- filltab(ap10, ap10List[[listName]], model)
}
rm(maxiList, miniList, cohort, cohortNames, componentScore, componentScores, model, models, pft, results, score, listName)
#formatting table
ap10 <- rowNames(ap10, "PFT")
ap10 <- ap10[grepl("*TIME$", ap10$PFT) | #only display longitudinal results
                grepl("^(FEV1|FVC|RATIO)$", ap10$PFT) | 
                grepl("^(HABC|RAS|META)$", ap10$PFT), ]






###APPENDIX 11 - dietary pattern & COPD###
#set COPD datasets (exclude all observations where FEV1 & FVC are missing or not meeting quality control criteria)
HABCcopd <- HABCfinal %>% filter((!is.na(FEV1) & !is.na(FVC)) & (QCFEV1 %in% c(3, 4) & QCFVC %in% c(3, 4)))
RAScopd <- RASfinal %>% filter((!is.na(FEV1) & !is.na(FVC)) & (QCFEV1 %in% c(3, 4) & QCFVC %in% c(3, 4)))
HABCcopdp <- HABCcopd %>% group_by(ID)
RAScopdp <- RAScopd %>% group_by(ID)
HABCcopdi <- HABCcopd %>% #exclude all baseline COPD = 1 (use this dataset for survival analysis)
  group_by(ID) %>% filter(COPDprev == 0) %>% ungroup()
RAScopdi <- RAScopd %>%
  group_by(ID) %>% filter(COPDprev == 0) %>% ungroup()
#define functinos for COPD analyses
COPDmodel <- function(model = "mini", diet, df, analysis, race = FALSE) {
  dat <- get(df, envir = .GlobalEnv) #follows similar function to 'runModel' function, but updated for COPD analyses
  if (model == "mini") {
    nParticipants <- length(unique(dat$ID)) # for the mini model
    nEve <- if (analysis == "prev") length(unique(dat$ID[dat$COPDany == 1])) else NA
  } else if (model == "maxi") {
    # Subset the data where SMKCAT & PACKYR1 are not missing
    datMaxi <- subset(dat, !is.na(SMKCAT) & !is.na(PACKYR1))
    nParticipants <- length(unique(datMaxi$ID)) # for the maxi model
    nEve <- if (analysis == "prev") length(unique(datMaxi$ID[datMaxi$COPDany == 1])) else NA
  } 
  cohort <- if (startsWith(df, "RAS")) {
    "RAS"
  } else if (startsWith(df, "HABC")) {
    "HABC"
  }
  base <- paste("AGE + AGE2 + HTCM + HT2",
                if(!race) {"+ RACE"} else {""},
                if(!(cohort == "RAS" || grepl("Male|Female", df))) {"+ GENDER"}, "+", diet)
  full <- if(model == "maxi" && !grepl("Never", df)) {"+ SITE + SMKCAT + PACKYR1"} else if(model == "maxi" && grepl("Never", df)) {"+ SITE"} else {""}
  if (analysis == "prev") {
    form <- paste("COPDany ~", base, full)
    mp <- glm(as.formula(form), data = dat, family = binomial)
    prevrez <- c(Name = diet, OR = exp(summary(mp)$coef[diet, "Estimate"]), CI_L = exp(confint(mp)[diet, 1]), CI_H = exp(confint(mp)[diet, 2]), pval = summary(mp)$coef[diet, "Pr(>|z|)"],  N = nParticipants, Events = nEve)
    save <- paste0(cohort, "_COPDp_", diet)
    results <- list(save = save, prev = prevrez)
  } else if (analysis == "inc") {
    form <- paste("Surv(COPDTIME, COPD) ~", base, full)
    mi <- coxph(as.formula(form), data = dat)
    nEve <- summary(mi)$nevent
    increz <- c(Name = diet, HR = summary(mi)$coefficients[diet, "exp(coef)"], CI_L = exp(confint(mi)[diet, 1]), CI_H = exp(confint(mi)[diet, 2]), pval = summary(mi)$coefficients[diet, "Pr(>|z|)"],  N = nParticipants, Events = nEve)
    save <- paste0(cohort, "_COPDi_", diet)
    results <- list(save = save, inc = increz)
  }
  return(results)
}
COPDmodel <- function(model = "mini", diet, df, analysis, race = FALSE) {
  dat <- get(df, envir = .GlobalEnv) #follows similar function to 'runModel' function, but updated for COPD analyses
  if (model == "mini") {
    nParticipants <- length(unique(dat$ID)) # for the mini model
    nEve <- if (analysis == "prev") length(unique(dat$ID[dat$COPDany == 1])) else NA
  } else if (model == "maxi") {
    # Subset the data where SMKCAT & PACKYR1 are not missing
    datMaxi <- subset(dat, !is.na(SMKCAT) & !is.na(PACKYR1))
    nParticipants <- length(unique(datMaxi$ID)) # for the maxi model
    nEve <- if (analysis == "prev") length(unique(datMaxi$ID[datMaxi$COPDany == 1])) else NA
  } 
  cohort <- if (startsWith(df, "RAS")) {
    "RAS"
  } else if (startsWith(df, "HABC")) {
    "HABC"
  }
  base <- paste("AGE + AGE2 + HTCM + HT2",
                if(!race) {"+ RACE"} else {""},
                if(!(cohort == "RAS" || grepl("Male|Female", df))) {"+ GENDER"}, "+", diet)
  full <- if(model == "maxi" && !grepl("Never", df)) {"+ SITE + SMKCAT + PACKYR1"} else if(model == "maxi" && grepl("Never", df)) {"+ SITE"} else {""}
  basei <- paste("AGE + tt(AGE) + AGE2 + tt(AGE2) + HTCM + tt(HTCM) + HT2 + tt(HT2)",
                if(!race) {"+ RACE"} else {""},
                if(!(cohort == "RAS" || grepl("Male|Female", df))) {"+ strata(GENDER)"}, "+", diet)
  fulli <- if(model == "maxi" && !grepl("Never", df)) {"+ strata(SITE) + strata(SMKCAT) + tt(PACKYR1)"} else if(model == "maxi" && grepl("Never", df)) {"+ strata(SITE)"} else {""}
  if (analysis == "prev") {
    form <- paste("COPDany ~", base, full)
    mp <- glm(as.formula(form), data = dat, family = binomial)
    prevrez <- c(Name = diet, OR = exp(summary(mp)$coef[diet, "Estimate"]), CI_L = exp(confint(mp)[diet, 1]), CI_H = exp(confint(mp)[diet, 2]), pval = summary(mp)$coef[diet, "Pr(>|z|)"],  N = nParticipants, Events = nEve)
    save <- paste0(cohort, "_COPDp_", diet)
    results <- list(save = save, prev = prevrez)
  } else if (analysis == "inc") {
    form <- paste("Surv(COPDTIME, COPD) ~", basei, fulli)
    mi <- coxph(as.formula(form), data = dat, tt = function(x,t,...) x*t)
    nEve <- summary(mi)$nevent
    increz <- c(Name = diet, HR = summary(mi)$coefficients[diet, "exp(coef)"], CI_L = exp(confint(mi)[diet, 1]), CI_H = exp(confint(mi)[diet, 2]), pval = summary(mi)$coefficients[diet, "Pr(>|z|)"],  N = nParticipants, Events = nEve)
    save <- paste0(cohort, "_COPDi_", diet)
    results <- list(save = save, inc = increz)
  }
  return(results)
}
copdMeta <- function(ListX) {
  for(analysisType in names(ListX)) { #'prev' or 'inc'
    for(modelType in names(ListX[[analysisType]])) { #'mini' or 'maxi'
      currentList <- ListX[[analysisType]][[modelType]]
      metaResults <- list()
      for (name in names(currentList)) {
        if (grepl("HABC", name)) {
          counterName <- gsub("HABC", "RAS", name)
          if (counterName %in% names(currentList)) { #extract data
            effectSizeKey <- ifelse(analysisType == "prev", "OR", "HR")
            effectSize <- as.numeric(currentList[[name]][[effectSizeKey]])
            CI_L <- as.numeric(currentList[[name]][["CI_L"]])
            CI_H <- as.numeric(currentList[[name]][["CI_H"]])
            effectSizeCounter <- as.numeric(currentList[[counterName]][[effectSizeKey]])
            CI_L_Counter <- as.numeric(currentList[[counterName]][["CI_L"]])
            CI_H_Counter <- as.numeric(currentList[[counterName]][["CI_H"]])
            combinedEffectSizes <- c(effectSize, effectSizeCounter) #combine HABC/RAS 
            combinedCIs_L <- c(CI_L, CI_L_Counter)
            combinedCIs_H <- c(CI_H, CI_H_Counter)
            combinedNs <- sum(as.numeric(currentList[[name]][["N"]]), as.numeric(currentList[[counterName]][["N"]]))
            combinedEvents <- sum(as.numeric(currentList[[name]][["Events"]]), as.numeric(currentList[[counterName]][["Events"]]))
            metaResults[[name]] <- copdMetaFunk(combinedEffectSizes, combinedCIs_L, combinedCIs_H, analysisType, combinedNs, combinedEvents)
          }
        }
      } #rename META resutls
      metaResults <- setNames(metaResults, gsub("HABC", "META", names(metaResults)))
      ListX[[analysisType]][[modelType]] <- c(currentList, metaResults)
    }
  }
  return(ListX)
}
copdMetaFunk <- function(effectSizes, CIs_L, CIs_H, analysisType, N, Events) {
  log_effectSizes <- log(effectSizes) # calculate the log of the ORs/HRs and their variance
  var_log_effectSizes <- ((log(CIs_H) - log(CIs_L)) / (2 * qnorm(0.975)))^2
  meta_result <- rma(yi = log_effectSizes, vi = var_log_effectSizes, method = "REML")
  meta_effectSize <- exp(meta_result$b[1]) # convert results back to OR/HR scale and extract value
  lower_CI <- exp(meta_result$ci.lb[1])
  upper_CI <- exp(meta_result$ci.ub[1])
  pval <- meta_result$pval
  i2 <- meta_result$I2
  Hetpval = meta_result$QEp
    effectSizeLabel <- ifelse(analysisType == "prev", "OR", "HR") #"OR" or "HR" depending on analysis type
  meta_rez <- setNames(as.character(c(meta_effectSize, lower_CI, upper_CI, pval, i2, Hetpval, N, Events)), 
                       c(effectSizeLabel, "CI_L", "CI_H", "pval", "i2", "Hetpval", "N", "Events")) #format as characters
  return(meta_rez)
}
roundingCOPD <- function(tlist) {
  for (model in c("prev", "inc")) {
    for (minmax in c("mini", "maxi")) {
      for (sub in names(tlist[[model]][[minmax]])) {
        itemList <- tlist[[model]][[minmax]][[sub]]  #get the sub-list/named character vector
        keys <- c("OR", "HR", "CI_L", "CI_H", "pval")
        for (key in keys) {
          if (key %in% names(itemList)) {  #check if 'key' exists (will skip OR/HR for analyses without one of those - prev/inc)
            tlist[[model]][[minmax]][[sub]][[key]] <- digi(itemList[key], digits = 3, p_value = (key == "pval"))
          }
        }
        if (startsWith(sub, "META")) { #only round i2 & Hetpval for meta-analyses
            tlist[[model]][[minmax]][[sub]][["i2"]] <- digi(itemList["i2"], is_i2 = TRUE)
            tlist[[model]][[minmax]][[sub]][["Hetpval"]] <- digi(itemList["Hetpval"], p_value = TRUE)
        }
      }
    }
  }
  return(tlist)
}
COPDrez <- function(cohorts, race = FALSE) {
  models <- c("mini", "maxi")
  dietScores <- c("HEI2020score", "DASHscore")
  resultsList <- list(prev = list(), inc = list()) # Initialize lists for 'prev' and 'inc' analyses
  for(cohort in cohorts) {
    if (grepl("copdp", cohort)) {
      analysis <- "prev"
    } else if (grepl("copdi", cohort)) {
      analysis <- "inc"
    }
    for(dietScore in dietScores) {
      for(model in models) {
        cat("Running", analysis, "analysis with", model, "model:", cohort, ",", dietScore, "\n")
        result <- COPDmodel(model, dietScore, cohort, analysis, race)
        key <- paste(cohort, dietScore, model, sep = "_")
        if (analysis == "prev") {
          if (!is.list(resultsList$prev[[model]])) {
            resultsList$prev[[model]] <- list()
          }
          resultsList$prev[[model]][[key]] <- result$prev
        } else if (analysis == "inc") {
          if (!is.list(resultsList$inc[[model]])) {
            resultsList$inc[[model]] <- list()
          }
          resultsList$inc[[model]][[key]] <- result$inc
        }
      }
    }
  }
  resultsList <- copdMeta(resultsList)
  resultsList <- roundingCOPD(resultsList)
  return(resultsList)
}
cohorts <- c("HABCcopdp", "HABCcopdi", "RAScopdp", "RAScopdi")
ap11List <- COPDrez(cohorts)
ap11 <- data.frame(COPD = character(), miniB = character(), miniSE = character(), miniP = character(),
                     maxiB = character(), maxiSE = character(), maxiP = character(), stringsAsFactors = FALSE) #empty table
cohortNames <- c("HABC", "RAS", "META")
previncs <- c("prev", "inc")
for (previnc in previncs) { #add row for PFT
  ap11 <- rbind(ap11, data.frame(COPD = previnc, mini = NA, maxi = NA))
  for (cohort in cohortNames) { #add row for cohort
    ap11 <- rbind(ap11, data.frame(COPD = paste0(cohort, "_", previnc), mini = NA, maxi = NA))
    for (score in dietScores) { #add each diet score (+longitudinal) row
      ap11 <- rbind(ap11, data.frame(COPD = paste0(cohort, "_", previnc, "_", score), mini = NA, maxi = NA))
    }
  }
}
extract_result <- function(cohort, score, type, period) { #function to find matching name from list & extract results
  base_name <- paste0(cohort, "copd", ifelse(period == "prev", "p_", "i_"), score, "_", type)
  result_list <- ap11List[[period]][[type]][[base_name]] #extract matching part of list
  OR_HR_value <- ifelse(period == "prev", result_list["OR"], result_list["HR"])
  return(paste(OR_HR_value, " (", result_list["CI_L"], " - ", result_list["CI_H"], ")", sep = ""))   
}
for (i in 1:nrow(ap11)) { #fill table
  elements <- unlist(strsplit(as.character(ap11$COPD[i]), "_"))
  cohort <- elements[1] #find matching name in table
  score <- ifelse(length(elements) >= 3, elements[3], NA)
  period <- ifelse(grepl("inc", ap11$COPD[i]), "inc", "prev")
  if (!(score %in% c("HEI2020score", "DASHscore"))) {next} #leave blank if not score (empty row)
  ap11$mini[i] <- extract_result(cohort, score, "mini", period) #fill values
  ap11$maxi[i] <- extract_result(cohort, score, "maxi", period)
}
rm(result, cohort, cohortNames, cohorts, elements, extract_result, i, period, previnc, previncs, score)
#formatting table
ap11 <- rowNames(ap11, "COPD")



integrateNEventsCorrectly <- function(resultsList, df) {
  # Helper function to extract N and Events
  getNEvents <- function(cohort, analysisType, modelType) {
    # Construct the key to get the right value
    key <- paste0(cohort, "copd", substr(analysisType, 1, 1), "_HEI2020score_", modelType)
    n <- resultsList[[analysisType]][[modelType]][[key]]["N"]
    events <- resultsList[[analysisType]][[modelType]][[key]]["Events"]
    return(paste("Events =", events, "(N =", n,")"))
  }
  
  # Define cohorts and their positions in the dataframe based on the screenshot
  cohorts <- c("HABC", "RAS", "META")
  positions <- c("prev" = 2, "inc" = 3) # Rows where the cohorts' values should be placed
  
  # Loop through the dataframe rows
  for (cohort in cohorts) {
    for (analysisType in c("prev", "inc")) {
      for (modelType in c("mini", "maxi")) {
        # Get the correct N and Events values
        nEventsValue <- getNEvents(cohort, analysisType, modelType)
        
        # Calculate the correct row based on cohort and analysisType
        rowName <- paste(cohort, analysisType, sep = "_")
        rowIndex <- which(df$COPD == rowName)
        
        # Check and insert the values into the correct cell if the row exists
        if (length(rowIndex) == 1) {
          df[rowIndex, modelType] <- nEventsValue
        }
      }
    }
  }
  
  return(df)
}
ap11 <- integrateNEventsCorrectly(ap11List, ap11)









###APPENDIX 4 - heterogeneity considerations for meta-analyses###
#needs to be run last so that it can incorporate all of the values from the analyses
#create blank table for appendix 4 with rows for all analyses
ap4 <- data.frame(Analysis = character(),
                  miniFEV1 = character(), miniFVC = character(), miniRATIO = character(),
                  maxiFEV1 = character(), maxiFVC = character(), maxiRATIO = character(), stringsAsFactors = FALSE) #empty table
analyses <- c("unstrat", "sex", "race", "smoking status", "diet quality", "copd")
stratGroups <- list(unstrat = NULL, sex = NULL, race = c("white", "black"),
                    `smoking status` = c("never", "ever"), `diet quality` = c("healthy", "unhealthy"), copd = c("prevalence", "incidence"))
for (analysis in analyses) { #make rows to fill with heterogeneity values
  if (analysis == "unstrat" | analysis == "sex") { #if unstratified, only add diet scores (not row for strat)
    ap4 <- rbind(ap4, data.frame(Analysis = analysis,
                                 miniFEV1 = NA, miniFVC = NA, miniRATIO = NA, 
                                 maxiFEV1 = NA, maxiFVC = NA, maxiRATIO = NA))
    for (score in dietScores) {
      ap4 <- rbind(ap4, data.frame(Analysis = paste0(analysis, "_", score), 
                                   miniFEV1 = NA, miniFVC = NA, miniRATIO = NA, 
                                   maxiFEV1 = NA, maxiFVC = NA, maxiRATIO = NA))
      ap4 <- rbind(ap4, data.frame(Analysis = paste0(analysis, "_", paste0(score, "*TIME")), 
                                   miniFEV1 = NA, miniFVC = NA, miniRATIO = NA, 
                                   maxiFEV1 = NA, maxiFVC = NA, maxiRATIO = NA))
    }
  } else {
    ap4 <- rbind(ap4, data.frame(Analysis = analysis,
                                 miniFEV1 = NA, miniFVC = NA, miniRATIO = NA, 
                                 maxiFEV1 = NA, maxiFVC = NA, maxiRATIO = NA))
    if (!is.null(stratGroups[[analysis]])) {
      for (strat in stratGroups[[analysis]]) { #add row for stratification
        ap4 <- rbind(ap4, data.frame(Analysis = paste0(analysis, "_", strat),
                                     miniFEV1 = NA, miniFVC = NA, miniRATIO = NA, 
                                     maxiFEV1 = NA, maxiFVC = NA, maxiRATIO = NA))
        for (score in dietScores) { #add each diet score (+longitudinal) row
          ap4 <- rbind(ap4, data.frame(Analysis = paste0(analysis, "_", strat, "_", score), 
                                       miniFEV1 = NA, miniFVC = NA, miniRATIO = NA, 
                                       maxiFEV1 = NA, maxiFVC = NA, maxiRATIO = NA))
          ap4 <- rbind(ap4, data.frame(Analysis = paste0(analysis, "_", strat, "_", paste0(score, "*TIME")), 
                                       miniFEV1 = NA, miniFVC = NA, miniRATIO = NA, 
                                       maxiFEV1 = NA, maxiFVC = NA, maxiRATIO = NA))
        }
      }
    }
  }
}
populate <- function(ap4, miniList, maxiList, analysis, score, pft) {
  #access the i2 and i2pval values directly from the vectors
  mini_cross_i2 <- miniList[[paste0("META_", pft, "_", score)]]$cross['i2']
  mini_cross_Hetpval <- miniList[[paste0("META_", pft, "_", score)]]$cross['Hetpval']
  mini_long_i2 <- miniList[[paste0("META_", pft, "_", score)]]$long['i2']
  mini_long_Hetpval <- miniList[[paste0("META_", pft, "_", score)]]$long['Hetpval']
  maxi_cross_i2 <- maxiList[[paste0("META_", pft, "_", score)]]$cross['i2']
  maxi_cross_Hetpval <- maxiList[[paste0("META_", pft, "_", score)]]$cross['Hetpval']
  maxi_long_i2 <- maxiList[[paste0("META_", pft, "_", score)]]$long['i2']
  maxi_long_Hetpval <- maxiList[[paste0("META_", pft, "_", score)]]$long['Hetpval']
  #format the values: "i2 (Hetpval)"
  mini_cross_value <- paste0(mini_cross_i2, " (", mini_cross_Hetpval, ")")
  mini_long_value <- paste0(mini_long_i2, " (", mini_long_Hetpval, ")")
  maxi_cross_value <- paste0(maxi_cross_i2, " (", maxi_cross_Hetpval, ")")
  maxi_long_value <- paste0(maxi_long_i2, " (", maxi_long_Hetpval, ")")
  #match to corresponding row in ap4 & fill
  ap4[ap4$Analysis == paste0(analysis, "_", score), paste0("mini", pft)] <- mini_cross_value
  ap4[ap4$Analysis == paste0(analysis, "_", score, "*TIME"), paste0("mini", pft)] <- mini_long_value
  ap4[ap4$Analysis == paste0(analysis, "_", score), paste0("maxi", pft)] <- maxi_cross_value
  ap4[ap4$Analysis == paste0(analysis, "_", score, "*TIME"), paste0("maxi", pft)] <- maxi_long_value
  return(ap4)
}
ap4strats <- function(ap4, analysis) {
  listName <- switch(analysis,
                     "unstrat" = list(miniList = t3List$miniList, maxiList = t3List$maxiList),
                     "sex" = list(miniList = ap7List$MminiList, maxiList = ap7List$MmaxiList),
                     "race_white" = list(miniList = ap5List$WminiList, maxiList = ap5List$WmaxiList),
                     "race_black" = list(miniList = ap5List$BminiList, maxiList = ap5List$BmaxiList),
                     "smoking status_never" = list(miniList = ap8List$NminiList, maxiList = ap8List$NmaxiList),
                     "smoking status_ever" = list(miniList = ap8List$EminiList, maxiList = ap8List$EmaxiList),
                     "diet quality_healthy" = list(miniList = ap9List$HminiList, maxiList = ap9List$HmaxiList),
                     "diet quality_unhealthy" = list(miniList = ap9List$UminiList, maxiList = ap9List$UmaxiList),
                     "copd_prevalence" = list(miniList = ap11List$pminiList, maxiList = ap11List$pmaxiList),
                     "copd_incidence" = list(miniList = ap11List$iminiList, maxiList = ap11List$imaxiList)
  )
  for (score in dietScores) {
    for (pft in PFTs) {
      ap4 <- populate(ap4, listName$miniList, listName$maxiList, analysis, score, pft)
    }
  }
  
  return(ap4)
}
stratifiedAnalyses <- c("unstrat", "sex", "race_white", "race_black", "smoking status_never", "smoking status_ever",
                        "diet quality_healthy", "diet quality_unhealthy")
for (analysis in stratifiedAnalyses) {
  ap4 <- ap4strats(ap4, analysis)
}
#formatting table
ap4 <- rowNames(ap4, "Analysis")



########################EXPORT TABLES TO EXCEL#####################

tables <- list("Table 1" = table1, "Table 2" = table2, "Table 3" = table3,
               "Appendix 3" = ap3, "Appendix 4" = ap4, "Appendix 5" = ap5, 
               "Appendix 6" = ap6, "Appendix 7" = ap7, "Appendix 8" = ap8, 
               "Appendix 9" = ap9, "Appendix 10" = ap10, "Appendix 11" = ap11)
wb <- createWorkbook() #create workbook
for (sheet_name in names(tables)) { #loop through tables and add each as its own sheet
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet = sheet_name, x = tables[[sheet_name]])
}
filename <- "//rschfs1x/userrs/R-Z/wav25_RS/Desktop/updated code/tables_workbook.xlsx"
saveWorkbook(wb, file = filename, overwrite = TRUE)




###FIGURES; Forest plots###

plot3 <- table3 %>% select(-miniB, -miniSE, -miniP, -maxiP)
plot3cross <- plot3[c(12, 13, 15, 28, 29, 31), ]
plot6 <- ap6 %>% select(-NminiB, -NminiSE, -NminiP, -NmaxiP, -EminiB, -EminiSE, -EminiP, -EmaxiP)
plot6cross <- plot6[c(12, 13, 15, 28, 29, 31), ]
nevercross <- plot6cross[, c(1:3)]
colnames(nevercross) <- c("PFT", "maxiB", "maxiSE")
evercross <- plot6cross[, c(1, 4, 5)]
colnames(evercross) <- c("PFT", "maxiB", "maxiSE")


fev1cross <- plot3cross[c(1:3), ]
fev1cross <- rbind(fev1cross, evercross[c(1:3), ])
fev1cross <- rbind(fev1cross, nevercross[c(1:3), ])

fvccross <- plot3cross[c(4:6), ]
fvccross <- rbind(fvccross, evercross[c(4:6), ])
fvccross <- rbind(fvccross, nevercross[c(4:6), ])

plot2 <- rbind(fev1cross, fvccross)


plot3long <- plot3[c(2,4,6,7,9,11,12,14,16),]
plot5 <- ap5 %>% select(-WminiB, -WminiSE, -WminiP, -WmaxiP, -BminiB, -BminiSE, -BminiP, -BmaxiP)
plot5long <- plot5[c(2,4,6,7,9,11,12,14,16),]
whitelong <- plot5long[, c(1:3)]
colnames(whitelong) <- c("PFT", "maxiB", "maxiSE")
blacklong <- plot5long[, c(1,4,5)]
colnames(blacklong) <- c("PFT", "maxiB", "maxiSE")

fev1long <- plot3long
fev1long <- rbind(fev1long, whitelong)
fev1long <- rbind(fev1long, blacklong)
plot3 <- fev1long

plot3 <- plot3 %>%
  mutate(maxiText = NA_character_)


plot3 <- plot3 %>%
  mutate(PFT = case_when(
    PFT == "HABC" ~ "Health ABC",
    PFT == "META" ~ "Meta-analysis",
    TRUE ~ PFT  # Keep all other values unchanged
  ))

rows_to_update <- plot3$PFT %in% c("Health ABC", "RAS", "Meta-analysis")
plot3$PFT[rows_to_update] <- paste0(plot3$PFT[rows_to_update], " (", plot3$maxiSE[rows_to_update], ")")
plot3$maxiSE[rows_to_update] <- NA

plot3 <- plot3 %>%
  mutate(PFT = case_when(
    PFT == "DASH * Time" ~ "DASH",
    PFT == "HEI2020 * Time" ~ "HEI-2020",
    TRUE ~ PFT  # Keep all other values unchanged
  ))


cols_to_convert <- c("maxiB", "maxiSE")
plot3[cols_to_convert] <- lapply(plot3[cols_to_convert], as.numeric)

std_devs <- c("HEI-2020" = 10, "DASH" = 5)


# Function to scale the values
scale_values <- function(value, pattern, std_devs) {
  if (str_detect(pattern, "HEI-2020")) {
    return(value * std_devs["HEI-2020"])
  } else if (str_detect(pattern, "DASH")) {
    return(value * std_devs["DASH"])
  } else {
    return(value)
  }
}

plot3 <- plot3 %>%
  rowwise() %>%
  mutate(
    maxiB = scale_values(maxiB, PFT, std_devs),
    maxiSE = scale_values(maxiSE, PFT, std_devs)
  ) %>%
  ungroup()


rows_to_convert <- plot3$PFT %in% c("HEI-2020", "DASH")

plot3 <- plot3 %>%
  mutate(
    maxiCIL = if_else(rows_to_convert, maxiB - 1.96 * maxiSE, NA_real_),
    maxiCIH = if_else(rows_to_convert, maxiB + 1.96 * maxiSE, NA_real_)
  )

plot3[] <- lapply(plot3, function(x) { #round to 1 dp
  if (is.numeric(x)) {
    round(x, digits = 1)
  } else {
    x  # Skip non-numeric values
  }
})

# Step 4: Prepare text labels for the forest plot
plot3 <- plot3 %>%
  mutate(
    significant = ifelse((maxiCIL > 0 & maxiCIH > 0) | (maxiCIL < 0 & maxiCIH < 0), "*", ""),
    maxiText = if_else(rows_to_convert,
                       paste0(maxiB, " (", maxiCIL, ", ", maxiCIH, ")", significant),
                       maxiText)
  )


# Insert a blank row in the middle to create space between cohorts
plot3 <- rbind(plot3[1:9, ], NA, plot3[10:18, ], NA, plot3[19:nrow(plot3), ])



# Calculate a proportional box size based on the standard error or other measure
proportional_boxsize <- 0.5 / plot3$maxiSE  # Example: inverse of standard error
proportional_boxsize[is.na(proportional_boxsize)] <- 0.1 # Handle NA values or set a minimum size

# Remove NA values and adjust formatting
plot3$maxiText <- ifelse(is.na(plot3$maxiB) | is.na(plot3$maxiCIL) | is.na(plot3$maxiCIH),
                         "",
                         sprintf("%.1f (%.1f, %.1f)%s", plot3$maxiB, plot3$maxiCIL, plot3$maxiCIH, plot3$significant))

plot3$PFT <- ifelse(plot3$PFT %in% c("HEI-2020", "DASH"), 
                    paste0("\t\t", plot3$PFT), # Use a long string of spaces
                    plot3$PFT)


forestplot(
  labeltext = list(plot3$PFT, plot3$maxiText),
  mean = cbind(plot3$maxiB),
  lower = cbind(plot3$maxiCIL),
  upper = cbind(plot3$maxiCIH),
  xlab = expression("Longitudinal change in " * FEV[1] * " (mL/year)"),
  graph.pos = 2,
  is.summary = c(rep(FALSE, nrow(plot2))),
  clip = c(-5, 10), # Adjust the clipping based on your actual data
  col = fpColors(box = "black"),
  lineheight = grid::unit(1.2, "lines"), # Adjust line height for better spacing
  xticks = c(-5, 0, 5, 10), # Explicitly set the tick marks to include 0
  boxsize = proportional_boxsize, # Use proportional box sizes
  align = c("l", "l"), # Align both the text labels and numbers to the left
  txt_gp = fpTxtGp(
    label = gpar(cex = 1, fontfamily = "Times"),   # Apply font to labels
    ticks = gpar(cex = 1, fontfamily = "Times"),   # Apply font to tick labels
    xlab = gpar(cex = 1, fontfamily = "Times"),    # Apply font to xlab
    title = gpar(cex = 1, fontfamily = "Times")    # Apply font to title (if any)
  )
)



# Adjust FEV1 and FVC labels
grid.text("All participants", x = unit(0.16, "npc"), y = unit(0.975, "npc"), just = "left", gp = gpar(cex = 1.1, fontfamily = "Times"))
grid.text("White participants", x = unit(0.16, "npc"), y = unit(0.685, "npc"), just = "left", gp = gpar(cex = 1.1, fontfamily = "Times"))
grid.text("Black participants", x = unit(0.16, "npc"), y = unit(0.395, "npc"), just = "left", gp = gpar(cex = 1.1, fontfamily = "Times"))

grid.text(expression(beta * " (95% CI)"), x = unit(0.897, "npc"), y = unit(0.97, "npc"), just = "center", gp = gpar(cex = 1.1, fontfamily = "Times"))


# Capture the current plot
p <- recordPlot()

# Save the plot as a PDF with high resolution
png("plot3.png", width = 10, height = 8, units = "in", res = 300)
replayPlot(p) # Replay the captured plot
dev.off()























# Assuming df is your data frame
new_names <- c("All participants", "Ever smokers", "Never smokers", 
               "All participants", "Ever smokers", "Never smokers")
plot2$PFT[plot2$PFT == "META"] <- new_names
rows_to_update <- plot2$PFT %in% c("All participants", "Ever smokers", "Never smokers")
plot2$PFT[rows_to_update] <- paste0(plot2$PFT[rows_to_update], " (", plot2$maxiSE[rows_to_update], ")")
plot2$maxiSE[rows_to_update] <- NA





plot2 <- plot2 %>%
  mutate(
    maxiText = if_else(PFT %in% c("HABC", "RAS", "META"), as.character(maxiSE), NA_character_),
    maxiSE = if_else(PFT %in% c("HABC", "RAS", "META"), NA_character_, maxiSE)
  )


cols_to_convert <- c("maxiB", "maxiSE")
plot2[cols_to_convert] <- lapply(plot2[cols_to_convert], as.numeric)

std_devs <- c("HEI2020" = 10, "DASH" = 5)

# Function to scale the values
scale_values <- function(value, pattern, std_devs) {
  if (str_detect(pattern, "HEI2020")) {
    return(value * std_devs["HEI2020"])
  } else if (str_detect(pattern, "DASH")) {
    return(value * std_devs["DASH"])
  } else {
    return(value)
  }
}

plot2 <- plot2 %>%
  rowwise() %>%
  mutate(
    maxiB = scale_values(maxiB, PFT, std_devs),
    maxiSE = scale_values(maxiSE, PFT, std_devs)
  ) %>%
  ungroup()

# Identify rows where PFT matches the criteria
rows_to_convert <- plot2$PFT %in% c("HEI2020", "HEI2020 * Time", "DASH", "DASH * Time")

plot2 <- plot2 %>%
  mutate(
    maxiCIL = if_else(rows_to_convert, maxiB - 1.96 * maxiSE, NA_real_),
    maxiCIH = if_else(rows_to_convert, maxiB + 1.96 * maxiSE, NA_real_)
  )

plot2[] <- lapply(plot2, function(x) { #round to 1 dp
  if (is.numeric(x)) {
    round(x, digits = 1)
  } else {
    x  # Skip non-numeric values
  }
})

# Step 4: Prepare text labels for the forest plot
plot2 <- plot2 %>%
  mutate(
    significant = ifelse((maxiCIL > 0 & maxiCIH > 0) | (maxiCIL < 0 & maxiCIH < 0), "*", ""),
    maxiText = if_else(rows_to_convert,
                       paste0(maxiB, " (", maxiCIL, ", ", maxiCIH, ")", significant),
                       maxiText)
  )





# Insert a blank row in the middle to create space between FEV1 and FVC
plot2 <- rbind(plot2[1:9, ], NA, plot2[10:nrow(plot2), ])

# Calculate a proportional box size based on the standard error or other measure
proportional_boxsize <- 4 / plot2$maxiSE  # Example: inverse of standard error
proportional_boxsize[is.na(proportional_boxsize)] <- 0.1 # Handle NA values or set a minimum size

plot2$PFT <- gsub("HEI2020", "HEI-2020", plot2$PFT)

# Remove NA values and adjust formatting
plot2$maxiText <- ifelse(is.na(plot2$maxiB) | is.na(plot2$maxiCIL) | is.na(plot2$maxiCIH),
                         "",
                         sprintf("%.1f (%.1f, %.1f)%s", plot2$maxiB, plot2$maxiCIL, plot2$maxiCIH, plot2$significant))

plot2$PFT <- ifelse(plot2$PFT %in% c("HEI-2020", "DASH"), 
                    paste0("\t\t", plot2$PFT), # Use a long string of spaces
                    plot2$PFT)


forestplot(
  labeltext = list(plot2$PFT, plot2$maxiText),
  mean = cbind(plot2$maxiB),
  lower = cbind(plot2$maxiCIL),
  upper = cbind(plot2$maxiCIH),
  xlab = expression("Cross-sectional difference in " * FEV[1] * " or FVC (mL)"),
  graph.pos = 2,
  is.summary = c(rep(FALSE, nrow(plot2))),
  clip = c(-15, 60), # Adjust the clipping based on your actual data
  col = fpColors(box = "black"),
  lineheight = grid::unit(1.2, "lines"), # Adjust line height for better spacing
  xticks = c(-15, 0, 15, 30, 45, 60), # Explicitly set the tick marks to include 0
  boxsize = proportional_boxsize, # Use proportional box sizes
  align = c("l", "l"), # Align both the text labels and numbers to the left
  txt_gp = fpTxtGp(
    label = gpar(cex = 1, fontfamily = "Times"),   # Apply font to labels
    ticks = gpar(cex = 1, fontfamily = "Times"),   # Apply font to tick labels
    xlab = gpar(cex = 1, fontfamily = "Times"),    # Apply font to xlab
    title = gpar(cex = 1, fontfamily = "Times")    # Apply font to title (if any)
  )
)



# Adjust FEV1 and FVC labels
grid.text("FEV1", x = unit(0.16, "npc"), y = unit(0.975, "npc"), just = "left", gp = gpar(cex = 1.1, fontfamily = "Times"))
grid.text("FVC", x = unit(0.16, "npc"), y = unit(0.535, "npc"), just = "left", gp = gpar(cex = 1.1, fontfamily = "Times"))

grid.text(expression(beta * " (95% CI)"), x = unit(0.876, "npc"), y = unit(0.95, "npc"), just = "center", gp = gpar(cex = 1.1, fontfamily = "Times"))


# Capture the current plot
p <- recordPlot()

# Save the plot as a PDF with high resolution
png("plot2.png", width = 10, height = 8, units = "in", res = 300)
replayPlot(p) # Replay the captured plot
dev.off()