# set working directory and load packages
setwd("//rschfs1x/userrs/R-Z/wav25_RS/Desktop/updated code")
required_packages <- c("tidyverse", "dplyr", "haven", "purrr", "rlang", "foreign", "zoo")
# Install packages that are not already installed
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
# Load the required packages
lapply(required_packages, library, character.only = TRUE)
#load in data - set paths first
#HABCDATA
habcdatpath <- "//rschfs1x/userrs/R-Z/wav25_RS/Desktop/updated code/HABCdata"
y2clnvis <- read_sas(file.path(habcdatpath, "y2clnvis.sas7bdat"))
y2read <- read_sas(file.path(habcdatpath, "y2read.sas7bdat"))
ph <- read_sas(file.path(habcdatpath, "ph.sas7bdat"))
y1read <- read_sas(file.path(habcdatpath, "y1read.sas7bdat"))
y1calc <- read_sas(file.path(habcdatpath, "y1calc.sas7bdat"))
y5read <- read_sas(file.path(habcdatpath, "y5read.sas7bdat"))
y5calc <- read_sas(file.path(habcdatpath, "y5calc.sas7bdat"))
y8read <- read_sas(file.path(habcdatpath, "y8read.sas7bdat"))
y8calc <- read_sas(file.path(habcdatpath, "y8calc.sas7bdat"))
y10read <- read_sas(file.path(habcdatpath, "y10read.sas7bdat"))
y10calc <- read_sas(file.path(habcdatpath, "y10calc.sas7bdat"))
y4calc <- read_sas(file.path(habcdatpath, "y4calc.sas7bdat"))
#RAS DATA
rasdatpath <- "//rschfs1x/userrs/R-Z/wav25_RS/Desktop/updated code/RASdata"
ras <- read_sas(file.path(rasdatpath, "ras.sas7bdat"))
####HABC COHORT####
#merge data
merged <- merge(merge(y2clnvis, y2read, by = "HABCID"), ph, by = "HABCID")
df <- merged
#only keep variables needed
df <- df[, c(
  "HABCID",
  "B2BACNSZ", "B2OATMSZ", "B2FIBRSZ", "B2TOTLSZ", "B2SPEKSZ", "B2COTT", "B2FAT1", "B2FAT2", "B2SAUSSZ",
  "B2PANCSZ", "B2COTTSZ", "B2CHESSZ", "B2YOGRSZ", "B2EGGS", "B2BACN", "B2SAUS", "B2PANC", "B2OATM",
  "B2CERE", "B2FIBR", "B2TOTL", "B2SPEK", "B2CHES", "B2BROC", "B2FAT4", "B2FRYS", "B2POTA", "B2YAMS",
  "B2RICE", "B2STUF", "B2BEAN", "B2CORN", "B2PEAS", "B2MVEG", "B2SPIN", "B2GRNS", "B2SLAW", "B2GSAL",
  "B2TOMA", "B2DRES", "B2FRYSSZ", "B2OTHV", "B2SOUP", "B2OTHS", "B2DRESSZ", "B2SOUPSZ", "B2OTHSSZ",
  "B2STEW", "B2BURG", "B2BEEF", "B2LIVR", "B2PORK", "B2MIXD", "B2FRCH", "B2CRAB", "B2TUNA", "B2FISH",
  "B2OTHF", "B2HDOG", "B2MEAT", "B2FAT8", "B2FISHSZ", "B2OTHFSZ", "B2HDOGSZ", "B2MEATSZ", "B2MACA",
  "B2PIZZ", "B2PIZZSZ", "B2FAT5", "B2PAST", "B2FAT7", "B2FRCHSZ", "B2CHIK", "B2CHIKSZ", "B2MIXDSZ",
  "B2BURGSZ", "B2BEEFSZ", "B2LIVRSZ", "B2PORKSZ", "B2MUFF", "B2ROLL", "B2WBRD", "B2DBRD", "B2CBRD",
  "B2MARG", "B2MAYO", "B2PNUT", "B2KTCH", "B2SNCK", "B2NUTS", "B2CRCK", "B2GRAV", "B2KTCHSZ",
  "B2PNUTSZ", "B2MAYOSZ", "B2MARGSZ", "B2DBRDSZ", "B2WBRDSZ", "B2ROLLSZ", "B2DONU", "B2CAKE", "B2COOK",
  "B2ICEC", "B2PIES", "B2OTHP", "B2PUDD", "B2CHOC", "B2FAT9", "B2FAT10", "B2FAT11", "B2OJ", "B2HIC",
  "B2V8", "B2FRUT", "B2MLK2", "B2SOFT", "B2BEER", "B2WINE", "B2COFF", "B2TEA", "B2SHOT", "B2OJC",
  "B2MILK3", "B2SOFTSZ", "B2BEERSZ", "B2COFFSZ", "B2TEASZ", "B2CRMSZ", "B2MLK4SZ", "B2SUGRSZ", "B2SUPL",
  "B2CRM", "B2MLK4", "B2SUGR", "B2POTASZ", "B2YAMSSZ", "B2RICESZ", "B2STUFSZ", "B2BEANSZ", "B2CORNSZ",
  "B2PEASSZ", "B2BROCSZ", "B2MVEGSZ", "B2SPINSZ", "B2GRNSSZ", "B2SLAWSZ", "B2GSALSZ", "B2TOMASZ",
  "B2OTHVSZ", "B2PASTSZ", "B2MACASZ", "B2STEWSZ", "B2CRABSZ", "B2TUNASZ", "B2GRAVSZ", "B2SNCKSZ",
  "B2NUTSSZ", "B2CRCKSZ", "B2CBRDSZ", "B2MUFFSZ", "B2DONUSZ", "B2CAKESZ", "B2COOKSZ", "B2ICECSZ",
  "B2PIESSZ", "B2OTHPSZ", "B2PUDDSZ", "B2CHOCSZ", "B2OJSZ", "B2HICSZ", "B2V8SZ", "B2FRUTSZ", "B2MLK2SZ",
  "B2WINESZ", "B2MILK", "B2YOGR", "B2FAT13", "B2VEG", "B2FRUIT", "B2MARGST", "B2BUTTER", "B2MARGTB",
  "B2OLVOIL", "B2LARD", "B2LOWFAT", "B2VEGOIL", "B2CRISCO", "B2DK", "B2BLEND", "B2PAM", "B2APASZ",
  "B2CNPCSZ", "B2OTHSZ", "B2BANA", "B2APPL", "B2ORAN", "B2GFRU", "B2CANT", "B2PEAC", "B2APA", "B2CNPC",
  "B2OTH", "B2CANTSZ", "B2PEACSZ", "B2BANASZ", "B2APPLSZ", "B2ORANSZ", "B2GFRUSZ", "B2EGGSSZ", "B2FAT3",
  "B2FAT6", "B2WBUTTR", "B2SPRAY", "B2FAT12", "B2FRIED", "B2DATE", "FFQCALS", "RACE", "GENDER",
  "FFQOLEA", "FFQLINA", "FFQSFAT", "FFQSOD", "FFQPSWT"
)]
df <- df %>%
  filter(!is.na(FFQCALS) & FFQCALS != 0) %>%
  {print(paste("Count after excluding missing/0 calorie data:", nrow(.))); .} # Print count and pass the data frame to the next step
#function to standardize items as amount per 1,000kcal, standardizing across frequency of consumption and portion size
#fcolor (frequency color), refers to the frequency of consumption, pcolor refers to the portion size of the item (both are in categories that i've color-coded)
#ASSUMPTIONS:
# 1) used gregorian calendar mean value of 30.436875 days per month for 'month' values
# 2) when a range is given for frequency on an FFQ item (e.g., item consumed 2-3 x/mo),
# median consumption frequency will be used (thus, e.g., 2.5 x/mo)
#'frequency' will correpond to the frequency variable in the HABC dataset
#'portsize' will correspond to the portion size variable in the HABC dataset
#'portionquant' is set manually set for each item, and is just the reference portion size that is pulled directly from the FFQ for each item
#'percent' is set manually, depending on if the item is a part of a mix across multiple food categories
#'FFQCALS' corresponds to the total amount of calories for the item --this is pre-determined--
per1000kcal <- function(fcolor, pcolor, frequency, porsize, portionquant, percent, FFQCALS) {
  # frequency constants
  fqval <- switch(
    fcolor,
    "lightblue" = c(365/((1+11)/2), 30.436875, 30.436875/((2+3)/2), 7, 7/2, 7/((3+4)/2), 7/((5+6)/2), 1),
    "blue" = c(30.436875, 30.436875/((2+3)/2), 7, 7/2, 7/((3+4)/2), 7/((5+6)/2), 1, 2),
    "darkblue" = c(30.436875/((1+3)/2), 7, 7/((2+4)/2), 7/((5+6)/2), 1, (2+3)/2, 4, 5)
  )
  # Frequency codes
  fq <- switch(
    as.character(frequency),
    `1` = 0,
    `2` = fqval[1],
    `3` = fqval[2],
    `4` = fqval[3],
    `5` = fqval[4],
    `6` = fqval[5],
    `7` = fqval[6],
    `8` = fqval[7],
    `9` = fqval[8],
    NA
  )
  # Adjust missing values for porsize
  if (is.na(porsize)) {
    porsizeX <- switch(
      pcolor,
      "red" = 2,
      "lightgreen" = 2,
      "green" = 2,
      "darkgreen" = 1,
      "yellow" = 2
    )
    porsize <- porsizeX
  }
  # Portion size codes
  multiplier <- switch(
    pcolor,
    "red" = switch(as.character(porsize), `1` = 0.5, `2` = 1, `3` = 2, `4` = 4, NA),
    "lightgreen" = switch(as.character(porsize), `1` = 0.5, `2` = 1, `3` = 2, `4` = 3, NA),
    "green" = switch(as.character(porsize), `1` = 1/8, `2` = 1/4, `3` = 1/2, `4` = 1, NA),
    "darkgreen" = as.numeric(porsize),
    "yellow" = switch(as.character(porsize), `1` = 2/3, `2` = 1, `3` = 4/3, `4` = 5/3, NA),
    NA
  )
  # Calculate var
  var <- ifelse(frequency == 1, 0, (((portionquant / fq) / FFQCALS) * 1000) * multiplier * percent)
  return(var)
}
####HEI-2020####
#TOTAL FRUITS
df$b2apapk <- mapply(per1000kcal, "lightblue", "red", df$B2APA, df$B2APASZ, 1, 1, df$FFQCALS)
df$b2applpk <- mapply(per1000kcal, "lightblue", "lightgreen", df$B2APPL, df$B2APPLSZ, 1.25, 1, df$FFQCALS)
df$b2banapk <- mapply(per1000kcal, "lightblue", "lightgreen", df$B2BANA, df$B2BANASZ, 3/4, 1, df$FFQCALS)
df$b2cantpk <- mapply(per1000kcal, "lightblue", "green", df$B2CANT, df$B2CANTSZ, 4, 1, df$FFQCALS)
df$b2cnpcpk <- mapply(per1000kcal, "lightblue", "red", df$B2CNPC, df$B2CNPCSZ, 1, 1, df$FFQCALS)
df$b2frutpk <- mapply(per1000kcal, "darkblue", "red", df$B2FRUT, df$B2FRUTSZ, 1, 1, df$FFQCALS)
df$b2gfrupk <- mapply(per1000kcal, "lightblue", "lightgreen", df$B2GFRU, df$B2GFRUSZ, 1, 1, df$FFQCALS)
df$b2hicpkf <- mapply(per1000kcal, "darkblue", "red", df$B2HIC, df$B2HICSZ, 1, 1, df$FFQCALS)
df$b2ojpk <- mapply(per1000kcal, "darkblue", "red", df$B2OJ, df$B2OJSZ, 1, 1, df$FFQCALS)
df$b2oranpk <- mapply(per1000kcal, "lightblue", "lightgreen", df$B2ORAN, df$B2ORANSZ, 1, 1, df$FFQCALS)
df$b2othpk <- mapply(per1000kcal, "lightblue", "red", df$B2OTH, df$B2OTHSZ, 1, 1, df$FFQCALS)
df$b2peacpk <- mapply(per1000kcal, "lightblue", "lightgreen", df$B2PEAC, df$B2PEACSZ, 1, 1, df$FFQCALS)
df$sumtotalfruit <- rowSums(df[, c("b2apapk", "b2applpk", "b2banapk", "b2cantpk", "b2cnpcpk", "b2frutpk",
                                   "b2gfrupk", "b2hicpkf", "b2ojpk", "b2oranpk", "b2othpk", "b2peacpk")], na.rm = TRUE)
#WHOLE FRUIT
df$sumwholefruit <- rowSums(df[, c("b2applpk", "b2banapk", "b2cantpk", 
                                   "b2gfrupk", "b2oranpk", "b2othpk", "b2peacpk")], na.rm = TRUE)
#TOTAL VEGETABLES
df$b2beanpkv <- mapply(per1000kcal, "lightblue", "red", df$B2BEAN, df$B2BEANSZ, 1/2, 1, df$FFQCALS)
df$b2brocpk <- mapply(per1000kcal, "lightblue", "red", df$B2BROC, df$B2BROCSZ, 1, 1, df$FFQCALS)
df$b2cornpk <- mapply(per1000kcal, "lightblue", "red", df$B2CORN, df$B2CORNSZ, 1, 1, df$FFQCALS)
df$b2grnspk <- mapply(per1000kcal, "lightblue", "red", df$B2GRNS, df$B2GRNSSZ, 1, 1, df$FFQCALS)
df$b2gsalpk <- mapply(per1000kcal, "lightblue", "red", df$B2GSAL, df$B2GSALSZ, 1/2, 1, df$FFQCALS)
df$b2mixdpkv <- mapply(per1000kcal, "lightblue", "red", df$B2MIXD, df$B2MIXDSZ, 1, 0.10, df$FFQCALS)
df$b2mvegpk <- mapply(per1000kcal, "lightblue", "red", df$B2MVEG, df$B2MVEGSZ, 1, 1, df$FFQCALS)
df$b2othspkv <- mapply(per1000kcal, "lightblue", "red", df$B2OTHS, df$B2OTHSSZ, 1, 0.05, df$FFQCALS)
df$b2othvpk <- mapply(per1000kcal, "lightblue", "red", df$B2OTHV, df$B2OTHVSZ, 1, 1, df$FFQCALS)
df$b2pastpkv <- mapply(per1000kcal, "lightblue", "red", df$B2PAST, df$B2PASTSZ, 1, 0.05, df$FFQCALS)
df$b2peaspkv <- mapply(per1000kcal, "lightblue", "red", df$B2PEAS, df$B2PEASSZ, 1/2, 1, df$FFQCALS)
df$b2pizzpkv <- mapply(per1000kcal, "lightblue", "darkgreen", df$B2PIZZ, df$B2PIZZSZ, 3/4, 0.05, df$FFQCALS)
df$b2potapk <- mapply(per1000kcal, "lightblue", "red", df$B2POTA, df$B2POTASZ, 1, 1, df$FFQCALS)
df$b2slawpk <- mapply(per1000kcal, "lightblue", "red", df$B2SLAW, df$B2SLAWSZ, 1, 0.5, df$FFQCALS)
df$b2souppkv <- mapply(per1000kcal, "lightblue", "red", df$B2SOUP, df$B2SOUPSZ, 3/4, 0.75, df$FFQCALS)
df$b2spinpk <- mapply(per1000kcal, "lightblue", "red", df$B2SPIN, df$B2SPINSZ, 3/4, 1, df$FFQCALS)
df$b2stewpkv <- mapply(per1000kcal, "lightblue", "red", df$B2STEW, df$B2STEWSZ, 3/4, 0.3, df$FFQCALS)
df$b2tomapk <- mapply(per1000kcal, "lightblue", "red", df$B2TOMA, df$B2TOMASZ, 1, 1, df$FFQCALS)
df$b2v8pk <- mapply(per1000kcal, "darkblue", "red", df$B2V8, df$B2V8SZ, 1, 1, df$FFQCALS)
df$b2yamspk <- mapply(per1000kcal, "lightblue", "red", df$B2YAMS, df$B2YAMSSZ, 3/4, 1, df$FFQCALS)
df$sumtotalveg <- rowSums(df[, c("b2beanpkv", "b2brocpk", "b2cornpk", "b2grnspk", 
                                 "b2gsalpk", "b2mixdpkv", "b2mvegpk", "b2othspkv", "b2othvpk", 
                                 "b2pastpkv", "b2peaspkv", "b2pizzpkv", "b2potapk", "b2slawpk", 
                                 "b2souppkv", "b2spinpk", "b2stewpkv", "b2tomapk", "b2v8pk", "b2yamspk")], na.rm = TRUE)
#BEANS & GREENS
df$sumgreensandbeans <- rowSums(df[, c("b2beanpkv", "b2brocpk", "b2grnspk", "b2gsalpk", "b2peaspkv")], na.rm = TRUE)
#WHOLE GRAINS
df$b2dbrdpk <- mapply(per1000kcal, "blue", "darkgreen", df$B2DBRD, df$B2DBRDSZ, 1, 1, df$FFQCALS)
df$b2fibrpk <- mapply(per1000kcal, "lightblue", "red", df$B2FIBR, df$B2FIBRSZ, 1, 1, df$FFQCALS)
df$b2oatmpk <- mapply(per1000kcal, "lightblue", "red", df$B2OATM, df$B2OATMSZ, 1, 1, df$FFQCALS)
df$b2muffpkw <- mapply(per1000kcal, "blue", "red", df$B2MUFF, df$B2MUFFSZ, 1.5, 0.5, df$FFQCALS)
df$b2ricepkw <- mapply(per1000kcal, "lightblue", "red", df$B2RICE, df$B2RICESZ, 2, 0.05, df$FFQCALS)
df$b2snckpk <- mapply(per1000kcal, "blue", "red", df$B2SNCK, df$B2SNCKSZ, 1, 1, df$FFQCALS)
df$b2spekpk <- mapply(per1000kcal, "lightblue", "red", df$B2SPEK, df$B2SPEKSZ, 1, 1, df$FFQCALS)
df$b2totlpk <- mapply(per1000kcal, "lightblue", "red", df$B2TOTL, df$B2TOTLSZ, 1, 1, df$FFQCALS)
df$sumwholegrains <- rowSums(df[, c("b2dbrdpk", "b2fibrpk", "b2oatmpk", "b2muffpkw", "b2ricepkw", "b2snckpk", "b2spekpk", "b2totlpk")], na.rm = TRUE)
#DAIRY
df$b2chespk <- mapply(per1000kcal, "lightblue", "darkgreen", df$B2CHES, df$B2CHESSZ, 1/4, 1, df$FFQCALS)
df$b2cottpkd <- mapply(per1000kcal, "lightblue", "red", df$B2COTT, df$B2COTTSZ, 1/2, 1, df$FFQCALS)
df$b2macapkd <- mapply(per1000kcal, "lightblue", "red", df$B2MACA, df$B2MACASZ, 1, 0.4, df$FFQCALS)
df$b2milkpk <- mapply(per1000kcal, "lightblue", "red", df$B2MILK, 2, 1/2, 1, df$FFQCALS) #note, no portion size variable, standard portion applied
df$b2mixdpkd <- mapply(per1000kcal, "lightblue", "red", df$B2MIXD, df$B2MIXDSZ, 1, 0.05, df$FFQCALS)
df$b2mlk4pk <- mapply(per1000kcal, "darkblue", "red", df$B2MLK4, df$B2MLK4SZ, 1/16, 1, df$FFQCALS)
df$b2othspkd <- mapply(per1000kcal, "lightblue", "red", df$B2OTHS, df$B2OTHSSZ, 3/4, 0.10, df$FFQCALS)
df$b2pastpkd <- mapply(per1000kcal, "lightblue", "red", df$B2PAST, df$B2PASTSZ, 1, 0.10, df$FFQCALS)
df$b2pizzpkd <- mapply(per1000kcal, "lightblue", "darkgreen", df$B2PIZZ, df$B2PIZZSZ, 3/4, 0.20, df$FFQCALS)
df$b2yogrpk <- mapply(per1000kcal, "lightblue", "red", df$B2YOGR, df$B2YOGRSZ, 1, 1, df$FFQCALS)
df$b2suplpk <- mapply(per1000kcal, "lightblue", "red", df$B2SUPL, 2, 1, 1, df$FFQCALS) #note, no portion size variable, standard portion applied
df$sumdairy <- rowSums(df[, c("b2chespk", "b2cottpkd", "b2macapkd", "b2milkpk", "b2mixdpkd",
                              "b2mlk4pk", "b2othspkd", "b2pastpkd", "b2pizzpkd", "b2suplpk", "b2yogrpk")], na.rm = TRUE)
#TOTAL PROTEIN
df$b2bacnpk <- mapply(per1000kcal, "lightblue", "darkgreen", df$B2BACN, df$B2BACNSZ, 0.5, 1, df$FFQCALS)
df$b2beanpkp <- mapply(per1000kcal, "lightblue", "red", df$B2BEAN, df$B2BEANSZ, 2, 1, df$FFQCALS)
df$b2beefpk <- mapply(per1000kcal, "lightblue", "red", df$B2BEEF, df$B2BEEFSZ, 3, 1, df$FFQCALS)
df$b2burgpk <- mapply(per1000kcal, "lightblue", "red", df$B2BURG, df$B2BURGSZ, 4, 1, df$FFQCALS)
df$b2chikpk <- mapply(per1000kcal, "lightblue", "red", df$B2CHIK, df$B2CHIKSZ, 3, 1, df$FFQCALS)
df$b2crabpk <- mapply(per1000kcal, "lightblue", "red", df$B2CRAB, df$B2CRABSZ, 3, 1, df$FFQCALS)
df$b2eggspk <- mapply(per1000kcal, "lightblue", "darkgreen", df$B2EGGS, df$B2EGGSSZ, 1, 1, df$FFQCALS)
df$b2fishpk <- mapply(per1000kcal, "lightblue", "red", df$B2FISH, df$B2FISHSZ, 8, 1, df$FFQCALS)
df$b2frchpk <- mapply(per1000kcal, "lightblue", "red", df$B2FRCH, df$B2FRCHSZ, 7, 1, df$FFQCALS)
df$b2hdogpk <- mapply(per1000kcal, "lightblue", "darkgreen", df$B2HDOG, df$B2HDOGSZ, 1.5, 1, df$FFQCALS)
df$b2livrpk <- mapply(per1000kcal, "lightblue", "red", df$B2LIVR, df$B2LIVRSZ, 3, 1, df$FFQCALS)
df$b2meatpk <- mapply(per1000kcal, "lightblue", "darkgreen", df$B2MEAT, df$B2MEATSZ, 1, 1, df$FFQCALS)
df$b2mixdpkp <- mapply(per1000kcal, "lightblue", "red", df$B2MIXD, df$B2MIXDSZ, 8, 0.40, df$FFQCALS)
df$b2nutspk <- mapply(per1000kcal, "blue", "red", df$B2NUTS, df$B2NUTSSZ, 2, 1, df$FFQCALS)
df$b2othfpk <- mapply(per1000kcal, "lightblue", "red", df$B2OTHF, df$B2OTHFSZ, 3, 1, df$FFQCALS)
df$b2othspkp <- mapply(per1000kcal, "lightblue", "red", df$B2OTHS, df$B2OTHSSZ, 6, 0.30, df$FFQCALS)
df$b2pastpkp <- mapply(per1000kcal, "lightblue", "red", df$B2PAST, df$B2PASTSZ, 8, 0.25, df$FFQCALS)
df$b2peaspkp <- mapply(per1000kcal, "lightblue", "red", df$B2PEAS, df$B2PEASSZ, 2, 1, df$FFQCALS)
df$b2pizzpkp <- mapply(per1000kcal, "lightblue", "darkgreen", df$B2PIZZ, df$B2PIZZSZ, 6, 0.15, df$FFQCALS)
df$b2pnutpk <- mapply(per1000kcal, "blue", "darkgreen", df$B2PNUT, df$B2PNUTSZ, 1, 1, df$FFQCALS)
df$b2porkpk <- mapply(per1000kcal, "lightblue", "red", df$B2PORK, df$B2PORKSZ, 3, 1, df$FFQCALS)
df$b2sauspk <- mapply(per1000kcal, "lightblue", "darkgreen", df$B2SAUS, df$B2SAUSSZ, 1, 1, df$FFQCALS)
df$b2souppkp <- mapply(per1000kcal, "lightblue", "red", df$B2SOUP, df$B2SOUPSZ, 6, 0.25, df$FFQCALS)
df$b2stewpkp <- mapply(per1000kcal, "lightblue", "red", df$B2STEW, df$B2STEWSZ, 6, 0.70, df$FFQCALS)
df$b2tunapk <- mapply(per1000kcal, "lightblue", "red", df$B2TUNA, df$B2TUNASZ, 3, 1, df$FFQCALS)
df$sumtotalprotein <- rowSums(df[, c("b2bacnpk", "b2beanpkp", "b2beefpk", "b2burgpk", "b2chikpk", "b2crabpk", 
                                     "b2eggspk", "b2fishpk", "b2frchpk", "b2hdogpk", "b2livrpk", "b2meatpk", 
                                     "b2mixdpkp", "b2nutspk", "b2othfpk", "b2othspkp", "b2peaspkp", "b2pastpkp", "b2pizzpkp", 
                                     "b2pnutpk", "b2porkpk", "b2sauspk", "b2souppkp", "b2stewpkp", "b2tunapk")], na.rm = TRUE)
#SEAFOOD AND PLANT PROTEINS
df$sumseafoodandplant <- rowSums(df[, c("b2beanpkp", "b2crabpk", "b2fishpk", "b2nutspk", 
                                        "b2othfpk", "b2peaspkp", "b2pnutpk", "b2tunapk")], na.rm = TRUE)
#REFINED GRAINS
df$b2cakepk <- mapply(per1000kcal, "blue", "red", df$B2CAKE, df$B2CAKESZ, 1.75, 1, df$FFQCALS)
df$b2cookpk <- mapply(per1000kcal, "blue", "red", df$B2COOK, df$B2COOKSZ, 0.25, 1, df$FFQCALS)
df$b2crckpk <- mapply(per1000kcal, "blue", "red", df$B2CRCK, df$B2CRCKSZ, 1.5, 1, df$FFQCALS)
df$b2donupk <- mapply(per1000kcal, "blue", "red", df$B2DONU, df$B2DONUSZ, 1.25, 1, df$FFQCALS)
df$b2fryspk <- mapply(per1000kcal, "lightblue", "red", df$B2FRYS, df$B2FRYSSZ, 4, 1, df$FFQCALS)
df$b2gravpk <- mapply(per1000kcal, "blue", "darkgreen", df$B2GRAV, df$B2GRAVSZ, 0.5, 1, df$FFQCALS)
df$b2macapkr <- mapply(per1000kcal, "lightblue", "red", df$B2MACA, df$B2MACASZ, 8, 0.66, df$FFQCALS)
df$b2mixdpkr <- mapply(per1000kcal, "lightblue", "red", df$B2MIXD, df$B2MIXDSZ, 8, 0.45, df$FFQCALS)
df$b2muffpkr <- mapply(per1000kcal, "blue", "red", df$B2MUFF, df$B2MUFFSZ, 1.5, 0.5, df$FFQCALS)
df$b2othppk <- mapply(per1000kcal, "blue", "red", df$B2OTHP, df$B2OTHPSZ, 3, 1, df$FFQCALS)
df$b2othspkr <- mapply(per1000kcal, "lightblue", "red", df$B2OTHS, df$B2OTHSSZ, 6, 0.40, df$FFQCALS)
df$b2pancpk <- mapply(per1000kcal, "lightblue", "darkgreen", df$B2PANC, df$B2PANCSZ, 1, 1, df$FFQCALS)
df$b2pastpkr <- mapply(per1000kcal, "lightblue", "red", df$B2PAST, df$B2PASTSZ, 8, 0.60, df$FFQCALS)
df$b2piespk <- mapply(per1000kcal, "blue", "red", df$B2PIES, df$B2PIESSZ, 2.5, 1, df$FFQCALS)
df$b2pizzpkr <- mapply(per1000kcal, "lightblue", "darkgreen", df$B2PIZZ, df$B2PIZZSZ, 6, 0.60, df$FFQCALS)
df$b2ricepkr <- mapply(per1000kcal, "lightblue", "red", df$B2RICE, df$B2RICESZ, 2, 0.95, df$FFQCALS)
df$b2rollpk <- mapply(per1000kcal, "blue", "lightgreen", df$B2ROLL, df$B2ROLLSZ, 1.75, 1, df$FFQCALS)
df$b2stufpk <- mapply(per1000kcal, "lightblue", "red", df$B2STUF, df$B2STUFSZ, 1, 1, df$FFQCALS)
df$b2wbrdpk <- mapply(per1000kcal, "blue", "darkgreen", df$B2WBRD, df$B2WBRDSZ, 1, 1, df$FFQCALS)
df$sumrefinedgrains <- rowSums(df[, c("b2cakepk", "b2cookpk", "b2crckpk", "b2donupk", "b2fryspk",
                                      "b2gravpk", "b2macapkr", "b2mixdpkr", "b2muffpkr", "b2othppk", "b2othspkr", "b2pancpk", 
                                      "b2pastpkr", "b2piespk", "b2pizzpkr", "b2ricepkr", "b2rollpk", "b2snckpk", 
                                      "b2stufpk", "b2wbrdpk")], na.rm = TRUE)
#FATTY ACIDS
df$faratio <- (df$FFQOLEA + df$FFQLINA) / df$FFQSFAT
#SODIUM
df$totalsodium <- (df$FFQSOD / 1000) / df$FFQCALS * 1000
#ADDED SUGARS
df$addsug <- df$FFQPSWT
#SATURATED FAT
df$satfat <- ((df$FFQSFAT * 9) / df$FFQCALS) * 100
#HEI-2020 SCORE CALCULATIONS
df$totalfruitscore <- pmin((df$sumtotalfruit * 5) / 0.8, 5)
df$wholefruitscore <- pmin((df$sumwholefruit * 5) / 0.4, 5)
df$totalvegscore <- pmin((df$sumtotalveg * 5) / 1.1, 5)
df$greensandbeansscore <- pmin((df$sumgreensandbeans * 5) / 0.2, 5)
df$wholegrainsscore <- pmin((df$sumwholegrains * 10) / 1.5, 10)
df$dairyscore <- pmin((df$sumdairy * 10) / 1.3, 10)
df$totalproteinscore <- pmin((df$sumtotalprotein * 5) / 2.5, 5)
df$seafoodandplantscore <- pmin((df$sumseafoodandplant * 5) / 0.8, 5)
df$refinedgrainsscore <- pmax(pmin(10 - (10 * (df$sumrefinedgrains - 1.8) / (4.3 - 1.8)), 10), 0)
df$fattyacidscore <- ifelse(df$FFQSFAT == 0 & (df$FFQOLEA + df$FFQLINA) == 0, 0, pmin(pmax(10 * ((df$faratio - 1.2) / (2.5 - 1.2)), 0), 10))
df$sodiumscore <- pmax(pmin(10 - (10 * (df$totalsodium - 1.1) / (2.0 - 1.1)), 10), 0)
df$addsugarscore <- pmax(pmin(10 - (10 * (df$addsug - 6.5) / (26 - 6.5)), 10), 0)
df$satfatscore <- pmax(pmin(10 - (10 * (df$satfat - 8) / (16 - 8)), 10), 0)
HEIscores <- c("totalfruitscore", "wholefruitscore", "totalvegscore", "greensandbeansscore",
                   "wholegrainsscore", "dairyscore", "totalproteinscore", "seafoodandplantscore", 
                   "refinedgrainsscore", "fattyacidscore", "sodiumscore", "addsugarscore", "satfatscore")
df[HEIscores] <- lapply(df[HEIscores], function(x) ifelse(df$FFQCALS == 0 | is.na(df$FFQCALS), NA, x))
df$HEI2020score <- ifelse(df$FFQCALS == 0 | is.na(df$FFQCALS), NA, 
                          rowSums(df[HEIscores], na.rm = TRUE))

####DASH####
#FRUITS
df$sumfruitdash <- df$sumtotalfruit

#VEGETABLES
df$sumvegdash <- rowSums(df[, c("b2brocpk", "b2cornpk", "b2grnspk", "b2gsalpk", "b2mixdpkv",
                                "b2mvegpk", "b2othvpk", "b2pastpkv", "b2pizzpkv", "b2slawpk",
                                "b2spinpk", "b2stewpkv", "b2tomapk", "b2v8pk", "b2yamspk")], na.rm = TRUE)
#NUTS & LEGUMES
df$sumnutslegdash <- rowSums(df[, c("b2beanpkv", "b2nutspk", "b2peaspkv", "b2pnutpk")], na.rm = TRUE)
#WHOLE GRAINS
df$sumwholegrainsdash <- df$sumwholegrains
#LOW-FAT DAIRY
df$b2cottpklf <- mapply(per1000kcal, "lightblue", "red", df$B2COTT, df$B2COTTSZ, 1/2, 0.5, df$FFQCALS)
df$b2milkpklf <- ifelse(df$B2MILK3 == 1, 0, #conditional check for low-fat status
                      mapply(per1000kcal, "lightblue", "red", df$B2MILK, 2, 1/2, 1, df$FFQCALS)) #note, no portion size variable, standard portion applied
df$b2yogrpklf <- ifelse(df$B2FAT2 %in% c(2, 3), 0, #conditional check for low-fat status
                      mapply(per1000kcal, "lightblue", "red", df$B2YOGR, df$B2YOGRSZ, 1, 1, df$FFQCALS))
df$sumdairydash <- rowSums(df[, c("b2cottpklf", "b2milkpklf", "b2yogrpklf")], na.rm = TRUE)
#RED AND PROCESSED MEAT
df$sumredprocdash <- rowSums(df[, c("b2bacnpk", "b2beefpk", "b2burgpk", "b2hdogpk", "b2meatpk", 
                                    "b2mixdpkp", "b2pastpkp", "b2pizzpkp", "b2porkpk", "b2sauspk", "b2souppkp")], na.rm = TRUE)
#SUGAR-SWEETENED BEVERAGES
df$b2hicpks <- mapply(per1000kcal, "darkblue", "red", df$B2HIC, df$B2HICSZ, 6, 1, df$FFQCALS)
df$b2mlk2pk <- mapply(per1000kcal, "darkblue", "red", df$B2MLK2, df$B2MLK2SZ, 8, 1, df$FFQCALS)
df$b2softpk <- mapply(per1000kcal, "darkblue", "yellow", df$B2SOFT, df$B2SOFTSZ, 12, 1, df$FFQCALS)
df$b2sugrpk <- mapply(per1000kcal, "darkblue", "darkgreen", df$B2SUGR, df$B2SUGRSZ, 1, 1, df$FFQCALS)
df$sumssb <- rowSums(df[, c("b2hicpks", "b2mlk2pk", "b2softpk", "b2sugrpk")], na.rm = TRUE)
#SODIUM
df$totalsodiumdash <- ((df$FFQSOD/1000)/df$FFQCALS)*1000
#DASH SCORE CALCULATIONS
DASHvars <- c("sumfruitdash", "sumvegdash", "sumnutslegdash", "sumwholegrainsdash", 
              "sumdairydash", "sumredprocdash", "sumssb", "totalsodiumdash")
df <- df %>% mutate_at(vars(all_of(DASHvars)), ~ifelse(is.na(.), 0, .))
quint <- lapply(DASHvars, function(x) {
  quantile(df[[x]], probsprobs = seq(0, 1, 0.2), na.rm = TRUE)
})
DASHscores <- function(x, q, reverse = FALSE) {
  score_val <- case_when(
    x >= q[4] ~ ifelse(reverse, 1, 5),
    x >= q[3] ~ ifelse(reverse, 2, 4),
    x >= q[2] ~ ifelse(reverse, 3, 3),
    x >= q[1] ~ ifelse(reverse, 4, 2),
    TRUE ~ ifelse(reverse, 5, 1))
  return(score_val)
}
ScoreNames <- c("dashfruitscore", "dashvegscore", "dashnutslegscore", "dashwholegrainsscore", 
                "dashdairyscore", "dashredprocscore", "dashssbscore", "dashsodscore")
for (i in seq_along(DASHvars)) {
  var <- DASHvars[i]
  quintile <- quint[[i]] #set quintiles for each component
  rev <- var %in% c("sumredprocdash", "sumssb", "totalsodiumdash") #set reverse for components scored in the opposite direction
  df[[ScoreNames[i]]] <- mapply(function(value, FFQCALS) {
    if (is.na(FFQCALS) | FFQCALS == 0) return(NA) #set score to NA if FFQCALS is NA or 0
    DASHscores(value, quintile, rev) #apply scoring function
  }, df[[var]], df$FFQCALS)
}
df$DASHscore <- rowSums(df[, ScoreNames], na.rm = FALSE)
#final scores dataset
habcscores <- df %>%
  select(HABCID, HEI2020score, totalfruitscore, 
         wholefruitscore, totalvegscore, greensandbeansscore, wholegrainsscore, 
         dairyscore, totalproteinscore, seafoodandplantscore, refinedgrainsscore, 
         fattyacidscore, sodiumscore, addsugarscore, satfatscore, FFQCALS, DASHscore, dashfruitscore, dashvegscore, 
         dashnutslegscore, dashwholegrainsscore, dashdairyscore, 
         dashredprocscore, dashssbscore, dashsodscore)

###CREATING LONG FORM DATASET FOR ANALYSIS###
cov <- ph %>% select(HABCID, RACE, GENDER, SITE, DOB, CV1AGE, CV1DATE, CV2DATE, CV5DATE, CV8DATE, CV10DATE)
QC <- function(data, PFTvar, QCvar, year) { #function to create QC-PFT variables for FEV1 & FVC
  QCsuffix <- if(QCvar == "QCFEV1") {
    "QCFEV1" } else if (QCvar == "QCFVC") {
      "QCFVC" } else { stop("VONDY ERROR: VARIABLE MISSING") }
  QCFx <- paste0("Y", year, QCsuffix) #these first two steps make a QC variable based on FEV1 or FVC, depending on input
  data %>%
    mutate(SPIROAVA1 = ifelse({{ PFTvar }} > 0, 1, 0),
           QCACCEPT1 = ifelse(SPIROAVA1 == 0, NA, 
                              ifelse(SPIROAVA1 == 1 & .data[[QCvar]] > 0, 1, 0)),
           !!sym(QCFx) := .data[[QCvar]]) %>%
    filter(SPIROAVA1 == 1 & QCACCEPT1 == 1) #filter to only include acceptable PFTs
}
process <- function(read, calc, year) { #function to process PFT data and create resultant datasets
  #define variable names for FEV1, FVC
  fev1 <- paste0("Y", year, "FEV1")
  fvc <- paste0("Y", year, "FVC")
  #fev1
  fev1dat <- read %>%
    mutate(!!sym(fev1) := pmax(BES_FEV1, SCN_FEV1, TRD_FEV1, na.rm = TRUE)) %>%
    QC(!!sym(fev1), "QCFEV1", year) %>%
    select(HABCID, !!sym(fev1), !!sym(paste0("Y", year, "QCFEV1")))
  #fvc
  fvcdat <- read %>%
    mutate(!!sym(fvc) := pmax(BES_FVC, SCN_FVC, TRD_FVC, na.rm = TRUE)) %>%
    QC(!!sym(fvc), "QCFVC", year) %>%
    select(HABCID, !!sym(fvc), !!sym(paste0("Y", year, "QCFVC")))
  #merge fev1 and fvc datasets
  pftdat <- full_join(fev1dat, fvcdat, by = "HABCID") %>%
    mutate(!!sym(paste0("Y", year, "RATIO")) := !!sym(fev1) / !!sym(fvc),
           !!sym(paste0("Y", year, "AIROBS")) := case_when(
             is.na(!!sym(paste0("Y", year, "RATIO"))) ~ NA_real_,
             !!sym(paste0("Y", year, "RATIO")) < 0.7 ~ 1,
             !!sym(paste0("Y", year, "RATIO")) >= 0.7 ~ 0
           ))
  #smoking status and packyears
  if (year == 5) {calc$P2SH <- y4calc$D2SH} #year 5 height variable missing, use year 4 instead
  htvar <- ifelse(year %in% c(8, 10), "S4SH", "P2SH") #year 8 and 10 use S4SH, year 1 uses P2SH for height
  htvarname <- if(year == 5) "Y4HT" else paste0("Y", year, "HT") #creates custom name scheme for y4's HT variable within y5 dataset
  calcdat <- calc %>%
    mutate(!!sym(htvarname) := .data[[htvar]]/1000,
           !!sym(paste0("Y", year, "LLNFVC")) := LLNFVC,
           !!sym(paste0("Y", year, "LLNFEV1")) := LLNFEV1,
           !!sym(paste0("Y", year, "LLNRATIO")) := LLNFEV1R,
           !!sym(paste0("Y", year, "PREDFVC")) := PREDFVC,
           !!sym(paste0("Y", year, "PREDFEV1")) := PREDFEV1,
           !!sym(paste0("Y", year, "PREDRATIO")) := PRDFEV1R) %>%
    select(HABCID,
           !!sym(paste0("SMK", year)),
           !!sym(htvarname), 
           !!sym(paste0("Y", year, "LLNFVC")), 
           !!sym(paste0("Y", year, "LLNFEV1")), 
           !!sym(paste0("Y", year, "LLNRATIO")), 
           !!sym(paste0("Y", year, "PREDFVC")), 
           !!sym(paste0("Y", year, "PREDFEV1")), 
           !!sym(paste0("Y", year, "PREDRATIO")),
           !!if (year == 1) sym("PACKYR1") else NULL,
           !!if (year == 1) sym("EDUC") else NULL)
  #list of resultant datasets
  dats <- list(fev1 = fev1dat, fvc = fvcdat, pft = pftdat, calc = calcdat)
  #naming datasets based on year
  names(dats) <- c(paste0("y", year, "fev1"), paste0("y", year, "fvc"), 
                   paste0("y", year, "pft"), paste0("y", year, "calc"))
  return(dats)
}
y1results <- process(y1read, y1calc, 1)
y5results <- process(y5read, y5calc, 5)
y8results <- process(y8read, y8calc, 8)
y10results <- process(y10read, y10calc, 10)
#print counts of results
rezlist <- list(y1 = y1results, y5 = y5results, y8 = y8results, y10 = y10results)
sapply(rezlist, function(year) sapply(c("fev1", "fvc", "pft"), function(pft) dim(year[[paste0("y", sub(".*y(\\d+).*", "\\1", names(year)[1]), pft)]])[1]))
merged <- Reduce(function(x, y) merge(x, y, by = "HABCID", all = TRUE), 
                      list(cov, habcscores, y1results$y1pft, y1results$y1calc, 
                           y5results$y5pft, y5results$y5calc, 
                           y8results$y8pft, y8results$y8calc, 
                           y10results$y10pft, y10results$y10calc))
#create dummy variables for race, sex, site location
merged$RACEDUMMY <- ifelse(merged$RACE == 1, 0, ifelse(merged$RACE == 2, 1, NA)) #0 = white, 1 = black
merged$GENDERDUMMY <- ifelse(merged$GENDER == 1, 1, ifelse(merged$GENDER == 2, 0, NA)) #0 = female, 1 = male
merged$SITEDUMMY <- ifelse(merged$SITE == 1, 0, ifelse(merged$SITE == 2, 1, NA)) #0 = memphis, 1 = pittsburgh
#create smoking status variables; first, carry over smoking status from prior visit where missing
merged <- merged %>%
  mutate(
    SMKSTAT1 = SMK1,
    SMKSTAT5 = ifelse(is.na(SMK5), SMK1, SMK5),
    SMKSTAT8 = ifelse(is.na(SMK8) & is.na(SMK5), SMK1,
                      ifelse(is.na(SMK8), SMK5, SMK8)),
    SMKSTAT10 = ifelse(is.na(SMK10) & is.na(SMK8) & is.na(SMK5), SMK1,
                       ifelse(is.na(SMK10) & is.na(SMK8), SMK5,
                              ifelse(is.na(SMK10), SMK8, SMK10)))
  )
#create longitudinal smoking status variable (0 = never, 1 = persistent, 2 = former, 3 = intermittent)
merged <- merged %>%
  mutate(
    SMKCAT = case_when(
      is.na(SMKSTAT1) & is.na(SMKSTAT5) & is.na(SMKSTAT8) & is.na(SMKSTAT10) ~ NA_real_,
      SMKSTAT1 == SMKSTAT5 & SMKSTAT1 == SMKSTAT8 & SMKSTAT1 == SMKSTAT10 ~ SMKSTAT1,
      SMKSTAT1 != SMKSTAT5 | SMKSTAT1 != SMKSTAT8 | SMKSTAT1 != SMKSTAT10 |
        SMKSTAT5 != SMKSTAT8 | SMKSTAT5 != SMKSTAT10 |
        SMKSTAT8 != SMKSTAT10 ~ 3,
      TRUE ~ NA_real_
    )
  )
table(merged$SMKCAT, useNA = "ifany") #5 participants missing all smoking data
#carrying over height data from prior visit where missing
summary(merged[, c("Y1HT", "Y4HT", "Y8HT", "Y10HT")]) #Y1HT = 0 miss, Y5HT = 682 miss, Y8HT = 1425 miss, Y10HT = 1621 miss
merged <- merged %>%
  mutate(
    Y4HT = case_when(
      !is.na(Y4HT) ~ Y4HT,
      is.na(Y1HT) & is.na(Y4HT) & is.na(Y8HT) & is.na(Y10HT) ~ NA_real_,
      !is.na(Y1HT) & is.na(Y4HT) & !is.na(Y8HT) ~ (Y1HT + Y8HT) / 2,
      !is.na(Y1HT) & is.na(Y4HT) & is.na(Y8HT) & !is.na(Y10HT) ~ (Y1HT + Y10HT) / 2,
      !is.na(Y1HT) & is.na(Y4HT) & is.na(Y8HT) & is.na(Y10HT) ~ Y1HT,
      TRUE ~ Y4HT
    )
  ) %>%
  mutate(
    Y8HT = case_when(
      !is.na(Y8HT) ~ Y8HT,
      is.na(Y1HT) & is.na(Y4HT) & is.na(Y8HT) & is.na(Y10HT) ~ NA_real_,
      !is.na(Y4HT) & is.na(Y8HT) & !is.na(Y10HT) ~ (Y4HT + Y10HT) / 2,
      is.na(Y4HT) & is.na(Y8HT) & !is.na(Y10HT) ~ (Y1HT + Y10HT) / 2,
      is.na(Y1HT) & is.na(Y4HT) & is.na(Y8HT) & !is.na(Y10HT) ~ Y10HT,
      is.na(Y10HT) ~ Y4HT,
      TRUE ~ Y8HT
    )
  ) %>%
  mutate(
    Y10HT = if_else(!is.na(Y10HT), Y10HT, Y8HT)
  )
merged <- merged %>%
  mutate( #creating 'PFT-count' variable; total number of PFT measurements for each participant
    N_FEV1 = rowSums(!is.na(select(., Y1FEV1, Y5FEV1, Y8FEV1, Y10FEV1)), na.rm = TRUE)
  )
#reshape dataset to long
reshape <- function(df) {
  longlist <- list()
  timepoints <- c(1, 5, 8, 10)
  years <- c("Y1", "Y5", "Y8", "Y10")
  for (i in seq_along(timepoints)) {
    timepoint <- timepoints[i]
    year <- years[i]
    htvar <- ifelse(year == "Y5", "Y4HT", paste0(year, "HT")) #year 5 uses Y4HT, other years use the same year's HT variable
    long <- df %>%
      mutate(
        ID = row_number(),
        FEV1 = as.numeric(.[[paste0(year, "FEV1")]]),
        QCFEV1 = factor(.[[paste0(year, "QCFEV1")]]),
        FVC = as.numeric(.[[paste0(year, "FVC")]]),
        QCFVC = factor(.[[paste0(year, "QCFVC")]]),
        RATIO = as.numeric(.[[paste0(year, "RATIO")]]),
        AIROBS = ifelse(.[[paste0(year, "RATIO")]] < 0.7, 1, 0),
        LLNFEV1 = as.numeric(.[[paste0(year, "LLNFEV1")]]),
        PREDFEV1 = as.numeric(.[[paste0(year, "PREDFEV1")]]),
        LLNFVC = as.numeric(.[[paste0(year, "LLNFVC")]]),
        PREDFVC = as.numeric(.[[paste0(year, "PREDFVC")]]),
        LLNRATIO = as.numeric(.[[paste0(year, "LLNRATIO")]]),
        PREDRATIO = as.numeric(.[[paste0(year, "PREDRATIO")]]),
        HT = as.numeric(.[[htvar]]),
        COUNT = i,
        CURRSTAT = coalesce(!!!syms(rev(paste0("SMK", timepoints)))) #current or most recent smoking status
      )
    longlist[[i]] <- long
  }
  longdat <- bind_rows(longlist) %>%
    arrange(HABCID, COUNT)
  return(longdat)
}
longpft <- reshape(merged)
longpft <- longpft %>%
  mutate(CVDATE = case_when(
    COUNT == 1 & is.na(CV1DATE) ~ NA,
    COUNT == 1 ~ CV1DATE,
    COUNT == 2 & is.na(CV5DATE) ~ NA,
    COUNT == 2 ~ CV5DATE,
    COUNT == 3 & is.na(CV8DATE) ~ NA,
    COUNT == 3 ~ CV8DATE,
    COUNT == 4 & is.na(CV10DATE) ~ NA,
    COUNT == 4 ~ CV10DATE,
    TRUE ~ NA
  ))
longpft <- longpft %>%
  group_by(HABCID) %>% #create 'TIME' variable; time since first clinic visit
  mutate(TIME = as.numeric(CVDATE - first(CVDATE), units = "days") / 365.25)
longpft <- longpft %>%
  mutate( #recreating 'duration' variable Bonnie created
    DURATION = case_when(
      !is.na(Y1AIROBS) ~ as.numeric(TIME),
      is.na(Y1AIROBS) & !is.na(Y5AIROBS) & COUNT == 2 ~ 0,
      is.na(Y1AIROBS) & !is.na(Y5AIROBS) & COUNT == 3 ~ as.numeric((CV8DATE - CV5DATE) / 365),
      is.na(Y1AIROBS) & !is.na(Y5AIROBS) & COUNT == 4 ~ as.numeric((CV10DATE - CV5DATE) / 365),
      is.na(Y1AIROBS) & is.na(Y5AIROBS) & !is.na(Y8AIROBS) & COUNT == 3 ~ 0,
      is.na(Y1AIROBS) & is.na(Y5AIROBS) & !is.na(Y8AIROBS) & COUNT == 4 ~ as.numeric((CV10DATE - CV8DATE) / 365),
      is.na(Y1AIROBS) & is.na(Y5AIROBS) & is.na(Y8AIROBS) & !is.na(Y10AIROBS) & COUNT == 4 ~ 0
    ))
#dummy variables for current smoking status
longpft$CSMK <- ifelse(is.na(longpft$CURRSTAT), NA, ifelse(longpft$CURRSTAT == 1, 1, 0))
longpft$FSMK <- ifelse(is.na(longpft$CURRSTAT), NA, ifelse(longpft$CURRSTAT == 2, 1, 0))
longpft <- longpft %>%
  mutate(
    SMKCAT2 = case_when(
      SMKCAT == 3 & CSMK == 0 ~ 3,
      SMKCAT == 3 & CSMK == 1 ~ 4,
      TRUE ~ SMKCAT
    ))
#secondary smoking category variable & duration
longpft <- longpft %>%
  group_by(HABCID) %>%
  mutate( #create '%of predicted PFT' variables
    RATIO = RATIO * 100,  #turn FEV1/FVC ratio into percentage (to match existing variables in the HABC dataset)
    PCTPREDFEV1 = FEV1 / PREDFEV1,
    PCTPREDFVC = FVC / PREDFVC,
    PCTPREDRATIO = (100 * RATIO) / PREDRATIO,
    LOWFEV1 = case_when(
      is.na(FEV1) ~ NA_real_,
      FEV1 < LLNFEV1 ~ 1,
      FEV1 >= LLNFEV1 ~ 0
    ),
    LOWFVC = case_when(
      is.na(FVC) ~ NA_real_,
      FVC < LLNFVC ~ 1,
      FVC >= LLNFVC ~ 0
    ),
    LOWRATIO = case_when(
      is.na(RATIO) ~ NA_real_,
      RATIO < LLNRATIO ~ 1,
      RATIO >= LLNRATIO ~ 0
    )
  ) %>%
  ungroup()
#final data transformations
longpft <- longpft %>%
  mutate(
    HTCM = HT * 100, #convert height to cm
    AGE = CV1AGE + TIME, #age at each visit
    HT2 = (HTCM - mean(HTCM, na.rm = T))^2, #height squared variable (non-linearity assessment variable)
    AGE2 = (AGE - mean(AGE, na.rm = T))^2, #age squared variable (non-linearity assessment variable))
  ) %>%
  mutate( #change smoking categories to words, add NA
    SMKCAT = case_when(
      SMKCAT == 0 ~ "Never",
      SMKCAT == 1 ~ "Persistent",
      SMKCAT == 2 ~ "Former",
      SMKCAT == 3 ~ "Intermittent",
      TRUE ~ "NA"
    ))
classifyCOPD <- function(val) {#function to classify COPD types
  if (all(val == 0)) {
    return("NoCOPD")
  } else if (any(val == 1)) {
    firsti <- which(val == 1)[1] #identify first incident COPD
    subsequent <- val[firsti:length(val)] #check subsequent values of COPD
    if (length(subsequent) > 1 && subsequent[2] == 0) {      
      return("Transient")
    } else {
      rleSubsequent <- rle(subsequent)
      runs <- rleSubsequent$lengths[rleSubsequent$values == 1] #vector within 'rle'
      if (any(runs >= 2)) {
        return("Chronic")
      } else {
        return("Acute")
      }
    }
  }
}
copdp <- longpft %>% #COPD variables
  filter((!is.na(FEV1) & !is.na(FVC))) %>%
  mutate(COPD = ifelse(RATIO < 70, 1, 0)) %>% #create *LONG FORM* COPD variable
  group_by(HABCID) %>% arrange(HABCID, CVDATE) %>%
  mutate(COPDTIME = as.numeric(CVDATE - first(CVDATE), units = "days") / 365.25) %>% #create COPD 'time' variable for years elapsed since first COPD visit (baseline)
  mutate(COPDprev = ifelse(first(COPD) == 1, 1, 0)) %>% #create COPD prevalence at baseline variable (binary: 0/1)
  mutate(COPDany = ifelse(any(COPD == 1), 1, 0))  # Create 'COPDany' variable to identify participants with any COPD diagnosis (this we will use for the COPD prevalence analysis)
copdi <- copdp %>% group_by(HABCID) %>% filter(COPDprev == 0 & !is.na(COPD) & (QCFEV1 %in% c(3, 4) & QCFVC %in% c(3, 4))) %>%
  mutate(COPDclass = classifyCOPD(COPD)) %>%
  select(HABCID, TIME, COPDclass)
copdp <- copdp %>% select(HABCID, TIME, COPD, COPDTIME, COPDprev, COPDany)
#merge datasets
longpft <- merge(longpft, copdp, by = c("HABCID", "TIME"), all.x = TRUE)
longpft <- merge(longpft, copdi, by = c("HABCID", "TIME"), all.x = TRUE)
final <- longpft %>% #cut unnecessary variables from final dataset
  select(-DOB, -CV1AGE, -CV1DATE, -CV2DATE, -CV5DATE, -CV8DATE, -CV10DATE,
         -Y1FEV1, -Y1QCFEV1, -Y1FVC, -Y1QCFVC, -Y1RATIO, -Y1AIROBS, -SMK1, -Y1HT, -Y1LLNFVC,
         -Y1LLNFEV1, -Y1LLNRATIO, -Y1PREDFVC, -Y1PREDFEV1, -Y1PREDRATIO, -Y5FEV1, -Y5QCFEV1,
         -Y5FVC, -Y5QCFVC, -Y5RATIO, -Y5AIROBS, -SMK5, -Y4HT, -Y5LLNFVC, -Y5LLNFEV1, -Y5LLNRATIO,
         -Y5PREDFVC, -Y5PREDFEV1, -Y5PREDRATIO, -Y8FEV1, -Y8QCFEV1, -Y8FVC, -Y8QCFVC, -Y8RATIO,
         -Y8AIROBS, -SMK8, -Y8HT, -Y8LLNFVC, -Y8LLNFEV1, -Y8LLNRATIO, -Y8PREDFVC, -Y8PREDFEV1,
         -Y8RATIO, -Y8PREDRATIO, -Y10FEV1, -Y10QCFEV1, -Y10FVC, -Y10QCFVC, -Y10RATIO, -Y10AIROBS, -SMK10,
         -Y10HT, -Y10LLNFVC, -Y10LLNFEV1, -Y10LLNRATIO, -Y10PREDFVC, -Y10PREDFEV1, -Y10PREDRATIO,
         -SMKSTAT1, -SMKSTAT5, -SMKSTAT8, -SMKSTAT10, -ID, -RACEDUMMY, -GENDERDUMMY, -SITEDUMMY,
         -CURRSTAT, -CSMK, -FSMK, -SMKCAT2, -PCTPREDFEV1, -PCTPREDFVC, -PCTPREDRATIO, -LOWFEV1,
         -LOWFVC, -LOWRATIO, -N_FEV1, -AIROBS, -LLNFEV1, -PREDFEV1, -LLNFVC, -PREDFVC, -LLNRATIO, -PREDRATIO, -COUNT)
#final cuts
print(length(unique(final$HABCID))) #full participant count; n = 3,075
final <- final %>% #remove participants with no diet data
  filter(!is.na(FFQCALS))
print(length(unique(final$HABCID))) #n = 2,713
final <- final %>% #remove participants with no PFT data (missing both FEV1 and FVC)
  filter(!is.na(FEV1) | !is.na(FVC))
print(length(unique(final$HABCID))) #n = 2,646
final <- final %>%
  group_by(HABCID) %>%
  mutate(N = row_number()) %>% #create sequential PFT number variable
  mutate(N_PFT = n()) %>% #count number of PFTs per participant
  ungroup()
HABCfinal <- final
counts <- function(data) { #function to print counts of participants
  FFQ <- !is.na(data$FFQCALS) & data$FFQCALS != 0
  vars <- c("FEV1", "FVC", "RATIO", "FEV1_QC", "FVC_QC", "RATIO_QC", "FEV1_SMK", "FVC_SMK", "RATIO_SMK", "FEV1_PACK", "FVC_PACK", "RATIO_PACK")
  listy <- list(
    FEV1 = !is.na(data$FEV1),
    FVC = !is.na(data$FVC),
    RATIO = !is.na(data$RATIO),
    FEV1_QC = data$QCFEV1 %in% c(3, 4),
    FVC_QC = data$QCFVC %in% c(3, 4),
    RATIO_QC = data$QCFEV1 %in% c(3, 4) & data$QCFVC %in% c(3, 4),
    SMK = !is.na(data$SMKCAT),
    PACK = !is.na(data$PACKYR1)
  )
  vartie <- function(base, ladder) {
    sum(!duplicated(data$HABCID[base & Reduce(`&`, ladder)]))
  }
  diet <- sum(!duplicated(data$HABCID[FFQ]))
  cat("Total participants:", diet, "\n")
  for (var in vars) {
    ladder <- switch(var,
                        FEV1 = list(listy$FEV1),
                        FVC = list(listy$FVC),
                        RATIO = list(listy$RATIO),
                        FEV1_QC = list(listy$FEV1, listy$FEV1_QC),
                        FVC_QC = list(listy$FVC, listy$FVC_QC),
                        RATIO_QC = list(listy$RATIO, listy$RATIO_QC),
                        FEV1_SMK = list(listy$FEV1, listy$FEV1_QC, listy$SMK),
                        FVC_SMK = list(listy$FVC, listy$FVC_QC, listy$SMK),
                        RATIO_SMK = list(listy$RATIO, listy$RATIO_QC, listy$SMK),
                        FEV1_PACK = list(listy$FEV1, listy$FEV1_QC, listy$SMK, listy$PACK),
                        FVC_PACK = list(listy$FVC, listy$FVC_QC, listy$SMK, listy$PACK),
                        RATIO_PACK = list(listy$RATIO, listy$RATIO_QC, listy$SMK, listy$PACK)
    )
    n <- vartie(FFQ, ladder)
    cat("Total participants with", var, "data:", n, "\n")
  }
}
counts(HABCfinal) #print counts
#save final dataset as SAS compatible file and as CSV
write.foreign(HABCfinal, datafile = "HABCfinal.dat",
              codefile="HABCfinal.sas", package="SAS")
write_csv(HABCfinal, "HABCfinal.csv")
rm(list = setdiff(ls(), "ras")) #remove all objects except 'ras'

###RAS COHORT###
df <- ras %>% #filter dataset to only be 1 observation and no missing diet data to calculate diet scores
  arrange (participantid, pftdate) %>%
  group_by(participantid) %>%
  slice(1) %>%
  ungroup() %>%
  filter(!is.na(calories) & calories != 0)
#function to calculate dietary pattern scores
pk <- function(freq, por, kcal, porsz, mixed) {
  var <- ((((freq*por)/365)/kcal)*1000)*porsz*mixed
  return(var)
}
###HEI-2020###
#Total fruit
df$ffrq11pkf <- mapply(pk, df$ffrq11, df$fpor11, df$calories, 1/8, 0.05)
df$ffrq52pk <- mapply(pk, df$ffrq52, df$fpor52, df$calories, 1/2, 1)
df$ffrq53pk <- mapply(pk, df$ffrq53, df$fpor53, df$calories, 3/4, 1)
df$ffrq54pk <- mapply(pk, df$ffrq54, df$fpor54, df$calories, 1/2, 1)
df$ffrq55pk <- mapply(pk, df$ffrq55, df$fpor55, df$calories, 1/2, 1)
df$ffrq56pk <- mapply(pk, df$ffrq56, df$fpor56, df$calories, 1/2, 1)
df$ffrq57pk <- mapply(pk, df$ffrq57, df$fpor57, df$calories, 3/4, 1)
df$ffrq58pk <- mapply(pk, df$ffrq58, df$fpor58, df$calories, 1, 1)
df$ffrq59pk <- mapply(pk, df$ffrq59, df$fpor59, df$calories, 1, 1)
df$ffrq60pk <- mapply(pk, df$ffrq60, df$fpor60, df$calories, 2, 1)
df$ffrq61pk <- mapply(pk, df$ffrq61, df$fpor61, df$calories, 1, 1)
df$bfrq7pk <- mapply(pk, df$bfrq7, df$bpor7, df$calories, 1, 1)
df$bfrq8pk <- mapply(pk, df$bfrq8, df$bpor8, df$calories, 1, 1)
df$bfrq9pk <- mapply(pk, df$bfrq9, df$bpor9, df$calories, 1, 1)
df$sumtotalfruit <- rowSums(df[,c("ffrq11pkf", "ffrq52pk", "ffrq53pk", "ffrq54pk",
                                  "ffrq55pk", "ffrq56pk", "ffrq57pk", "ffrq58pk",
                                  "ffrq59pk", "ffrq60pk", "ffrq61pk", "bfrq7pk",
                                  "bfrq8pk", "bfrq9pk")], na.rm = TRUE)
#Whole fruit
df$sumwholefruit <- rowSums(df[,c("ffrq52pk", "ffrq53pk", "ffrq54pk", "ffrq55pk",
                                  "ffrq56pk", "ffrq57pk", "ffrq58pk", "ffrq59pk",
                                  "ffrq60pk", "ffrq61pk")], na.rm = TRUE)
#Total vegetables
df$ffrq26pkv <- mapply(pk, df$ffrq26, df$fpor26, df$calories, 2, 0.25)
df$ffrq36pkv <- mapply(pk, df$ffrq36, df$fpor36, df$calories, 2, 0.10)
df$ffrq37pkv <- mapply(pk, df$ffrq37, df$fpor37, df$calories, 2, 0.25)
df$ffrq38pkv <- mapply(pk, df$ffrq38, df$fpor38, df$calories, 2, 0.05)
df$ffrq39pkv <- mapply(pk, df$ffrq39, df$fpor39, df$calories, 2, 0.5)
df$ffrq41pkv <- mapply(pk, df$ffrq41, df$fpor41, df$calories, 1.5, 0.2)
df$ffrq42pkv <- mapply(pk, df$ffrq42, df$fpor42, df$calories, 2, 0.05)
df$ffrq43pkv <- mapply(pk, df$ffrq43, df$fpor43, df$calories, 0.5, 1)
df$ffrq44pkv <- mapply(pk, df$ffrq44, df$fpor44, df$calories, 2, 0.05)
df$ffrq45pkv <- mapply(pk, df$ffrq45, df$fpor45, df$calories, 2, 0.15)
df$ffrq46pk <- mapply(pk, df$ffrq46, df$fpor46, df$calories, 2, 1)
df$ffrq47pkv <- mapply(pk, df$ffrq47, df$fpor47, df$calories, 2, 0.15)
df$ffrq48pkv <- mapply(pk, df$ffrq48, df$fpor48, df$calories, 2, 1)
df$ffrq49pkv <- mapply(pk, df$ffrq49, df$fpor49, df$calories, 2, 1)
df$ffrq62pk <- mapply(pk, df$ffrq62, df$fpor62, df$calories, 2, 1)
df$ffrq64pk <- mapply(pk, df$ffrq64, df$fpor64, df$calories, 3/4, 1)
df$ffrq65pk <- mapply(pk, df$ffrq65, df$fpor65, df$calories, 1, 1)
df$ffrq66pk <- mapply(pk, df$ffrq66, df$fpor66, df$calories, 0.5, 1)
df$ffrq67pk <- mapply(pk, df$ffrq67, df$fpor67, df$calories, 0.5, 1)
df$ffrq68pk <- mapply(pk, df$ffrq68, df$fpor68, df$calories, 1, 1)
df$ffrq69pk <- mapply(pk, df$ffrq69, df$fpor69, df$calories, 1, 1)
df$ffrq70pkv <- mapply(pk, df$ffrq70, df$fpor70, df$calories, 1, 1)
df$ffrq71pkv <- mapply(pk, df$ffrq71, df$fpor71, df$calories, 1, 1)
df$ffrq72pk <- mapply(pk, df$ffrq72, df$fpor72, df$calories, 1, 1)
df$ffrq73pk <- mapply(pk, df$ffrq73, df$fpor73, df$calories, 3/4, 1)
df$ffrq74pk <- mapply(pk, df$ffrq74, df$fpor74, df$calories, 1, 1)
df$ffrq75pk <- mapply(pk, df$ffrq75, df$fpor75, df$calories, 3/4, 1)
df$ffrq76pk <- mapply(pk, df$ffrq76, df$fpor76, df$calories, 1, 1)
df$ffrq77pk <- mapply(pk, df$ffrq77, df$fpor77, df$calories, 1/2, 1)
df$ffrq78pk <- mapply(pk, df$ffrq78, df$fpor78, df$calories, 1/48, 1)
df$ffrq79pk <- mapply(pk, df$ffrq79, df$fpor79, df$calories, 1/2, 1)
df$ffrq81pk <- mapply(pk, df$ffrq81, df$fpor81, df$calories, 1, 1)
df$ffrq82pkv <- mapply(pk, df$ffrq82, df$fpor82, df$calories, 1, 1)
df$ffrq83pkv <- mapply(pk, df$ffrq83, df$fpor83, df$calories, 1, 1)
df$ffrq84pk <- mapply(pk, df$ffrq84, df$fpor84, df$calories, 1, 0.5)
df$ffrq85pkv <- mapply(pk, df$ffrq85, df$fpor85, df$calories, 1, 0.5)
df$ffrq95pk <- mapply(pk, df$ffrq95, df$fpor95, df$calories, 1/2, 1)
df$bfrq6pk <- mapply(pk, df$bfrq6, df$bpor6, df$calories, 1, 1)
df$sumtotalveg <- rowSums(df[,c("ffrq26pkv", "ffrq36pkv", "ffrq37pkv", "ffrq38pkv", "ffrq39pkv",
                                "ffrq41pkv", "ffrq42pkv", "ffrq43pkv", "ffrq44pkv", "ffrq45pkv",
                                "ffrq46pk", "ffrq47pkv", "ffrq48pkv", "ffrq49pkv", "ffrq62pk",
                                "ffrq64pk", "ffrq65pk", "ffrq66pk", "ffrq67pk",
                                "ffrq68pk", "ffrq69pk", "ffrq70pkv", "ffrq71pkv",
                                "ffrq72pk", "ffrq73pk", "ffrq74pk", "ffrq75pk",
                                "ffrq76pk", "ffrq77pk", "ffrq78pk", "ffrq79pk",
                                "ffrq81pk", "ffrq82pkv", "ffrq83pkv", "ffrq84pk",
                                "ffrq85pkv", "ffrq95pk", "bfrq6pk")], na.rm = TRUE)
#Greens and beans
df$sumgreensandbeans <- rowSums(df[,c("ffrq43pkv", "ffrq48pkv", "ffrq49pkv", "ffrq62pk", "ffrq68pk",
                                      "ffrq69pk", "ffrq70pkv", "ffrq71pkv", "ffrq76pk", "ffrq82pkv",
                                      "ffrq83pkv", "ffrq84pk")], na.rm = TRUE)
#Whole grains
df$ffrq1pk <- mapply(pk, df$ffrq1, df$fpor1, df$calories, 1.5, 1)
df$ffrq2pk <- mapply(pk, df$ffrq2, df$fpor2, df$calories, 3, 1)
df$ffrq7pk <- mapply(pk, df$ffrq7, df$fpor7, df$calories, 1, 1)
df$ffrq11pkw <- mapply(pk, df$ffrq11, df$fpor11, df$calories, 1, 0.2)
df$ffrq12pkw <- mapply(pk, df$ffrq12, df$fpor12, df$calories, 1, 0.25)
df$ffrq15pk <- mapply(pk, df$ffrq15, df$fpor15, df$calories, 1, 1)
df$ffrq16pk <- mapply(pk, df$ffrq16, df$fpor16, df$calories, 1, 1)
df$ffrq86pkw <- mapply(pk, df$ffrq86, df$fpor86, df$calories, 2, 0.05)
df$sumwholegrains <- rowSums(df[,c("ffrq1pk", "ffrq2pk", "ffrq7pk", "ffrq11pkw",
                                   "ffrq12pkw", "ffrq15pk", "ffrq16pk", "ffrq86pkw")], na.rm = TRUE)
#Dairy
df$ffrq3pk <- mapply(pk, df$ffrq3, df$fpor3, df$calories, 3/4, 1)
df$ffrq36pkd <- mapply(pk, df$ffrq36, df$fpor36, df$calories, 2, 0.05)
df$ffrq38pkd <- mapply(pk, df$ffrq38, df$fpor38, df$calories, 2, 0.05)
df$ffrq40pkd <- mapply(pk, df$ffrq40, df$fpor40, df$calories, 2, 0.33)
df$ffrq42pkd <- mapply(pk, df$ffrq42, df$fpor42, df$calories, 2, 0.20)
df$ffrq44pkd <- mapply(pk, df$ffrq44, df$fpor44, df$calories, 2, 0.15)
df$ffrq45pkd <- mapply(pk, df$ffrq45, df$fpor45, df$calories, 2, 0.15)
df$ffrq47pkd <- mapply(pk, df$ffrq47, df$fpor47, df$calories, 2, 0.5)
df$ffrq85pkd <- mapply(pk, df$ffrq85, df$fpor85, df$calories, 1, 0.15)
df$ffrq88pk <- mapply(pk, df$ffrq88, df$fpor88, df$calories, 1, 1)
df$ffrq89pk <- mapply(pk, df$ffrq89, df$fpor89, df$calories, 1, 1)
df$ffrq90pk <- mapply(pk, df$ffrq90, df$fpor90, df$calories, 1, 1)
df$ffrq91pk <- mapply(pk, df$ffrq91, df$fpor91, df$calories, 1, 1)
df$ffrq92pk <- mapply(pk, df$ffrq92, df$fpor92, df$calories, 1/2, 1)
df$bfrq1pk <- mapply(pk, df$bfrq1, df$bpor1, df$calories, 1, 1)
df$bfrq2pk <- mapply(pk, df$bfrq2, df$bpor2, df$calories, 1, 1)
df$bfrq5pk <- mapply(pk, df$bfrq5, df$bpor5, df$calories, 1/16, 0.5)
df$bfrq10pk <- mapply(pk, df$bfrq10, df$bpor10, df$calories, 1, 1)
df$sumdairy <- rowSums(df[,c("ffrq3pk", "ffrq36pkd", "ffrq38pkd", "ffrq40pkd", "ffrq42pkd",
                             "ffrq44pkd", "ffrq45pkd", "ffrq47pkd", "ffrq85pkd",
                             "ffrq88pk", "ffrq89pk", "ffrq90pk", "ffrq91pk",
                             "ffrq92pk", "bfrq1pk", "bfrq2pk", "bfrq5pk", "bfrq10pk")], na.rm = TRUE)
#Total protein
df$ffrq12pkp <- mapply(pk, df$ffrq12, df$fpor12, df$calories, 1, 0.25)
df$ffrq19pk <- mapply(pk, df$ffrq19, df$fpor19, df$calories, 4, 1)
df$ffrq20pk <- mapply(pk, df$ffrq20, df$fpor20, df$calories, 2, 1)
df$ffrq21pk <- mapply(pk, df$ffrq21, df$fpor21, df$calories, 2, 1)
df$ffrq22pk <- mapply(pk, df$ffrq22, df$fpor22, df$calories, 4, 1)
df$ffrq23pk <- mapply(pk, df$ffrq23, df$fpor23, df$calories, 4, 1)
df$ffrq24pk <- mapply(pk, df$ffrq24, df$fpor24, df$calories, 2, 1)
df$ffrq25pk <- mapply(pk, df$ffrq25, df$fpor25, df$calories, 2, 1)
df$ffrq26pkp <- mapply(pk, df$ffrq26, df$fpor26, df$calories, 16, 0.75)
df$ffrq27pk <- mapply(pk, df$ffrq27, df$fpor27, df$calories, 6, 1)
df$ffrq28pk <- mapply(pk, df$ffrq28, df$fpor28, df$calories, 3, 1)
df$ffrq29pk <- mapply(pk, df$ffrq29, df$fpor29, df$calories, 6, 1)
df$ffrq30pk <- mapply(pk, df$ffrq30, df$fpor30, df$calories, 8, 1)
df$ffrq31pk <- mapply(pk, df$ffrq31, df$fpor31, df$calories, 8, 1)
df$ffrq32pkp <- mapply(pk, df$ffrq32, df$fpor32, df$calories, 4, 0.5)
df$ffrq33pk <- mapply(pk, df$ffrq33, df$fpor33, df$calories, 4, 1)
df$ffrq34pk <- mapply(pk, df$ffrq34, df$fpor34, df$calories, 6, 1)
df$ffrq35pk <- mapply(pk, df$ffrq35, df$fpor35, df$calories, 6, 1)
df$ffrq36pkp <- mapply(pk, df$ffrq36, df$fpor36, df$calories, 16, 0.45)
df$ffrq37pkp <- mapply(pk, df$ffrq37, df$fpor37, df$calories, 16, 0.75)
df$ffrq38pkp <- mapply(pk, df$ffrq38, df$fpor38, df$calories, 16, 0.2)
df$ffrq41pkp <- mapply(pk, df$ffrq41, df$fpor41, df$calories, 12, 0.50)
df$ffrq42pkp <- mapply(pk, df$ffrq42, df$fpor42, df$calories, 16, 0.15)
df$ffrq43pkp <- mapply(pk, df$ffrq43, df$fpor43, df$calories, 4, 1)
df$ffrq44pkp <- mapply(pk, df$ffrq44, df$fpor44, df$calories, 10, 0.2)
df$ffrq45pkp <- mapply(pk, df$ffrq45, df$fpor45, df$calories, 15, 0.35)
df$ffrq48pkp <- mapply(pk, df$ffrq48, df$fpor48, df$calories, 16, 1)
df$ffrq49pkp <- mapply(pk, df$ffrq49, df$fpor49, df$calories, 16, 1)
df$ffrq51pkp <- mapply(pk, df$ffrq51, df$fpor51, df$calories, 16, 0.33)
df$ffrq70pkp <- mapply(pk, df$ffrq70, df$fpor70, df$calories, 4, 1)
df$ffrq71pkp <- mapply(pk, df$ffrq71, df$fpor71, df$calories, 4, 1)
df$ffrq82pkp <- mapply(pk, df$ffrq82, df$fpor82, df$calories, 4, 1)
df$ffrq83pkp <- mapply(pk, df$ffrq83, df$fpor83, df$calories, 4, 1)
df$ffrq85pkp <- mapply(pk, df$ffrq85, df$fpor85, df$calories, 8, 0.25)
df$ffrq93pkp <- mapply(pk, df$ffrq93, df$fpor93, df$calories, 4, 0.5)
df$sumtotalprotein <- rowSums(df[,c("ffrq12pkp", "ffrq19pk", "ffrq20pk", "ffrq21pk",
                                    "ffrq22pk", "ffrq23pk", "ffrq24pk", "ffrq25pk",
                                    "ffrq26pkp", "ffrq27pk", "ffrq28pk", "ffrq29pk",
                                    "ffrq30pk", "ffrq31pk", "ffrq32pkp", "ffrq33pk",
                                    "ffrq34pk", "ffrq35pk", "ffrq36pkp", "ffrq37pkp",
                                    "ffrq38pkp", "ffrq41pkp", "ffrq42pkp", "ffrq43pkp",
                                    "ffrq44pkp", "ffrq45pkp", "ffrq48pkp", "ffrq49pkp",
                                    "ffrq51pkp", "ffrq82pkp", "ffrq83pkp", "ffrq85pkp", "ffrq93pkp")], na.rm = TRUE)
#Seafood and plant proteins
df$sumseafoodandplant <- rowSums(df[,c("ffrq19pk", "ffrq26pkp", "ffrq32pkp", "ffrq33pk",
                                       "ffrq34pk", "ffrq35pk", "ffrq43pkp", "ffrq48pkp",
                                       "ffrq49pkp", "ffrq70pkp", "ffrq71pkp", "ffrq82pkp", "ffrq83pkp")], na.rm = TRUE)
#Refined grains
df$ffrq4pk <- mapply(pk, df$ffrq4, df$fpor4, df$calories, 2, 1)
df$ffrq5pk <- mapply(pk, df$ffrq5, df$fpor5, df$calories, 1, 1)
df$ffrq6pk <- mapply(pk, df$ffrq6, df$fpor6, df$calories, 1, 1)
df$ffrq8pk <- mapply(pk, df$ffrq8, df$fpor8, df$calories, 1, 1)
df$ffrq11pkr <- mapply(pk, df$ffrq11, df$fpor11, df$calories, 1, 0.75)
df$ffrq12pkr <- mapply(pk, df$ffrq12, df$fpor12, df$calories, 1, 0.5)
df$ffrq17pk <- mapply(pk, df$ffrq17, df$fpor17, df$calories, 1.75, 1)
df$ffrq18pk <- mapply(pk, df$ffrq18, df$fpor18, df$calories, 2, 1)
df$ffrq32pkr <- mapply(pk, df$ffrq32, df$fpor32, df$calories, 4, 0.5)
df$ffrq36pkr <- mapply(pk, df$ffrq36, df$fpor36, df$calories, 16, 0.4)
df$ffrq38pkr <- mapply(pk, df$ffrq38, df$fpor38, df$calories, 16, 0.7)
df$ffrq39pkr <- mapply(pk, df$ffrq39, df$fpor39, df$calories, 16, 0.5)
df$ffrq40pkr <- mapply(pk, df$ffrq40, df$fpor40, df$calories, 16, 0.66)
df$ffrq41pkr <- mapply(pk, df$ffrq41, df$fpor41, df$calories, 12, 0.4)
df$ffrq42pkr <- mapply(pk, df$ffrq42, df$fpor42, df$calories, 16, 0.6)
df$ffrq44pkr <- mapply(pk, df$ffrq44, df$fpor44, df$calories, 10, 0.6)
df$ffrq45pkr <- mapply(pk, df$ffrq45, df$fpor45, df$calories, 16, 0.35)
df$ffrq47pkr <- mapply(pk, df$ffrq47, df$fpor47, df$calories, 2, 0.10)
df$ffrq50pk <- mapply(pk, df$ffrq50, df$fpor50, df$calories, 12, 1)
df$ffrq51pkr <- mapply(pk, df$ffrq51, df$fpor51, df$calories, 16, 0.66)
df$ffrq80pk <- mapply(pk, df$ffrq80, df$fpor80, df$calories, 4, 1)
df$ffrq85pkr <- mapply(pk, df$ffrq85, df$fpor85, df$calories, 4, 0.15)
df$ffrq86pkr <- mapply(pk, df$ffrq86, df$fpor86, df$calories, 2, 0.95)
df$ffrq93pkr <- mapply(pk, df$ffrq93, df$fpor93, df$calories, 4, 0.5)
df$ffrq100pk <- mapply(pk, df$ffrq100, df$fpor100, df$calories, 0.5, 1)
df$ffrq101pk <- mapply(pk, df$ffrq101, df$fpor101, df$calories, 1.5, 1)
df$sumrefinedgrains <- rowSums(df[,c("ffrq4pk", "ffrq5pk", "ffrq6pk", "ffrq8pk",
                                     "ffrq11pkr", "ffrq12pkr", "ffrq17pk", "ffrq18pk",
                                     "ffrq32pkr", "ffrq36pkr", "ffrq38pkr",
                                     "ffrq39pkr", "ffrq40pkr", "ffrq41pkr", "ffrq42pkr",
                                     "ffrq44pkr", "ffrq45pkr", "ffrq47pkr", "ffrq50pk", "ffrq51pkr", "ffrq80pk",
                                     "ffrq85pkr", "ffrq86pkr", "ffrq93pkr", "ffrq100pk", "ffrq101pk")], na.rm = TRUE)
#Fatty acid ratio
df$faratio <- (df$mfatot+df$pfatot)/df$sfatot
#Sodium
df$totalsodium <- ((df$sodium/1000)/df$calories)*1000
#Added sugar
df$addsug <- ((df$addsugar*3.4)/df$calories)*100
#Saturated fat
df$satfat <- ((df$sfatot*9)/df$calories)*100
#HEI-2020 score calculations
df$totalfruitscore <- pmin((df$sumtotalfruit * 5) / 0.8, 5)
df$wholefruitscore <- pmin((df$sumwholefruit * 5) / 0.4, 5)
df$totalvegscore <- pmin((df$sumtotalveg * 5) / 1.1, 5)
df$greensandbeansscore <- pmin((df$sumgreensandbeans * 5) / 0.2, 5)
df$wholegrainsscore <- pmin((df$sumwholegrains * 10) / 1.5, 10)
df$dairyscore <- pmin((df$sumdairy * 10) / 1.3, 10)
df$totalproteinscore <- pmin((df$sumtotalprotein * 5) / 2.5, 5)
df$seafoodandplantscore <- pmin((df$sumseafoodandplant * 5) / 0.8, 5)
df$refinedgrainsscore <- pmax(pmin(10 - (10 * (df$sumrefinedgrains - 1.8) / (4.3 - 1.8)), 10), 0)
df$fattyacidscore <- ifelse(df$sfatot == 0 & (df$mfatot + df$pfatot) == 0, 0, pmin(pmax(10 * ((df$faratio - 1.2) / (2.5 - 1.2)), 0), 10))
df$sodiumscore <- pmax(pmin(10 - (10 * (df$totalsodium - 1.1) / (2.0 - 1.1)), 10), 0)
df$addsugarscore <- pmax(pmin(10 - (10 * (df$addsug - 6.5) / (26 - 6.5)), 10), 0)
df$satfatscore <- pmax(pmin(10 - (10 * (df$satfat - 8) / (16 - 8)), 10), 0)
HEIscores <- c("totalfruitscore", "wholefruitscore", "totalvegscore", "greensandbeansscore",
               "wholegrainsscore", "dairyscore", "totalproteinscore", "seafoodandplantscore", 
               "refinedgrainsscore", "fattyacidscore", "sodiumscore", "addsugarscore", "satfatscore")
df[HEIscores] <- lapply(df[HEIscores], function(x) ifelse(df$calories == 0 | is.na(df$calories), NA, x))
df$HEI2020score <- ifelse(df$calories == 0 | is.na(df$calories), NA, 
                          rowSums(df[HEIscores], na.rm = TRUE))
mean(df$HEI2020score, na.rm = TRUE)

###DASH SCORE###
#Fruits
df$sumfruitdash <- df$sumtotalfruit
#Vegetables
df$sumvegdash <- rowSums(df[,c("ffrq26pkv", "ffrq36pkv", "ffrq38pkv", "ffrq39pkv", "ffrq41pkv", 
                               "ffrq42pkv", "ffrq46pk", "ffrq49pkv", "ffrq62pk",
                               "ffrq64pk", "ffrq65pk", "ffrq66pk", "ffrq67pk",
                               "ffrq68pk", "ffrq69pk", "ffrq72pk", "ffrq73pk",
                               "ffrq74pk", "ffrq75pk", "ffrq76pk", "ffrq77pk",
                               "ffrq78pk", "ffrq79pk", "ffrq84pk", "ffrq95pk", "bfrq6pk")], na.rm = TRUE)
#Nuts and legumes
df$ffrq19pkl <- mapply(pk, df$ffrq19, df$fpor19, df$calories, 1/2, 1)
df$sumnutslegdash <- rowSums(df[,c("ffrq19pkl", "ffrq37pkv", "ffrq43pkv", "ffrq44pkv",
                               "ffrq45pkv", "ffrq48pkv", "ffrq49pkv", "ffrq70pkv", "ffrq71pkv",
                               "ffrq82pkv", "ffrq83pkv")], na.rm = TRUE)
#Whole grains
df$sumwholegrainsdash <- df$sumwholegrains
#Low-fat dairy
df <- df %>%
  mutate(
    bfrq1_L = case_when( #low-fat dairy version of milk as a beverage
      bfrq1 == 0 ~ 0, #if no overall milk consumption, no low-fat milk consumption
      bfrq1 > 0 & (adj8a_3 == 1 | adj8a_4 == 1 | adj8a_5 == 1) & (is.na(adj8a_1) & is.na(adj8a_2)) ~ bfrq1, #if milk consumption is only low-fat, then low-fat milk consumption is set to overall milk consumption
      bfrq1 > 0 & (adj8a_1 == 1 | adj8a_2 == 1) & (is.na(adj8a_3) & is.na(adj8a_4) & is.na(adj8a_5)) ~ 0, #if milk consumption is only high-fat, then low-fat milk consumption is set to 0
      bfrq1 > 0 & (adj8a_3 == 1 | adj8a_4 == 1 | adj8a_5 == 1) & (adj8a_1 == 1 | adj8a_2 == 1) ~ bfrq1 / 2, #if milk consumption is low-fat AND high-fat, then low-fat milk consumption is set to half of overall milk consumption
      is.na(bfrq1) & (adj8a_1 == 1 | adj8a_2 == 1) ~ 0, #if overall milk consumption data is missing, but they marked only consuming high-fat milk, then low-fat milk consumption is set to 0
      bfrq1 > 0 & (is.na(adj8a_1) & is.na(adj8a_2) & is.na(adj8a_3) & is.na(adj8a_4) & is.na(adj8a_5) & is.na(adj8a_6) | adj8a_6 == 1) #if no data on type of milk consumption for beverage, check type of milk consumption for cereal use;
        & (adj6a_4 == 1 | adj6a_5 == 1 | adj6a_6 == 1) & (is.na(adj6a_1) & is.na(adj6a_2) & is.na(adj6a_3) & is.na(adj6a_7)) ~ bfrq1, #if cereal use milk is only low-fat, set to low-fat
      bfrq1 > 0 & (is.na(adj8a_1) & is.na(adj8a_2) & is.na(adj8a_3) & is.na(adj8a_4) & is.na(adj8a_5) & is.na(adj8a_6) | adj8a_6 == 1) #if no data on type of milk consumption for beverage, check type of milk consumption for cereal use;
        & (adj6a_1 == 1 | adj6a_2 == 1 | adj6a_3 == 1 | adj6a_7 == 1) & (is.na(adj6a_4) & is.na(adj6a_5) & is.na(adj6a_6)) ~ 0, #if cereal use milk is only high-fat, set to 0
      bfrq1 > 0 & (is.na(adj8a_1) & is.na(adj8a_2) & is.na(adj8a_3) & is.na(adj8a_4) & is.na(adj8a_5) & is.na(adj8a_6) & is.na(adj6a_7) | adj8a_6 == 1) #if no data on type of milk consumption for beverage, check type of milk consumption for cereal use;
        & ((adj6a_1 == 1 | adj6a_2 == 1 | adj6a_3 == 1 | adj6a_7 == 1) & (adj6a_4 == 1 | adj6a_5 == 1 | adj6a_6 == 1)) ~ bfrq1/2, #if mixed-use (high & low fat), divide by 2
      bfrq1 > 0 & (is.na(adj8a_1) & is.na(adj8a_2) & is.na(adj8a_3) & is.na(adj8a_4) & is.na(adj8a_5) & is.na(adj8a_6) | adj8a_6 == 1) #if no data on type of milk consumption for beverage
        & (is.na(adj6a_1) & is.na(adj6a_2) & is.na(adj6a_3) & is.na(adj6a_4) & is.na(adj6a_5) & is.na(adj6a_6) & is.na(adj6a_7) | adj6a_8 == 1) #OR cereal, check type of milk consumption for coffee & tea use;
          & (adj7a_4 == 1 | adj7a_5 == 1 | adj7a_6 == 1) & (is.na(adj7a_1) & is.na(adj7a_2) & is.na(adj7a_3) & is.na(adj7a_7)) ~ bfrq1, #if coffee & tea use milk is only low-fat, set to low-fat
      bfrq1 > 0 & (is.na(adj8a_1) & is.na(adj8a_2) & is.na(adj8a_3) & is.na(adj8a_4) & is.na(adj8a_5) & is.na(adj8a_6) | adj8a_6 == 1) #if no data on type of milk consumption for beverage
      & (is.na(adj6a_1) & is.na(adj6a_2) & is.na(adj6a_3) & is.na(adj6a_4) & is.na(adj6a_5) & is.na(adj6a_6) & is.na(adj6a_7) | adj6a_8 == 1) #OR cereal, check type of milk consumption for coffee & tea use;
         & (adj7a_1 == 1 | adj7a_2 == 1 | adj7a_3 == 1 | adj7a_7 == 1) & (is.na(adj7a_4) & is.na(adj7a_5) & is.na(adj7a_6)) ~ 0, #if coffee & tea use milk is only high-fat, set to 0
      bfrq1 > 0 & (is.na(adj8a_1) & is.na(adj8a_2) & is.na(adj8a_3) & is.na(adj8a_4) & is.na(adj8a_5) & is.na(adj8a_6) | adj8a_6 == 1) #if no data on type of milk consumption for beverage
       & (is.na(adj6a_1) & is.na(adj6a_2) & is.na(adj6a_3) & is.na(adj6a_4) & is.na(adj6a_5) & is.na(adj6a_6) & is.na(adj6a_7) | adj6a_8 == 1) #OR cereal, check type of milk consumption for coffee & tea use;
         & ((adj7a_1 == 1 | adj7a_2 == 1 | adj7a_3 == 1 | adj7a_7 == 1) & (adj7a_4 == 1 | adj7a_5 == 1 | adj7a_6 == 1)) ~ bfrq1/2, #if mixed-use (high & low fat), divide by 2
      TRUE ~ NA_real_  #all other results NA (NAs for bfrq1 - overall milk consumption, or all milk types are NA)
    ),
    ffrq3_L = case_when( #low-fat dairy version of milk consumed on cereals
      ffrq3 == 0 ~ 0, #similar logic flow to bfrq1_L
      ffrq3 > 0 & (adj6a_4 == 1 | adj6a_5 == 1 | adj6a_6 == 1) & (is.na(adj6a_1) & is.na(adj6a_2) & is.na(adj6a_3) & is.na(adj6a_7)) ~ ffrq3, 
      ffrq3 > 0 & (adj6a_1 == 1 | adj6a_2 == 1 | adj6a_3 == 1 | adj6a_7 == 1) & (is.na(adj6a_4) & is.na(adj6a_5) & is.na(adj6a_6)) ~ 0,
      ffrq3 > 0 & (adj6a_1 == 1 | adj6a_2 == 1 | adj6a_3 == 1 | adj6a_7 == 1) & (adj6a_4 == 1 | adj6a_5 == 1 | adj6a_6 == 1) ~ ffrq3 / 2, 
      is.na(ffrq3) & (adj6a_1 == 1 | adj6a_2 == 1 | adj6a_3 == 1 | adj6a_7 == 1) ~ 0,
      ffrq3 > 0 & (is.na(adj6a_1) & is.na(adj6a_2) & is.na(adj6a_3) & is.na(adj6a_4) & is.na(adj6a_5) & is.na(adj6a_6) & is.na(adj6a_7) | adj6a_8 == 1)
        & (adj8a_3 == 1 | adj8a_4 == 1 | adj8a_5 == 1) & (is.na(adj8a_1) & is.na(adj8a_2)) ~ ffrq3,
      ffrq3 > 0 & (is.na(adj6a_1) & is.na(adj6a_2) & is.na(adj6a_3) & is.na(adj6a_4) & is.na(adj6a_5) & is.na(adj6a_6) & is.na(adj6a_7) | adj6a_8 == 1)
        & (adj8a_1 == 1 | adj8a_2 == 1) & (is.na(adj8a_3) & is.na(adj8a_4) & is.na(adj8a_5)) ~ 0,
      ffrq3 > 0 & (is.na(adj6a_1) & is.na(adj6a_2) & is.na(adj6a_3) & is.na(adj6a_4) & is.na(adj6a_5) & is.na(adj6a_6) & is.na(adj6a_7) | adj6a_8 == 1)
        & ((adj8a_3 == 1 | adj8a_4 == 1 | adj8a_5 == 1) & (adj8a_1 == 1 | adj8a_2 == 1)) ~ ffrq3/2,
      ffrq3 > 0 & (is.na(adj6a_1) & is.na(adj6a_2) & is.na(adj6a_3) & is.na(adj6a_4) & is.na(adj6a_5) & is.na(adj6a_6) & is.na(adj6a_7) | adj6a_8 == 1)
       & (is.na(adj8a_1) & is.na(adj8a_2) & is.na(adj8a_3) & is.na(adj8a_4) & is.na(adj8a_5) & is.na(adj8a_6) | adj8a_6 == 1)
          & (adj7a_4 == 1 | adj7a_5 == 1 | adj7a_6 == 1) & (is.na(adj7a_1) & is.na(adj7a_2) & is.na(adj7a_3) & is.na(adj7a_7)) ~ ffrq3,
      ffrq3 > 0 & (is.na(adj6a_1) & is.na(adj6a_2) & is.na(adj6a_3) & is.na(adj6a_4) & is.na(adj6a_5) & is.na(adj6a_6) & is.na(adj6a_7) | adj6a_8 == 1)
       & (is.na(adj8a_1) & is.na(adj8a_2) & is.na(adj8a_3) & is.na(adj8a_4) & is.na(adj8a_5) & is.na(adj8a_6) | adj8a_6 == 1)
          & (adj7a_1 == 1 | adj7a_2 == 1 | adj7a_3 == 1 | adj7a_7 == 1) & (is.na(adj7a_4) & is.na(adj7a_5) & is.na(adj7a_6)) ~ 0,
      ffrq3 > 0 & (is.na(adj6a_1) & is.na(adj6a_2) & is.na(adj6a_3) & is.na(adj6a_4) & is.na(adj6a_5) & is.na(adj6a_6) & is.na(adj6a_7) | adj6a_8 == 1)
        & (is.na(adj8a_1) & is.na(adj8a_2) & is.na(adj8a_3) & is.na(adj8a_4) & is.na(adj8a_5) & is.na(adj8a_6) | adj8a_6 == 1)
          & ((adj7a_1 == 1 | adj7a_2 == 1 | adj7a_3 == 1 | adj7a_7 == 1) & (adj7a_4 == 1 | adj7a_5 == 1 | adj7a_6 == 1)) ~ ffrq3/2,
      TRUE ~ NA_real_ 
    ),
    bfrq5_L = case_when( #low-fat dairy version of milk consumed with coffee/tea
      bfrq5 == 0 ~ 0, #similar logic flow to bfrq1_L
      bfrq5 > 0 & (adj7a_4 == 1 | adj7a_5 == 1 | adj7a_6 == 1) & (is.na(adj7a_1) & is.na(adj7a_2) & is.na(adj7a_3) & is.na(adj7a_7)) ~ bfrq5, 
      bfrq5 > 0 & (adj7a_1 == 1 | adj7a_2 == 1 | adj7a_3 == 1 | adj7a_7 == 1) & (is.na(adj7a_4) & is.na(adj7a_5) & is.na(adj7a_6)) ~ 0,
      bfrq5 > 0 & (adj7a_1 == 1 | adj7a_2 == 1 | adj7a_3 == 1 | adj7a_7 == 1) & (adj7a_4 == 1 | adj7a_5 == 1 | adj7a_6 == 1) ~ bfrq5 / 2, 
      is.na(bfrq5) & (adj7a_1 == 1 | adj7a_2 == 1 | adj7a_3 == 1 | adj7a_7 == 1) ~ 0,
      bfrq5 > 0 & (is.na(adj7a_1) & is.na(adj7a_2) & is.na(adj7a_3) & is.na(adj7a_4) & is.na(adj7a_5) & is.na(adj7a_6) & is.na(adj7a_7) | adj7a_8 == 1)
        & (adj8a_3 == 1 | adj8a_4 == 1 | adj8a_5 == 1) & (is.na(adj8a_1) & is.na(adj8a_2)) ~ bfrq5,
      bfrq5 > 0 & (is.na(adj7a_1) & is.na(adj7a_2) & is.na(adj7a_3) & is.na(adj7a_4) & is.na(adj7a_5) & is.na(adj7a_6) & is.na(adj7a_7) | adj7a_8 == 1)
        & (adj8a_1 == 1 | adj8a_2 == 1) & (is.na(adj8a_3) & is.na(adj8a_4) & is.na(adj8a_5)) ~ 0,
      bfrq5 > 0 & (is.na(adj7a_1) & is.na(adj7a_2) & is.na(adj7a_3) & is.na(adj7a_4) & is.na(adj7a_5) & is.na(adj7a_6) & is.na(adj7a_7) | adj7a_8 == 1)
        & ((adj8a_3 == 1 | adj8a_4 == 1 | adj8a_5 == 1) & (adj8a_1 == 1 | adj8a_2 == 1)) ~ bfrq5/2,
      bfrq5 > 0 & (is.na(adj7a_1) & is.na(adj7a_2) & is.na(adj7a_3) & is.na(adj7a_4) & is.na(adj7a_5) & is.na(adj7a_6) & is.na(adj7a_7) | adj7a_8 == 1)
        & (is.na(adj8a_1) & is.na(adj8a_2) & is.na(adj8a_3) & is.na(adj8a_4) & is.na(adj8a_5) & is.na(adj8a_6) | adj8a_6 == 1)
          & (adj6a_4 == 1 | adj6a_5 == 1 | adj6a_6 == 1) & (is.na(adj6a_1) & is.na(adj6a_2) & is.na(adj6a_3) & is.na(adj6a_7)) ~ bfrq5,
      bfrq5 > 0 & (is.na(adj7a_1) & is.na(adj7a_2) & is.na(adj7a_3) & is.na(adj7a_4) & is.na(adj7a_5) & is.na(adj7a_6) & is.na(adj7a_7) | adj7a_8 == 1)
       & (is.na(adj8a_1) & is.na(adj8a_2) & is.na(adj8a_3) & is.na(adj8a_4) & is.na(adj8a_5) & is.na(adj8a_6) | adj8a_6 == 1)
         & (adj6a_1 == 1 | adj6a_2 == 1 | adj6a_3 == 1 | adj6a_7 == 1) & (is.na(adj6a_4) & is.na(adj6a_5) & is.na(adj6a_6)) ~ 0,
      bfrq5 > 0 & (is.na(adj7a_1) & is.na(adj7a_2) & is.na(adj7a_3) & is.na(adj7a_4) & is.na(adj7a_5) & is.na(adj7a_6) & is.na(adj7a_7) | adj7a_8 == 1)
       & (is.na(adj8a_1) & is.na(adj8a_2) & is.na(adj8a_3) & is.na(adj8a_4) & is.na(adj8a_5) & is.na(adj8a_6) | adj8a_6 == 1)
         & ((adj6a_1 == 1 | adj6a_2 == 1 | adj6a_3 == 1 | adj6a_7 == 1) & (adj6a_4 == 1 | adj6a_5 == 1 | adj6a_6 == 1)) ~ bfrq5/2,
      TRUE ~ NA_real_ 
    )
  )
df$bfrq1_Lpk <- mapply(pk, df$bfrq1_L, df$bpor1, df$calories, 1, 1)
df$bfrq5_Lpk <- mapply(pk, df$bfrq5_L, df$bpor5, df$calories, 1/16, 1)
df$ffrq3_Lpk <- mapply(pk, df$ffrq3_L, df$fpor3, df$calories, 3/4, 1)
df$sumdairydash <- rowSums(df[,c("bfrq1_L", "bfrq5_L", "ffrq3_L",
                                 "ffrq88pk", "ffrq89pk", "ffrq91pk")], na.rm = TRUE)
#Red and processed meat
df$sumredprocdash <- rowSums(df[,c("ffrq21pk", "ffrq22pk", "ffrq23pk", "ffrq24pk",
                                   "ffrq25pk", "ffrq27pk", "ffrq28pk", "ffrq29pk",
                                   "ffrq36pkp", "ffrq37pkp", "ffrq38pkp", "ffrq41pkp",
                                   "ffrq42pkp", "ffrq44pkp", "ffrq45pkp", "ffrq93pkp")], na.rm = TRUE)

#Sugar-sweetened beverages
df$bfrq11pk <- mapply(pk, df$bfrq11, df$bpor11, df$calories, 12, 1)
df$bfrq12pk <- mapply(pk, df$bfrq12, df$bpor12, df$calories, 12, 1)
df$sumssb <- rowSums(df[,c("bfrq2pk", "bfrq9pk", "bfrq11pk", "bfrq12pk")], na.rm = TRUE)
#Sodium
df$totalsodiumdash <- ((df$sodium/1000)/df$calories)*1000
#DASH SCORE CALCULATIONS
DASHvars <- c("sumfruitdash", "sumvegdash", "sumnutslegdash", "sumwholegrainsdash", 
              "sumdairydash", "sumredprocdash", "sumssb", "totalsodiumdash")
df <- df %>% mutate_at(vars(all_of(DASHvars)), ~ifelse(is.na(.), 0, .))
quint <- lapply(DASHvars, function(x) {
  quantile(df[[x]], probs = c(0.2, 0.4, 0.6, 0.8), na.rm = TRUE)
})
DASHscores <- function(x, q, reverse = FALSE) {
  score_val <- case_when(
    x >= q[4] ~ ifelse(reverse, 1, 5),
    x >= q[3] ~ ifelse(reverse, 2, 4),
    x >= q[2] ~ ifelse(reverse, 3, 3),
    x >= q[1] ~ ifelse(reverse, 4, 2),
    TRUE ~ ifelse(reverse, 5, 1))
  return(score_val)
}
ScoreNames <- c("dashfruitscore", "dashvegscore", "dashnutslegscore", "dashwholegrainsscore", 
                "dashdairyscore", "dashredprocscore", "dashssbscore", "dashsodscore")
for (i in seq_along(DASHvars)) {
  var <- DASHvars[i]
  quintile <- quint[[i]] #set quintiles for each component
  rev <- var %in% c("sumredprocdash", "sumssb", "totalsodiumdash") #set reverse for components scored in the opposite direction
  df[[ScoreNames[i]]] <- mapply(function(value, calories) {
    if (is.na(calories) | calories == 0) return(NA) #set score to NA if FFQCALS is NA or 0
    DASHscores(value, quintile, rev) #apply scoring function
  }, df[[var]], df$calories)
}
df$DASHscore <- rowSums(df[, ScoreNames], na.rm = FALSE)
#final scores dataset to merge with main dataset later
rasscores <- df %>%
  select(participantid, HEI2020score, totalfruitscore, 
         wholefruitscore, totalvegscore, greensandbeansscore, wholegrainsscore, 
         dairyscore, totalproteinscore, seafoodandplantscore, refinedgrainsscore, 
         fattyacidscore, sodiumscore, addsugarscore, satfatscore, calories, DASHscore, dashfruitscore, dashvegscore, 
         dashnutslegscore, dashwholegrainsscore, dashdairyscore, 
         dashredprocscore, dashssbscore, dashsodscore)
###PFT data transformation###
#cut to variables needed
rascut <- ras[, c(
  "participantid", "csmoker", "eversmoker", "calories", "educat", "smoknum", "smokyr", "SMOKCIG", "birthdt", "Smoker", "CIGPDAY", "Gender",
  "race", "fev1", "fvc", "FVC_score", "FEV1_score", "fev1_2", "fvc_2", "fev1_3", "fvc_3", "InstitutionName", "pftdate", "cleanedheight_330",
  "BJ_Best_fev1", "BJ_best_fvc"
)]
#data transformation
df <- rascut
df <- df %>%
  group_by(participantid) %>% arrange(participantid, pftdate) %>%  #arrange by participant and date
  mutate(time = as.numeric(pftdate - first(pftdate), units = "days") / 365.25) %>% #create 'time' variable for years elapsed since first visit (baseline)
  ungroup()
df <- df %>% distinct(participantid, time, .keep_all = TRUE) #found a duplicate row (participantid/time pair), and removed it
df$CIGPDAY <- as.numeric(df$CIGPDAY)
for(i in 1:2) {
  df <- df %>% #harmonize smoking variables (they're kind of a mess/all over the place)
    arrange(participantid, pftdate) %>%
    mutate(SMOKCIG = ifelse(SMOKCIG == "", NA, SMOKCIG)) %>%
    group_by(participantid) %>%
    mutate(eversmoker = ifelse((first(eversmoker) == 1 & (first(csmoker) == 1 | first(Smoker) == 0 | first(SMOKCIG) == "Y" | first(CIGPDAY) >0)), 1, eversmoker)) %>% #if the first instance of a participant is marked as eversmoker, all instances set to eversmoker
    ungroup() %>%
    mutate(
      csmoker = case_when(
        (is.na(csmoker) | csmoker == 0) & eversmoker == 1 & (Smoker == 0 | SMOKCIG == "Y" | CIGPDAY >0) ~ 1,
        is.na(csmoker) & eversmoker == 0 ~ 0,
        is.na(csmoker) & SMOKCIG == "N" ~ 0,
        is.na(csmoker) & (eversmoker == 1 & Smoker == 2) ~ 0,
        csmoker == 1 & Smoker == 2 & SMOKCIG == "N" ~ 0,
        TRUE ~ csmoker
      ),
      eversmoker = case_when(
        is.na(eversmoker) & (csmoker == 1 | Smoker == 0 | Smoker == 2) ~ 1,
        is.na(eversmoker) & (csmoker == 0 | Smoker == 1) ~ 0,
        SMOKCIG == "Y" | CIGPDAY >0 | !is.na(smoknum) | !is.na(smokyr) ~ 1, #if they have any smoking history, not a 'never' smoker
        TRUE ~ eversmoker
      ),
      SMOKCIG = case_when(
        (is.na(SMOKCIG) & eversmoker == 0) ~ "N",
        (is.na(SMOKCIG) & csmoker == 0) ~ "N",
        is.na(SMOKCIG) & (Smoker == 0 | csmoker == 1)  ~ "Y",
        is.na(SMOKCIG) & (Smoker == 2 | Smoker == 1) ~ "N",
        TRUE ~ SMOKCIG
      )) %>%
    group_by(participantid) %>%
    mutate(SMOKCIG = na.locf(SMOKCIG, na.rm = FALSE)) %>% #fill in missing SMOKCIG values with participant's previous non-missing value
    mutate(eversmoker = cummax(eversmoker)) %>%
    ungroup()
}
df <- df %>% #adding/transforming variables for analysis
  group_by(participantid) %>% arrange(participantid, pftdate) %>%
  rename(fev1_1 = fev1, 
         fvc_1 = fvc) %>%
  mutate(fev1 = pmax(fev1_1, fev1_2, fev1_3, na.rm = TRUE), #select best FEV1 & FVC measurements
         fvc = pmax(fvc_1, fvc_2, fvc_3, na.rm = TRUE)) %>%
  mutate(fev1 = if_else(is.na(fev1), BJ_Best_fev1, fev1)) %>%
  select(-fev1_1, -fev1_2, -fev1_3, -fvc_1, -fvc_2, -fvc_3, -BJ_Best_fev1, -BJ_best_fvc) %>%
  mutate(ratio = fev1/fvc*100) %>% #create FEV1/FVC ratio variable (turned into a %)
  mutate(SmokingStatus = case_when( #create 5-level smoking status variable
    any(csmoker == 1) & !all(csmoker == 1) ~ "Intermittent",
    all(is.na(csmoker)) & all(is.na(eversmoker)) ~ "Unknown",
    all(csmoker == 0) & all(eversmoker == 0) ~ "Never",
    all(csmoker == 0) & any(eversmoker == 1) ~ "Former",
    all(csmoker == 1) ~ "Persistent",
    TRUE ~ NA_character_
  )) %>%
  mutate( #create baseline packyears of smoking variable (and convert smoknum and smokyr to numeric values from categories)
    smoknum = case_when(
      smoknum == 1 ~ 2.5, #<5 cig/day
      smoknum == 2 ~ 9.5, #5-14 cig/day
      smoknum == 3 ~ 19.5, #15-24 cig/day
      smoknum == 4 ~ 29.5, #25-34 cig/day
      smoknum == 5 ~ 39.5, #35-44 cig/day
      smoknum == 6 ~ 49.5, #45-54 cig/day
      smoknum == 7 ~ 60, #>55 cig/day
      is.na(smoknum) & !is.na(first(CIGPDAY)) ~ mean(CIGPDAY, na.rm = TRUE), #if smoknum is NA, use CIGPDAY as a proxy for smoknum
      TRUE ~ smoknum
    ),
    smokyr = case_when(
      smokyr == 1 ~ 2.5, #1-5y
      smokyr == 2 ~ 7, #6-8y
      smokyr == 3 ~ 14.5, #9-20y
      smokyr == 4 ~ 24.5, #21-28y
      smokyr == 5 ~ 34.5, #29-40y
      smokyr == 6 ~ 44.5, #41-48y
      smokyr == 7 ~ 60, #49-71+y
      TRUE ~ smokyr
    ),
    packyears = ifelse(eversmoker == 0, 0, smoknum/20 * smokyr)
  ) %>%
  mutate(packyears = na.locf(packyears, na.rm = FALSE)) %>% #sets packyears variable across all observations
  rename(height = cleanedheight_330) %>% #rename height variable
  mutate(age = as.numeric(difftime(pftdate, birthdt, units = "days"))/365.25) %>% #create age variable (at each pft measurement)
  mutate(education = case_when(
    educat %in% c(1, 2) ~ "< High School",
    educat %in% c(3, 4) ~ "High School/GED",
    educat %in% c(5, 6, 7, 8) ~ "Postsecondary",
    TRUE ~ NA_character_)) %>% #create 3-category education variable
  ungroup() %>% #ungroup to do row-wise variables
  mutate(ht2 = (height - mean(height, na.rm = T))^2, #centered, height squared variable (non-linearity assessment variable)
         age2 = (age - mean(age, na.rm = T))^2) %>% #centered, age squared variable (non-linearity assessment variable)
  select(-birthdt, -educat)
classifyCOPD <- function(val) {
  if (all(val == 0)) {
    return("NoCOPD")
  } else if (any(val == 1)) {
    firsti <- which(val == 1)[1] # First incident COPD
    subsequent <- val[firsti:length(val)] #check subsequent values of COPD
    if (length(subsequent) > 1 && subsequent[2] == 0) {      
      return("Transient")
    } else {
      rleSubsequent <- rle(subsequent)
      runs <- rleSubsequent$lengths[rleSubsequent$values == 1] #vector within 'rle'
      if (any(runs >= 2)) {
        return("Chronic")
      } else {
        return("Acute")
      }
    }
  }
}
copdp <- df %>% #COPD variables
  filter((!is.na(fev1) & !is.na(fvc))) %>%
  mutate(COPD = ifelse(ratio < 70, 1, 0)) %>% #create *LONG FORM* COPD variable
  group_by(participantid) %>% arrange(participantid, pftdate) %>%
  mutate(COPDTIME = as.numeric(pftdate - first(pftdate), units = "days") / 365.25) %>% #create COPD 'time' variable for years elapsed since first COPD visit (baseline)
  mutate(COPDprev = ifelse(first(COPD) == 1, 1, 0)) %>% #create COPD prevalence at baseline variable (binary: 0/1)
  mutate(COPDany = ifelse(any(COPD == 1), 1, 0))  # Create 'COPDany' variable to identify participants with any COPD diagnosis (this we will use for the COPD prevalence analysis)
copdi <- copdp %>% group_by(participantid) %>% filter(COPDprev == 0 & !is.na(COPD) & (FEV1_score %in% c(3, 4) & FVC_score %in% c(3, 4))) %>%
  mutate(COPDclass = classifyCOPD(COPD)) %>%
  select(participantid, time, COPDclass)
copdp <- copdp %>% select(participantid, time, COPD, COPDTIME, COPDprev, COPDany)
#merge datasets
df <- merge(df, copdp, by = c("participantid", "time"), all.x = TRUE)
df <- merge(df, copdi, by = c("participantid", "time"), all.x = TRUE)
merged <- merge(df, rasscores, by="participantid", all.x = TRUE)
final <- merged %>% select(-calories.x) %>% rename(calories = calories.y) %>% #remove duplicate calories variable and rename
  select(-csmoker, -eversmoker, -smoknum, -smokyr, -SMOKCIG, -Smoker) #remove not-needed variables
#REMOVE PARTICIPANTS
print(length(unique(final$participantid))) #full participant count; n = 2,921
final <- final %>% #remove participants with no diet data
  filter(!is.na(calories))
print(length(unique(final$participantid))) #n = 2,890
final <- final %>% #remove participants with no PFT data (missing both FEV1 and FVC)
  filter(!is.na(fev1) | !is.na(fvc))
print(length(unique(final$participantid))) #n = 2,884
final <- final %>% group_by(participantid) %>%
  mutate(count = row_number()) %>% #sequential variable PFT measurements
  mutate(n_pft = n()) #count number of PFTs per participant
RASfinal <- final
counts <- function(data) { #function to count participants
  kcal <- !is.na(data$calories) & data$calories != 0
  vars <- c("fev1", "fvc", "ratio", "FEV1_QC", "FVC_QC", "RATIO_QC", "FEV1_SMK", "FVC_SMK", "RATIO_SMK", "FEV1_PACK", "FVC_PACK", "RATIO_PACK")
  listy <- list(
    fev1 = !is.na(data$fev1),
    fvc = !is.na(data$fvc),
    ratio = !is.na(data$ratio),
    FEV1_QC = data$FEV1_score %in% c(3, 4),
    FVC_QC = data$FVC_score %in% c(3, 4),
    RATIO_QC = data$FEV1_score %in% c(3, 4) & data$FVC_score %in% c(3, 4),
    SMK = !is.na(data$SmokingStatus) & data$SmokingStatus !='Unknown',
    PACK = !is.na(data$packyears)
  )
  diet <- sum(!duplicated(data$participantid[kcal]))
  cat("Total participants:", diet, "\n")
  for (var in vars) {
    specific_conds <- switch(var,
                             fev1 = listy$fev1,
                             fvc = listy$fvc,
                             ratio = listy$ratio,
                             FEV1_QC = listy$fev1 & listy$FEV1_QC,
                             FVC_QC = listy$fvc & listy$FVC_QC,
                             RATIO_QC = listy$ratio & listy$RATIO_QC,
                             FEV1_SMK = listy$fev1 & listy$FEV1_QC & listy$SMK,
                             FVC_SMK = listy$fvc & listy$FVC_QC & listy$SMK,
                             RATIO_SMK = listy$ratio & listy$RATIO_QC & listy$SMK,
                             FEV1_PACK = listy$fev1 & listy$FEV1_QC & listy$SMK & listy$PACK,
                             FVC_PACK = listy$fvc & listy$FVC_QC & listy$SMK & listy$PACK,
                             RATIO_PACK = listy$ratio & listy$RATIO_QC & listy$SMK & listy$PACK)
    
    if (!is.null(specific_conds)) {
      count <- sum(!duplicated(data$participantid[specific_conds & kcal]))
      cat("Total participants with", var, "data:", count, "\n")
    } else {
      cat("Variable", var, "not found or not handled\n")
    }
  }
}
counts(RASfinal) #print final counts
#write final dataset to file
write.foreign(RASfinal, datafile = "RASfinal.dat",
              codefile="RASfinal.sas", package="SAS")
write_csv(RASfinal, "RASfinal.csv")
rm(list = ls())


