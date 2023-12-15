library(tidyverse)
library(skimr)
library(maps)
library(plotly)

########## 1. Import data ##########
#---- Set working directory ----
getwd()
setwd("/Users/jackward/Documents/01_Projects/01_PhD/02_Indonesia_SGAM/02_Data/01_Geochemistry/01_R_Import")
getwd()

#---- Import data ----
## Manually compiled data 
Indonesia_SGAM_2023_Carn_and_Pyle_2001 <- read_csv("Indonesia_SGAM_2023_Carn_and_Pyle_2001_Import.csv")
Indonesia_SGAM_2023_Claproth_1989 <- read_csv("Indonesia_SGAM_2023_Claproth_1989_Import.csv")
Indonesia_SGAM_2023_Cooke_2017 <- read_csv("Indonesia_SGAM_2023_Cooke_2017_Import.csv")
Indonesia_SGAM_2023_Edwards_et_al_1991 <- read_csv("Indonesia_SGAM_2023_Edwards_et_al_1991_Import.csv")
Indonesia_SGAM_2023_Foden_1986 <- read_csv("Indonesia_SGAM_2023_Foden_1986_Import.csv")
Indonesia_SGAM_2023_Gertisser_et_al_2012 <- read_csv("Indonesia_SGAM_2023_Gertisser_et_al_2012_Import.csv")
Indonesia_SGAM_2023_Kirchenbaur_et_al_2022 <- read_csv("Indonesia_SGAM_2023_Kirchenbaur_et_al_2022_HFSEs_Import.csv")
Indonesia_SGAM_2023_Nicholls_and_Whitford_1983 <- read_csv("Indonesia_SGAM_2023_Nicholls_and_Whitford_1983_Import.csv")
Indonesia_SGAM_2023_Sendjaja_et_al_2009 <- read_csv("Indonesia_SGAM_2023_Sendjaja_et_al_2009_Import.csv")
Indonesia_SGAM_2023_Setijadji_et_al_2006 <- read_csv("Indonesia_SGAM_2023_Setijadji_et_al_2006_Import.csv")
Indonesia_SGAM_2023_Stolz_et_al_1988 <- read_csv("Indonesia_SGAM_2023_Stolz_et_al_1988_Import.csv")
Indonesia_SGAM_2023_Turner_and_Foden_2001 <- read_csv("Indonesia_SGAM_2023_Turner_and_Foden_2001_Import.csv")
Indonesia_SGAM_2023_Turner_et_al_2003 <- read_csv("Indonesia_SGAM_2023_Turner_et_al_2003_Import.csv")
Indonesia_SGAM_2023_Van_Bergen_et_al_1992 <- read_csv("Indonesia_SGAM_2023_Van_Bergen_et_al_1992_Import.csv")
Indonesia_SGAM_2023_Yu_et_al_2022 <- read_csv("Indonesia_SGAM_2023_Yu_et_al_2022_Import.csv")
Indonesia_SGAM_2023_Wheller_et_al_1987 <- read_csv("Indonesia_SGAM_2023_Wheller_et_al_1987_Import.csv")
Indonesia_SGAM_2023_Foden_1983 <- read_csv("Indonesia_SGAM_2023_Foden_1983_Import.csv")
Indonesia_SGAM_2023_Chesner_and_Rose_1993 <- read_csv("Indonesia_SGAM_2023_Chesner_and_Rose_Import.csv")
Indonesia_SGAM_2023_Gardner_et_al_2012 <- read_csv("Indonesia_SGAM_2023_Gardner_et_al_2013_Import.csv") 
Indonesia_SGAM_2023_Gertisser_et_al_2023 <- read_csv("Indonesia_SGAM_2023_Gertisser_et_al_2022_Import.csv")

## GEOROC Data
Indonesia_SGAM_2023_GEOROC_Sunda <- read_csv("Indonesia_SGAM_2023_GEOROC_Sunda.csv")
Indonesia_SGAM_2023_GEOROC_Banda <- read_csv("Indonesia_SGAM_2023_GEOROC_Banda.csv")

########## 2. Format manually compiled data ##########
#---- Check column names ----
colnames(Indonesia_SGAM_2023_Carn_and_Pyle_2001)
colnames(Indonesia_SGAM_2023_Claproth_1989)
colnames(Indonesia_SGAM_2023_Cooke_2017)
colnames(Indonesia_SGAM_2023_Edwards_et_al_1991)
colnames(Indonesia_SGAM_2023_Foden_1986)
colnames(Indonesia_SGAM_2023_Gertisser_et_al_2012)
colnames(Indonesia_SGAM_2023_Kirchenbaur_et_al_2022)
colnames(Indonesia_SGAM_2023_Nicholls_and_Whitford_1983)
colnames(Indonesia_SGAM_2023_Sendjaja_et_al_2009)
colnames(Indonesia_SGAM_2023_Setijadji_et_al_2006)
colnames(Indonesia_SGAM_2023_Stolz_et_al_1988)
colnames(Indonesia_SGAM_2023_Turner_and_Foden_2001)
colnames(Indonesia_SGAM_2023_Turner_et_al_2003)
colnames(Indonesia_SGAM_2023_Van_Bergen_et_al_1992)
colnames(Indonesia_SGAM_2023_Yu_et_al_2022)
colnames(Indonesia_SGAM_2023_Wheller_et_al_1987)
colnames(Indonesia_SGAM_2023_Foden_1983)
colnames(Indonesia_SGAM_2023_Chesner_and_Rose_1993)
colnames(Indonesia_SGAM_2023_Gardner_et_al_2012)
colnames(Indonesia_SGAM_2023_Gertisser_et_al_2023)

#---- Convert relevant sample names to character ----
## Foden (1983)
Indonesia_SGAM_2023_Foden_1983$SampleName <- as.character(Indonesia_SGAM_2023_Foden_1983$SampleName)

## Chesner and Rose (1993)
Indonesia_SGAM_2023_Chesner_and_Rose_1993$SampleName <- as.character(Indonesia_SGAM_2023_Chesner_and_Rose_1993$SampleName)

#---- Bind rows ----
Indonesia_SGAM_2023_Manual_Merged <- bind_rows(Indonesia_SGAM_2023_Carn_and_Pyle_2001, Indonesia_SGAM_2023_Claproth_1989, Indonesia_SGAM_2023_Cooke_2017, 
                                               Indonesia_SGAM_2023_Edwards_et_al_1991, Indonesia_SGAM_2023_Foden_1986, Indonesia_SGAM_2023_Gertisser_et_al_2012,
                                               Indonesia_SGAM_2023_Kirchenbaur_et_al_2022, Indonesia_SGAM_2023_Nicholls_and_Whitford_1983, Indonesia_SGAM_2023_Sendjaja_et_al_2009,
                                               Indonesia_SGAM_2023_Setijadji_et_al_2006, Indonesia_SGAM_2023_Stolz_et_al_1988, Indonesia_SGAM_2023_Turner_and_Foden_2001,
                                               Indonesia_SGAM_2023_Turner_et_al_2003, Indonesia_SGAM_2023_Van_Bergen_et_al_1992, Indonesia_SGAM_2023_Yu_et_al_2022,
                                               Indonesia_SGAM_2023_Wheller_et_al_1987, Indonesia_SGAM_2023_Foden_1983, Indonesia_SGAM_2023_Chesner_and_Rose_1993,
                                               Indonesia_SGAM_2023_Gardner_et_al_2012, Indonesia_SGAM_2023_Gertisser_et_al_2023)

#---- Check column names ----
colnames(Indonesia_SGAM_2023_Manual_Merged)

colnames(Indonesia_SGAM_2023_Manual_Merged %>% 
  select(sort(names(.))))

#---- Correct column classes ----
Indonesia_SGAM_2023_Manual_Merged[, c(1:5, 9:11, 75:76)] <- sapply(Indonesia_SGAM_2023_Manual_Merged[, c(1:5, 9:11, 75:76)], as.character)
Indonesia_SGAM_2023_Manual_Merged[, c(6:8, 12:74, 77:99)] <- sapply(Indonesia_SGAM_2023_Manual_Merged[, c(6:8, 12:74, 77:99)], as.numeric)

#---- Add data compilation source column ----
Indonesia_SGAM_2023_Manual_Merged$CompilationSource <- "Manual"

#---- Remove Pleistocene samples ----
Indonesia_SGAM_2023_Manual_Merged <- Indonesia_SGAM_2023_Manual_Merged %>% 
  filter(!SampleName %in% c("TT2 A", "TT12", "I25Pe1"))

Indonesia_SGAM_2023_Manual_Merged <- Indonesia_SGAM_2023_Manual_Merged %>% 
  filter(Volcano != "Cikuray")

#---- Skim data ----
skim(Indonesia_SGAM_2023_Manual_Merged)

#---- Plot data on map ----
Indonesia_SGAM_2023_Manual_Merged %>% 
  ggplot(aes(x = Long, y = Lat)) + 
  geom_point() +
  borders("world", colour = "black", fill = NA, size = 0.25) +
  xlim(90, 135) +
  ylim(-12, 15) +
  coord_fixed() +
  theme_bw()
  
########## 3. Format GEOROC data ##########
#---- Check column names ----
Sunda_Colnames <- colnames(Indonesia_SGAM_2023_GEOROC_Sunda)
Banda_Colnames <- colnames(Indonesia_SGAM_2023_GEOROC_Banda)

identical(Sunda_Colnames, Banda_Colnames)

#---- Remove unrelated papers ----
## Sunda 
Indonesia_SGAM_2023_GEOROC_Sunda <- Indonesia_SGAM_2023_GEOROC_Sunda %>% 
  filter(!CITATIONS %in% c("[2800]", "[3059]", "[3118]", "[3139]", "[3147]", "[3324]", "[3500]", "[3602]", "[3734]", "[3735]", "[3741]", "[3743]", "[3744]", "[4012]", "[4098]", "[4469]",
                           "[4651]", "[4652]", "[4674]", "[4776]", "[4777]", "[4929]", "[4938]", "[4947]", "[5580]", "[5868]", "[6639]", "[7940]", "[8721]", "[8794]",
                           "[9256]", "[9261]", "[10312]", "[10712]", "[11042]", "[11443]", "[12453]", "[14434]", "[14459]", "[14626]", "[14827]", "[14956]", "[14971]",
                           "[16040]", "[17006]", "[17163]", "[17558]", "[19058]", "[19571]", "[20041]", "[20049]", "[20507]", "[20823]", "[22262]", "[22596]", "[22961]",
                           "[23498]", "[23697]", "[23708]", "[24006]", "[24168]", "[24953]", "[25131]", "[25819]", "[26032]", "[26209]", "[26352]", "[26412]", "[26169]"))

## Banda
Indonesia_SGAM_2023_GEOROC_Banda <- Indonesia_SGAM_2023_GEOROC_Banda %>% 
  filter(!CITATIONS %in% c("[3070]", "[3083]", "[3119]", "[3128]", "[3131]", "[3139]", "[3147]", "[3733]", "[3910]", "[5178]", "[7176]", "[8599]", "[12598]", "[15265]",
                           "[15465]", "[19040]"))

#---- Merge GEOROC data frames ----
Indonesia_SGAM_2023_GEOROC_Merged <- bind_rows(Indonesia_SGAM_2023_GEOROC_Sunda, Indonesia_SGAM_2023_GEOROC_Banda)

#---- Remove manually compiled papers ----
Indonesia_SGAM_2023_GEOROC_Merged <- Indonesia_SGAM_2023_GEOROC_Merged %>% 
  filter(!CITATIONS %in% c("[3146]", "[3436]", "[3732]", "[3737]", "[3885]", "[3968]", "[4583]", "[4588]", "[4647]", "[6098]", "[15296]", "[15822]", "[25541]", "[25995]", "[16924]"))

#---- Correct column names ----
## Check column names
colnames(Indonesia_SGAM_2023_GEOROC_Merged)

## Correct column names
colnames(Indonesia_SGAM_2023_GEOROC_Merged) <- c("Citations", "TectonicSetting", "Location", "LocationComment", "LatitudeMin", "LatitudeMax", "LongitudeMin", 
                  "LongitudeMax", "LandorSea", "ElevationMin", "ElevationMax", "SampleName", "RockName", "MinAge",
                  "MaxAge", "Geol", "Age", "EruptionDay", "EruptionMonth", "EruptionYear", "RockTexture",
                  "RockType", "DrillDepthMin", "DrillDepthMax", "Alteration", "Mineral", "Material", "SiO2",
                  "TiO2", "B2O3", "Al2O3", "Cr2O3", "Fe2O3", "FeO", "FeOT", "CaO", "MgO", "MnO", "NiO", "K2O", "Na2O", "P2O5", "H2O", "H2OP", "H2OM",
                  "H2OT", "CO2_WT", "CO1", "F_WT", "Cl_WT", "Cl2", "OH", "CH4", "SO2", "SO3", "SO4","S_WT", "LOI", "Volatiles", "O", "Others", "He_CCM_G", "He_CCMSTP_G",
                  "He3_CCMSTP_G", "He3_AT_G", "He4_CCM_G", "He4_CCMSTP_G", "He4_AT_G", "He4_Mole_G", "He4_NCC_G", "He_NCC_G", "Li", "Be", "B", "C","CO2_PPM",
                  "F_PPM", "Na", "Mg", "Al", "P", "S_PPM", "Cl_PPM", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Rb", 
                  "Sr", "Y", "Zr", "Nb", "Mo", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Cs", "Ba", "La","Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb",
                  "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf","Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Th", "U", "143Nd144Nd", "Nd143_Nd144_INI",
                  "Epsilon_Nd", "87Sr86Sr", "Sr87_Sr86_INI", "206Pb204Pb", "Pb206_Pb204_INI", "207Pb204Pb", "Pb207_Pb204_INI", "208Pb204Pb", "Pb208_Pb204_INI",
                  "Os184_Os188", "Os186_Os188", "Os187_Os186", "Os187_Os188", "Re187_Os186", "Re187_Os188", "Hf176_Hf177", "He3_He4", "He3_He4_R_RA", "He4_He3", 
                  "He4_He3_R_RA", "K40_Ar40", "Ar40_K40", "UniqueID")

## Check column names
colnames(Indonesia_SGAM_2023_GEOROC_Merged)

#---- Select relevant columns ----
Indonesia_SGAM_2023_GEOROC_Merged <- Indonesia_SGAM_2023_GEOROC_Merged %>% 
  select(Citations, SampleName, LatitudeMin, LongitudeMin, RockType, Material, SiO2, TiO2, Al2O3, Cr2O3, Fe2O3, FeO, FeOT, CaO, MgO, MnO, NiO, K2O, Na2O, P2O5, LOI, Rb,
         Ba, Th, U, Nb, Ta, La, Ce, Pb, Pr, Sr, Nd, Sm, Zr, Hf, Eu, Tb, Dy, Er, Y, Yb, Lu, V, Sc, `143Nd144Nd`, `87Sr86Sr`, `206Pb204Pb`, `207Pb204Pb`, `208Pb204Pb`)

#---- Rename Citations, LatitudeMin, and LongitudeMin ----
## Rename Citations
Indonesia_SGAM_2023_GEOROC_Merged <- Indonesia_SGAM_2023_GEOROC_Merged %>% 
  rename_at("Citations", ~"Reference")

## Rename LatitudeMin
Indonesia_SGAM_2023_GEOROC_Merged <- Indonesia_SGAM_2023_GEOROC_Merged %>% 
  rename_at("LatitudeMin", ~"Lat")

## Rename LongitudeMin
Indonesia_SGAM_2023_GEOROC_Merged <- Indonesia_SGAM_2023_GEOROC_Merged %>% 
  rename_at("LongitudeMin", ~"Long")

## Check column names
colnames(Indonesia_SGAM_2023_GEOROC_Merged)

#---- Remove remaining Pleistocene samples ----
Indonesia_SGAM_2023_GEOROC_Merged <- Indonesia_SGAM_2023_GEOROC_Merged %>% 
  filter(!SampleName %in% c("s_LAVA 1 [25424]", "s_LAVA 2 [25424]", "s_67150 [3113]", "s_67154 [3113]"))

#---- Split "Material" column ----
Indonesia_SGAM_2023_GEOROC_Merged <- Indonesia_SGAM_2023_GEOROC_Merged %>% 
  separate(Material, c("Material", "MaterialRef"), sep = " ")

#---- Select whole rock volcanic analyses ----
Indonesia_SGAM_2023_GEOROC_Merged <- Indonesia_SGAM_2023_GEOROC_Merged %>% 
  filter(RockType == "VOL" & Material == "WR")

#---- Plot data on map ----
GEOROC_SampleMap <- Indonesia_SGAM_2023_GEOROC_Merged %>% 
  ggplot(aes(x = Long, y = Lat, colour = Reference)) + 
  geom_point() +
  borders("world", colour = "black", fill = NA, size = 0.25) +
  xlim(90, 135) +
  ylim(-12, 15) +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "none")

ggplotly(GEOROC_SampleMap)

#---- Remove rows with data from two or more references (references that should have been removed or include data from manually compiled papers or in wrong spot) ----
Indonesia_SGAM_2023_GEOROC_Merged <- Indonesia_SGAM_2023_GEOROC_Merged %>% 
  filter(!Reference %in% c("[2800][4938]", "[4647][4590]", "[4647][2608]", "[3122][3144][21977]", "[3122][21977][15235]", "[25541][15259]", "[3122][21977]",
                           "[3436][3113]", "[3885][3058]", "[3436][3146]", "[3885][3145]", "[3500][3734]", "[3737][3735]", "[4647][2713]", "[6098][7092]",
                           "[3885][3145][3085][3058]", "[3885][3436][3113]", "[3741][4938]", "[3059][3436][3146][3113]", "[3059][3436][3113]", "[3500][3734][22052]",
                           "[3500][3732]", "[14971][22961]", "[14956][22961]", "[14956][4652][22961]", "[14971][4776]", "[6098][3734][3500]", "[6098][22052][7092]",
                           "[3500][3740]", "[6098][3740][3500]", "[3741][13600][4938]", "[3735][22052]", "[12598][7176]"))

#---- Remove point in Philippine Sea ----
Indonesia_SGAM_2023_GEOROC_Merged <- Indonesia_SGAM_2023_GEOROC_Merged %>% 
  filter(!Lat == 4.525830)

#---- Add "Tambora" as volcano name to reference 25424 ----
Indonesia_SGAM_2023_GEOROC_Merged$Volcano <- ifelse(Indonesia_SGAM_2023_GEOROC_Merged$Reference == "[25424]", "Tambora", NA)

#---- Add data compilation source column ----
Indonesia_SGAM_2023_GEOROC_Merged$CompilationSource <- "GEOROC"

#---- Plot data on map ----
GEOROC_SampleMap <- Indonesia_SGAM_2023_GEOROC_Merged %>% 
  ggplot(aes(x = Long, y = Lat, colour = Reference)) + 
  geom_point() +
  borders("world", colour = "black", fill = NA, size = 0.25) +
  xlim(90, 135) +
  ylim(-12, 15) +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "none")

ggplotly(GEOROC_SampleMap)

#---- Check data classes ----
skim(Indonesia_SGAM_2023_GEOROC_Merged)

########## 4. Merge data sets ##########
#---- Check column names ----
## Manual compilation
colnames(Indonesia_SGAM_2023_Manual_Merged %>% 
           select(sort(names(.))))

## GEOROC compilation
colnames(Indonesia_SGAM_2023_GEOROC_Merged %>% 
           select(sort(names(.))))

#---- Merge data frames ----
Indonesia_SGAM_2023_Final <- bind_rows(Indonesia_SGAM_2023_Manual_Merged, Indonesia_SGAM_2023_GEOROC_Merged)

#---- Check data frame ----
skim(Indonesia_SGAM_2023_Final)

########## 5. Export data sets ##########
#---- Manual compilation formatted ----
write_csv(Indonesia_SGAM_2023_Manual_Merged, 
          "/Users/jackward/Documents/01_Projects/01_PhD/02_Indonesia_SGAM/02_Data/01_Geochemistry/02_R_Export/Indonesia_SGAM_EPSL_Final_Manual_Merged.csv")

#---- GEOROC formatted ----
write_csv(Indonesia_SGAM_2023_GEOROC_Merged, 
          "/Users/jackward/Documents/01_Projects/01_PhD/02_Indonesia_SGAM/02_Data/01_Geochemistry/02_R_Export/Indonesia_SGAM_EPSL_Final_GEOROC_Merged.csv")

#---- Combined formatted ----
write_csv(Indonesia_SGAM_2023_Final, 
          "/Users/jackward/Documents/01_Projects/01_PhD/02_Indonesia_SGAM/02_Data/01_Geochemistry/02_R_Export/Indonesia_SGAM_EPSL_Final.csv")

########## 6. Export EPSL supplementary data set ########## 
#---- Create supplementary data file ----
Indonesia_SGAM_Supplementary_Data <- Indonesia_SGAM_2023_Manual_Merged %>% 
  select(SampleName, Volcano, Lat, Long, Reference, SiO2, TiO2, Al2O3, Fe2O3, MnO, MgO, CaO, Na2O, K2O, P2O5, LOI, FeO, Rb, Ba, Th, U,
         Nb, Ta, La, Ce, Pb, Pr, Sr, Nd, Sm, Zr, Hf, Eu, Tb, Dy, Er, Y, Yb, Lu)

#---- Export supplementary data file ----
write_csv(Indonesia_SGAM_Supplementary_Data, "/Users/jackward/Documents/01_Projects/01_PhD/02_Indonesia_SGAM/10_Revisions/04_Supplementary_Data/Indonesia_SGAM_EPSL_Supplementary_Geochem_Final.csv")























