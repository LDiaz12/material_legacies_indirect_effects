# Respo Code for Material Legacy Indirect Effects P/R #####
# Created by Nyssa Silbiger
# Modified by Laurel Diaz
# Modified on 05/15/2024


#############################
#Install Packages 
#############################

if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented')
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix')
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra')
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') devtools::install_github('colin-olito/LoLinR')
if ("lubridate" %in% rownames(installed.packages()) == 'FALSE') install.packages('lubridate')
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron')
if ("tidyverse" %in% rownames(installed.packages()) == 'FALSE') install.packages('tidyverse')
if ("here" %in% rownames(installed.packages()) == 'FALSE') install.packages('here')
if ("patchwork" %in% rownames(installed.packages()) == 'FALSE') install.packages('patchwork')
if ("PNWColors" %in% rownames(installed.packages()) == 'FALSE') install.packages('PNWColors')

#rm(list=ls())

#############################
#Read in required libraries
#############################
##### Include Versions of libraries
library(segmented)
library(plotrix)
library(gridExtra)
library(LoLinR)
library(lubridate)
library(chron)
library(patchwork)
library(tidyverse)
library(here)
library(PNWColors)
library(ggrepel)


#############################
# get the file path
#############################

#set the path to all of the raw oxygen datasheets from presens software
here()
path.p<-here("Data","RespoFiles", "Initial") #the location of all your respirometry files

# bring in all of the individual files
file.names<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders

#basename above removes the subdirectory name from the file, re-name as file.names.full
file.names.full<-list.files(path = path.p, pattern = "csv$", recursive = TRUE) 

#generate a 4 column dataframe with specific column names
# data is in umol.L.sec
RespoR <- data.frame(matrix(NA, nrow=length(file.names.full), ncol=4)) # use instead of tidyverse tibble
colnames(RespoR) <- c("FileName", "Intercept", "umol.L.sec","Temp.C")

#Load your respiration data file, with all the times, water volumes(mL), surface area
# starting with initial respo files to start since there's less
Sample.Info <- read_csv(here("Data","RespoFiles","Initial", "InitialRespoMetaData.csv"))
surface_area <- read_csv(here("Data", "Data_Raw", "Growth", "SA", "MO24BEAST_SA_calculated.csv"))

surface_area <- surface_area %>%
  select(CORAL_NUM, GENOTYPE, SURFACE_AREA = SA_cm_2) %>%
  mutate(CORAL_NUM = as.character(CORAL_NUM))


Sample.Info <- Sample.Info %>%
  left_join(surface_area)

#Sample.Info$START <- as.POSIXct(Sample.Info$START,format="%H:%M:%S", tz = "") #convert time from character to time
#Sample.Info$END <- as.POSIXct(Sample.Info$END,format="%H:%M:%S", tz = "") #convert time from character to time

#empty chamber volume
ch.vol <- 450 #mL


#############################
# PROCESS O2 DATA
#############################


###forloop##### 
for(i in 1:length(file.names.full)) {
  FRow <- which(Sample.Info$FileName==strsplit(file.names[i], '.csv')) # stringsplit renames our file
  Respo.Data1 <- read_csv(file.path(path.p, file.names.full[i]), skip = 1) %>%
    select(Date, Time, Value, Temp) %>% # keep only what we need: Time stamp per 1sec, Raw O2 value per 1sec, in situ temp per 1sec
    unite(Date,Time,col="Time",remove=T, sep = " ") %>%
    drop_na()

  Respo.Data1 <-  Respo.Data1[-c(1:60),] %>% #starting 1 minute in
    mutate(sec = 1:n())  # create a new column for every second for the regression
  
  # Get the filename without the .csv
  rename<- sub(".csv","", file.names[i])
  
  
  ### plot and export the raw data ####
  p1<- ggplot(Respo.Data1, aes(x = sec, y = Value)) +
    geom_point(color = "dodgerblue") +
    labs(
      x = 'Time (seconds)',
      y = expression(paste(' O'[2],' (',mu,'mol/L)')),
      title = "original"
    )
  
  # thin the data by every 20 seconds to speed it up
  Respo.Data.orig<-Respo.Data1 # save original unthinned data
  Respo.Data1 <- thinData(Respo.Data1 ,by=5)$newData1 #thin data by every 5
  Respo.Data1$sec <- as.numeric(rownames(Respo.Data1)) #maintain numeric values for time
  Respo.Data1$Temp<-NA # add a new column to fill with the thinned data
  Temp_thinned <- thinData(Respo.Data.orig, xy = c(1, 3), by = 5)$newData1[, 2] #5
  Respo.Data1$Temp <-  Temp_thinned #thin data by every 20 points for the temp values
  
  p2 <- ggplot(Respo.Data1, aes(x = sec, y = Value))+
    geom_point(color = "dodgerblue")+
    labs(
      x = 'Time (seconds)',
      y = expression(paste(' O'[2],' (',mu,'mol/L)')),
      title = "thinned"
    )
  
  ##Olito et al. 2017: It is running a bootstrapping technique and calculating the rate based on density
  #option to add multiple outputs method= c("z", "e "pc")
  Regs  <-  rankLocReg(xall=Respo.Data1$sec, yall=Respo.Data1$Value, alpha=0.5, method="pc", verbose=TRUE)  
  
  # Print across two pages so use baseplot to create the pdf
  pdf(paste0(here("Output","RespoOutput"),"/", rename,"thinning.pdf"))
  
  p1+p2
  plot(Regs) # plot the results of Regs
  dev.off()
  
  # fill in all the O2 consumption and rate data
  RespoR[i,2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
  RespoR[i,1] <- rename #stores the file name in the Date column
  RespoR[i,4] <- mean(Respo.Data1$Temp, na.rm=T)  #stores the Temperature from the incubation in the Temp.C column
}  


#############################
# save file
#############################

#export raw data and read back in as a failsafe 
#this allows me to not have to run the for loop again 

write_csv(RespoR, here("Data","RespoFiles", "Initial", "Initial_Respo_R.csv"))  


#############################
# post-processing: normalize rates
#############################
RespoR <- read_csv(here("Data","RespoFiles","Initial", "Initial_Respo_R.csv"))

RespoR2 <- RespoR %>%
  right_join(Sample.Info, by = c("GENOTYPE", "CORAL_NUM", "LIGHT_DARK", "RUN_NUM", "DATE")) %>% # Join the raw respo calculations with the metadata
  mutate(umol.sec = umol.L.sec*ch.vol/1000) %>% #Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
  mutate_if(sapply(., is.character), as.factor) %>% #convert character columns to factors
  mutate(BLANK = as.factor(BLANK)) #make blank column a factor

#RespoR2$SURFACE_AREA <- as.numeric(RespoR2$SURFACE_AREA)

RespoR_Normalized <- RespoR2 %>%
  group_by(LIGHT_DARK, BLANK, RUN_NUM) %>% #  add run num here if one blank per run
  summarise(umol.sec = mean(umol.sec, na.rm=TRUE)) %>% # get mean value of blanks per run and remove NAs
  filter(BLANK == 1) %>% # only keep the actual blanks
  select(LIGHT_DARK, RUN_NUM, BLANK, blank.rate = umol.sec) %>%
  right_join(RespoR2, RespoR, by = c("LIGHT_DARK", "RUN_NUM")) %>% 
  mutate(umol.sec.corr = umol.sec - blank.rate, # subtract the blank rates from the raw rates
         umol.cm2.hr = (umol.sec.corr*3600)/SURFACE_AREA,
         umol.cm2.hr_uncorr = (umol.sec*3600)/SURFACE_AREA) %>% 
  select(DATE, CORAL_NUM, GENOTYPE, LIGHT_DARK, RUN_NUM, umol.cm2.hr, CHAMBER, 
                Temp.C, umol.cm2.hr_uncorr) %>%
  drop_na()

#View(RespoR_Normalized)


#############################
#Account for blank rate by light/Dark and Block
#############################

#View(RespoR)


# CALCULATING R AND GP

# make the respiration values positive (pull out data for dark treatments)
RespoR_Normalized_DARK <- RespoR_Normalized %>% 
  ungroup() %>%
  filter(LIGHT_DARK == "DARK") %>% 
  mutate(R = umol.cm2.hr*-1, ##mmol cm2 hr
         R_uncorr = umol.cm2.hr_uncorr*-1) %>% 
  mutate(R = ifelse(R < 0, 0, R), # for any values below 0, make 0
         R_uncorr = ifelse(R_uncorr < 0, 0, R_uncorr)) %>% 
  select(-c(umol.cm2.hr, umol.cm2.hr_uncorr, LIGHT_DARK)) # all dark run rates get R for respiration

# all light run rates get NP for net photosynthesis
RespoR_Normalized_LIGHT <- RespoR_Normalized %>% 
  ungroup() %>%
  filter(LIGHT_DARK == "LIGHT") %>% 
  #mutate(mmol.gram.hr = ifelse(mmol.gram.hr < 0, 0, mmol.gram.hr), # for any values below 0, make 0
         #mmol.gram.hr_uncorr = ifelse(mmol.gram.hr_uncorr < 0, 0, mmol.gram.hr_uncorr)) %>% 
  select(DATE, CORAL_NUM, GENOTYPE, NP = umol.cm2.hr, NP_uncorr = umol.cm2.hr_uncorr)

# rejoin data into single df
RespoR_Normalized2 <- full_join(RespoR_Normalized_DARK, RespoR_Normalized_LIGHT) %>%
  mutate(GP = NP + R, 
         GP_uncorr = NP_uncorr + R_uncorr)



#############################
# save file
#############################

write_csv(RespoR_Normalized2, here("Data","RespoFiles","Initial", "Respo_RNormalized_InitialRates.csv"))  




###############################################################################################
##### END OF CALCULATIONS
###############################################################################################
