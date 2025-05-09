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
path.p<-here("Data","RespoFiles") #the location of all your respirometry files

# bring in all of the individual files
file.names<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders

#basename above removes the subdirectory name from the file, re-name as file.names.full
file.names.full<-list.files(path = path.p, pattern = "csv$", recursive = TRUE) 

#empty chamber volume
ch.vol <- 450 #mL 

#Load your respiration data file, with all the times, water volumes(mL), surface area
# starting with initial respo files to start since there's less
InitialRespoMeta <- read_csv(here("Data","RespoFiles","Initial", "InitialRespoMetaData.csv"))

Metadata <- read_csv(here("Data","MO24BEAST_Metadata.csv"))

# join the data together
Metadata$CORAL_NUN <- as.numeric(Metadata$CORAL_NUM)
InitialRespoMeta$CORAL_NUM <- as.numeric(InitialRespoMeta$CORAL_NUM)
Sample.Info <- left_join(InitialRespoMeta, Metadata) %>%
  select(DATE, CORAL_NUM, GENOTYPE, LIGHT_DARK, RUN_NUM, CHAMBER, VOLUME, BLANK, START, END,
         TIME_IN, TIME_OUT, SURFACE_AREA)

##### Make sure times are consistent ####

# make start and stop times real times, so that we can join the data frames
Sample.Info <- Sample.Info %>%
  unite(DATE,START,col="START_DATE_TIME",remove=F, sep=" ") %>% 
  unite(DATE,END,col="END_DATE_TIME",remove=F, sep=" ") %>% 
  mutate(START_TIME = mdy_hms(START_DATE_TIME)) %>% 
  mutate(END_TIME = mdy_hms(END_DATE_TIME)) %>% 
  mutate(DATE = mdy(DATE))

Sample.Info <- Sample.Info %>%
  mutate(
    START_DATE_TIME = mdy_hms(START_DATE_TIME),
    END_DATE_TIME = mdy_hms(END_DATE_TIME)
  )


#############################
# PROCESS O2 DATA
#############################

#generate a 4 column dataframe with specific column names
# data is in umol.L.sec
RespoR <- data.frame(matrix(NA, nrow=length(file.names.full), ncol=4)) # use instead of tidyverse tibble
colnames(RespoR) <- c("FileID","Intercept", "umol.L.sec","Temp.C")

###forloop##### 
for(i in 1:length(file.names.full)) {
  FRow <- as.numeric(which(Sample.Info$FileID==file.names.full[i])) # stringsplit this renames our file
  Respo.Data1 <- read_csv(file.path(path.p, paste0(file.names.full[i])), # reads in each file in list
    skip = 1) %>%
    dplyr::select(Date, Time, Value, Temp) %>% # keep only what we need: Time stamp per 1sec, Raw O2 value per 1sec, in situ temp per 1sec
    unite(Date,Time,col="Time",remove=T, sep = " ") %>% 
    mutate(Time = mdy_hms(Time)) %>% 
    drop_na()
  
  Respo.Data1 <- Respo.Data1 %>%
    filter(between(Time, Sample.Info$START_DATE_TIME[FRow], Sample.Info$END_DATE_TIME[FRow])) # select only data between start and stop time
  
  
  Respo.Data1 <-  Respo.Data1[-c(1:60),] %>% 
    mutate(sec = 1:n())  # create a new column for every second for the regression
  
  # Get the filename without the .csv
  rename<- sub(".csv","", filenames_final[i])
  
  
  ### plot and export the raw data ####
  p1<- ggplot(Respo.Data1, aes(x = sec, y = Value)) +
    geom_point(color = "dodgerblue") +
    labs(
      x = 'Time (seconds)',
      y = expression(paste(' O'[2],' (',mu,'mol/L)')),
      title = "original"
    )
  
  # thin the data by every 20 seconds to speed it up
  Respo.Data.orig<-Respo.Data1 # save original unthinned data #there is no thin() anymore, created alternative 
  newRespo<-tibble(
    Time=as.numeric(),
    Value=as.numeric(),
    Temp=as.numeric(),
    sec=as.numeric()
  )
  for(j in 1:nrow(Respo.Data.orig)) { # alternative thinning strategy
    if(j%%20==0){
      newRespo<-rbind(newRespo,Respo.Data1[j,])
    }
  }
  Respo.Data1<-newRespo # assign thinned data to previous df
  
  # Respo.Data1 <- Thin(Respo.Data1 ,By=20)$newData1 #thin data by every 20 points for all the O2 values
  # Respo.Data1$sec <- as.numeric(rownames(Respo.Data1 )) #maintain numeric values for time
  # Respo.Data1$Temp<-NA # add a new column to fill with the thinned data
  # Respo.Data1$Temp <-  Thin(Respo.Data.orig,xy = c(1,3),by=20)#$newData1[,2] #thin data by every 20 points for the temp values
  
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
  pdf(paste0(here("Output","RespoOutput","ThinningPlots"),"/", rename,"thinning.pdf"))
  
  plot(Regs) # plot the results of Regs
  plot(p2+p1) # use patchwork to bring the raw and thinned data together
  dev.off()
  
  # fill in all the O2 consumption and rate data
  # need clarity on what this is
  RespoR[i,2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
  RespoR[i,1] <- paste0(rename,".csv") #stores the file name in the Date column
  RespoR[i,4] <- mean(Respo.Data1$Temp, na.rm=T)  #stores the Temperature from the incubation in the Temp.C column
}  



#############################
# save file
#############################

#export raw data and read back in as a failsafe 
#this allows me to not have to run the for loop again 

write_csv(RespoR, here("Data","RespoFiles","Respo_R.csv"))  



#############################
# post-processing: noramlize rates
#############################
RespoR <- read_csv(here("Data","RespoFiles","Respo_R.csv"))

# adjust temporarily labeled sample ID's
Sample.Info <- Sample.Info %>% # needed to rename for above file processing due to duplicate runs with same names
  select() %>% ##NEED SURFACE AREA, VOLUME
  full_join(BioMeta) # need to rejoin biometadata because of issues with run 6 dark run joining volume and afdw

# Calculate Respiration rate
RespoR2 <- RespoR %>%
  drop_na(FileID) %>% # drop NAs
  left_join(Sample.Info) %>% # Join the raw respo calculations with the metadata
  left_join(BioMeta) %>% 
  arrange(FileID) %>% 
  mutate(Ch.Volume.L = Ch.Volume.ml * 0.001) %>% # mL to L conversion
  mutate(umol.sec = umol.L.sec*Ch.Volume.L) %>% #Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
  mutate_if(sapply(., is.character), as.factor) %>% #convert character columns to factors
  mutate(BLANK = as.factor(BLANK)) #make the blank column a factor

# Remove duplicates when assemblages were run multiple times
# anti <- RespoR2 %>% 
#   filter(Date == "2022-07-14") %>% 
#   distinct(SampleID) %>% 
#   left_join(RespoR2) %>% 
#   filter(Date == "2022-07-13")
# 
# RespoR2 <- RespoR2 %>% 
#   anti_join(anti)


#View(RespoR2)



#############################
#Account for blank rate by light/Dark and Block
#############################

#View(RespoR)

RespoR_Normalized <- RespoR2 %>%
  arrange(FileID) %>% 
  group_by(run_block, BLANK,light_dark) %>% # also add block here if one blank per block
  # group_by(BLANK)%>% # also add block here if one blank per block
  summarise(umol.sec = mean(umol.sec, na.rm=TRUE)) %>% # get mean value of blanks per run
  filter(BLANK == 1) %>% # only keep the actual blanks
  group_by(run_block, light_dark) %>% 
  dplyr::select(blank.rate = umol.sec) %>% # rename the blank rate column
  right_join(RespoR2) %>% # join with the respo data
  arrange(FileID) %>% 
  mutate(umol.sec.corr = umol.sec - blank.rate, # subtract the blank rates from the raw rates
         mmol.gram.hr = 0.001*(umol.sec.corr*3600)/AFDW.g, # convert to mmol g-1 hr-1
         mmol.gram.hr_uncorr = 0.001*(umol.sec*3600)/AFDW.g) %>% 
  filter(is.na(BLANK)) %>% # remove the Blank data
  ungroup() %>% 
  dplyr::select(Date, SampleID, AT, ET, light_dark, run_block, AFDW.g, Ch.Volume.ml, mmol.gram.hr, chamber_channel, 
                Temp.C, mmol.gram.hr_uncorr)  #keep only what we need

# CALCULATING R AND GP

# make the respiration values positive (pull out data for dark treatments)
RespoR_Normalized_dark <- RespoR_Normalized %>% 
  filter(light_dark == "DARK") %>% 
  mutate(mmol.gram.hr = mmol.gram.hr*-1, ##mmol cm2 hr
         mmol.gram.hr_uncorr = mmol.gram.hr_uncorr*-1) %>% 
  mutate(mmol.gram.hr = ifelse(mmol.gram.hr < 0, 0, mmol.gram.hr), # for any values below 0, make 0
         mmol.gram.hr_uncorr = ifelse(mmol.gram.hr_uncorr < 0, 0, mmol.gram.hr_uncorr)) %>% 
  mutate(P_R = "R") # all dark run rates get R for respiration

# all light run rates get NP for net photosynthesis
RespoR_Normalized_light <- RespoR_Normalized %>% 
  filter(light_dark == "LIGHT") %>% 
  mutate(mmol.gram.hr = ifelse(mmol.gram.hr < 0, 0, mmol.gram.hr), # for any values below 0, make 0
         mmol.gram.hr_uncorr = ifelse(mmol.gram.hr_uncorr < 0, 0, mmol.gram.hr_uncorr)) %>% 
  mutate(P_R = "NP")

# rejoin data into single df
RespoR_Normalized2 <- full_join(RespoR_Normalized_light, RespoR_Normalized_dark) %>% 
  drop_na(mmol.gram.hr) # removes anticipated sampleID's that were not actually run


#make column for GP and group by fragment ID and temp to keep R and NP together
RespoR_NormalizedGP <- RespoR_Normalized2 %>% 
  group_by(SampleID, AT, ET) %>% 
  summarize(mmol.gram.hr = sum(mmol.gram.hr),
            mmol.gram.hr_uncorr = sum(mmol.gram.hr_uncorr), # NP + R = GP
            Temp.C = mean(Temp.C)) %>% # get mean temperature of light and dark runs
  mutate(P_R="GP") %>% # Label for Gross Photosynthesis
  mutate(light_dark = "LIGHT") %>% 
  mutate(mmol.gram.hr = ifelse(mmol.gram.hr < 0, 0, mmol.gram.hr), # for any values below 0, make 0
         mmol.gram.hr_uncorr = ifelse(mmol.gram.hr_uncorr < 0, 0, mmol.gram.hr_uncorr))

# rejoin for full df with NP, R, and GP rates
RespoR_Normalized_Full <- RespoR_Normalized2 %>% 
  dplyr::select(SampleID, AT, ET, light_dark, P_R, mmol.gram.hr, mmol.gram.hr_uncorr, Temp.C) %>% 
  full_join(RespoR_NormalizedGP)


#############################
# save file
#############################

#write_csv(RespoR_Normalized_Full , here("Data","RespoFiles","Respo_RNormalized_AllRates.csv"))  




###############################################################################################
##### END OF CALCULATIONS
###############################################################################################



### PLOTS AND STATS IN NEC_Assemblage.R