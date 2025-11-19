### PI curve #####################################################
# Created by Jamie Kerlin
##################################################################

### Load libraries ################################################
library(tidyverse)
library(here)
library(patchwork)
library(segmented)
library(plotrix)
library(gridExtra)
library(LoLinR)
library(lubridate)
library(chron)

### Load data ###################################################
#set the path to raw oxygen datasheets
path.p<-here("Data",
             "RespoFiles",
             "PRUS_PI_CURVE_O2.csv") #location of PI curve

#file.names<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders


#basename above removes the subdirectory name from the file
#add file names that include the subdirectory name (note, these are the same for this example, but I often have lots of subfolders for different Runs)
#file.names.full<-list.files(path = path.p, pattern = "csv$", recursive = TRUE) #list all csv file names in the folder and subfolders

#generate an empty 3 column dataframe with specific column names
Photo.R <- data.frame(matrix(NA, nrow=length(file.names), ncol=3))
colnames(Photo.R) <- c("Intercept", "umol.L.sec","Temp.C") # name the columns


#Load Sample Info
Sample.Info <- read_csv(here("Data", "RespoFiles", "PRUS_PI_CURVE_O2.csv")) %>%
  mutate(new = "PI_Ch") %>%
  unite("fragment.ID.full", c(new, chamber.channel), sep = "", remove = FALSE) %>%
  unite("fragment.ID.full", c(fragment.ID.full, run, fragment.ID), sep = "_", remove = FALSE)

sa_volume <- read_csv(here("Short_term", "Data_Raw", "PR", "PI_curve", 
                           "PI_curve_2021_08_07", "Metadata", "sample_sa_volume.csv"))


# make start and stop times real times
Sample.Info$start.time <- as.POSIXct(Sample.Info$start.time,format="%H:%M:%S", tz = "") #convert time from character to time
Sample.Info$stop.time <- as.POSIXct(Sample.Info$stop.time,format="%H:%M:%S", tz = "") #convert time from character to time


PR<-c('Photo','Resp')

# for every file in list calculate O2 uptake or release rate and add the data to the Photo.R dataframe
for(i in 1:length(file.names.full)) { # for every file in list calculate O2 uptake or release rate and add the data to the Photo.R dataframe
  
  #find the lines in sample info that have the same file name that is being brought it
  FRow<-which(Sample.Info$fragment.ID==strsplit(file.names[i],'.csv'))
  
  # read in the O2 data one by one
  Photo.Data1 <-read.csv(file.path(path.p,file.names.full[i]), skip = 1, header=T) # skips the first line
  Photo.Data1  <- Photo.Data1[,c("Time","Value","Temp")] #subset columns of interest
  Photo.Data1$Time <- as.POSIXct(Photo.Data1$Time,format="%H:%M:%S", tz = "") #convert time from character to time
  Photo.Data1 <- na.omit(Photo.Data1)
  
  
  # clean up some of the data
  n<-dim(Photo.Data1)[1] # length of full data
  Photo.Data1 <-Photo.Data1[(n-120):(n-3),] #start at data point ~2 minute in to avoid excess noise from start of run and remove last 3 lines containing text
  n<-dim(Photo.Data1)[1] #list length of trimmed data
  Photo.Data1$sec <- (1:n) #set seconds by one from start to finish of run in a new column
  
  
  #Save plot prior to and after data thinning to make sure thinning is not too extreme
  rename <- sub(".csv","", file.names[i]) # remove all the extra stuff in the file name
  
  pdf(paste0("Short_term/Data_Raw/PR/PI_curve/PI_curve_2021_08_07/Output/",rename,"thinning.pdf")) # open the graphics device
  
  par(omi=rep(0.3, 4)) #set size of the outer margins in inches
  par(mfrow=c(1,2)) #set number of rows and columns in multi plot graphic
  plot(Value ~ sec, data=Photo.Data1 , xlab='Time (seconds)', ylab=expression(paste(' O'[2],' (',mu,'mol/L)')), type='n', axes=FALSE) #plot (empty plot to fill) data as a function of time
  usr  <-  par('usr') # extract the size of the figure margins
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA) # put a grey background on the plot
  whiteGrid() # make a grid
  box() # add a box around the plot
  points(Photo.Data1 $Value ~ Photo.Data1 $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
  axis(1) # add the x axis
  axis(2, las=1) # add the y-axis
  
  # This the data to make the code run faster
  Photo.Data.orig<-Photo.Data1#save original unthinned data
  Photo.Data1 <-  thinData(Photo.Data1 ,by=20)$newData1 #thin data by every 20 points for all the O2 values
  Photo.Data1$sec <- as.numeric(rownames(Photo.Data1 )) #maintain numeric values for time
  Photo.Data1$Temp<-NA # add a new column to fill with the thinned data
  Photo.Data1$Temp <-  thinData(Photo.Data.orig,xy = c(1,3),by=20)$newData1[,2] #thin data by every 20 points for the temp values
  
  # plot the thinned data
  plot(Value ~ sec, data=Photo.Data1 , xlab='Time (seconds)', ylab=expression(paste(' O'[2],' (',mu,'mol/L)')), type='n', axes=FALSE) #plot thinned data
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  box()
  points(Photo.Data1 $Value ~ Photo.Data1 $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
  axis(1)
  axis(2, las=1)
  ##Olito et al. 2017: It is running a bootstrapping technique and calculating the rate based on density
  #option to add multiple outputs method= c("z", "eq", "pc")
  Regs  <-  rankLocReg(xall=Photo.Data1$sec, yall=Photo.Data1$Value, alpha=0.5, method="pc", verbose=TRUE)  
  
  # add the regression data
  plot(Regs)
  dev.off()
  
  # fill in all the O2 consumption and rate data
  Photo.R[i,2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
  Photo.R[i,1] <- rename #stores the file name in the Date column
  Photo.R[i,4] <- mean(Photo.Data1$Temp, na.rm=T)  #stores the Temperature in the Temp.C column
  #Photo.R[i,5] <- PR[j] #stores whether it is photosynthesis or respiration
  
  
  # rewrite the file everytime... I know this is slow, but it will save the data that is already run
}


write.csv(Photo.R, 'Short_term/Data_Raw/PR/PI_curve/PI_curve_2021_08_07/Output/Photo.R.csv')  
View(Photo.R)

# Calculate P and R rate

fragIDs <- Sample.Info %>%
  dplyr::select(fragment.ID.full, fragment.ID, run) 

Photo.R1 <- left_join(Photo.R, fragIDs)
Photo.R2 <- left_join(Photo.R1, sa_volume)
View(Photo.R2)


#Convert sample volume to L
Photo.R2$volume <- Photo.R2$volume/1000 #calculate volume

#Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
Photo.R2$umol.sec <- Photo.R2$umol.L.sec*Photo.R2$volume

#Account for blank rate by temperature
#convert character columns to factors
Photo.R3 <- Photo.R2 %>%
  mutate_if(sapply(., is.character), as.factor)
View(Photo.R3)

#make the blank column a factor
Photo.R4 <- Photo.R3 %>%
  mutate(BLANK = if_else(fragment.ID == "BLANK1" | fragment.ID == "BLANK2", 1, 0))

photo.blnk <- aggregate(umol.sec ~ run*BLANK, data=Photo.R4, mean)
# pull out only the blanks

photo.blnk2 <- photo.blnk %>%
  filter(BLANK == 1) %>%
  dplyr::select(!BLANK) %>%
  rename(blank.rate = umol.sec) 

# join the blank data with the rest of the data
Photo.R5 <- left_join(Photo.R4, photo.blnk2)

# subtract the blanks######################
Photo.R6 <- Photo.R5 %>%
  mutate(umol.sec.corr = umol.sec - blank.rate)

View(Photo.R6)

#### Normalize to surface area #####

Photo.R7 <- Photo.R6 %>%
  mutate(umol.cm2.hr = (umol.sec.corr*3600)/surf.area.cm2) %>% #convert to hour and normalize to surface area
  filter(BLANK == 0) %>% # remove blanks
  dplyr::select(!Light_Dark)  #dont need this column
  # mutate(rate.ln = log(umol.cm2.hr + 0.1)) #log the rates #produces NAs because of negative rates

Photo.R7 %>%
  write_csv(here("Short_term", "Data_Raw", "PR", "PI_curve", 
                 "PI_curve_2021_08_07", "Output", "PI_rates.csv"))

PhotoMeans<- Photo.R7 %>%
  group_by(run)%>%
  summarise(rates.mean = mean(umol.cm2.hr), se = sd(umol.cm2.hr)/sqrt(n()))


# plot the raw data with the means on top
ggplot()+
  theme_bw()+  
  #geom_point(data=Photo.R, aes(x=run, y=umol.cm2.hr, alpha = 0.05), position = position_dodge(width = 0.2), size=4)+
  geom_point(data=PhotoMeans, aes(x=run, y=rates.mean),  size=1)+
  geom_line(data = PhotoMeans,  aes(x=run, y=rates.mean), size=1)+
  geom_errorbar(data = PhotoMeans, aes(x = run, ymin=rates.mean-se, ymax=rates.mean+se, width=.2))
#facet_wrap(~ Species, labeller = labeller(.multi_line = FALSE))+
ggsave('PI_curve_sites1/Output/RespirationRates.png')

#Mo'orea PI curve fit
#pulling out numeric for everything, pull put for high and pull out for low and do a curve 
PAR <- as.numeric(Photo.R7$run)

Pc <- as.numeric(Photo.R7$umol.cm2.hr)

plot(PAR,Pc,xlab="",
     ylab="",
     xlim=c(0,max(PAR)),
     ylim=c(-1,1.2),
     cex.lab=0.8,cex.axis=0.8,cex=1, main="",
     adj=0.05)

#set plot info

mtext(expression("Photon Flux Density ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1)

#add labels

mtext(expression(Photosynthetic~Rate* " ("*mu*"mol "*O[2]*" "*cm^-2*h^-1*")"),side=2,line=2,cex=1)

#add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)

curve.nlslrc= nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,start=list(Am=(max(Pc)-min(Pc)),AQY=0.001,Rd=-min(Pc),theta=0.6)) 

my.fit <- summary(curve.nlslrc)
#summary of model fit


#draw the curve using the model fit
#hyperbolic tangent, more common way to fit 

mor.curve.fitting <- curve((1/(2*summary(curve.nlslrc)$coef[4,1]))*(summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1]-sqrt((summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1])^2-4*summary(curve.nlslrc)$coef[2,1]*summary(curve.nlslrc)$coef[4,1]*summary(curve.nlslrc)$coef[1,1]*x))-summary(curve.nlslrc)$coef[3,1],lwd=2,col="blue",add=T)




#Amax (max gross photosytnthetic rate)

Pmax.gross <- my.fit$parameters[1]


#AQY (apparent quantum yield) alpha

AQY <- my.fit$parameters[2]


#Rd (dark respiration)

Rd <- my.fit$parameters[3]


# Ik light saturation point

Ik <- Pmax.gross/AQY

# add a line to figure for the Ik or saturation point 

abline(v=Ik, col="red", lty=3, lwd = 3)
text(x = 300, y = 1.0, label = "Ik", srt = 0)

dev.copy(png,'Short_term/Data_Raw/PR/PI_curve/PI_curve_2021_08_07/Output/PIcurve.png', width = 5, height = 4, units = "in", res = 1200)
dev.off()


# Ic light compensation point

Ic <- Rd/AQY


# Net photosynthetic rates

Pmax.net <- Pmax.gross-Rd


#output parameters into a table

Mor.PI.Output <- as.data.frame(cbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic)) %>%
  write_csv(here("Short_term", "Data_Raw", "PR", "PI_curve",
                 "PI_curve_2021_08_07", "Output", "PIcurve_values.csv"))



pdf("Short_term/Data_Raw/PR/PI_curve/PI_curve_2021_08_07/Output/PIcurve.pdf")

getwd()

dev.off()

#randomize my tanks and samples
x3 <- sample (1:10, 10)
x3

###############################################################################################################
#load three PI curves and light spectrum to then make into one plot
library(png)
library(grid)
library(gridExtra)

plot1 <- readPNG('PI_curve_sites1/Output/Moorea_Sum19_PI_Sites1.png')
plot2 <- readPNG('PI_curve_sites2/Output/Moorea_Sum19_PI_Sites2.png')
plot3 <- readPNG('PI_curve_sites3/Output/Moorea_Sum19_PI_Sites3.png')
plot4 <- readPNG('../../Repositories/Light_Spectrum/Output/blueandmultilight.png')

arrange <- grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3), rasterGrob(plot4), ncol=2, nrow=2)


ggsave("PI_curve_sites1/Output/PI_curve_all.pdf", device = "pdf", arrange, width = 8, height = 6)
