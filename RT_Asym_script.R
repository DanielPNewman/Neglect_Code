# setwd(("S:\R-MNHS-SPP\Bellgrove-data\11.Megan_ONeill\Analysis Scripts_R"))

setwd(("S:/R-MNHS-SPP/Bellgrove-data/11.Megan_ONeill/Analysis Scripts_R"))

### download relevant libraries
# install.packages("car")
# install.packages("nlme")
# install.packages("reshape")
# install.packages("pastecs")
# install.packages("psych")
# install.packages("multcomp")
# install.packages("ggplot2")
# install.packages("compute.es")
# install.packages("ez")
# install.packages("tidyr")

## Install relevant libraries 
library(foreign)
library(car)
library(ggplot2)
library(reshape)
library(pastecs)
library(psych)
library(plyr)
library(multcomp)
library(reshape2)
library(compute.es)
library(ez)
library(lme4)
library(effects)
# library(tidyr)

#Import data:
# data <- read.csv("D:/11.Megan_ONeill/Analysis Scripts_Matlab/master_matrix_R.csv", header=FALSE)
data <- read.csv("S:/R-MNHS-SPP/Bellgrove-data/11.Megan_ONeill/Analysis Scripts_Matlab/master_matrix_R_Healthy_Younger.csv", header=FALSE)
#Import IDs:
# ID <- read.table("D:/11.Megan_ONeill/Analysis Scripts_Matlab/ID_vector.csv", quote="\"")
ID <- read.table("S:/R-MNHS-SPP/Bellgrove-data/11.Megan_ONeill/Analysis Scripts_Matlab/ID_vector_Healthy_Younger.csv", quote="\"")
data$IDnum<-data[,1]
#Replace the participant numbers with IDs:
data[,1]<-ID[,1]
#Remove the seperate ID vector now it has been included into data dataframe
rm(ID)

#Rename data columns:

data<-rename(data, c("V1"="ID", "V2"="ResponseMapping","V3"="TotalTrialNumber","V4"="Trial","V5"="ITI",
                     "V6"="Hemifield","V7"="MotionDirection","V8"="RespLR",
                     "V9"="Accuracy","V10"="FixationBreak","V11"="Artefact500msPre",
                     "V12"="ArtefactPost","V13"="RejectedTrial","V14"="Threshold_cohLevel",
                     "V15"="RT","V16"="PreAlphaPower","V17"="PreAlphaPower_LeftHemi",
                     "V18"="PreAlphaPower_RightHemi","V19"="PreAlphaAsym","V20"="PrePupilDiameter",
                     "V21"="Sex","V22"="Age","V23"="TimeOfDay",
                     "V24"="NumberRepeatedTrials"))


#Make the required columns into factors:
data$ITI <- factor(data$ITI)
data$Hemifield <- factor(data$Hemifield)
data$MotionDirection <- factor(data$MotionDirection)
data$Sex <- factor(data$Sex)
data$Accuracy <- factor(data$Accuracy)

#Rename factor Levels:
data$ITI <- revalue(data$ITI, c("1"="1800ms", "2"="2800ms", "3"="3800ms"))
data$Hemifield <- revalue(data$Hemifield, c("1"="Left", "2"="Right"))
data$MotionDirection <- revalue(data$MotionDirection, c("1"="Uppward", "2"="Downward"))
data$Sex <- revalue(data$Sex, c("1"="Male", "2"="Female"))
data$Accuracy <- revalue(data$Accuracy, c("1"="Hit", "0"="Miss", "2"="WrongButton"))

#Re-class required vectors into Logicals:
data$FixationBreak<-as.logical(data$FixationBreak)
data$Artefact500msPre<-as.logical(data$Artefact500msPre)
data$ArtefactPost<-as.logical(data$ArtefactPost)
data$RejectedTrial<-as.logical(data$RejectedTrial)

###############Data Cleaning on Single Trial Level######################

#Remove trials with fixation-breaks and rejected trials
data<-data[!data$FixationBreak,]
data<-data[!data$RejectedTrial,]

#Check number of Trials for each participant by running the function 'length', 
#on "data$RT" for each group, broken down by ID 
num_trials1 <- ddply(data, c("ID"), summarise,
                     Trials    = length(RT))


#Create a variable numbering all the valid trials for each ID
data$ValidTrialNum<-(matrix(unlist(tapply(data[,1], data$ID, seq_along)), nrow=length(data[,1]), byrow=T))
data$ValidTrialNum<-c(data$ValidTrialNum) #change if from a one dimensional materix to a vector


#Remove trials where RT=0 (i.e. they did not respond)
data<-data[data$RT!=0,]
#Remove trials with missing values :
data<-data[complete.cases(data),] 

data$RT_log<-log(data$RT)

#####Z-score each participant's RT_log data 
#calculate mean and sd 
m <- tapply(data$RT_log,data$ID,mean)
s <- sqrt(tapply(data$RT_log,data$ID,var))
#calculate RT_log.Z and save it inside data.frame
data$RT_log.Z <- (data$RT_log-m[data$ID])/s[data$ID]
#check that Z scores have mean=0 and std=1 
RT_log.Z_checker <- ddply(data, c("ID"), summarise,
                      N    = length(RT_log.Z ),
                      mean = round(mean(RT_log.Z )),
                      sd   = sd(RT_log.Z ),
                      se   = sd / sqrt(N) )
summary(RT_log.Z_checker)
##Remove trials where absolute RT.Z>2.5 (i.e. remove outlier RTs)
data<-data[!abs(data$RT_log.Z)>3,]

#Histogram trial by trial RT:
hist_RT <- ggplot(data, aes(RT))  + geom_histogram(aes(y=..count..), colour="black", fill="white") #look at RT collapsed accross light condition
hist_RT + stat_function(fun = dnorm, args = list(mean = mean(data$RT, na.rm = TRUE), sd = sd(data$RT, na.rm = TRUE)), colour = "black", size = 1) + geom_density()
hist_RT + stat_function(fun = dnorm, args = list(mean = mean(data$RT, na.rm = TRUE), sd = sd(data$RT, na.rm = TRUE)), colour = "black", size = 1) + geom_density()


########################
########################
RT_random_intercepts_only<-lmer(log(RT) ~ 1 + (1|ID) +(1|ITI) +(1|MotionDirection) + (1|Hemifield), data = data, REML=FALSE, na.action = na.omit)
RT_Hemifield<-update(RT_random_intercepts_only, .~. + Hemifield)
anova(RT_random_intercepts_only, RT_Hemifield)

require(effects)
plot(allEffects(RT_Hemifield))







source("summarySE.R") 
source("summarySEwithin.R") #function to calculate Std.Er of mean
source("normDataWithin.R")
plotdata <- summarySEwithin(data, measurevar="RT", withinvars=c( "Hemifield"), idvar="ID")
ggplot(plotdata, aes(x=Hemifield, y=RT)) +
    geom_bar(position=position_dodge(.9), colour="Black", stat="identity") + 
    geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=RT-ci, ymax=RT+ci)) + #can change "se" to "ci" if I want to use 95%ci instead
    geom_hline(yintercept=0) +  coord_cartesian(ylim = c(400, 700)) +
    xlab("Hemifield") + ylab("RT (ms)") +
    theme(axis.title.x = element_text(face="bold", size=14),
          axis.text.x  = element_text(face="bold", angle=0,  size=12)) +
    theme(axis.title.y = element_text(face="bold", size=14),
          axis.text.y  = element_text(angle=0, vjust=0.5, size=12)) +
    theme(legend.title = element_text(size=14, face="bold")) +
    theme(legend.text = element_text(size = 12, face = "bold")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
     ggtitle("Healthy Younger")

plotdata

