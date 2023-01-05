## Script to create the compound scores for the baseline data and to then plot them against the CTG-repeat length, for 2 different scaling approaches ##

## Libraries
library("readxl")
library(dplyr)
library(mice)

## Functions
IQR_filter <- function(x, n){
  #This function removes outliers from a vector based on the Interquartile Range Rule
  #The upper and lower boundaries are respectively Q1-n*IQR and Q3+n*IQR
  Q <- quantile(x, probs=c(0.25, 0.75), na.rm = T)
  iqr <- IQR(x, na.rm = T)
  low <- unname(Q[1] - n*iqr)
  up <- unname(Q[2] + n*iqr)
  x[x < low | x > up] <- NA
  return(x)
}
assign_q <- function(quantiles, value){
  ## Scaling approach 2: quantile distribution based
  #quantiles: questionnaire specific quantile boundaries
  #value: specific questionnaire value
  #returns matching quantile
  if(is.na(value)){
    return(NA)
  }
  else{
    vec <- c(quantiles, value)
    vec <- sort(vec)
    return((min(which(value==vec)))/length(quantiles))
  }
}

## Load data
pdf <- read_excel("Z:/Documents/Outcome measures overview (V5)_corrected.xlsx", sheet=1)
df <- as.data.frame(pdf)

## Calculate stroop interference score at V2
stroop2V2 <- ((60 - df$StroopCardIIErrorsV2) / 60 ) / df$StroopCardIITimeV2
stroop3V2 <- ((60 - df$StroopCardIIIErrorsV2) / 60) / df$StroopCardIIITimeV2
StroopInterferenceV2 <- stroop3V2 / stroop2V2

## Calculate V2 TMT score
TMTV2 <- df$TMTBV2 / df$TMTAV2

## Select the desired variables from the data
bdf <- data.frame(
  df$DM1ActivCV2,
  df$SMWTV2, df$AEScScoreV2, df$FDSSV2, df$MDHIV2,
  df$CISFatigueV2, df$BDIFsV2, df$INQOLQolScoreV2, StroopInterferenceV2, df$MeanENMOV2, df$M5ENMOV2,
  TMTV2, df$McGillPainV2, df$ASBQV2, df$SSLDScoreV2,  
  df$SSLNScoreV2, df$SSLIScoreV2, df$JFCSV2, df$ICQV2, df$IMQV2, df$SES28V2, df$CISactivityV2, df$V2Mode
)

## Clean up variable names
colnames(bdf) <- gsub(x=colnames(bdf), pattern = "df.", replacement = "")
colnames(bdf) <- gsub(x=colnames(bdf), pattern = "V2", replacement = "")

## Correct specific names
bdf$INQoL <- bdf$INQOLQolScore
bdf$INQOLQolScore <- NULL

## Recode outcome measures: higher scores beneficial
bdf$MDHI <- bdf$MDHI*(-1)
bdf$FDSS <- bdf$FDSS*(-1)
bdf$CISFatigue <- bdf$CISFatigue*(-1)
bdf$INQoL <- bdf$INQoL*(-1)
bdf$BDIFs <- bdf$BDIFs*(-1)
bdf$AEScScore <- bdf$AEScScore*(-1)
bdf$TMT <- bdf$TMT*(-1)
bdf$McGillPain <- bdf$McGillPain*(-1)
bdf$ASBQ <- bdf$ASBQ*(-1)
bdf$SSLDScore <- bdf$SSLDScore*(-1)
bdf$SSLNScore <- bdf$SSLNScore*(-1)
bdf$JFCS <- bdf$JFCS*(-1)
bdf$IMQ <- bdf$IMQ*(-1)
bdf$CISactivity <- bdf$CISactivity*(-1)

## Change variable type to numeric so that NA values are recognized as such
bdf$Mode<-as.numeric(bdf$Mode)

## Filter outliers
#Set Mode as the last column not to take it when filtering the outliers
bdf<-bdf[,c(1:21,23,22)]
#Generate backup
cdf <- bdf[,1:22]
sum(is.na(cdf))
for(OM in colnames(cdf)){
  cdf[,OM] <- IQR_filter(cdf[,OM], 3)
}
sum(is.na(cdf))
#Accept all changes
bdf <- data.frame(cdf,bdf$Mode)

## Exclude patients with more than 25% missing values
pat_NA <- apply(bdf[1:255,1:22], 1, function(x) sum(is.na(x)))
pat_NA_ID <- c()
for (i in 1:255){ 
  if (pat_NA[i] > 5.5) {
    pat_NA_ID <- c(pat_NA_ID,i)}
}
bdf <- bdf[!row.names(bdf) %in% pat_NA_ID,]

## Impute NA values
imp<-mice(bdf[,1:23], maxit = 60)
k_bdf<-complete(imp)
plot(imp)

## Scale with the R base function
kq_bdf<-as.data.frame(scale(k_bdf[1:22]),center=F,scale=T)
kq_bdf<-data.frame(kq_bdf,bdf$bdf.Mode)

## Divide the clinical outcome measures into the clusters
physical<-data.frame(kq_bdf$DM1ActivC,kq_bdf$SMWT,kq_bdf$M5ENMO,kq_bdf$MeanENMO)
cognitive<-data.frame(kq_bdf$TMT,kq_bdf$StroopInterference,kq_bdf$ASBQ,kq_bdf$SSLDScore
                      ,kq_bdf$SSLIScore,kq_bdf$SSLNScore,kq_bdf$AEScScore)
fatigue<-data.frame(kq_bdf$FDSS,kq_bdf$MDHI,kq_bdf$CISFatigue,kq_bdf$BDIFs,kq_bdf$McGillPain
                    ,kq_bdf$INQoL,kq_bdf$CISactivity,kq_bdf$SES28,kq_bdf$JFCS,kq_bdf$ICQ
                    ,kq_bdf$IMQ)

## Clean up variable names
colnames(physical) <- gsub(x=colnames(physical), pattern = "kq_bdf.", replacement = "")
colnames(cognitive) <- gsub(x=colnames(cognitive), pattern = "kq_bdf.", replacement = "")
colnames(fatigue) <- gsub(x=colnames(fatigue), pattern = "kq_bdf.", replacement = "")

## Create the compound scores
physical_score<-rowMeans(physical)
cognitive_score<-rowMeans(cognitive)
fatigue_score<-rowMeans(fatigue)

## Plot the compound scores against the CTG
plot(kq_bdf$bdf.bdf.Mode,physical_score,ylab = 'Physical Score',xlab ='CTG Repeat')
plot(kq_bdf$bdf.bdf.Mode,cognitive_score,ylab = 'Cognitive Score',xlab ='CTG Repeat')
plot(kq_bdf$bdf.bdf.Mode,fatigue_score,ylab = 'Fatigue Score',xlab ='CTG Repeat')

## Scaling approach 2: quantile distribution based
#Change all scores in the baseline dataframe to decile scores
#outer loop: apply, calculates quantile boundaries per questionnaire
#inner loop: sapply, uses these boundaries to change every value
probs <- seq(0.1, 1, 0.1)
qm <- apply(k_bdf[1:22], 2, function(y){
  vec <- y
  quant <- unname(quantile(vec, probs = probs, na.rm=T))
  
  sapply(vec, function(x){
    assign_q(quant, x)
  })
})
kq_bdf <- as.data.frame(qm)
kq_bdf<-data.frame(kq_bdf,bdf$bdf.Mode)

## Divide the clinical outcome measures into the clusters
physical<-data.frame(kq_bdf$DM1ActivC,kq_bdf$SMWT,kq_bdf$M5ENMO,kq_bdf$MeanENMO)
cognitive<-data.frame(kq_bdf$TMT,kq_bdf$StroopInterference,kq_bdf$ASBQ,kq_bdf$SSLDScore
                      ,kq_bdf$SSLIScore,kq_bdf$SSLNScore,kq_bdf$AEScScore)
fatigue<-data.frame(kq_bdf$FDSS,kq_bdf$MDHI,kq_bdf$CISFatigue,kq_bdf$BDIFs,kq_bdf$McGillPain
                    ,kq_bdf$INQoL,kq_bdf$CISactivity,kq_bdf$SES28,kq_bdf$JFCS,kq_bdf$ICQ
                    ,kq_bdf$IMQ)

## Clean up variable names
colnames(physical) <- gsub(x=colnames(physical), pattern = "kq_bdf.", replacement = "")
colnames(cognitive) <- gsub(x=colnames(cognitive), pattern = "kq_bdf.", replacement = "")
colnames(fatigue) <- gsub(x=colnames(fatigue), pattern = "kq_bdf.", replacement = "")

## Create the compound scores
physical_score<-rowMeans(physical)
cognitive_score<-rowMeans(cognitive)
fatigue_score<-rowMeans(fatigue)

## Plot the compound scores against the CTG
plot(kq_bdf$bdf.bdf.Mode,physical_score,ylab = 'Physical Score',xlab ='CTG Repeat')
plot(kq_bdf$bdf.bdf.Mode,cognitive_score,ylab = 'Cognitive Score',xlab ='CTG Repeat')
plot(kq_bdf$bdf.bdf.Mode,fatigue_score,ylab = 'Fatigue Score',xlab ='CTG Repeat')

