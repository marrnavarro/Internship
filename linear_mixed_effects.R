## Script to fit the linear mixed-effects models for each scaling approach and each compound score ##

## Libraries
library(readxl)
library(dplyr)
library(mice)
library(wrapr)

## Functions
IQR_filter <- function(x, n){
  # This function removes outliers from a vector based on the Interquartile Range Rule
  # The upper and lower boundaries are respectively Q1-n*IQR and Q3+n*IQR
  Q <- quantile(x, probs=c(0.25, 0.75), na.rm = T)
  iqr <- IQR(x, na.rm = T)
  low <- unname(Q[1] - n*iqr)
  up <- unname(Q[2] + n*iqr)
  x[x < low | x > up] <- NA
  return(x)
}
assign_q <- function(quantiles, value){
  #Scaling approach 2: quantile distribution based
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
load(file = "Z:/Analysis_projects/RvC_recognition_rnaseq_diff_expr/for_paper_counts/metadata/samples.RDATA")
#Result of voom package: Transformed count data to log2-counts per million (logCPM)
hgnc_symbol <- read.table("Z:/Analysis_projects/RvC_recognition_rnaseq_diff_expr/Analyses_rnaseq_recognition/metadata/ENSG_geneSymbol2.txt", sep =",", header=TRUE)
v<-get(load(file="Z:/Analysis_projects/RvC_recognition_rnaseq_diff_expr/for_paper_counts/fits_paper/v_visit.RDATA"))
rownames(v$weights) <- rownames(v$E)
colnames(v$weights) <- colnames(v$E)

## Calculate stroop interference score at V2 and V4
stroop2V2 <- ((60 - df$StroopCardIIErrorsV2) / 60 ) / df$StroopCardIITimeV2
stroop3V2 <- ((60 - df$StroopCardIIIErrorsV2) / 60) / df$StroopCardIIITimeV2
StroopInterferenceV2 <- stroop3V2 / stroop2V2
stroop2V4 <- ((60 - df$StroopCardIIErrorsV4) / 60 ) / df$StroopCardIITimeV4
stroop3V4 <- ((60 - df$StroopCardIIIErrorsV4) / 60) / df$StroopCardIIITimeV4
StroopInterferenceV4 <- stroop3V4 / stroop2V4

## Calculate V2 and V4 TMT score
TMTV2 <- df$TMTBV2 / df$TMTAV2
TMTV4 <- df$TMTBV4 / df$TMTAV4

## Select the desired variables from the data
bdf <- data.frame(
  df$DM1ActivCV2,
  df$SMWTV2, df$AEScScoreV2, df$FDSSV2, df$MDHIV2,
  df$CISFatigueV2, df$BDIFsV2, df$INQOLQolScoreV2, StroopInterferenceV2, df$MeanENMOV2, df$M5ENMOV2,
  TMTV2, df$McGillPainV2, df$ASBQV2, df$SSLDScoreV2,  
  df$SSLNScoreV2, df$SSLIScoreV2, df$JFCSV2, df$ICQV2, df$IMQV2, df$SES28V2, df$CISactivityV2, df$PatientID...1
)
bdf4 <- data.frame(
  df$DM1ActivCV4, 
  df$SMWTV4, df$AEScScoreV4, df$FDSSV4, df$MDHIV4,
  df$CISFatigueV4, df$BDIFsV4, df$INQOLQolScoreV4, StroopInterferenceV4, df$MeanENMOV4, df$M5ENMOV4,
  TMTV4, df$McGillPainV4, df$ASBQV4, df$SSLDScoreV4,  
  df$SSLNScoreV4, df$SSLIScoreV4, df$JFCSV4, df$ICQV4, df$IMQV4, df$SES28V4, df$CISactivityV4, df$PatientID...1
)

## Clean up variable names
colnames(bdf) <- gsub(x=colnames(bdf), pattern = "df.", replacement = "")
colnames(bdf) <- gsub(x=colnames(bdf), pattern = "V2", replacement = "")
colnames(bdf4) <- gsub(x=colnames(bdf4), pattern = "df.", replacement = "")

## Correct specific names
bdf$INQoL <- bdf$INQOLQolScore
bdf$INQOLQolScore <- NULL
bdf4$INQoLV4 <- bdf4$INQOLQolScoreV4
bdf4$INQOLQolScoreV4 <- NULL

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

bdf4$MDHIV4 <- bdf4$MDHIV4*(-1)
bdf4$FDSSV4 <- bdf4$FDSSV4*(-1)
bdf4$CISFatigueV4 <- bdf4$CISFatigueV4*(-1)
bdf4$INQoLV4 <- bdf4$INQoLV4*(-1)
bdf4$BDIFsV4 <- bdf4$BDIFsV4*(-1)
bdf4$AEScScoreV4 <- bdf4$AEScScoreV4*(-1)
bdf4$TMTV4 <- bdf4$TMTV4*(-1)
bdf4$McGillPainV4 <- bdf4$McGillPainV4*(-1)
bdf4$ASBQV4 <- bdf4$ASBQV4*(-1)
bdf4$SSLDScoreV4 <- bdf4$SSLDScoreV4*(-1)
bdf4$SSLNScoreV4 <- bdf4$SSLNScoreV4*(-1)
bdf4$JFCSV4 <- bdf4$JFCSV4*(-1)
bdf4$IMQV4 <- bdf4$IMQV4*(-1)
bdf4$CISactivityV4 <- bdf4$CISactivityV4*(-1)
colnames(bdf4) <- gsub(x=colnames(bdf4), pattern = "V4", replacement = "")

## Combine V2 and V4 data frames
bdf<-bind_rows(bdf,bdf4)                                                        

## Add timepoint identifier
bdf$PatientID[1:255]<-paste(bdf$PatientID[1:255],'V2',sep='_')
bdf$PatientID[256:510]<-paste(bdf$PatientID[256:510],'V4',sep='_')

## Filter outliers
#Set PatientID column as row names because we will not filter outliers there
row.names(bdf)<-bdf$PatientID
bdf <- bdf[,c(1:21,23,22)] 
#Generate backup
cdf <- bdf[,1:22]                                             
sum(is.na(cdf))                                                                 
for(OM in colnames(cdf)){                                                       
  cdf[,OM] <- IQR_filter(cdf[,OM], 3)
}
sum(is.na(cdf))
#Accept all changes
bdf <- cdf

## Impute NA values
imp<-mice(bdf, maxit = 100)
k_bdf<-complete(imp)
plot(imp) 

## Select only the patients with gene expression data
row.names(k_bdf) <- gsub(x=row.names(k_bdf), pattern = "P", replacement = "")
k_bdf<-subset(k_bdf, row.names(k_bdf) %in% colnames(v$E))

## Scale with the R base function
kq_bdf<-as.data.frame(scale(k_bdf),center=F,scale=T)           

## Divide the clinical outcome measures into the clusters
physical<-data.frame(kq_bdf$DM1ActivC,kq_bdf$SMWT,kq_bdf$M5ENMO,kq_bdf$MeanENMO)
cognitive<-data.frame(kq_bdf$TMT,kq_bdf$StroopInterference,kq_bdf$ASBQ,kq_bdf$SSLDScore
                      ,kq_bdf$SSLIScore,kq_bdf$AEScScore,kq_bdf$SSLNScore)
fatigue<-data.frame(kq_bdf$FDSS,kq_bdf$MDHI,kq_bdf$CISFatigue,kq_bdf$BDIFs
                    ,kq_bdf$INQoL,kq_bdf$CISactivity,kq_bdf$SES28,kq_bdf$JFCS,kq_bdf$ICQ
                    ,kq_bdf$IMQ,kq_bdf$McGillPain)

## Clean up variable names
colnames(physical) <- gsub(x=colnames(physical), pattern = "kq_bdf.", replacement = "")
colnames(cognitive) <- gsub(x=colnames(cognitive), pattern = "kq_bdf.", replacement = "")
colnames(fatigue) <- gsub(x=colnames(fatigue), pattern = "kq_bdf.", replacement = "")

## Generate compound scores
physical_score<-data.frame(rowMeans(physical),row.names(kq_bdf))
cognitive_score<-data.frame(rowMeans(cognitive),row.names(kq_bdf))
fatigue_score<-data.frame(rowMeans(fatigue),row.names(kq_bdf))

## Align compound scores and samples data frame based on patient ID and visit 
samples$PatientID<-paste(samples$PatientID,samples$Visit,sep='_')
matchorder<-match_order(physical_score$row.names.kq_bdf.,samples$PatientID)
physical_score<-physical_score[matchorder,]
matchorder<-match_order(fatigue_score$row.names.kq_bdf.,samples$PatientID)
fatigue_score<-fatigue_score[matchorder,]
matchorder<-match_order(cognitive_score$row.names.kq_bdf.,samples$PatientID)
cognitive_score<-cognitive_score[matchorder,]
samples$PatientID <- gsub(x=samples$PatientID, pattern = "_V2", replacement = "")
samples$PatientID <- gsub(x=samples$PatientID, pattern = "_V4", replacement = "")

## Create design matrix
design <- model.matrix(~ samples[,"Visit"] + physical_score$rowMeans.physical. + fatigue_score$rowMeans.fatigue. + cognitive_score$rowMeans.cognitive.)
colnames(design) <- c("(Intercept", "CBT",'Physical','Fatigue','Cognitive')

## Fit mixed-effects model for each gene                                         
lmer_fit <- list()
for (gene in rownames(v$E)){
  lmer_fit[[gene]] <- suppressMessages(
    lmerTest::lmer(v$E[gene,] ~ design[,"CBT"] + design[,'Physical'] +
                     (1|samples[,"PatientID"]), weights = v$weights[gene,]))
}

## Extract regression coefficients and p-values
lmer_fit_coefficients <- lapply(names(lmer_fit), function(x) summary(lmer_fit[[x]])[["coefficients"]])
names(lmer_fit_coefficients) <- names(lmer_fit)

lmer_fit_values <- list()
lmer_fit_values[["CBT"]] <- do.call(rbind, lapply(names(lmer_fit_coefficients), function(x){lmer_fit_coefficients[[x]][2,]}))
lmer_fit_values[["compound"]] <- do.call(rbind, lapply(names(lmer_fit_coefficients), function(x){lmer_fit_coefficients[[x]][3,]}))

## Store relevant results per predictor in list, add gene names and FDR correction
for (fit in names(lmer_fit_values)){
  lmer_fit_values[[fit]] <- cbind(names(lmer_fit_coefficients),
                                  hgnc_symbol$hgnc_symbol[hgnc_symbol$ensembl_gene_id %in% names(lmer_fit_coefficients)],
                                  lmer_fit_values[[fit]], 
                                  p.adjust(lmer_fit_values[[fit]][, "Pr(>|t|)"], method = "fdr"))
  colnames(lmer_fit_values[[fit]]) <-c("ENSG", "hgnc_symbol", "Estimate", "Std.Error", "df", "t_value", "p.val", "FDR")
}

## Convert certain rows to numeric
lmer_fit_values[["compound"]] <- data.frame(lmer_fit_values[["compound"]])
lmer_fit_values[["compound"]][,3:8] <- as.numeric(unlist(lmer_fit_values[["compound"]][,3:8]))

## Temporarily store and save results
compound_fit_new <- lmer_fit_values$compound
write.table(lmer_fit_values[["compound"]], file = paste("Z:/Analysis_projects/MN_analyses/mixed_effects/rbase_scaling/new/physical_lmer_fit_compound_values.csv", sep=""))
save(lmer_fit_values, file = "Z:/Analysis_projects/MN_analyses/mixed_effects/rbase_scaling/new/physical_compound_coef.RDATA")
save(lmer_fit, file = "Z:/Analysis_projects/MN_analyses/mixed_effects/rbase_scaling/new/physical_compound_fits.RDATA")

## Scale with the min max method
probs <- seq(0.1, 1, 0.1)
qm <- apply(k_bdf, 2, function(y){
  vec <- y
  quant <- unname(quantile(vec, probs = probs, na.rm=T))
  
  sapply(vec, function(x){
    assign_q(quant, x)
  })
})
kq_bdf <- as.data.frame(qm)

## Divide clinical outcome measures into the clusters
physical<-data.frame(kq_bdf$DM1ActivC,kq_bdf$SMWT,kq_bdf$M5ENMO,kq_bdf$MeanENMO)
cognitive<-data.frame(kq_bdf$TMT,kq_bdf$StroopInterference,kq_bdf$ASBQ,kq_bdf$SSLDScore
                      ,kq_bdf$SSLIScore,kq_bdf$AEScScore,kq_bdf$SSLNScore)
fatigue<-data.frame(kq_bdf$FDSS,kq_bdf$MDHI,kq_bdf$CISFatigue,kq_bdf$BDIFs
                    ,kq_bdf$INQoL,kq_bdf$CISactivity,kq_bdf$SES28,kq_bdf$JFCS,kq_bdf$ICQ
                    ,kq_bdf$IMQ,kq_bdf$McGillPain)

## Clean up variable names
colnames(physical) <- gsub(x=colnames(physical), pattern = "kq_bdf.", replacement = "")
colnames(cognitive) <- gsub(x=colnames(cognitive), pattern = "kq_bdf.", replacement = "")
colnames(fatigue) <- gsub(x=colnames(fatigue), pattern = "kq_bdf.", replacement = "")

## Generate compound scores
physical_score<-data.frame(rowMeans(physical),row.names(kq_bdf))
cognitive_score<-data.frame(rowMeans(cognitive),row.names(kq_bdf))
fatigue_score<-data.frame(rowMeans(fatigue),row.names(kq_bdf))

## Align compound scores and samples data frame based on patient ID and visit 
samples$PatientID<-paste(samples$PatientID,samples$Visit,sep='_')
matchorder<-match_order(physical_score$row.names.kq_bdf.,samples$PatientID)
physical_score<-physical_score[matchorder,]
matchorder<-match_order(fatigue_score$row.names.kq_bdf.,samples$PatientID)
fatigue_score<-fatigue_score[matchorder,]
matchorder<-match_order(cognitive_score$row.names.kq_bdf.,samples$PatientID)
cognitive_score<-cognitive_score[matchorder,]
samples$PatientID <- gsub(x=samples$PatientID, pattern = "_V2", replacement = "")
samples$PatientID <- gsub(x=samples$PatientID, pattern = "_V4", replacement = "")

## Create design matrix
design <- model.matrix(~ samples[,"Visit"] + physical_score$rowMeans.physical. + fatigue_score$rowMeans.fatigue. + cognitive_score$rowMeans.cognitive.)
colnames(design) <- c("(Intercept", "CBT",'Physical','Fatigue','Cognitive')

## Fit mixed effects model for each gene                                         
lmer_fit <- list()
for (gene in rownames(v$E)){
  lmer_fit[[gene]] <- suppressMessages(
    lmerTest::lmer(v$E[gene,] ~ design[,"CBT"] + design[,'Fatigue'] +
                     (1|samples[,"PatientID"]), weights = v$weights[gene,]))
}

## Extract regression coefficients and p-values
lmer_fit_coefficients <- lapply(names(lmer_fit), function(x) summary(lmer_fit[[x]])[["coefficients"]])
names(lmer_fit_coefficients) <- names(lmer_fit)

lmer_fit_values <- list()
lmer_fit_values[["CBT"]] <- do.call(rbind, lapply(names(lmer_fit_coefficients), function(x){lmer_fit_coefficients[[x]][2,]}))
lmer_fit_values[["compound"]] <- do.call(rbind, lapply(names(lmer_fit_coefficients), function(x){lmer_fit_coefficients[[x]][3,]}))

## Store relevant results per predictor in list, add gene names and FDR correction
for (fit in names(lmer_fit_values)){
  lmer_fit_values[[fit]] <- cbind(names(lmer_fit_coefficients),
                                  hgnc_symbol$hgnc_symbol[hgnc_symbol$ensembl_gene_id %in% names(lmer_fit_coefficients)],
                                  lmer_fit_values[[fit]], 
                                  p.adjust(lmer_fit_values[[fit]][, "Pr(>|t|)"], method = "fdr"))
  colnames(lmer_fit_values[[fit]]) <-c("ENSG", "hgnc_symbol", "Estimate", "Std.Error", "df", "t_value", "p.val", "FDR")
}

## Convert certain rows to numeric
lmer_fit_values[["compound"]] <- data.frame(lmer_fit_values[["compound"]])
lmer_fit_values[["compound"]][,3:8] <- as.numeric(unlist(lmer_fit_values[["compound"]][,3:8]))

## Temporarily store and save results
compound_fit_new <- lmer_fit_values$compound
write.table(lmer_fit_values[["compound"]], file = paste("Z:/Analysis_projects/MN_analyses/mixed_effects/minmax_scaling/new/physical_lmer_fit_compound_values.csv", sep=""))
save(lmer_fit_values, file = "Z:/Analysis_projects/MN_analyses/mixed_effects/minmax_scaling/new/physical_compound_coef.RDATA")
save(lmer_fit, file = "Z:/Analysis_projects/MN_analyses/mixed_effects/minmax_scaling/new/physical_compound_fits.RDATA")
