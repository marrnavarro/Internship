## Script to perform a PCA on gene expression data to see whether the sex and the age of the patients are correlated enough with the first PCs

## Libraries
library(readxl)
library(factoextra)
library(ggplot2)
library(corrplot)

## Load data
pdf <- read_excel("Z:/Documents/Outcome measures overview (V5)_corrected.xlsx", sheet=1)
df <- as.data.frame(pdf)
#Result of voom package: Transformed RNA-seq count data to log2-counts per million (logCPM)
load("Z:/Analysis_projects/RvC_recognition_rnaseq_diff_expr/for_paper_counts/fits_paper/v_visit.RDATA")
rna_counts <- data.frame(v$E)
rna_counts<-rna_counts[,order(colnames(rna_counts))]
rna_counts<-rna_counts[,c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,
                          29,31,33,35,37,39,41,43,45,47,49,51,53,
                          2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,
                          46,48,50,52,54)]

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
df2 <- data.frame(
  df$PatientID...1, df$Sex, df$AgeBaseline, df$DM1ActivCV2, 
  df$SMWTV2, df$AEScScoreV2, df$FDSSV2, df$MDHIV2, 
  df$CISFatigueV2, df$BDIFsV2, df$INQOLQolScoreV2, df$MeanENMOV2, df$M5ENMOV2,
  df$McGillPainV2, df$ASBQV2, df$SSLDScoreV2, StroopInterferenceV2, TMTV2,
  df$SSLNScoreV2, df$SSLIScoreV2, df$JFCSV2, df$ICQV2, df$IMQV2, df$SES28V2, df$CISactivityV2
)
df4<-data.frame(
  df$PatientID...1, df$Sex, df$AgeBaseline,
  df$DM1ActivCV4,
  df$SMWTV4, df$AEScScoreV4, df$FDSSV4, df$MDHIV4,
  df$CISFatigueV4, df$BDIFsV4, df$INQOLQolScoreV4, StroopInterferenceV4, df$MeanENMOV4, df$M5ENMOV4,
  TMTV4, df$McGillPainV4, df$ASBQV4, df$SSLDScoreV4,  
  df$SSLNScoreV4, df$SSLIScoreV4, df$JFCSV4, df$ICQV4, df$IMQV4, df$SES28V4, df$CISactivityV4
)

## Clean up variable names
colnames(df2) <- gsub(x=colnames(df2), pattern = "df.", replacement = "")
colnames(df2) <- gsub(x=colnames(df2), pattern = "V2", replacement = "")
colnames(df4) <- gsub(x=colnames(df4), pattern = "df.", replacement = "")
colnames(df4) <- gsub(x=colnames(df4), pattern = "V4", replacement = "")
colnames(rna_counts) <- gsub(x=colnames(rna_counts), pattern = "rna_counts.", replacement = "")

## Select only the patients with RNA-seq data
patients <- gsub(x=colnames(rna_counts), pattern = "_V4", replacement = "")
patients <- gsub(x=patients, pattern = "_V2", replacement = "")
df2$PatientID...1 <- gsub(x=df2$PatientID...1, pattern = "P", replacement = "")
df2_27<-subset(df2, df2$PatientID...1 %in% patients)
df4$PatientID...1 <- gsub(x=df4$PatientID...1, pattern = "P", replacement = "")
df4_27<-subset(df4, df4$PatientID...1 %in% patients)

## Combine V2 and V4 data frames
#Make the Patient column be the rows
df227<-df2_27[,2:length(colnames(df2_27))]
rownames(df227)<-df2_27[,1]
df427<-df4_27[,2:length(colnames(df4_27))]
rownames(df427)<-df4_27[,1]
df27<-rbind(df227,df427)
rownames(df27)<-colnames(rna_counts)

## PCA
pca<-prcomp(t(rna_counts), scale=T, center=T)

#Scree plot
var<-pca$sdev**2 / sum(pca$sdev**2)
qplot(c(1:54),var)+geom_line()+xlab('Principal Component')+ylab('Variance Explained')+
  ggtitle('Scree Plot')

#Correlation between each PC and sex, correlation between each PC and age
pca_df <- as.data.frame(pca$x)
pca_df <- cbind(pca_df, df27)
pca_tests <- data.frame() 
for (OM in colnames(pca_df)[55:78]){
  for (PCA in colnames(pca_df)[1:54]){  
    formula <- paste(PCA, OM, sep=" ~ ")
    fit <- lm(formula = formula, data=pca_df)
    pca_tests[PCA,OM] <- summary(fit)$r.squared
  }
}
pca_tests <- as.matrix(pca_tests)
corrplot(pca_tests, tl.cex = 0.65,
         is.corr = F
)
#The same but only for the first 17 PCs (for visualization purposes)
pca_tests <- data.frame() 
for (OM in colnames(pca_df)[55:78]){
  for (PCA in colnames(pca_df)[1:17]){  
    formula <- paste(PCA, OM, sep=" ~ ")
    fit <- lm(formula = formula, data=pca_df)
    pca_tests[PCA,OM] <- summary(fit)$r.squared
  }
}
pca_tests <- as.matrix(pca_tests)
corrplot(pca_tests, tl.cex = 0.95,
         is.corr = F
)

## Plot the PCs most correlated with sex and age
p<-ggplot(pca_df,aes(PC9,PC8,color=Sex))+geom_point()
p+labs(x='PC9',y='PC8',colour='Sex')
p<-ggplot(pca_df,aes(PC8,PC10,color=AgeBaseline))+geom_point()
p+labs(x='PC8',y='PC10',colour='Age')

