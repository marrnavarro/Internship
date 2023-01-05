## Script to cluster the clinical outcome measures and to validate the resulting clusters ##

## Libraries
library("readxl")
library(dplyr)
library(mice)
library('ClusterR')
library(cluster)
library(factoextra)
library(corrplot)
library(e1071)
library(fpc)
library(umap)
library(dbscan)
library(clValid)

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
  #Scaling function: quantile distribution based
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
  df$DM1ActivCV2,  #df$L5ENMOV2, df$AESIV2, df$CSIV2,
  df$SMWTV2, df$AEScScoreV2, df$FDSSV2, df$MDHIV2,
  df$CISFatigueV2, df$BDIFsV2, df$INQOLQolScoreV2, StroopInterferenceV2, df$MeanENMOV2, df$M5ENMOV2,
  TMTV2, df$McGillPainV2, df$ASBQV2, df$SSLDScoreV2,  
  df$SSLNScoreV2, df$SSLIScoreV2, df$JFCSV2, df$ICQV2, df$IMQV2, df$SES28V2, df$CISactivityV2
)
bdf4 <- data.frame(
  df$DM1ActivCV4,  #df$L5ENMOV2, df$AESIV2, df$CSIV2,
  df$SMWTV4, df$AEScScoreV4, df$FDSSV4, df$MDHIV4,
  df$CISFatigueV4, df$BDIFsV4, df$INQOLQolScoreV4, StroopInterferenceV4, df$MeanENMOV4, df$M5ENMOV4,
  TMTV4, df$McGillPainV4, df$ASBQV4, df$SSLDScoreV4,  
  df$SSLNScoreV4, df$SSLIScoreV4, df$JFCSV4, df$ICQV4, df$IMQV4, df$SES28V4, df$CISactivityV4
)

## Clean up variable names and change all variable types to numeric
colnames(bdf) <- gsub(x=colnames(bdf), pattern = "df.", replacement = "")
colnames(bdf) <- gsub(x=colnames(bdf), pattern = "V2", replacement = "")
bdf <- sapply(bdf[,1:length(colnames(bdf))], function(x){as.numeric(as.character(x))})
bdf <- as.data.frame(bdf)
colnames(bdf4) <- gsub(x=colnames(bdf4), pattern = "df.", replacement = "")
bdf4 <- sapply(bdf4[,1:length(colnames(bdf4))], function(x){as.numeric(as.character(x))})
bdf4 <- as.data.frame(bdf4)

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

## Filter outliers
cdf <- bdf
sum(is.na(cdf))
for(OM in colnames(cdf)){
  cdf[,OM] <- IQR_filter(cdf[,OM], 3)
}
sum(is.na(cdf))
#Manually check flagged outliers
apply(bdf, 2, function(x) length(which(!is.na(x)))) - apply(cdf, 2, function(x) length(which(!is.na(x))))
boxplot(bdf$M5ENMO)
boxplot(cdf$M5ENMO)
boxplot(bdf$MeanENMO)
boxplot(cdf$MeanENMO)
boxplot(bdf$TMT)
boxplot(cdf$TMT)
boxplot(bdf$BDIFs)
boxplot(cdf$BDIFs)
#Accept all changes
bdf <- cdf

## Exclude patients with more than 25% missing values
pat_NA <- apply(bdf[1:510,1:22], 1, function(x) sum(is.na(x)))
pat_NA_ID <- c()
for (i in 1:510){ 
  if (pat_NA[i] > 5.5) {
    pat_NA_ID <- c(pat_NA_ID,i)}
}
bdf <- bdf[!row.names(bdf) %in% pat_NA_ID,]

## Impute NA values (for PCA, HDBSCAN, and C-means)
imp<-mice(bdf, maxit = 35)
k_bdf<-complete(imp)
plot(imp)

## Scaling
#Change all scores in the baseline dataframe to decile scores
#outer loop: apply, calculates quantile boundaries per questionnaire
#inner loop: sapply, uses these boundaries to change every value
probs <- seq(0.1, 1, 0.1)
qm <- apply(k_bdf, 2, function(y){
  vec <- y
  quant <- unname(quantile(vec, probs = probs, na.rm=T))
  
  sapply(vec, function(x){
    assign_q(quant, x)
  })
})
kq_bdf <- as.data.frame(qm)
#do the same for the non-imputed data
qm_noimp <- apply(bdf, 2, function(y){
  vec <- y
  quant <- unname(quantile(vec, probs = probs, na.rm=T))
  
  sapply(vec, function(x){
    assign_q(quant, x)
  })
})
q_bdf <- as.data.frame(qm_noimp)

## PCA
pca<-prcomp((k_bdf), scale=T, center=T)

#Scree plot
var<-pca$sdev**2 / sum(pca$sdev**2)
qplot(c(1:22),var)+geom_line()+xlab('Principal Component')+ylab('Variance Explained')+
  ggtitle('Scree Plot')

#Correlation between each PC and each clinical outcome
pca_df <- as.data.frame(pca$x)
pca_df <- cbind(pca_df, k_bdf)
pca_tests <- data.frame() 
for (OM in colnames(pca_df)[23:44]){
  for (PCA in colnames(pca_df)[1:22]){  
    formula <- paste(PCA, OM, sep=" ~ ")
    fit <- lm(formula = formula, data=pca_df)
    pca_tests[PCA,OM] <- summary(fit)$r.squared
  }
}
pca_tests <- as.matrix(pca_tests)
corrplot(pca_tests,
         is.corr = F
         )

## Agglomerative Hierarchical Clustering
q_d <- dist(t(q_bdf), method = "euclidian")
x <- dist(t(q_bdf))

q_hc_w <- hclust(q_d,method= "ward.D")
plot(q_hc_w,cex=0.75)

q_hc_c <- hclust(q_d,method= "complete") 
plot(q_hc_c,cex=0.75)

q_hc_a <- hclust(q_d, method= "average")
plot(q_hc_a,cex=0.75)

## Fuzzy C-means Clustering
eu<-cmeans(t(kq_bdf),centers=3,dist='euclidean')
man<-cmeans(t(kq_bdf),centers=3,dist='manhattan')
#Plots
corrplot(eu$membership,is.corr=FALSE)
fviz_cluster(list(data=t(kq_bdf),cluster=eu$cluster))
corrplot(man$membership,is.corr=FALSE)
fviz_cluster(list(data=t(kq_bdf),cluster=man$cluster))

#The same but with 4 centers
eu4<-cmeans(t(kq_bdf),centers=4,dist='euclidean')
man4<-cmeans(t(kq_bdf),centers=4,dist='manhattan')
#Plots
corrplot(eu4$membership,is.corr=FALSE)
fviz_cluster(list(data=t(kq_bdf),cluster=eu4$cluster))
corrplot(man4$membership,is.corr=FALSE)
fviz_cluster(list(data=t(kq_bdf),cluster=man4$cluster))

## HDBSCAN with UMAP
#Check which is the best number of neighbors (the one with the highest mean cluster-membership probability)
nei<-3
graph<-c()
for (i in 1:17){
  nei<-nei+1
  u<-umap(t(kq_bdf),min_dist=0.1,n_neighbors=nei)
  hb<-hdbscan(u$layout,minPts = 2)
  graph<-c(graph,mean(hb$membership_prob))
}
plot(c(4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),graph,xlab = 'Number of Neighbors',ylab='Mean Membership')
lines(4:20,graph)

u<-umap(t(kq_bdf),min_dist=0.1,n_neighbors=5)
hb<-hdbscan(u$layout,minPts = 2)
plot(hb$hc)
fviz_cluster(list(data=t(kq_bdf),cluster=hb$cluster))

## Evaluate performance of the best clustering results
#Agglomerative hierarchical clustering
cluster.stats(d=q_d,cutree(q_hc_a,k=4))$dunn
cluster.stats(d=q_d,cutree(q_hc_w,k=3))$dunn
#HDBSCAN
mean(hb$membership_prob)
hb$cluster_scores
