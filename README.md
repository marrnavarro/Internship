# Internship
This repository contains all the scripts I created and used to carry out my internship. Next, there is a description of each script.

**-cluster_outcomes.R:** script to cluster the clinical outcome measures using agglomerative hierarchical clustering, fuzzy c-means, HDBSCAN, and PCA. Previous to this, there are data preprocessing steps (recoding, filtering outliers, excluding patients with more than 25% missing data, imputing NA values, and scaling). Once the cluster results are generated, there is an evaluation of the best clustering results.

**-combined_scores_vs_ctg.R:** script to first create the three compound scores and then plot them against the CTG-repeat length. At the beginning, there are data preprocessing steps on the clinical outcome measures (recoding, filtering outliers, excluding patients with more than 25% missing data, imputing NA values, and scaling).

**-pca_age_sex.R:** script to perform a PCA on gene expression data and visualize the results to know whether the sex and age of the patients are enough correlated with the first PCs and, hence, are covariates that should be added in the mixed-effects models. At the end, the PCs that resulted the most correlated with the sex and age variables are plotted againsts sex and age to further check whether they contribute to the division of the patients.

**-linear_mixed_effects.R:** script to fit linear mixed-effects models on gene expression data. It does so for the Physical compound score, to fit the models of the Cognitive and Fatigue compound scores, the word 'Physical' in 'design[,'Physical']' has to be replaced by either 'Cognitive' or 'Fatigue'. The scripts starts with data preprocessing steps on the clinical outcome measures (recoding, filtering outliers, imputing NA values, selecting the observations with available gene expression data, and scaling), followed by the generation of the three compound scores. Then, the models are fitted and there is a FDR correction on the p-values.

**-volcano_plots_histograms.R:** script to.
