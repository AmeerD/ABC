# Adjusted Bayesian Completion Rates (ABC)

Sustainable Development Goal 4 Indicator 4.1.2, completion rates for primary, lower secondary, and upper secondary education, is defined as the "percentage of a cohort of children or young people aged 3-5 years above the intended age for the last grade of each level of education who have completed that grade." ([UIS](https://unstats.un.org/sdgs/metadata/files/Metadata-04-01-02.pdf)). Unlike enrolment, estimates of school completion rates require the use of household surveys, such as those perfomed by UNICEF Multiple Indicator Cluster Surveys ([MICS](https://mics.unicef.org/)) and The Demographic and Health Surveys Program ([DHS](https://dhsprogram.com/)). In our 2021 paper, [Adjusted Bayesian Completion Rates (ABC) Estimation](https://github.com/AmeerD/ABC) **INSERT SocArXiV link**, we introduce the ABC model to produce the first complete se of estimates for Indicator 4.1.2 for 153 countries. The model estimates a latent random walk corresponding to the true completion rates while adjusting for a suite of topical data challenges including gaps between survey waves, conflicting estimates, age-misreporting, and delayed completion. Further model details can be found in the paper or by directly viewing the [model](https://github.com/AmeerD/ABC/tree/main/models/ABC_indep.stan).

## Package Installation

The `abcR` package provides the necessary tools to execute the ABC model, including data pre-processing, modelling, post-processing, testing, and visualisation. 

To install the package, execute the following in R:
```r
devtools::install_github("AmeerD/ABC/abcR", ref="main")
library(abcR)
```

## Sample Data

We provide a sample dataset in the [data](https://github.com/AmeerD/ABC/tree/main/data) folder to illustrate the mechanisms of the model. Note that due to data sharing policies, the presented data contains ten anonymised countries and random noise has been added to all observed values. The magnitude of this error exceeds the sampling variation for a great majority of the data points. Readers may experiment with this data using the ABCtest.R script in the [scripts](https://github.com/AmeerD/ABC/tree/main/scripts) folder. This folder also contains the [drake](https://github.com/ropensci/drake) plans used in the actual modelling pipeline.
