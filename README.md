# Adjusted Bayesian Completion Rates (ABC)

Sustainable Development Goal 4 Indicator 4.1.2, completion rates for primary, lower secondary, and upper secondary education, is defined as the "percentage of a cohort of children or young people aged 3-5 years above the intended age for the last grade of each level of education who have completed that grade." ([UIS](https://unstats.un.org/sdgs/metadata/files/Metadata-04-01-02.pdf)). Unlike enrolment, estimates of school completion rates require the use of household surveys, such as those perfomed by UNICEF Multiple Indicator Cluster Surveys ([MICS](https://mics.unicef.org/)) and The Demographic and Health Surveys Program ([DHS](https://dhsprogram.com/)). In our 2021 paper, [Adjusted Bayesian Completion Rates (ABC) Estimation](https://github.com/AmeerD/ABC) **INSERT SocArXiV link**, we introduce the ABC model to produce the first complete set of estimates for Indicator 4.1.2 for 153 countries. The model estimates a latent random walk corresponding to the true completion rates while adjusting for a suite of topical data challenges including gaps between survey waves, conflicting estimates, age-misreporting, and delayed completion. Further model details can be found in the paper or by directly viewing the [model](https://github.com/AmeerD/ABC/tree/main/models/ABC_indep.stan). 

The purpose of this repository is to consolidate the school completion rate analysis pipeline and make it accessible to the public. We have collected all relevant functions and prepared the `abcR` package. It provides the necessary tools to excute the ABC model, and is divided into the following sections:
* Auxiliary functions used either throughout the pipeline, or in the raw data formatting and pre-processing.
* Modelling functions that prepare the cleaned raw data for model execution as well as model execution functions.
* Post-processing functions to digest the Stan output.
* Extraction functions to retrieve quantities of interest.
* Testing functions to assess performance.
* Aggregation functions to combine country level output into region or income group level output.
* Plotting functions for various visualisations.

In addition to the `abcR` package, the ABC Stan file is provided along with a sample dataset and test script to guide users in how to use the provided functions. For a more involved example, the [drake](https://github.com/ropensci/drake) plans used in the actual modelling pipeline are provided [here](https://github.com/AmeerD/ABC/tree/main/scripts/sample_plans.R).

## Package Installation

To install the `abcR` package, execute the following in R:
```r
devtools::install_github("AmeerD/ABC/abcR", ref="main")
library(abcR)
```

## Sample Data

We provide a sample dataset in the [data](https://github.com/AmeerD/ABC/tree/main/data) folder to illustrate the mechanisms of the model. Note that due to data sharing policies, the presented data contains ten anonymised countries and random noise has been added to all observed values. The magnitude of this error exceeds the sampling variation for a great majority of the data points. Readers may experiment with this data using the ABCtest.R script in the [scripts](https://github.com/AmeerD/ABC/tree/main/scripts) folder. 
