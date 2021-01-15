# Adjusted Bayesian Completion Rates (ABC)

## Package

The `abcR` package provides the necessary tools to execute the ABC model, including data pre-processing, modelling, post-processing, testing, and visualisation. 

To install the package, execute the following in R:
```r
devtools::install_github("AmeerD/ABC/abcR", ref="main")
```

## Sample Data

We provide a sample dataset in the [data](https://github.com/AmeerD/ABC/tree/main/data) folder to illustrate the mechanisms of the model. Note that due to data sharing policies, the presented data contains ten anonymised countries and random noise has been added to all observed values. The magnitude of this error exceeds the sampling variation for a great majority of the data points. Readers may experiment with this data using the ABCtest.R script in the [scripts](https://github.com/AmeerD/ABC/tree/main/scripts) folder. This folder also contains the [drake](https://github.com/ropensci/drake) plans used in the actual modelling pipeline.
