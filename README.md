# Drug Efficacy Analysis

This project aims to assess the risk of event in a sample population prescribed with either Drug A or B.

As the sample patient population was monitored for the occurrence of a specific event after start of treatment, 
a survival analysis approach is deemed appropriate to analyse this dataset.

In order to reduce confounding, especially in an unbalanced real-world patient dataset, several measures were taken 
to reduce bias as much as possible. For instance, missing data and standardized mean differences were examined, 
of which the result was used to fit a propensity model.

Matching covariates produces a more well-balanced dataset, that is consequently fitted to estimate treatment effect, 
and used to generate survival curves and hazard ratios. 

## Building and Running

To execute the code in the repository, open RStudio, select all code and run
