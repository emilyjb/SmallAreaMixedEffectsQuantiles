# Code for small area estimation based on the mixed effects linearly interpolated generalized Pareto distribution
This repository contains code to implement the small area estimation methodology described in the paper "Prediction of Small Area Quantiles for the Conservation Effects Assessment Project Using a Mixed Effects Quantile Regression Model" by Emily Berg & Danhyang Lee.
### Main simulation programs
The following two programs produce output similar to the simulation output in the manuscript.
* No transformation (Section 4.1): 
    * SimulationNoTransformation.R -- run simulation
    * OutputNoTransformation.R -- tabulate output
* With transformation (Section 4.2):
    * SimulationWithTransformation.R -- run simulation
    * OutputWithTransformation.R -- tabulate output 
### Verification code and output
To allow the user to verify that the code is working as expected, we saved a test data set in the Rdata file "TestDataSetNoTrans20Oct2018.Rdata." Output from applying the procedure to the test data set is saved in the spreadsheet "SampleOutput20Oct2018.csv." The program "CodeToCompareToTestDataSets.R" reproduces the output in "SampleOutput20Oct2018.csv." For the ALD and LIGPD procedures, the output from "CodeToCompareToTestDataSets.R" will be identical to the outupt in the spreadsheet "SampleOutput20Oct2018.csv." For the NEB predictor, the output from "CodeToCompareToTestDataSets.R" will approximate the outupt in the spreadsheet "SampleOutput20Oct2018.csv" because the NEB predictor uses Monte Carlo to predict the population quantile.
