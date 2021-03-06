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
To allow the user to verify that the code is working as expected, we saved a test data set in the Rdata file "TestDataSetNoTrans28March2019.Rdata." Output from applying the procedure to the test data set is saved in the spreadsheet "SampleOutput28March2019.csv." The program "CodeToCompareToTestDataSets.R" reproduces the output in "SampleOutput28March2019.csv."
