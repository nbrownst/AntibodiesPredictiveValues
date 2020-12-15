# AntibodiesPredictiveValues

#These files include code and tables required to reproduce the calculations in our accompanying paper.



#The file "Code3Dplots.R" produces Figures 1 and 2, general plots of PPV and NPV for reasonable ranges of antibodies tests

#The resulting plots can be viewed in interactive mode. The axes can be rotated to better understand relationships between variables.



#The file "FDAtestsandLocations.R" calculates the PPV and false positive rates for tests operating under an Emergency Use Authorization

#This is the information calculated for Table 1 and Figures 3 and 4

#The file also calculates PPV plots and calculations for the case studies



#The file "formattedTestsMay23.csv" contains input information on sensitivity/specificity for EUA tests.
#An updated version for the revision is "formattedTestsNov23.csv".
#The 89 total tests were then saved with the new name: formattedTestsFull.csv
#If a particular company reported IgG, IgM and combined, we kept the combined only, which left 61 tests, saved in: formattedTestsCombOnly.csv

#The full and combined only files are needed to run "FDAtestsandLocations.R"



#The file "SensitivityAnalysis.R" accompanies the supplementary material, which includes a sensitivity analysis calculating 
effects on results using estimated prevalence from the positive testing rate in the case studies in New York and Chelsea, MA

