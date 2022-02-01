# nectar_growth_assay
This repository contains all the data and R code for the 2021 growth assay experiments




directory format is as follows

### code
contains r code

* code labelled after a treatment (e.g 4mM_h2o2) takes raw platereader csv files and runs grofit curve fitting then outputs the results into the "output" folder with the ending ...summary.csv
  + these files are all created and can be found in the output folder
* "analyze_all.R" takes all the output files created by each treatments code, aggregates them, and analyzes them
* "cogrowth_cfu_anaylsis.R" analyzes cogrowth assays
* "old_control" folder contains code analyzing data from prelinary data runs that had a different control treatment formulation. (not used in analysis)

### input
contains raw data files

* the files ending in cleaned.csv are raw files from plate reader that were cleaned
* the files ending in formatted.csv are outputs from r that formates and adjusts cleaned.csv files
    
### output
contains output csv files from grofit r code

### final_graphs
contains the output graphs created in "analyze_all.R" and "cogrowth_cfu_analysis.R"





