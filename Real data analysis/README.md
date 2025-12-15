## Real Data Analysis

This folder contains scripts to reproduce the real data analysis on pathway-informed correlation estimation for pan-gynecologic proteomics data.

**Step 1**  
Run `1_Data_extraction.R` to extract the required dataset. Ensure that the working directory is set to the location of this file before execution.

**Step 2**  
Run `SpCov_SCAD_fitting.m` to fit the SCAD-penalized Frobenius loss model with a pathway-wise block-diagonal structure, which enforces joint penalization across all within-pathway protein pairs.  
Execute this script for `case_id = 1` (BRCA), `case_id = 2` (CESC), `case_id = 3` (OV), `case_id = 4` (UCEC), and `case_id = 5` (UCS).

**Step 3**  
Run `3_post_analysis_plot.R`, ensuring that the working directory is set to the location of this file.  
Execute the script for `case_id = 1` (BRCA), `case_id = 2` (CESC), `case_id = 3` (OV), `case_id = 4` (UCEC), and `case_id = 5` (UCS).  

This step generates all figures and plots reported in the main manuscript.

