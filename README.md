## Table of Contents
*  [Overview](#overview)
*  [Details](#details)
*  [Requirements](#requirements)
*  [Usage](#usage)
*  [Troubleshooting](#troubleshooting)

## Overview
This software package comprises scripts that include the main functions for reproducing the results presented in the paper titled "Deconvolving mutational patterns of poliovirus outbreaks reveals its intrinsic fitness landscape".


## Details
#### Title of paper
Deconvolving mutational patterns of poliovirus outbreaks reveals its intrinsic fitness landscape
#### Authors
Ahmed A. Quadeer, John P. Barton, Arup K. Chakraborty, and Matthew R. McKay

## Requirements
1.  A PC with MATLAB (preferrably v2017a or later) installed on it with the following additional toolboxes:
    * Bioinformatics Toolbox
    * Statistics and Machine Learning Toolbox
    * Curve Fitting Toolbox
    * Parallel Computing Toolbox
    * MATLAB Distributed Computing Server
 
2.  For inferring maximum-entropy model parameters:
    * Adaptive Cluster Expansion (ACE) package, available at https://github.com/johnbarton/ACE 
 
3.  For visualizing antigenic sites on the vp1 crystal structure
    * Pymol, available at https://pymol.org/ 

4.  For determining and visualizing pathways connecting the peaks:
    * MATLAB supported compiler installed for compiling the C file used in implementing the zero temperature Monte Carlo computational method
    * Web version of Circos, available at http://mkweb.bcgsc.ca/tableviewer/visualize/ 

5.  For constructing and visualizing phylogenetic trees:
    * PASTA v1.6.4, available at https://github.com/smirarab/pasta 
    * Dendroscope v3.5.9, available at http://dendroscope.org/ 

6.  For contrasting the output of standard clustering methods with the peaks
    * Python v2.7
    * Python package Seaborn v0.8.1, available at https://seaborn.pydata.org/ 


## Usage
1.  Download the repository
2.  Open MATLAB
3.  Change directory to the downloaded repository 
4.  Run the script main.m 


For visualizing the step-by-step procedure and the corresponding output
1. Download the html folder
2. Open the "main.html" file in your browser

## Troubleshooting
For any questions or comments, please email at ahmedaq@gmail.com. 
