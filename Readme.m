%% Software Package

% This software package comprises scripts that include the main
% functions for reproducing the results presented in the paper titled 
% "Deconvolving mutational patterns of poliovirus outbreaks reveals
% its intrinsic fitness landscape".


%% Details

% Title of paper:
% Deconvolving mutational patterns of poliovirus outbreaks reveals its intrinsic fitness landscape
%
% Authors:
% Ahmed A. Quadeer, John P. Barton, Arup K. Chakraborty, and Matthew R. McKay


%% Requirements

% A PC with MATLAB (preferrably v2017a or later) installed on it with the following additional toolboxes:
%   - Bioinformatics Toolbox
%   - Statistics and Machine Learning Toolbox
%   - Curve Fitting Toolbox
%   - Parallel Computing Toolbox
%   - MATLAB Distributed Computing Server

% For inferring maximum-entropy model parameters:
%   - Adaptive Cluster Expansion (ACE) package, available at https://github.com/johnbarton/ACE

% For visualizing antigenic sites on the vp1 crystal structure
%   - Pymol, available at https://pymol.org/

% For identifying and visualizing pathways connecting the peaks:
%   - MATLAB supported compiler installed for compiling the provided C file used in implementing the 
%     zero temperature Monte Carlo computational method.
%   - Web version of Circos, available at http://mkweb.bcgsc.ca/tableviewer/visualize/

% For constructing and visualizing phylogenetic trees:
%   - PASTA v1.6.4, available at https://github.com/smirarab/pasta
%   - Dendroscope v3.5.9, available at http://dendroscope.org/

% For contrasting the output of standard clustering methods with the peaks
%   - Python v2.7
%   - Seaborn v0.8.1

%% Usage: 

% Run the script main.m 

% To visualize the step-by-step procedure and the corresponding output, 
% open the "main.html" file in the "html" folder.
