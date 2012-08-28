% Add paths that are needed to run the VB3 analysis 
% F.P. and M.L. 2012-08-08 

dir0=pwd;
addpath(genpath([dir0 filesep '.' filesep 'VB3']))
addpath(genpath([dir0 filesep '.' filesep 'HMMcore']))
addpath(genpath([dir0 filesep '.' filesep 'Tools']))
addpath([dir0 filesep '.'])
disp('Added local vbSPT paths')
disp('---------------------')