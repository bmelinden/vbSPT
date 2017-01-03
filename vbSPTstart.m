% Add paths that are needed to run the VB3 analysis 


dir0=fileparts(mfilename('fullpath'));% path to this file, even if called from other folder
addpath(genpath([dir0 filesep '.' filesep 'VB3']))
addpath(genpath([dir0 filesep '.' filesep 'HMMcore']))
addpath(genpath([dir0 filesep '.' filesep 'Tools']))
addpath(genpath([dir0 filesep '.' filesep 'external']))
addpath([dir0 filesep '.'])
disp('Added local vbSPT paths')
disp('---------------------')

clear dir0
