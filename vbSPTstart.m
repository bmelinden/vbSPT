% Add paths that are needed to run the VB3 analysis 
%
% F.P. 2012-07-04
% M.L. 2012-08-08 : slight changes to reflect new folder names

dir0=pwd;
addpath(genpath([dir0 filesep '.' filesep 'HMMcode' filesep 'VB3']))
addpath(genpath([dir0 filesep '.' filesep 'HMMcode' filesep 'HMMcore']))
addpath(genpath([dir0 filesep '.' filesep 'HMMcode' filesep 'Tools']))
addpath([dir0 filesep '.'])
disp('Added local vbSPT paths')
disp('---------------------')