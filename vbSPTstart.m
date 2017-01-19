function vbSPTstart()
% Add paths that are needed to run the VB3 analysis 
dir0=fileparts(mfilename('fullpath'));% path to this file, even if called from other folder
addpath(dir0)
addpath(fullfile(dir0,'VB3'))
addpath(fullfile(dir0,'HMMcore'))
addpath(fullfile(dir0,'Tools'))
addpath(fullfile(dir0,'external'))
disp('Added local vbSPT paths from')
disp(dir0)
disp('-----------------------')
