% Add paths that are needed to run the VB3 analysis 


dir0=pwd;
addpath(genpath([dir0 filesep '.' filesep 'VB3']))
addpath(genpath([dir0 filesep '.' filesep 'HMMcore']))
addpath(genpath([dir0 filesep '.' filesep 'Tools']))
addpath(genpath([dir0 filesep '.' filesep 'external']))
addpath([dir0 filesep '.'])
disp('Added local vbSPT paths')
disp('---------------------')

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if (isOctave)
  disp('Loading Octave packages')
  pkg load specfun
  pkg load struct
 end

clear dir0
