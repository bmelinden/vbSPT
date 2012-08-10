%%%% simulation parameter file %%%%
% this is an example of an HMM analysis runinput file, which specifies
% everything the code needs to know about how to analyz a particular data
% set. To run the HMM analysis, copy this file and edit the vaiables as
% needed, and start the analysis with 
% >> result=VB3_HMManalysis('runinputfilename')
%
% M.L. F.P. 2012-04-17




% Where to load the data and save the result
inputfile='./InputData/testdata_VB3_HMM.mat';
trajectoryfield='finalTraj';  % field name of data
% e.g., a=load(inputfile,trajectoryfield)
% X=a.(trajectoryfield){k} is the position matrix of trajectory k

% This computation can be sped up by running the different optimization
% attempts in parallell (recommended).
parallelize_config=true; % if true; use the settings below.
parallel_start='matlabpool open';  % executed before the parallelizable loop
parallel_end  ='matlabpool close'; % executed after the parallelizable loop
% To change the parallelization settings, either edit the parallel_start
% and parallel_end command strings, or set parallelize_config=false and set
% the configuration outside of the analysis program.

% Where to save the result
outputfile='./Results/testresult_VB3_HMM'; % NOTE: This file will be overwritten !!!
% String to describe the job for own use
jobID='a short string to describe this job';

% Data properties
timestep=3e-3;    % sample time: this sets the time unit in seconds of the analysis
dim=1;      % dimensionality of the data (the first dim columns of the 
            % coordinate matrices will be used)
trjLmin=2;  % trajectories with fewer than trjLmin positions are excluded from 
            % the analysis. Default=2 (recommended; this variable is
            % included for testing purposes).

% Convergence and things to compute:
runs=20;    % number of attempts at each model size (make a multiple of 
            % number of cores when running in parallel)
maxHidden=10;    % maximum number of hidden states to consider. 
            % Aim for ~2 times the number of hidden states.

maxIter=[]; % maximum number of VB iterations ([]: use default values)
relTolF=[]; % convergence criterion for relative change in likelihood bound.
tolPar =[]; % convergence criterion for M-step parameters (leave non-strict).
% For details on the convergence criteria: see help text of VB3_VBEMiterator
% These fields are here mostly for completeness: There should in general be
% no need to change the default values. 

% Compute Viterbi path etc for the best model (sends the 'estimate' argument to
% VB3_VBEMiterator).
stateEstimate=false;

% Run bootstrapping fits in the end. Can be done manually afterwards as well. 
% Set bootstrapNum=0 to disable bootstrapping.
bootstrapNum=100; % number of bootstrap resamplings
fullBootstrap=false; % if true, the best model for all model sizes are bootstrapped, 
                    % and the win probability is estimated for each size 
                    % (computationally costly).

% Parameters for generating initial conditions
init_D=[0.1 10]*1e6; % interval for diffusion constant initial guess given in 
                     % [L^2/s], where L is the length unit ofthe input data.
init_tD=[0.1 1];     % interval for mean dwell time initial guess given in [s].


% Prior distributions
% units: same time units as timestep ([s]), same length unit as input data.
prior_piStrength=5;  % prior strength of initial state distribution (assumed uniform).
prior_D=1e6;         % prior diffusion constant.
prior_Dstrength=5;   % strength of diffusion constant prior, number of pseudocounts (positive).
prior_tD=10*timestep;    % mean dwell time in [s].
prior_tDstrength=2*prior_tD/timestep;  % transition rate strength, number of pseudocounts.

% Our tests suggest that prior_tDstrength=2*prior_tD (plus unit conversion)
% is a good rule of thumb and that it is more important to have a weak prior
% than to get the prior dwell times right. Hence, we believe this prior is
% a good first choice for most people, irrespective of their choice for
% init_tD.

