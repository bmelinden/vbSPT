% An example script of putting together indata in different ways that should 
% be passed on as arguments to 'VB3_generateSynthData'.
%
%
% F.P. 2012-07-04

%% Define misc parameters

% Choose if to use parallel computing
do_parallel = true;

% Define timestep
timestep = 0.003  % s

% Define localization accuracy
locAccuracy = 0.020 % um

% Define cell geometry parameters
cylL = 2000 % nm
cylRadius = 400  % nm

% Define trajectory parameters
numTraj = 500
avTrajLength = 10
stdTrajLength = 5
shortestTraj = 3

% generate trajectory lengths
trajLengths = ceil(ones(1,numTraj) .* avTrajLength + stdTrajLength .* randn(1, numTraj));
trajLengths(find(trajLengths <= shortestTraj)) = [];

% one can also insert an arbitrary vector of trajectory lengths here as well.

%% Define states (max 5)

% Occupancy:
occProb = [];  % Leaving this setting as '[]' will ensure a steady state occupation.

% Diffusion coeffs:
D1 = 1.0      % um2/s
D2 = 3.0    % um2/s
D3 = 0      % um2/s
D4 = 0      % um2/s
D5 = 0      % um2/s
D = [D1 D2 D3 D4 D5];  % um2/s
clear D1 D2 D3 D4 D5
D(find(D==0)) = [];

if ~isempty(occProb) && length(D)~=length(occProb)
   error('Wrong number of states defined in some step. Leave the states not to be used as 0 in all parameters') 
end

%% Define the transfer matrix

% The transfer matrix A contains the probabilities of transitioning between 
% states during one timestep such that Aij is the probability of transitioning 
% between state i and state j during one timestep.
% The transition matrix can of course be defined manually as for example:
% transMat = [0.9998 0.0001 0.0001;
%            0.0001 0.9998 0.0001;
%            0.0001 0.0001 0.9998];
% or by writing up transition rates as below.

% Transition rates:
k12 = 15     % s^-1 
k13 = 0     % s^-1
k14 = 0     % s^-1 
k15 = 0     % s^-1 
k21 = 30     % s^-1 
k23 = 0     % s^-1 
k24 = 0     % s^-1 
k25 = 0     % s^-1
k31 = 0     % s^-1 
k32 = 0     % s^-1
k34 = 0     % s^-1 
k35 = 0     % s^-1
k41 = 0     % s^-1 
k42 = 0     % s^-1
k43 = 0     % s^-1 
k45 = 0     % s^-1
k51 = 0     % s^-1 
k52 = 0     % s^-1
k53 = 0     % s^-1 
k54 = 0     % s^-1


k = [k12 k13 k14 k15; k21 k23 k24 k25; k31 k32 k34 k35; k41 k42 k43 k45; k51 k52 k53 k54]; % s^-1 
% Make into transitionprobability / timestep
k = k.*timestep;

%make transition matrix
transMat = [1-sum(k(1,:)), k(1, 1:4);
            k(2, 1), 1-sum(k(2,:)), k(2, 2:4);
            k(3, 1:2), 1-sum(k(3,:)), k(3, 3:4);
            k(4, 1:3), 1-sum(k(4,:)), k(4, 4);
            k(4, 1:4), 1-sum(k(5,:))];
           
clear k1* k2* k3* k4* k5*

% Remove the tranistions belonging to states that should not exist
transMat(length(D)+1:end, :) = [];
transMat(:, length(D)+1:end) = [];


%% Call the main function which generates the trajectories

if do_parallel
finalTraj = VB3_generateSynthData('timestep', timestep, 'cylinderLength', cylL, 'cylinderRadius', cylRadius,...
    'trajLengths', trajLengths, 'Dapp', D, 'transMat', transMat, 'occProb', occProb, 'locAccuracy', locAccuracy, 'parallel');
else
    finalTraj = VB3_generateSynthData('timestep', timestep, 'cylinderLength', cylL, 'cylinderRadius', cylRadius,...
    'trajLengths', trajLengths, 'Dapp', D, 'transMat', transMat, 'occProb', occProb, 'locAccuracy', locAccuracy);

end


