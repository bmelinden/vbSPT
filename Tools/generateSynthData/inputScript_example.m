% An example script of putting together indata in different ways that should 
% be passed on as arguments to 'VB3_generateSynthData'.
%



%% Define misc parameters

% Choose if to use parallel computing
do_parallel = true;

% Choose how many runs to do
runs = 1;

% Define timestep
timestep = 0.003;  % s

% Define spatial discretization step
stepSize = 5; % nm

% Define localization accuracy
locAccuracy = 20; % nm

% Define cell geometry parameters
cylL = 2000; % nm
cylRadius = 400;  % nm

% Define trajectory parameters
do_single = false;

if ~do_single
    numTraj = 500;
    avTrajLength = 10;
    shortestTraj = 2;
    
    % generate trajectories with exponential length distribution and minimum
    % length shortestTraj
    trajLengths = zeros(1,numTraj);
    for k=1:numTraj
        while(trajLengths(k)<shortestTraj)
            trajLengths(k)=ceil(-avTrajLength*log(1-rand));
        end
    end
    clear k
    
    % one can also insert an arbitrary vector of trajectory lengths here as well.
    %%%%%%%%%%%%%%%%
else
    trajLengths = [500000];
end

%% Define states (max 6)

% Occupancy:
occProb = [];  % Leaving this setting as '[]' will ensure a steady state occupation.

% Diffusion coeffs:
D1 = 1.0;      % um2/s
D2 = 3.0;    % um2/s
D3 = 0;      % um2/s
D4 = 0;      % um2/s
D5 = 0;      % um2/s
D6 = 0;      % um2/s
D = [D1 D2 D3 D4 D5 D6];  % um2/s
clear D1 D2 D3 D4 D5 D6
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
k12 = 15;     % s^-1 
k13 = 0;     % s^-1
k14 = 0;     % s^-1 
k15 = 0;     % s^-1 
k16 = 0;     % s^-1 
k21 = 30;     % s^-1 
k23 = 0;     % s^-1 
k24 = 0;     % s^-1 
k25 = 0;     % s^-1
k26 = 0;     % s^-1 
k31 = 0;     % s^-1 
k32 = 0;     % s^-1
k34 = 0;     % s^-1 
k35 = 0;     % s^-1
k36 = 0;     % s^-1 
k41 = 0;     % s^-1 
k42 = 0;     % s^-1
k43 = 0;     % s^-1 
k45 = 0;     % s^-1
k46 = 0;     % s^-1 
k51 = 0;     % s^-1 
k52 = 0;     % s^-1
k53 = 0;     % s^-1 
k54 = 0;     % s^-1
k56 = 0;     % s^-1
k61 = 0;     % s^-1 
k62 = 0;     % s^-1
k63 = 0;     % s^-1 
k64 = 0;     % s^-1
k65 = 0;     % s^-1


transRate = [0 k12 k13 k14 k15 k16;...
    k21 0 k23 k24 k25 k26;...
    k31 k32 0 k34 k35 k36;...
    k41 k42 k43 0 k45 k46;...
    k51 k52 k53 k54 0 k56;...
    k61 k62 k63 k64 k65 0]; % s^-1 

% fix diagonal so each state sums up to 0
transRate(~~eye(size(transRate))) = -sum(transRate, 2);

% Make into transitionprobability / timestep
transMat = expm(transRate.*timestep);
           
clear k1* k2* k3* k4* k5* k6*

% Remove the transitions belonging to states that should not exist
transRate(length(D)+1:end, :) = [];
transRate(:, length(D)+1:end) = [];
transMat(length(D)+1:end, :) = [];
transMat(:, length(D)+1:end) = [];

% disp('Transition probability matrix [/timestep]:');
% transMat
% disp('Transition rate matrix [/s]:');
% transRate


%% Call the main function which generates the trajectories

if do_single
    finalTraj = VB3_generateSynthData('runs', runs, 'timestep', timestep, 'stepSize', stepSize, 'cylinderLength', cylL, 'cylinderRadius', cylRadius,...
       'trajLengths', trajLengths, 'Dapp', D, 'transRate', transRate, 'occProb', occProb, 'locAccuracy', locAccuracy, 'singleTraj');
elseif do_parallel
    finalTraj = VB3_generateSynthData('runs', runs, 'timestep', timestep, 'stepSize', stepSize, 'cylinderLength', cylL, 'cylinderRadius', cylRadius,...
       'trajLengths', trajLengths, 'Dapp', D, 'transRate', transRate, 'occProb', occProb, 'locAccuracy', locAccuracy, 'parallel');
else
    finalTraj = VB3_generateSynthData('runs', runs, 'timestep', timestep, 'stepSize', stepSize, 'cylinderLength', cylL, 'cylinderRadius', cylRadius,...
       'trajLengths', trajLengths, 'Dapp', D, 'transRate', transRate, 'occProb', occProb, 'locAccuracy', locAccuracy);
end

