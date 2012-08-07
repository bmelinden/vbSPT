function finalTraj=VB3_generateSynthData(varargin)
%% finalTraj=VB3_generateSynthData(runinput, varargin)
%
% Generates synthetic trajectory data in an E. coli like geometry (can be 
% specified within this file)data partially based on a HMM model. All parameters
% from the HMM model except the ones given as options. This version accepts
% either a runinput file or a struct as the first argument or only options.
% This script uses parallel computing by default to disable it comment out the
% relevant rows (containing 'matlabpool') in this function.
% 
% Note: Currently this file saves the generated trajectories to files called 
% 'syntheticData_XX.mat', where XX is an autogenerated index to prevent overwriting,
% in the 'generated Data' subfolder. 
% 
% options:
% 'runs'       : if given, determine how many datasets to be generated with
%               the same input parameters.
% 'timestep'   : should be given in [s].
% 'locAccuracy': should be given in [nm]. Default = 0. It dose actually not
%               add to the diffusion constant that is set (which is the 
%               apperent diffusion) but merely removes its average
%               contribution from the apparent diffusion and then is added 
%               as a positioning uncertainty in every point. 
% 'transMat'   : the transition matrix. Should be an N*N matrix where N is 
%               the number of hidden states and the ij element gives the 
%               probability of transitioning from state i to state j during
%               one timestep. Each row should sum up to 1.
% 'occProb'    : a 1*N array with the occupancy of state 1 to N. Should sum
%               up to 1. If an empty array is given then use the steady state 
%               occupancy calculated from the transfer matrix. 
% 'Dapp'       : a 1*N array with the apparent diffusion constant in [um^/s].
% 'trajLengths': a 1*m array with the trajectory lengths for m trajectories. 
% 'cylinderLength'   : should be given in [nm]. Default value 2000 nm.
% 'cylinderRadius'   : should be given in [nm]. Default value 400 nm.
% 'parallel'    : if given use parallel computing.
%
% F.P. 2012-07-03
%
%
% Change log:
% F.P. 2012-07-04 : Added functionality that it saves in a subfolder to the
%                 script with index that avoids overwriting data.
% F.P. 2012-07-08 : Added so that a steady state occupancy is used if an
%                 empty argument is given for occProb.
% F.Pf 2012-07-10 : Added 'parallel' as an option.

tgenData=tic;



%% define default geometry

CylinderL = 2000 % nm (length of cylindrical part only)
Radius = 400  % nm    (spherical end caps, and cylinder radius)

%% read input

% if an existing file, generate options structure
if(ischar(varargin{1}) && exist(varargin{1}, 'file')==2)
    runinputfile = varargin{1};
    opt=VB3_getOptions(runinputfile);
    disp(['Read runinput file ' runinputfile])
    runinputExists = true;
    % if an option struct, read in the runinputfilename
elseif(isstruct(varargin{1}))
    opt=varargin{1};
    runinputfile=opt.runinputfile;
    disp(['Read options structure based on runinput file ' runinputfile ])
    runinputExists = true;
else
    runinputExists = false;
end
%% load data and read options
if runinputExists
    
    res=load(opt.outputfile, 'Wbest');
    
    Wbest = res.Wbest;
    clear res;
    
    % initiate options
    timestep = opt.timestep; % [s]
    locAccuracy = 0; %[nm]
    transMat = Wbest.est.Amean; % [/timestep]
    occProb = Wbest.est.Ptot;
    Dapp = Wbest.est.DdtMean./timeStep;
    if max(Dapp)<100    % if small assume its in um^2/s
        Dapp = Dapp*1e6;% convert to nm^2/s
    end
    trajLengths = Wbest.T;
else
    % initiate options
    timestep = 0; % [s]
    locAccuracy = 0; %[nm]
    transMat = 0; % [/timestep]
    occProb = 0;
    Dapp = 0;
    trajLengths = 0;
end
   
runs = 1;
do_steadystate = false;
do_parallel = false;

if(nargin>1)        % parse options
    % argument counter
    if runinputExists 
        k=2;
    else
        k=1;
    end
    kmax=nargin;  % stopping criterion
    while(k<=kmax)
        option=varargin{k};
        if(strcmpi(option,'parallel'))
            do_parallel = true;
            k=k+1;
        elseif(strcmpi(option,'timestep'))
            if(~isempty(varargin{k+1}))
                timestep=varargin{k+1};
                if(~isnumeric(timestep) || timestep<=0)
                    error('VB3_synthData: timestep option must be followed by a positive number.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'runs'))
            if(~isempty(varargin{k+1}))
                runs=varargin{k+1};
                if(~isnumeric(runs) || runs<=0 || runs~=round(runs))
                    error('VB3_synthData: runs option must be followed by a positive integer.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'locAccuracy'))
            if(~isempty(varargin{k+1}))
                locAccuracy=varargin{k+1};
                if(~isnumeric(locAccuracy) || locAccuracy<=0)
                    error('VB3_synthData: locAccuracy option must be followed by a positive number.')
                end
            end
            k=k+2;
       elseif(strcmpi(option,'cylinderLength'))
            if(~isempty(varargin{k+1}))
                CylinderL=varargin{k+1};
                if(~isnumeric(CylinderL) || CylinderL<=0)
                    error('VB3_synthData: CylinderL option must be followed by a positive number.')
                end
            end
            k=k+2;
       elseif(strcmpi(option,'cylinderRadius'))
            if(~isempty(varargin{k+1}))
                Radius=varargin{k+1};
                if(~isnumeric(Radius) || Radius<=0)
                    error('VB3_synthData: Radius option must be followed by a positive number.')
                end
            end
            k=k+2;
         elseif(strcmpi(option,'transMat'))
            if(~isempty(varargin{k+1}))
                transMat=varargin{k+1};
                [m, n] = size(transMat);
                temp = sum(transMat');
                if(~isnumeric(transMat) | m~=n | temp~=ones(1,m))
                    error('VB3_synthData: transMat option must be followed by a square numeric matrix where the rows sum up to 1.')
                end
            end
            k=k+2;
         elseif(strcmpi(option,'occProb'))
            if(~isempty(varargin{k+1}))
                occProb=varargin{k+1};
                [m, ~] = size(occProb);
                if(~isnumeric(occProb) | occProb<=0 | m ~= 1 | sum(occProb) ~= 1)
                    error('VB3_synthData: occProb option must be followed by a 1*N numeric array which sums up to 1 or [].')
                end
            else
                do_steadystate = true;
            end
            k=k+2;
        elseif(strcmpi(option,'Dapp'))
            if(~isempty(varargin{k+1}))
                Dapp=varargin{k+1};
                [m, ~] = size(occProb);
                if(~isnumeric(Dapp) | Dapp<=0 | m ~= 1 | max(Dapp) > 100)
                    error('VB3_synthData: Dapp option must be followed by a 1*N numeric array in units um^2/s.')
                end
                Dapp = Dapp*1e6;
            end
            k=k+2;
        elseif(strcmpi(option,'trajLengths'))
            if(~isempty(varargin{k+1}))
                trajLengths=varargin{k+1};
                [m, ~] = size(trajLengths);
                if(~isnumeric(trajLengths) | trajLengths<=0 | m ~= 1 | trajLengths~=round(trajLengths))
                    error('VB3_synthData: trajLengths option must be followed by a 1*N numeric integer array.')
                end
            end
            k=k+2;
        else
            error(['VB3_synthData: option ' option ' not recognized.'])
        end
    end
end

% Use steady state occupancy if desired
if do_steadystate
    occProb = transMat^1000;
    occProb = occProb(1,:)
end

if(timestep==0 | sum(sum(transMat))==0 | sum(occProb)==0 | sum(Dapp)==0 | sum(trajLengths)==0) 
    error('Not a valid input, either a runinputfile/struct has to be the first argument or all options must be specified.');
end


% Convert transMat to dwelltimes
dwellTimes=1./(1-diag(transMat))'

% Convert apparent diffusion to actual diffusion to go into the simulation
diffCoeff = Dapp-locAccuracy^2/timestep;

% Plot the trajectory length distribution
figure(1);
clf
hist(trajLengths,0:100);

%% Make the trajectories
% check if matlabpool is open, if yes close it
temp = matlabpool('size');
if temp~=0
matlabpool close
end

% start parallel computing if the option is given
if do_parallel
matlabpool open
end

for m=1:runs
tic    
[finalTraj, ~] = MakeTrajectories(timestep, CylinderL, Radius, trajLengths, diffCoeff, transMat, occProb, locAccuracy);
toc

p = mfilename('fullpath');
[p, ~, ~] = fileparts(p);
oldFold = cd([p filesep 'generated Data' filesep]);
cont = what;
ind = length(cont.mat)+1;
savename=sprintf('syntheticData_%02d.mat',ind);
save(savename)
cd(oldFold);
pause(1)
end

% close parallel computing
if do_parallel
matlabpool close
end

disp(['Finished generating synthetic data in ' num2str(toc(tgenData)/60) ' min.']);
end
