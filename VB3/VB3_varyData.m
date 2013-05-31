function modelRun = VB3_varyData(runinput, varargin)
%% wModels=VB3_varyData(runinput, varargin)
%
% Takes a previously used runinputfile or options struct as input. 
% Converges the best models for different number of states with an increasing amount of 
% input data, sampled from the original input data.
% This script uses parallel computing by default to disable it comment out the
% relevant rows (containing 'matlabpool') in this function.
%
% Options:
% 'noRandom'    : if given, the data subset isnt randomly generated but  truncates
%               from the data array.
% 'runs'        : an integer defining the number of runs you want. Default
%               is 1.
% 'numStates'    : an 1*m integer array with the statesizes that should be
%                   converged. Default is running the 3 state and 4 state
%                   models.
% 'save'          : If given, it saves the data with the ending '_varData'.
% 'sampleSize'   : an 1*N numeric array with how big portions of the data
%               should be analysed. Default is [0.05:0.05:1], i.e.
%               starting with 5% and increasing 5% until 100%.


%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_varyData, test inference robustness, in the vbSPT package
% =========================================================================
% 
% Copyright (C) 2013 Martin Lindén and Fredrik Persson
% 
% E-mail: bmelinden@gmail.com, freddie.persson@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.   
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
%  Additional permission under GNU GPL version 3 section 7
%  
%  If you modify this Program, or any covered work, by linking or combining it
%  with Matlab or any Matlab toolbox, the licensors of this Program grant you 
%  additional permission to convey the resulting work.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%% start of actual code


tvaryData=tic;
%% read input

% if an existing file, generate options structure
if(isstr(runinput) && exist(runinput)==2)
    runinputfile = runinput;
    opt=VB3_getOptions(runinputfile);
    disp(['Read runinput file ' runinputfile])
    % if an option struct, read in the runinputfilename
elseif(isstruct(runinput))
    opt=runinput;
    runinputfile=opt.runinputfile;
    disp(['Read options structure based on runinput file ' runinputfile ])
else
    error(['Not a valid input, aborting']);
end

% Get input data
X = VB3_readData(opt);
Ntraj = length(X); % number of trajectories

% Get results
res = load(opt.outputfile);

%% get options
do_save = false;
do_random = true;
numStates = [3 4];
sampleSize = [0.05:0.05:1];
runs = 1;

if(nargin>1)        % parse options
    
    k=1;        % argument counter
    kmax=nargin-1;  % stopping criterion
    while(k<=kmax)
        option=varargin{k};
        if(strcmpi(option,'save'))
           do_save = true;
           k=k+1;
        elseif(strcmpi(option,'noRandom'))
           do_random = false;
           k=k+1;
       elseif(strcmpi(option,'runs'))
            if(~isempty(varargin{k+1}))
                runs=varargin{k+1};
                if(~isnumeric(runs) | runs<=0 | runs~=round(runs))
                    error('VB3_varyData: runs option must be followed by a positive integer.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'numStates'))
            if(~isempty(varargin{k+1}))
                numStates=varargin{k+1};
                [m, ~] = size(numStates);
                if(~isnumeric(numStates) | m ~= 1 | numStates<=0 | numStates~=round(numStates))
                    error('VB3_varyData: numStates option must be followed by a positive integer 1*N array.')
                end
            end
            k=k+2;
         elseif(strcmpi(option,'sampleSize'))
            if(~isempty(varargin{k+1}))
                sampleSize=varargin{k+1};
                [m, ~] = size(sampleSize);
                if(~isnumeric(sampleSize) | m ~= 1 | max(sampleSize)>1 | sampleSize<=0)
                    error('VB3_varyData: sampleSize option must be followed by a positive numeric 1*N array, with no element larger than 1.')
                end
            end
            k=k+2;
        else
            error(['VB3_varyData: option ' option ' not recognized.'])
        end
    end
end


%% split up data and run

% If not randomly picking the data points then no need for multiple runs
if ~do_random
    runs = 1;
end

wModels = cell(1, length(sampleSize));
modelRun = cell(1, runs);

temp = matlabpool('size');
if temp~=0
matlabpool close
end

matlabpool open

for ind = 1:runs
    disp(['Run: ' num2str(ind) ' of ' num2str(runs)]);
    
parfor k=1:length(sampleSize)
    % Sample the data
    if do_random
        ind2 = randperm(Ntraj);
        ind2 = sort(ind2(1:ceil(Ntraj*sampleSize(k))));
        Xtemp = {X{ind2}};
    else
            Xtemp = {X{1:ceil(Ntraj*sampleSize(k))}};
    end
    
    for k2=1:length(numStates)
        w = VB3_VBEMiterator(res.WbestN{numStates(k2)},Xtemp,'outputLevel',1,'maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar,'slim');
        w.sampleSize = sampleSize(k);
        w.numStates = numStates(k2);
        wModels{k}{k2} = w;
    end
    
end
modelRun{ind} = wModels;
end

matlabpool close

if do_save
    if do_random
        saveName = [opt.outputfile, '_varDataRandom.mat'];
    else
        saveName = [opt.outputfile, '_varDataNotRandom.mat'];
    end
    if runs==1;
        save(saveName, 'wModels')
    else
        save(saveName, 'modelRun')
    end
end

if runs==1
    modelRun=wModels;
end

disp(['Finished ' opt.runinputfile ' varying amount of data in ' num2str(toc(tvaryData)/60) ' min.'])
end

