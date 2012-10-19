function sl = getSL(traj, varargin)
%% function sl = getSL(traj, varargin)
%
% For getting out the steplengths in X or XY for all the steps in 'traj' which
% can be a cell array of trajectories.
%
% Input:
% traj          : A cell array of trajectory arrays, where each line corresponds to a
%               position. Number of columns is arbitrary. 
% 
% Options:
% numDim        : An integer determine if to take the steplengths in only X
%               or XY. Default is 1, first column in the trajectory array
%               (often X).


%% Initiate options
numDim = 1;

%% Read options
if(nargin>1)        % parse options
    k = 1;
    kmax=nargin;  % stopping criterion
    while(k<kmax)
        option=varargin{k};
        if(strcmpi(option,'numDim'))
            if(~isempty(varargin{k+1}))
                numDim=varargin{k+1};
                if(~isnumeric(numDim) || numDim<=0 || numDim~=round(numDim) || numDim>2)
                    error('getSL: numDim option must be followed by a positive integer less than 3.')
                end
            end
            k=k+2;
        else
            error(['getSL: option ' option ' not recognized.'])
        end
    end
end

sl = [];

for ind = 1:length(traj)
    temp = traj{ind};
    
    if numDim == 1
        sl = [sl; temp(2:end, 1)-temp(1:end-1, 1)];
    elseif numDim == 2
        dCoords = temp(2:end, 1:2) - temp(1:end-1, 1:2);
        sqDispl = sum(dCoords.^2,2); % dx^2+dy^2
        sl = [sl; sqrt(sqDispl)];
    else
        error(['getSL only supports either 1 (X) or 2 dimensions (XY).'])
    end
    
end

sl = abs(sl);

figure; 
hist(sl, 100); 
title(['Histogram over steplengths in ' num2str(numDim) 'D']);
xlabel('Steplength [nm]');
ylabel('Counts');
end