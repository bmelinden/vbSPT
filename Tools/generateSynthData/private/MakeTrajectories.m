function [data, trajLengths] = MakeTrajectories(timeStep, CellL, CellR, trajLengths, diffCoeff, transMat, occProb, locAccuracy)

%% [data, trajLengths] = MakeTrajectories(timeStep, CellL, CellR, trajLengths, diffCoeff, transMat, occProb, locAccuracy)
%
% Starts trajectories in states according to the occupancy probability and calls
% 'simCellHMM' to generate the single trajectories.
%
% F.P. 2012-04-25
%  
% Change log:
% 

%% initiate

data = cell(1, length(trajLengths));

%% Start and generate trajectories
parfor j = 1:length(trajLengths)
    % Choose starting state for the trajectory
    state = find(rand<=cumsum(occProb),1);
    
    % Pass it on
    [TimePoints, Traj, n, state] = simCell(CellL, CellR, diffCoeff, transMat, trajLengths(j), timeStep, state, locAccuracy);
    j

    data{j} = Traj;
    trajLengths(j) = n;

end

end
