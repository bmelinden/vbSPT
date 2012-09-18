function [data, trajLengths] = MakeTrajectories(CellL, CellR, diffCoeff, transRate, trajLengths, timeStep, stepSize, locAccuracy, occProb)

%% [data, trajLengths] = MakeTrajectories(CellL, CellR, diffCoeff, transRate, trajLengths, timeStep, stepSize, locAccuracy, occProb)
%
% Starts trajectories in states according to the occupancy probability and calls
% 'simCellHMM' to generate the single trajectories.
%
% F.P. 2012-04-25

%% initiate

data = cell(1, length(trajLengths));


%% Start and generate trajectories
for trajNr = 1:length(trajLengths)
    % Choose starting state for the trajectory
    state = find(rand<=cumsum(occProb),1);
    
    % Pass it on
    [~, Traj, n, ~] = simCell(CellL, CellR, diffCoeff, transRate, trajLengths(trajNr), timeStep, stepSize, locAccuracy, state);
    
    data{trajNr} = Traj;
    trajLengths(trajNr) = n;

end

end
