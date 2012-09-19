function [finalTraj, trajL] = chopTrajExp(traj, trajLmean)

% For splitting one long simultaed trajectory into an exponential distribution 
% of shorter trajectories.
%
% Input:
%
% traj          : A trajectory array where each line corresponds to a
%               position. Number of columns is arbitrary. 
% trajLmean     : The mean of the exponential distribution that the 
%               trajectories should be split into.
%
%


finalTraj = cell(1, 2*length(traj)/trajLmean);
trajL = zeros(size(finalTraj));

ind = 1;
tind = 1;
tind2 = round(exprnd(trajLmean));
if tind2 <= 2
    tind2 = 2;
end

while tind2<=length(traj)+1
    
    finalTraj{ind} = traj(tind:tind2-1, :);
    trajL(ind) = tind2-tind;
    tind = tind2;
    tind2 = tind+round(exprnd(trajLmean));
    ind = ind+1;
end

% remove empty cells
finalTraj = finalTraj(~cellfun('isempty', finalTraj));
trajL(trajL==0) = [];

end