function [finalTraj, trajL] = chopTrajExp_constStep(traj,trajLmean,trajLmin)
% [finalTraj, trajL] = chopTrajExp_constStep(traj,trajLmean,trajLmin)
%
% For splitting one long simulated trajectory into an exponential distribution 
% of shorter trajectories. All steps are included, which means that the end
% point of one trajectory is the starting point of the next one.
%
% Input:
%
% traj          : A trajectory cell array where each line corresponds to a
%               position. Number of columns is arbitrary. 
% trajLmean     : The mean of the exponential distribution that the 
%               trajectories should be split into.
%
%


finalTraj = cell(1, 2*round(length(traj)/trajLmean));
trajL = zeros(size(finalTraj));

% generate trajectory lengths that sum to the right length
trajL=[];


i1=1;i2=1;
%i1=[1 cumsum(trajL(1:9)-1)'+1];
%i2=[cumsum(trajL(1:10)-1)'+1];
if(trajLmin==1)
    error('minimum trj length must exceed 1')
end
kk=0;
while(i2<=size(traj,1) && kk<size(traj,1))
    dL=round(exprnd(trajLmean)); % number of positions in this trajectory
    if(dL>=trajLmin)
        kk=kk+1;
        i1=i2;
        i2=i2+dL-1;
        if(i2<=size(traj,1))
            finalTraj{kk}=traj(i1:i2,:);
            trajL(kk)=dL;
        elseif(size(traj,1)-i1>=trajLmin)
            finalTraj{kk}=traj(i1:end,:);
            trajL(kk)=size(traj,1)-i1;
        end
    end
end
kk=kk-1;
finalTraj={finalTraj{1:kk}};
