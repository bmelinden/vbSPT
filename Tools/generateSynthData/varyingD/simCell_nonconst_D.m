
function [TimePoints, Traj, m] = simCell_nonconst_D(L, R, Dfun,  N, stepT,locAccuracy)
%% [TimePoints, Traj, m] = simCell_nonconst_D(L, R, Dfun,  N, stepT,locAccuracy)
% 
% Generates a trajectory within a confined E. coli like geometry. The
% starting state of the trajectory is given and the initial position is
% randomly positioned within the geometry.
%
% modified to handle space-dependent diffusion constant:
% L,R   standard cell length and radius (total cell length is L+2R)
% Dfun  function handle to diffusion constant, e.g.,
%       @(x,y,z)(50+(x-L/2)*0.01 for a linear gradient
%
% Traj(n,:) = [n y z Dfun(x,y,z) <D>], 
%           where <D> is the average diffusion constant between rwo n-1
%           and n.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simCell_nonconst_D, simulate difusion with spatially changing D.
% =========================================================================
% 
% Copyright (C) 2012 Martin Lind??n and Fredrik Persson
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


%% initiate with random positions
% Generate a random x-position
x = rand*(L+2*R)-R;

% Generate random y and z-position
y = 0;
ini_y = 0;
z = 0;
ini_z = 0;
while ini_y == 0
    y = rand*R*2-R;
    
    if x>0 && x<L && (z)^2+y^2<R^2
        ini_y = y;
    elseif x<0 && (z)^2+y^2+x^2<R^2
        ini_y = y;
    elseif x>L && (z)^2+y^2+(x-L)^2<R^2
        ini_y = y;
    end
end
while ini_z == 0
    z = rand*R*2-R;
    
    if x>0 && x<L && (z)^2+y^2<R^2
        ini_z = z;
    elseif x<0 && (z)^2+y^2+x^2<R^2
        ini_z = z;
    elseif x>L && (z)^2+y^2+(x-L)^2<R^2
        ini_z = z;
    end
end

%% Define variables
dt = 0;     % Timesteps (will follow a exp distribution)
dl = 5;     % Spatial discretisation of 5 nm
t_old = 0;
t = 0;
m = 1;

pos = [x y z];  % Start position
dpos = [0 0 0]; % Movement


Traj = zeros(N,5);
TimePoints = zeros(N,1);


%% Generate the trajectory
% Loop over all trajectory steps
Ddt_sum=0;
tfrac=0;dTfrac=0.05;
if(N<50)
    tfrac=2; % only write out progress for fairly long trajectories
end
for n = 1:N
    
    while t < (t_old+stepT)
        
        %% Calculate rates and waiting times
        % Diffusion
        x = pos(1);
        y = pos(2);
        z = pos(3);
        dpos = [0 0 0];
        
        % Diffusion rates
        pos_Dx=x+0.5*dl*[+1 -1 0 0 0 0];
        pos_Dy=y+0.5*dl*[0 0 +1 -1 0 0];
        pos_Dz=z+0.5*dl*[0 0 0 0 +1 -1];
        
        rate_D=zeros(1,6);
        for kd=1:6
            rate_D(kd)=Dfun(pos_Dx(kd),pos_Dy(kd),pos_Dz(kd))/dl^2;
        end
        
        % next, remove jumps out of the cell
        % Positive x direction
        if (x+dl)>0 && (x+dl)<L && (z)^2+y^2<R^2        % In the pipe section
        elseif (x+dl)<0 && (z)^2+y^2+(x+dl)^2<R^2       % in the left end cap
        elseif (x+dl)>L && (z)^2+y^2+((x+dl)-L)^2<R^2   % in the right end cap
        else
            rate_D(1) = 0;
        end
        % Negative x direction
        if (x-dl)>0 && (x-dl)<L && (z)^2+y^2<R^2        % In the pipe section
        elseif (x-dl)<0 && (z)^2+y^2+(x-dl)^2<R^2       % in the left end cap
        elseif (x-dl)>L && (z)^2+y^2+((x-dl)-L)^2<R^2   % in the right end cap
        else
            rate_D(2) = 0;
        end
        
        % Positive y direction
        if x>0 && x<L && (z)^2+(y+dl)^2<R^2
        elseif x<0 && (z)^2+(y+dl)^2+x^2<R^2
        elseif x>L && (z)^2+(y+dl)^2+(x-L)^2<R^2
        else
            rate_D(3) = 0;
        end
        % Negative y direction
        if x>0 && x<L && (z)^2+(y-dl)^2<R^2
        elseif x<0 && (z)^2+(y-dl)^2+x^2<R^2
        elseif x>L && (z)^2+(y-dl)^2+(x-L)^2<R^2
        else
            rate_D(4) = 0;
        end
        
        % Positive z direction
        if x>0 && x<L && (z+dl)^2+y^2<R^2
        elseif x<0 && (z+dl)^2+y^2+x^2<R^2
        elseif x>L && (z+dl)^2+y^2+(x-L)^2<R^2
        else
            rate_D(5) = 0;
        end
        % Negative z direction
        if x>0 && x<L && (z-dl)^2+y^2<R^2
        elseif x<0 && (z-dl)^2+y^2+x^2<R^2
        elseif x>L && (z-dl)^2+y^2+(x-L)^2<R^2
        else
            rate_D(6) = 0;
        end
        
        % All rates
        rate_tot = sum(rate_D);
        
        
        dt = -log(rand)/rate_tot;
        t = t+dt;
        
        %% Perform action
        action = find(rand<=cumsum(rate_D)/rate_tot, 1);
        
        switch action
            case 1
                dpos(1) = dl;
            case 2
                dpos(1) = -dl;
            case 3
                dpos(2) = dl;
            case 4
                dpos(2) = -dl;
            case 5
                dpos(3) = dl;
            case 6
                dpos(3) = -dl;
        end
        
        % update position
        pos = pos+dpos;
        Ddt_sum=Ddt_sum+rate_D(action)*dl^2*dt;
    end
    
    Ddt_sum=Ddt_sum-rate_D(action)*dl^2*(t-t_old-stepT); % remove overhanging time average
    meanD=Ddt_sum/stepT;                                 % time-averaged diffusion constant
    Ddt_sum=rate_D(action)*dl^2*(t-t_old-stepT);         % reinsert overhanging time average 
    
    %% Read out parameters
    Traj(n, :) = [pos Dfun(pos(1),pos(2),pos(3)) 0];
    if(n>1)
       Traj(n-1,5)=meanD; % so, Traj(n,5) is the average diffusion constant between Traj(n,:) and Traj(n+1,:)
    end
    
    TimePoints(n) = t_old+stepT;
    m = n;
    t_old = t_old+stepT;

    if(n/N-tfrac>dTfrac)
        disp(['finished : ' num2str(100*n/N) ' %'])
        tfrac=n/N;
    end
   
end

%% Add localization limitations
for j = 1:m
    newpos = Traj(j,1:3)+locAccuracy.*randn(1,3);
    Traj(j, 1:3) = newpos;
end

end



