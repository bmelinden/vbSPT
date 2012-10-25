
function [TimePoints, Traj, m, state] = simCell(L, R, diffCoeff, transRate, N, stepT, stepSize, locAccuracy, state)

%% [TimePoints, Traj, m, state] = simCell(L, R, diffCoeff, transRate, N, stepT, stepSize, locAccuracy, state)
% 
% Generates a trajectory within a confined E. coli like geometry. The
% starting state of the trajectory is given and the initial position is
% randomly positioned within the geometry.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_simCell, simulate a trajectory in a cell geometry, in the vbSPT package
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

%% Start of the actual code
% remove diagonal entries in transRate matrix
transRate=transRate-diag(diag(transRate));
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
dt = 0;         % Timesteps (will follow a exp distribution)
stepSize = stepSize;     % Spatial discretisation of stepSize nm
t_old = 0;
t = 0;
m = 1;


pos = [x y z];  % Start position
dpos = [0 0 0]; % Movement


Traj = zeros(N,4);
TimePoints = zeros(N,1);

transRate=transRate-diag(diag(transRate)); % remove negative diagonal elements, just to be sure

%% Generate the trajectory
% Loop over all trajectory steps
for n = 1:N
    
    while t < (t_old+stepT)
        
        %% Calculate rates and waiting times
        % Diffusion
        x = pos(1);
        y = pos(2);
        z = pos(3);
        dpos = [0 0 0];
        
        % Positive x direction
        if (x+stepSize)>0 && (x+stepSize)<L && (z)^2+y^2<R^2        % In the pipe section
            rate_Dx1 = diffCoeff(state)/stepSize^2;
        elseif (x+stepSize)<0 && (z)^2+y^2+(x+stepSize)^2<R^2       % in the left end cap
            rate_Dx1 = diffCoeff(state)/stepSize^2;
        elseif (x+stepSize)>L && (z)^2+y^2+((x+stepSize)-L)^2<R^2   % in the right end cap
            rate_Dx1 = diffCoeff(state)/stepSize^2;
        else
            rate_Dx1 = 0;
        end
        % Negative x direction
        if (x-stepSize)>0 && (x-stepSize)<L && (z)^2+y^2<R^2        % In the pipe section
            rate_Dx2 = diffCoeff(state)/stepSize^2;
        elseif (x-stepSize)<0 && (z)^2+y^2+(x-stepSize)^2<R^2       % in the left end cap
            rate_Dx2 = diffCoeff(state)/stepSize^2;
        elseif (x-stepSize)>L && (z)^2+y^2+((x-stepSize)-L)^2<R^2   % in the right end cap
            rate_Dx2 = diffCoeff(state)/stepSize^2;
        else
            rate_Dx2 = 0;
        end
        
        % Positive y direction
        if x>0 && x<L && (z)^2+(y+stepSize)^2<R^2
            rate_Dy1 = diffCoeff(state)/stepSize^2;
        elseif x<0 && (z)^2+(y+stepSize)^2+x^2<R^2
            rate_Dy1 = diffCoeff(state)/stepSize^2;
        elseif x>L && (z)^2+(y+stepSize)^2+(x-L)^2<R^2
            rate_Dy1 = diffCoeff(state)/stepSize^2;
        else
            rate_Dy1 = 0;
        end
        % Negative y direction
        if x>0 && x<L && (z)^2+(y-stepSize)^2<R^2
            rate_Dy2 = diffCoeff(state)/stepSize^2;
        elseif x<0 && (z)^2+(y-stepSize)^2+x^2<R^2
            rate_Dy2 = diffCoeff(state)/stepSize^2;
        elseif x>L && (z)^2+(y-stepSize)^2+(x-L)^2<R^2
            rate_Dy2 = diffCoeff(state)/stepSize^2;
        else
            rate_Dy2 = 0;
        end
        
        % Positive z direction
        if x>0 && x<L && (z+stepSize)^2+y^2<R^2
            rate_Dz1 = diffCoeff(state)/stepSize^2;
        elseif x<0 && (z+stepSize)^2+y^2+x^2<R^2
            rate_Dz1 = diffCoeff(state)/stepSize^2;
        elseif x>L && (z+stepSize)^2+y^2+(x-L)^2<R^2
            rate_Dz1 = diffCoeff(state)/stepSize^2;
        else
            rate_Dz1 = 0;
        end
        % Negative z direction
        if x>0 && x<L && (z-stepSize)^2+y^2<R^2
            rate_Dz2 = diffCoeff(state)/stepSize^2;
        elseif x<0 && (z-stepSize)^2+y^2+x^2<R^2
            rate_Dz2 = diffCoeff(state)/stepSize^2;
        elseif x>L && (z-stepSize)^2+y^2+(x-L)^2<R^2
            rate_Dz2 = diffCoeff(state)/stepSize^2;
        else
            rate_Dz2 = 0;
        end
        
        % Diffusion rates
        rate_D = [rate_Dx1, rate_Dx2, rate_Dy1, rate_Dy2, rate_Dz1, rate_Dz2];
        % Transitions between states
        rate_Trans = transRate(state, :);
        % All rates
        rates = [rate_D, rate_Trans];
        rate_tot = sum(rates);
        
        
        dt = -log(rand)/rate_tot;
        t = t+dt;
        
        %% Perform action
        action = find(rand<=cumsum(rates)/rate_tot, 1);
        
        if(sum(rate_D)==0)
            keyboard
        end
        
        switch action
            case 1
                dpos(1) = stepSize;
            case 2
                dpos(1) = -stepSize;
            case 3
                dpos(2) = stepSize;
            case 4
                dpos(2) = -stepSize;
            case 5
                dpos(3) = stepSize;
            case 6
                dpos(3) = -stepSize;
            otherwise
                new_state = find(rand<=(cumsum(rate_Trans./(sum(rate_Trans)))), 1);
                state = new_state;
%                  disp(['State changed at time ' num2str(t)]);
        end
        
        % update position
        pos = pos+dpos;
        
        
    end
    
    %% Read out parameters
    % Add localization error and possibly position error due to discretization
    % Positioning error not included since it would change the app Diff a
    % bit.
    
    Traj(n, :) = [pos+locAccuracy.*randn(1,3), state];
    
    TimePoints(n) = t_old+stepT;
    m = n;
    t_old = t_old+stepT;
end



end
