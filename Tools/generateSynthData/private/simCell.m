
function [TimePoints, Traj, m, state] = simCell(L, R, diffCoeff, transRate, N, stepT, state, locAccuracy)

%% [TimePoints, Traj, m, state] = simCell(L, R, diffCoeff, transRate, N, stepT, state, locAccuracy)
% 
% Generates a trajectory within a confined E. coli like geometry. The
% starting state of the trajectory is given and the initial position is
% randomly positioned within the geometry.

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
dl = 5;     % Spatial discretisation of 5 nm
t_old = 0;
t = 0;
m = 1;


pos = [x y z];  % Start position
dpos = [0 0 0]; % Movement


Traj = zeros(N,4);
TimePoints = zeros(N,1);


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
        if (x+dl)>0 && (x+dl)<L && (z)^2+y^2<R^2        % In the pipe section
            rate_Dx1 = diffCoeff(state)/dl^2;
        elseif (x+dl)<0 && (z)^2+y^2+(x+dl)^2<R^2       % in the left end cap
            rate_Dx1 = diffCoeff(state)/dl^2;
        elseif (x+dl)>L && (z)^2+y^2+((x+dl)-L)^2<R^2   % in the right end cap
            rate_Dx1 = diffCoeff(state)/dl^2;
        else
            rate_Dx1 = 0;
        end
        % Negative x direction
        if (x-dl)>0 && (x-dl)<L && (z)^2+y^2<R^2        % In the pipe section
            rate_Dx2 = diffCoeff(state)/dl^2;
        elseif (x-dl)<0 && (z)^2+y^2+(x-dl)^2<R^2       % in the left end cap
            rate_Dx2 = diffCoeff(state)/dl^2;
        elseif (x-dl)>L && (z)^2+y^2+((x-dl)-L)^2<R^2   % in the right end cap
            rate_Dx2 = diffCoeff(state)/dl^2;
        else
            rate_Dx2 = 0;
        end
        
        % Positive y direction
        if x>0 && x<L && (z)^2+(y+dl)^2<R^2
            rate_Dy1 = diffCoeff(state)/dl^2;
        elseif x<0 && (z)^2+(y+dl)^2+x^2<R^2
            rate_Dy1 = diffCoeff(state)/dl^2;
        elseif x>L && (z)^2+(y+dl)^2+(x-L)^2<R^2
            rate_Dy1 = diffCoeff(state)/dl^2;
        else
            rate_Dy1 = 0;
        end
        % Negative y direction
        if x>0 && x<L && (z)^2+(y-dl)^2<R^2
            rate_Dy2 = diffCoeff(state)/dl^2;
        elseif x<0 && (z)^2+(y-dl)^2+x^2<R^2
            rate_Dy2 = diffCoeff(state)/dl^2;
        elseif x>L && (z)^2+(y-dl)^2+(x-L)^2<R^2
            rate_Dy2 = diffCoeff(state)/dl^2;
        else
            rate_Dy2 = 0;
        end
        
        % Positive z direction
        if x>0 && x<L && (z+dl)^2+y^2<R^2
            rate_Dz1 = diffCoeff(state)/dl^2;
        elseif x<0 && (z+dl)^2+y^2+x^2<R^2
            rate_Dz1 = diffCoeff(state)/dl^2;
        elseif x>L && (z+dl)^2+y^2+(x-L)^2<R^2
            rate_Dz1 = diffCoeff(state)/dl^2;
        else
            rate_Dz1 = 0;
        end
        % Negative z direction
        if x>0 && x<L && (z-dl)^2+y^2<R^2
            rate_Dz2 = diffCoeff(state)/dl^2;
        elseif x<0 && (z-dl)^2+y^2+x^2<R^2
            rate_Dz2 = diffCoeff(state)/dl^2;
        elseif x>L && (z-dl)^2+y^2+(x-L)^2<R^2
            rate_Dz2 = diffCoeff(state)/dl^2;
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
            otherwise
                new_state = find(rand<=(cumsum(rate_Trans./(sum(rate_Trans)))), 1);
                state = new_state;
%                  disp(['State changed at time ' num2str(t)]);
        end
        
        % update position
        pos = pos+dpos;
        
        
    end
    
    %% Read out parameters
    
    Traj(n, :) = [pos, state];
    
    TimePoints(n) = t_old+stepT;
    m = n;
    t_old = t_old+stepT;
end

%% Add localization limitations
for j = 1:m
    newpos = Traj(j,1:3)+locAccuracy.*randn(1,3);
    Traj(j, 1:3) = newpos;
end

end



