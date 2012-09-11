
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
        rate_D = (6*diffCoeff(state)/dl^2);
        % Transitions between states
        rate_Trans = sum(transRate(state, :));
        % Additional misc rates
        rate_Misc = 0;
        rates = [rate_D, rate_Trans, rate_Misc];
        rate_tot = sum(rates);
        
        dt = -log(rand)/rate_tot;
        t = t+dt;
        
        %% Perform action
        
        action = find(rand<=cumsum(rates)/rate_tot, 1);
        
        switch action
            case 1  
                %% Diffusion
                
                % Generate random direction
                direction = ceil(rand*6);
                x = pos(1);
                y = pos(2);
                z = pos(3);
                dpos = [0 0 0];
                switch direction
                    case 1
                        if x>0 && x<L && (z+dl)^2+y^2<R^2
                            dpos = [0 0 dl];
                        elseif x<0 && (z+dl)^2+y^2+x^2<R^2
                            dpos = [0 0 dl];
                        elseif x>L && (z+dl)^2+y^2+(x-L)^2<R^2
                            dpos = [0 0 dl];
                        end
                    case 2
                        if x>0 && x<L && (z-dl)^2+y^2<R^2
                            dpos = [0 0 -dl];
                        elseif x<0 && (z-dl)^2+y^2+x^2<R^2
                            dpos = [0 0 -dl];
                        elseif x>L && (z-dl)^2+y^2+(x-L)^2<R^2
                            dpos = [0 0 -dl];
                        end
                    case 3
                        if x>0 && x<L && (z)^2+(y+dl)^2<R^2
                            dpos = [ 0 dl 0];
                        elseif x<0 && (z)^2+(y+dl)^2+x^2<R^2
                            dpos = [ 0 dl 0];
                        elseif x>L && (z)^2+(y+dl)^2+(x-L)^2<R^2
                            dpos = [ 0 dl 0];
                        end
                    case 4
                        if x>0 && x<L && (z)^2+(y-dl)^2<R^2
                            dpos = [ 0 -dl 0];
                        elseif x<0 && (z)^2+(y-dl)^2+x^2<R^2
                            dpos = [ 0 -dl 0];
                        elseif x>L && (z)^2+(y-dl)^2+(x-L)^2<R^2
                            dpos = [ 0 -dl 0];
                        end
                    case 5
                        if (x+dl)>0 && (x+dl)<L && (z)^2+y^2<R^2
                            dpos = [dl 0 0];
                        elseif (x+dl)<0 && (z)^2+y^2+(x+dl)^2<R^2
                            dpos = [dl 0 0];
                        elseif (x+dl)>L && (z)^2+y^2+((x+dl)-L)^2<R^2
                            dpos = [dl 0 0];
                        end
                    case 6
                        if (x-dl)>0 && (x-dl)<L && (z)^2+y^2<R^2
                            dpos = [-dl 0 0];
                        elseif (x-dl)<0 && (z)^2+y^2+(x-dl)^2<R^2
                            dpos = [-dl 0 0];
                        elseif (x-dl)>L && (z)^2+y^2+((x-dl)-L)^2<R^2
                            dpos = [-dl 0 0];
                        end
                end
                
                % Update position
                pos = pos+dpos;
                
            case 2
                %% Transition
                
                new_state = find(rand<=(cumsum(transRate(state,:))./(sum(transRate(state,:)))), 1);
                state = new_state;
                
            case 3
                %% Misc
                
                % Nothing so far
                
        end
        
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



