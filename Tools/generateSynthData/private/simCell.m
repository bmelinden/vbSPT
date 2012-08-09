
function [TimePoints, Traj, m, state] = simCell(L, R, diffCoeff, transMat, N, stepT, state, locAccuracy)

%% [TimePoints, Traj, m, state] = simCell(L, R, diffCoeff, transMat, N, stepT, state, locAccuracy)
% 
% Generates a trajectory within a confined E. coli like geometry. The
% starting state of teh trajectory is given and the initial position is
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
    
    if x>0 & x<L & (z)^2+y^2<R^2
        ini_y = y;
    elseif x<0 & (z)^2+y^2+x^2<R^2
        ini_y = y;
    elseif x>L & (z)^2+y^2+(x-L)^2<R^2
        ini_y = y;
    end
end
while ini_z == 0
    z = rand*R*2-R;
    
    if x>0 & x<L & (z)^2+y^2<R^2
        ini_z = z;
    elseif x<0 & (z)^2+y^2+x^2<R^2
        ini_z = z;
    elseif x>L & (z)^2+y^2+(x-L)^2<R^2
        ini_z = z;
    end
end

%% Define variables
dt = 0;         % Timesteps (will follow a exp distribution)
dx = 5;     % Spatial discretisation of 5 nm
t_old = 0;
t = 0;
n = 1;


pos = [x y z];  % Start position
dpos = [0 0 0]; % Movement


Traj = zeros(N,4);
TimePoints = zeros(N,1);
2*(x+R)/(L+2*R);

cumTransMat = cumsum(transMat, 2);

%% Generate the trajectory
% Loop over all trajectory steps
for n = 1:N
    new_state = find(rand<=cumTransMat(state,:),1);
    state = new_state;
    t;
    while t < (t_old+stepT)
               
        % Update position
        pos = pos+dpos;

        dt = -log(rand)/(6*diffCoeff(state)/dx^2);
        t = t+dt;
        % Generate random direction
        direction = ceil(rand*6);
        
        x = pos(1);
        y = pos(2);
        z = pos(3);
        dpos = [0 0 0];
        switch direction
            case 1
                if x>0 & x<L & (z+dx)^2+y^2<R^2
                    dpos = [0 0 dx];
                elseif x<0 & (z+dx)^2+y^2+x^2<R^2
                    dpos = [0 0 dx];
                elseif x>L & (z+dx)^2+y^2+(x-L)^2<R^2
                    dpos = [0 0 dx];
                end
            case 2
                if x>0 & x<L & (z-dx)^2+y^2<R^2
                    dpos = [0 0 -dx];
                elseif x<0 & (z-dx)^2+y^2+x^2<R^2
                    dpos = [0 0 -dx];
                elseif x>L & (z-dx)^2+y^2+(x-L)^2<R^2
                    dpos = [0 0 -dx];
                end
            case 3
                if x>0 & x<L & (z)^2+(y+dx)^2<R^2
                    dpos = [ 0 dx 0];
                elseif x<0 & (z)^2+(y+dx)^2+x^2<R^2
                    dpos = [ 0 dx 0];
                elseif x>L & (z)^2+(y+dx)^2+(x-L)^2<R^2
                    dpos = [ 0 dx 0];
                end
            case 4
                if x>0 & x<L & (z)^2+(y-dx)^2<R^2
                    dpos = [ 0 -dx 0];
                elseif x<0 & (z)^2+(y-dx)^2+x^2<R^2
                    dpos = [ 0 -dx 0];
                elseif x>L & (z)^2+(y-dx)^2+(x-L)^2<R^2
                    dpos = [ 0 -dx 0];
                end
            case 5
                if (x+dx)>0 & (x+dx)<L & (z)^2+y^2<R^2
                    dpos = [dx 0 0];
                elseif (x+dx)<0 & (z)^2+y^2+(x+dx)^2<R^2
                    dpos = [dx 0 0];
                elseif (x+dx)>L & (z)^2+y^2+((x+dx)-L)^2<R^2
                    dpos = [dx 0 0];
                end
            case 6
                if (x-dx)>0 & (x-dx)<L & (z)^2+y^2<R^2
                    dpos = [-dx 0 0];
                elseif (x-dx)<0 & (z)^2+y^2+(x-dx)^2<R^2
                    dpos = [-dx 0 0];
                elseif (x-dx)>L & (z)^2+y^2+((x-dx)-L)^2<R^2
                    dpos = [-dx 0 0];
                end
        end
    end
    Traj(n, :) = [pos, state];
    
    
    TimePoints(n) = t;
    m = n;
    t_old = t;
end

%% Add localization limitations
for j = 1:m
        newpos = Traj(j,1:3)+locAccuracy.*randn(1,3); 
        Traj(j, 1:3) = newpos;
end

end



