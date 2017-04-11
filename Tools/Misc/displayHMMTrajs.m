function displayHMMTrajs()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to visualize HMM model results in trajectory plots
% Place this script and the plotc.m file in the experiment folder (where the 
% runinputfiles are kept) and then modify parameters and choose run. 
% FP 120413
% ML 170411: script -> function, fix stateTrj length discrepancy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define variables
start = 0;      
stop =  0;     % stop=0 -> plot all trajectories
minTraj = 1;
trjColoring='viterbi';%'viterbi', 'sMaxP',, or '<D(t)>'

bgTraj = 0;    % 0/1
trueSimData = 0; % 0/1

lineWidth=1; % with of trajectory lines

% Get filename and path with "uigetfile"
[filename, pathname] = uigetfile({'*.mat'}, 'Select HMM .mat file');
if ( filename == 0 )
    disp('Error! No (or wrong) file selected!')
    return
end

% Load the mat file
full_filename = [ pathname, filename ];
a = load(full_filename);

% Read in D for the found states
D_states = a.Wbest.est.DdtMean/a.options.timestep;

% Read in the trajectory data
a.options.dim = 2;
X = VB3_readData(a.options);

hand = figure('Name',['HMM segmentation: ' trjColoring ]);
hold on;
axis equal
box on

if stop == 0
    start = 1;
    stop = length(X);
end
if bgTraj
    for i = start:stop
        plot(X{i}(:,1), X{i}(:,2), 'color', [0.8 0.8 0.8])
    end
end

longTrajs = [];
% Loop over trajectories
for i = start:stop
    if(size(X{i},1) > minTraj)
        % Note: state trajectories are always one shorter that position
        % trajectories, and needs extension to plot, but the last color
        % value is not actually used.
        switch lower(trjColoring)
            case 'viterbi'
                viter = double(a.Wbest.est2.viterbi{i}([1:end end]));
                plotc(X{i}(:,1), X{i}(:,2), D_states(viter), 'LineWidth',lineWidth);
            case 'smaxp'
                sMaxP=double(a.Wbest.est2.sMaxP{i}([1:end end]));
                plotc(X{i}(:,1), X{i}(:,2), double(sMaxP), 'LineWidth',lineWidth);
            case '<d(t)>'
                Dav = D_states*(a.Wbest.est2.pst{i}([1:end end],:))';
                plotc(X{i}(:,1), X{i}(:,2), Dav, 'LineWidth',lineWidth);
            otherwise
                error(['trjColoring : ' trjColoring ' not recognized'])
        end

        longTrajs = [longTrajs, i];
        colorbar;
    end
end
plotc([0 0], [0 1], [1 3])
hold off

if trueSimData
    hand2 = figure('Name','Simulated data');
    hold on
    if bgTraj
        for i = start:stop
            plot(X{i}(:,1), X{i}(:,2), 'color', [0.8 0.8 0.8])
        end
    end
    
    longTrajs = [];
    % Loop over trajectories
    for i = start:stop
        if(size(X{i},1) > minTraj)
            plotc(X{i}(:,1), X{i}(:,2), D_states(X{i}(:,4)), 'LineWidth',lineWidth);
            
            longTrajs = [longTrajs, i];
            colorbar;
        end
    end
    plotc([0 0], [0 1], [1 3])
    hold off
end


