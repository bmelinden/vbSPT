
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File for visualize HMM model results in trajectory plots
% FP 120413
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

%Define variables
start = 2000;
stop =  3000;
minTraj = 7;
discrete = 1;  % 0/1
viterbi = 0;   % 0/1
bgTraj = 1;    % 0/1

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
D_states = a.Wbest.est.DdtMean/a.options.dt;

% Read in the trajectory data
a.options.sourcefile = ['/Users/FredrikPersson/Desktop/WORK/UU/Projects/HMM/HMMcode/120408_Exp Data HFQ' '/' a.options.sourcefile];
a.options.dim = 4;
X = VB3_readData(a.options);

hand = figure;
hold on
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
       
        if discrete & viterbi
            viter = double(a.Wbest.est2.viterbi{i});
            Dvit = D_states(viter);
            plotc(X{i}(:,1), X{i}(:,2), D_states(viter), 'LineWidth',2);
        elseif discrete
            [val, ind] = max(a.Wbest.est2.pst{i}, [], 2);
            [I, J] = ind2sub(size(a.Wbest.est2.pst{i}), ind);
            I=I+1;
            plotc(X{i}(:,1), X{i}(:,2), double(I), 'LineWidth', 2);
        else
            Dav = D_states*a.Wbest.est2.pst{i}';
            plotc(X{i}(:,1), X{i}(:,2), Dav, 'LineWidth', 2);
        end
        
        longTrajs = [longTrajs, i];
        colorbar;
    end
end
plotc([0 0], [0 1], [1 3])
hold off

