
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File for visualize HMM model results in trajectory plots
% FP 120413
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all

%Define variables
start = 5000;
stop =  6000;
minTraj = 7;
discrete = 1;  % 0/1
viterbi = 1;   % 0/1
bgTraj = 1;    % 0/1

X = VB3_readData(options);
Wbest = VB3_VBEMiterator(Wbest, X, 'outputLevel', 2, 'slim', 'estimate');

% oldfold = cd('./120408_Exp Data HFQ');
options.dim = 2;
X = VB3_readData(options);%finalTraj;%
options.dim = 1;
% cd(oldfold);
D_states = Wbest.est.DdtMean;%Dapp;%

hand = figure;
hold on
if stop == 0
    start = 1;
    stop = length(X);
end
if bgTraj
for i = 2000:6000%start:stop
    plot(X{i}(:,1), X{i}(:,2), 'color', [0.8 0.8 0.8], 'Linewidth', 2.5)
end
end
%%
longTrajs = [];
% Loop over trajectories
for i = start:stop
    if(size(X{i},1) > minTraj)
       
        if discrete & viterbi
            viter = double(Wbest.est2.viterbi{i});%X{i}(:, 4);%
            Dvit = D_states(viter);
            plotc(X{i}(:,1), X{i}(:,2), D_states(viter), 'LineWidth',2.5);
        elseif discrete
            [val, ind] = max(Wbest.est2.pst{i}, [], 2);
            [I, J] = ind2sub(size(Wbest.est2.pst{i}), ind);
            I=I+1;
            plotc(X{i}(:,1), X{i}(:,2), double(I), 'LineWidth', 2.5);
        else
            Dav = D_states*Wbest.est2.pst{i}';
            plotc(X{i}(:,1), X{i}(:,2), Dav, 'LineWidth', 2.5);
        end
        
        longTrajs = [longTrajs, i];
        colorbar;
    end
end
plotc([0 0], [0 1], [1 3])
hold off

