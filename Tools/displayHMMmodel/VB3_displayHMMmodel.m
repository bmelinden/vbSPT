function VB3_displayHMMmodel(runinput)
% VB3_displayHMMmodel(runinput)
%
% File for getting out and graphically representing HMM model results from
% vbSPT. It can take either a runinputfile or options struct. If no
% argument is given then you get to choose the result file you want from an
% 'Open file' dialogue. It shows the result for the globally best model,
% however in the code you can set what modelsize should be shown.
%
% Updated to 1.1 model format 2013-11-07 ML.


%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_HMManalysis, runs data analysis in the vbSPT package
% =========================================================================
% 
% Copyright (C) 2013 Martin Lindén and Fredrik Persson
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

VB3_license('VB3_displayHMMmodel')


%% Options
modelSize = 0;



%% Load input

if(nargin==0)
%     clear all
    
    % Get filename and path with "uigetfile"
    [filename, pathname] = uigetfile({'*.mat'}, 'Select mat file');
    if ( filename == 0 )
        disp('Error! No (or wrong) file selected!')
        return
    end
    
    % Load the mat file
    full_filename = [ pathname, filename ];
    load(full_filename);
else
    % if an existing file, generate options structure
    if(isstr(runinput) && exist(runinput)==2)
        runinputfile = runinput;
        options=VB3_getOptions(runinputfile);
        disp(['Read runinput file ' runinputfile])
        % if an option struct, read in the runinputfilename
    elseif(isstruct(runinput))
        options=runinput;
        runinputfile=options.runinputfile;
        disp(['Read options structure based on runinput file ' runinputfile ])
    else
        error(['VB3_displayHMMmodel: input ' runinput ' not recognized.']);
    end
    
    filename = options.outputfile;
    load(filename);
end

%% Make backwards compatible
% temporary translations from old (pre 2012-06-11) parameter names
old={'Nmax','dt','pathestimate','BSnum','savefile','sourcefile'};
new={'maxHidden','timestep','stateEstimate','bootstrapNum','outputfile','inputfile'};
for kk=1:length(old)
    if(isfield(options,old{kk}))
        options.(new{kk})=options.(old{kk});
        options=rmfield(options,old{kk});
        warning(['This file uses an outdated parameter name: ' old{kk} '. Use ' new{kk} ' instead.'])
    end
end




%% Read in parameters given prior to HMM
timeStep = options.timestep        
dim = options.dim
if options.init_D(1)>100
    LscaleFactor = 1e6;
else
    LscaleFactor = 1;
end

%% Read in trajectory lengths
trajL = Wbest.T;
avTrajL = mean(trajL);

%% Read in HMM results
if modelSize == 0;
    numStates = Wbest.N     % Number of states
    
    [diffCoeff, ind] = sort(Wbest.est.DdtMean./timeStep/LscaleFactor); % Mean diff. coeff in um^2/s provided the data is in nm
    diffCoeffStd = Wbest.est.Ddtstd(ind)./timeStep/LscaleFactor;
    
    occTot = Wbest.est.Ptot(ind);       % Total occupation percentage
    dwellTime = Wbest.est.dwellMean(ind)'*timeStep;     % Mean dwelltime in each state

    % extract pure transition counts from 1.1 model
    wA=diag(Wbest.M.wa(:,2)-Wbest.PM.wa(:,2))+ Wbest.M.wB - Wbest.PM.wB;
    A=rowNormalize(wA);    
    %A = Wbest.M.wA - Wbest.PM.wA;   %The transition probability matrix with the prior values subtracted
    %A = spdiags (sum (A,2), 0, numStates, numStates) \ A; % Rownormalize the matrix
    transitions = A(ind, ind);
    
elseif exist('WbestN', 'var') && modelSize<=length(WbestN)
    numStates = WbestN{modelSize}.N;
    
    [diffCoeff, ind] = sort(WbestN{modelSize}.est.DdtMean./timeStep/LscaleFactor); % Mean diff. coeff in um^2/s provided the data is in nm
    diffCoeffStd = WbestN{modelSize}.est.Ddtstd(ind)./timeStep/LscaleFactor;
    
    occTot = WbestN{modelSize}.est.Ptot(ind);       % Total occupation percentage
    dwellTime = WbestN{modelSize}.est.dwellMean(ind)'*timeStep;     % Mean dwelltime in each state
    
    % extract pure transition counts from 1.1 model
    wA=diag(WbestN{modelSize}.M.wa(:,2)-WbestN{modelSize}.PM.wa(:,2)) ...
        + (WbestN{modelSize}.M.wB - WbestN{modelSize}.PM.wB);
    A=rowNormalize(wA);
    %A = WbestN{modelSize}.M.wA - WbestN{modelSize}.PM.wA;   %The transition probability matrix with the prior values subtracted
    %A = spdiags (sum (A,2), 0, numStates, numStates) \ A; % Rownormalize the matrix
    transitions = A(ind, ind);
else
    error(['VB3_displayHMMmodel: Could not find the best model of size ' modelSize '.']);
end

%% Print all the numbers in the Matlab prompt
disp('Num. Traj.:')
disp(num2str(length(trajL)));

disp('Av. Traj. length:')
disp(num2str(avTrajL));

disp('Diffusion Coeff [um^2/s]:')
disp(num2str(diffCoeff));
disp(num2str(diffCoeffStd));

disp('Total Occupation:')
disp(num2str(occTot)); 

disp('Dwell time [s]:')
disp(num2str(dwellTime)); 

disp('Transition probability per timestep:')
disp(num2str(transitions)); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for graphical representation

separation = 0.3;
textsep = 0.2;
textindent = 0.1;
precision = 3;

% Graphical representation for models with 2-4 states
figure;
hold on
text(2, 1.5, strcat('Num. Traj: ', num2str(length(trajL))));
text(2, 1.2, strcat('Av. Traj. length: ', num2str(avTrajL)));

if numStates == 2

    minX = 2; maxX = 8; minY = 2; maxY = 6;
    
    stateCoord1 = [4 4];
    stateCoord2 = [6 4];
    
    circle(stateCoord1, occTot(1), 20);
    text(stateCoord1(1), stateCoord1(2)+1, strcat('D1 = ', num2str(diffCoeff(1), precision)));
    circle(stateCoord2, occTot(2), 20);
    text(stateCoord2(1), stateCoord2(2)+1, strcat('D2 = ', num2str(diffCoeff(2), precision)));
    
    
    warning off;
    
    if transitions(1, 2) > 10^(-5)
    arrow(stateCoord1, stateCoord2, 'wid', round(transitions(1, 2)*100));
    end
    if transitions(2, 1) > 10^(-5)
    arrow(stateCoord2-[0, separation], stateCoord1-[0, separation], 'wid', round(transitions(2, 1)*100));
    end
    
    warning on;
    
    text((minX+textindent), (maxY-textindent)-textsep*0, strcat('p12 = ', num2str(transitions(1, 2), precision)));
    text((minX+textindent), (maxY-textindent)-textsep*1, strcat('p21 = ', num2str(transitions(2, 1), precision)));
    
    axis([2 8 2 6]);   
end

if numStates == 3
    
    minX = 0; maxX = 8; minY = 2; maxY = 8;
    
    stateCoord1 = [4 6];
    stateCoord2 = [2 4];
    stateCoord3 = [6 4];
    
    circle(stateCoord1, occTot(1), 20);
    text(stateCoord1(1), stateCoord1(2)+1, strcat('D1 = ', num2str(diffCoeff(1), precision)));
    circle(stateCoord2, occTot(2), 20);
    text(stateCoord2(1), stateCoord2(2)-1, strcat('D2 = ', num2str(diffCoeff(2), precision)));
    circle(stateCoord3, occTot(3), 20);
    text(stateCoord3(1), stateCoord3(2)-1, strcat('D3 = ', num2str(diffCoeff(3), precision)));
    
    
    
    warning off;
    
    if transitions(1, 2) > 10^(-5)
    arrow(stateCoord1, stateCoord2, 'wid', round(transitions(1, 2)*100));
    end
    if transitions(2, 1) > 10^(-5)
    [dx, dy] = pol2cart(3*pi/4, separation);
    arrow(stateCoord2+[dx, dy], stateCoord1+[dx, dy], 'wid', round(transitions(2, 1)*100));
    end
    
    if transitions(1, 3) > 10^(-5)
    arrow(stateCoord1, stateCoord3, 'wid', round(transitions(1, 3)*100));
    end
    if transitions(3, 1) > 10^(-5)
    [dx, dy] = pol2cart(pi/4, separation);
    arrow(stateCoord3+[dx, dy], stateCoord1+[dx, dy], 'wid', round(transitions(3, 1)*100));
    end
    
    if transitions(2, 3) > 10^(-5)
    arrow(stateCoord2, stateCoord3, 'wid', round(transitions(2, 3)*100));
    end
    if transitions(3, 2) > 10^(-5)
    arrow(stateCoord3-[0, separation], stateCoord2-[0, separation], 'wid', round(transitions(3, 2)*100));
    end
    
    warning on;
    
    text((minX+textindent), (maxY-textindent)-textsep*0, strcat('p12 = ', num2str(transitions(1, 2), precision)));
    text((minX+textindent), (maxY-textindent)-textsep*1, strcat('p21 = ', num2str(transitions(2, 1), precision)));
    text((minX+textindent), (maxY-textindent)-textsep*2, strcat('p13 = ', num2str(transitions(1, 3), precision)));
    text((minX+textindent), (maxY-textindent)-textsep*3, strcat('p31 = ', num2str(transitions(3, 1), precision)));
    text((minX+textindent), (maxY-textindent)-textsep*4, strcat('p23 = ', num2str(transitions(2, 3), precision)));
    text((minX+textindent), (maxY-textindent)-textsep*5, strcat('p32 = ', num2str(transitions(3, 2), precision)));
    
    axis([0 8 2 8]);
end

if numStates == 4
    
    minX = 0; maxX = 8; minY = 2; maxY = 8;
    
    stateCoord1 = [4 6];
    stateCoord2 = [6 6];
    stateCoord3 = [4 4];
    stateCoord4 = [6 4];
    
    circle(stateCoord1, occTot(1), 20);
    text(stateCoord1(1), stateCoord1(2)+1, strcat('D1 = ', num2str(diffCoeff(1), precision)));
    circle(stateCoord2, occTot(2), 20);
    text(stateCoord2(1), stateCoord2(2)+1, strcat('D2 = ', num2str(diffCoeff(2), precision)));
    circle(stateCoord3, occTot(3), 20);
    text(stateCoord3(1), stateCoord3(2)-1, strcat('D3 = ', num2str(diffCoeff(3), precision)));
    circle(stateCoord4, occTot(4), 20);
    text(stateCoord4(1), stateCoord4(2)-1, strcat('D4 = ', num2str(diffCoeff(4), precision)));
    
    
    
    warning off;
    
    if transitions(1, 2) > 10^(-5)
    arrow(stateCoord1, stateCoord2, 'wid', round(transitions(1, 2)*100));
    end
    if transitions(2, 1) > 10^(-5)
    arrow(stateCoord2+[0, separation], stateCoord1+[0, separation], 'wid', round(transitions(2, 1)*100));
    end
    
    if transitions(1, 3) > 10^(-5)
    arrow(stateCoord1, stateCoord3, 'wid', round(transitions(1, 3)*100));
    end
    if transitions(3, 1) > 10^(-5)
    arrow(stateCoord3-[separation, 0], stateCoord1-[separation, 0], 'wid', round(transitions(3, 1)*100));
    end
    
    if transitions(2, 3) > 10^(-5)
    arrow(stateCoord2, stateCoord3, 'wid', round(transitions(2, 3)*100));
    end
    if transitions(3, 2) > 10^(-5)
    [dx, dy] = pol2cart(3*pi/4, separation);
    arrow(stateCoord3+[dx, dy], stateCoord2+[dx, dy], 'wid', round(transitions(3, 2)*100));
    end
    
    if transitions(3, 4) > 10^(-5)
    arrow(stateCoord3, stateCoord4, 'wid', round(transitions(3, 4)*100));
    end
    if transitions(4, 3) > 10^(-5)
    arrow(stateCoord4-[0, separation], stateCoord3-[0, separation], 'wid', round(transitions(4, 3)*100));
    end
    
    if transitions(2, 4) > 10^(-5)
    arrow(stateCoord2, stateCoord4, 'wid', round(transitions(2, 4)*100));
    end
    if transitions(4, 2) > 10^(-5)
    arrow(stateCoord4+[separation, 0], stateCoord2+[separation, 0], 'wid', round(transitions(4, 2)*100));
    end
    
    if transitions(1, 4) > 10^(-5)
    arrow(stateCoord1, stateCoord4, 'wid', round(transitions(1, 4)*100));
    end
    if transitions(4, 1) > 10^(-5)
    [dx, dy] = pol2cart(pi/4, separation);
    arrow(stateCoord4+[dx, dy], stateCoord1+[dx, dy], 'wid', round(transitions(4, 1)*100));
    end
    
    warning on;
    
    text((minX+textindent), (maxY-textindent)-textsep*0, strcat('p12 = ', num2str(transitions(1, 2), precision)));
    text((minX+textindent), (maxY-textindent)-textsep*1, strcat('p21 = ', num2str(transitions(2, 1), precision)));
    text((minX+textindent), (maxY-textindent)-textsep*2, strcat('p13 = ', num2str(transitions(1, 3), precision)));
    text((minX+textindent), (maxY-textindent)-textsep*3, strcat('p31 = ', num2str(transitions(3, 1), precision)));
    text((minX+textindent), (maxY-textindent)-textsep*4, strcat('p14 = ', num2str(transitions(1, 4), precision)));
    text((minX+textindent), (maxY-textindent)-textsep*5, strcat('p41 = ', num2str(transitions(4, 1), precision)));
    text((minX+textindent), (maxY-textindent)-textsep*6, strcat('p23 = ', num2str(transitions(2, 3), precision)));
    text((minX+textindent), (maxY-textindent)-textsep*7, strcat('p32 = ', num2str(transitions(3, 2), precision)));
    text((minX+textindent), (maxY-textindent)-textsep*8, strcat('p24 = ', num2str(transitions(2, 4), precision)));
    text((minX+textindent), (maxY-textindent)-textsep*9, strcat('p42 = ', num2str(transitions(4, 2), precision)));
    text((minX+textindent), (maxY-textindent)-textsep*10, strcat('p34 = ', num2str(transitions(3, 4), precision)));
    text((minX+textindent), (maxY-textindent)-textsep*11, strcat('p43 = ', num2str(transitions(4, 3), precision)));
    
    axis([0 8 2 8]);
end


set(gcf, 'Name', strcat(filename, '_HMM'))
