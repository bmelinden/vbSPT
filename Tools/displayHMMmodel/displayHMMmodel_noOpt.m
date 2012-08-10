
function displayHMM_noOpt(model)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File for getting out and graphically representing HMM model results from
% VB3. Input model is a vbSPT model object.
% FP 2012-04-13

VB3_license('displayHMM_noOpt')

% Set parameters
timeStep = 3e-3; % s     
dim = 1;

% Read in HMM results
numStates = model.N     % Number of states
 
[diffCoeff, ind] = sort(model.est.DdtMean./timeStep/(10^6)); % Mean diff. coeff in um^2/s provided the data is in nm
diffCoeffStd = model.est.Ddtstd(ind)./timeStep/(10^6);

occTot = model.est.Ptot(ind);       % Total occupation percentage
dwellTime = model.est.dwellMean(ind)'*timeStep;     % Mean dwelltime in each state

A = model.M.wA - model.PM.wA;   %The transition probability matrix with the prior values subtracted
A = spdiags (sum (A,2), 0, numStates, numStates) \ A; % Rownormalize the matrix
transitions = A(ind, ind);

% Print all the numbers in the Matlab prompt
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
% Parameters for graphical representation

separation = 0.3;
textsep = 0.2;
textindent = 0.1;
precision = 3;

% Graphical representation for models with 2-4 states
if numStates == 2

    minX = 2; maxX = 8; minY = 2; maxY = 6;
    
    figure;
    hold on;
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
    
    figure;
    hold on;
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
    
    figure;
    hold on;
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


set(gcf, 'Name', 'HMM model')
