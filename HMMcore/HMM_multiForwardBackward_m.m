% [lnZ,wA,pst]=HMM_multiForwardBackward_m(Q,H,iEnd,doBackward) 
%
% performs forward-backward sweeps on HMM-type time series stacked on top
% of each other. 
%
% Q    : transition matrix (not necessarily normalized) transition matrix
% H    : emission likelihood, including initial state probabilities on
%        appropriate rows 
% iEnd : indices to ends of individual time series, e.g., individual
%        trajectories run from  1 -> iEnd(1), iEnd(1)+2 -> iEnd(2), etc.
%        This is so that the algorithm know at what lines the transition
%        matrix should be omitted from the calculations.
% doBackward: flag to activate a backward sweep and compute pst (default
%       false).
%        
% lnZ  : log of normalization constant for the forward sweep
% wA   : transition count matrix (requires backward sweep)
% pst  : pst(t,j)= P(s(t)==j), averaged over all state sequences (some
%        extra computational cost in addition to backward sweep).
%        
% This is the matlab implementation of HMM_multiForwardBackward.c, created
% for debugging purposes (much slower). 
%
% M.L. 2013-05-21

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rHMM_multiForwardBackward_m.m, part of HMMcore/
% =========================================================================
% 
% Copyright (C) 2013 Martin Lindén, E-mail: bmelinden@gmail.com
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
% Additional permission under GNU GPL version 3 section 7
% 
% If you modify this Program, or any covered work, by linking or combining it
% with Matlab or any Matlab toolbox, the licensors of this Program grant you 
% additional permission to convey the resulting work.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.

%% start of actual code
function [lnZ,wA,pst]=HMM_multiForwardBackward_m(Q,H,iEnd,doBackward) 

if(~exist('doBackward','var'))
    doBackward=false;
end

[T,N]=size(H);

% output
lnZ=0;
pst=[];
wA=[];


% forward sweep
Za=0;
alpha=zeros(T,N);

% backward sweep
if(doBackward)
    QT=Q';
    beta=zeros(T,N);
    Zb=0;
    wA=zeros(N,N);
    P=zeros(N,N);
    pst=zeros(T,N);
end

tStart=1;
for n=1:length(iEnd)
    tEnd=iEnd(n);
    
    % forward sweep
    alpha(tStart,:)=H(tStart,:);
    Za=sum(alpha(tStart,:));
    alpha(tStart,:)=alpha(tStart,:)/Za;    
    lnZ=lnZ+log(Za);
    for t=(tStart+1):tEnd
        alpha(t,:)=(alpha(t-1,:)*Q).*H(t,:);
        Za=sum(alpha(t,:));
        alpha(t,:)=alpha(t,:)/Za;
        lnZ=lnZ+log(Za);
    end % forward iterations    

    if(doBackward) % backward sweep
        beta(tEnd,:)=ones(1,N)/N;
        for t=tEnd-1:-1:tStart
            beta(t,:)=(beta(t+1,:).*H(t+1,:))*QT;
            Zb=sum(beta(t,:));
            beta(t,:)=beta(t,:)/Zb;
        end % backward iterations
        
        % transition counts
        for t=(tStart+1):tEnd
            ZW=0;
            for j=1:N
                for k=1:N
                    P(j,k)=alpha(t-1,j)*Q(j,k)*H(t,k)*beta(t,k);
                    ZW=ZW+P(j,k);
                end
            end
            wA=wA+P/ZW;
        end
    end
    tStart=tEnd+1; % start index for next trajectory
end
if(doBackward) 
        pst=rowNormalize(alpha.*beta);    
end
