function VB3_printModel(W,Wstd,dt,fstr,Dfactor)
% VB3_printfModel(W,Wstd,dt,fstr)
% W     : vbSPT model
% Wstd  : vbSPT bootstrap std error estimate
% dt    : time step [s]
% fstr  : fprintf floating point format string (default 6.3)
% Dfactor: rescaling factor for diffusion constant (default: 1e-6)
% Dunits : default 'um^2/s'
% A structured printout of the vbSPT model W and bootstrap errors Wstd. 
% ML 2017-01-26

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_printModel, runs data analysis in the vbSPT package
% =========================================================================
% 
% Copyright (C) 2017 Martin Lindén
% 
% E-mail: bmelinden@gmail.com
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

%VB3_license('VB3_displayHMMmodel')

%% start oc actual code

N=W.N;
if(~exist('fstr','var') || isempty(fstr))
    fstr='6.3';
end
if(~exist('Dfactor','var') || isempty(fstr))
    Dfactor=1e-6;
end
if(~exist('Dunit','var') || isempty(fstr))
    Dunit='um^2/s';
end


%%%%%%%%%%
disp( [int2str(N) '-state model, ' int2str(sum(W.E.wA(:))) ' steps in data.'])
fullLine(fstr,N);
%%%%%%%%%%
fprintf('D [%6s]     : ',Dunit)
for n=1:N
fprintf(['%' fstr 'f+-%-' fstr 'f '],W.est.DdtMean(n)/dt*Dfactor,Wstd.est.DdtMean(n)/dt*Dfactor)
end
fprintf('\n')
%%%%%%%%%%
fprintf('occupancy      : ')
for n=1:N
fprintf(['%' fstr 'f+-%-' fstr 'f '],W.est.Ptot(n),Wstd.est.Ptot(n))
end
fprintf('\n')
%%%%%%%%%%
fprintf('dwell times [s]: ')
for n=1:N
fprintf(['%' fstr 'f+-%-' fstr 'f '],W.est.dwellMean(n)*dt,Wstd.est.dwellMean(n)*dt)
end
fprintf('\n')
%%%%%%%%%%
fprintf('dwell times[dt]: ')
for n=1:N
fprintf(['%' fstr 'f+-%-' fstr 'f '],W.est.dwellMean(n),Wstd.est.dwellMean(n))
end
fprintf('\n')
%%%%%%%%%%
fullLine(fstr,N)
fprintf('trans. matrix  : ')
for m=1:N
    if(m>1)
        fprintf('               : ')
    end
    for n=1:N
        fprintf(['%' fstr 'f+-%-' fstr 'f '],W.est.Amean(m,n),Wstd.est.Amean(m,n))
    end
    fprintf('\n')
end
%%%%%%%%%%
fullLine(fstr,N)

rArr=sum(W.M.wB-W.PM.wB,1)/sum(W.T-1);
drArr=sum(Wstd.M.wB,1)/sum(W.T-1);
fprintf('arrival rates  : ')
for n=1:N
fprintf(['%' fstr 'f+-%-' fstr 'f '],rArr(n),drArr(n))
end
fprintf('\n')

rDep=sum(W.M.wB-W.PM.wB,2)'/sum(W.T-1);
drDep=sum(Wstd.M.wB,2)'/sum(W.T-1);
fprintf('arrival rates  : ')
for n=1:N
fprintf(['%' fstr 'f+-%-' fstr 'f '],rDep(n),drDep(n))
end
fprintf('\n')
fullLine(fstr,N)

rDep=sum(W.E.wA.*(1-eye(W.N)),2)'/sum(W.T-1);

%%fprintf('%6.1e ',)
%fprintf('\n')
%fprintf('depart. rates  : ')
%fprintf('%6.1e ',sum(W.E.wA.*(1-eye(W.N)),2)'/sum(W.T-1))
%fprintf('\n') 
end
function fullLine(fstr,N)
        fprintf('-----------------')
        for n=1:N
            for m=1:2*floor(str2double(fstr))+3
                fprintf('-')
            end
        end
        fprintf('\n')
    end
