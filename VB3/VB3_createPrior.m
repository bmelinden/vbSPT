function W=VB3_createPrior(runinput,N)
% W=V3_createPrior(runinput,N)
%
% Creates a model structure W with N states, and prior distributions
% according to runinput, either a runinputfile or an options struct.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_createPrior.m, model initialization in the vbSPT package
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

%% Parse input
% if an existing file, generate options structure
if(isstr(runinput) && exist(runinput)==2)
    runinputfile = runinput;
    opt=VB3_getOptions(runinputfile);
    % if an option struct, read in the runinputfilename
elseif(isstruct(runinput))
    opt=runinput;
else
    error(['Not a valid input, aborting VB3_createPrior']);
end

%% start of actual code
timestep=opt.timestep;              % sampling time step
D0=opt.prior_D;                     % prior diffusion constant
Dn=opt.prior_Dstrength;             % strength of diffusion constant prior

t0=opt.prior_tD/timestep;           % prior mean dwell time
t0Var=(opt.prior_tDstd/timestep)^2; % prior dwell time variance

% each emission variable gets same strength independent of model size
W.PM.n=Dn*ones(1,N);                
W.PM.c=4*D0*timestep*W.PM.n;

% initial states have a flat distribution
W.PM.wPi=ones(1,N); % flat initial state distribution

% conditional jump probabilities are also flat
W.PM.wB=ones(N,N)-eye(N);

% mean dwell times are Gamma-distributed
W.PM.wa=ones(N,2)+t0*(t0-1)/t0Var;
W.PM.wa(:,2)=W.PM.wa(:,1)*(t0-1);

% in case there is only one hidden states, B and a are not defined in the
% model:
if(N==1)
   W.PM.wB=0; 
   %W.PM.wa=[0 0]; % there really should not be an a-variable for N=1, but
   %this is taken care of elsewhere (mostly in VBEMiterator).
end

% misc
W.dim=opt.dim;
W.N=N;

if(t0<1)
    error('VB3_createPrior: prior dwell time must exceed sampling time step')
end
