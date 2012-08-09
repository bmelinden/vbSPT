function W=VB3_createPrior(opt,N)
% W=V3_createPrior(opt,N)
%
% creates a model structure W with N states, and prior distributions
% according to the options structure opt.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_createPrior.m, model initialization in the vbSPT package
% =========================================================================
% 
% Copyright (C) 2012 Martin Lindén and Fredrik Persson
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
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%% start of actual code


timestep=opt.timestep;              % sampling time step
At0=opt.prior_tD/timestep;          % prior mean dwell time
D0=opt.prior_D;                     % prior diffusion constant
Dn=opt.prior_Dstrength;             % strength of diffusion constant prior
An =opt.prior_tDstrength;           % transition probability strength (per row)
nPi=opt.prior_piStrength;           % strength of initial state prior

% each emission variable gets same strength independent of model size
W.PM.n=Dn*ones(1,N);                
W.PM.c=4*D0*timestep*W.PM.n;

% transition matrix 
A0=1;
if(N>1)
    A0=(1-1/At0)*eye(N)+1/At0/(N-1)*(ones(N,N)-eye(N));
end
W.PM.wA =An*A0;           % each row gets prior strength An

% initial state
W.PM.wPi=nPi*ones(1,N)/N; % flat initial state distribution

% misc
W.dim=opt.dim;
W.N=N;

if(At0<1)
    error('VB3_createPrior: prior dwell time must exceed sampling time step')
end
