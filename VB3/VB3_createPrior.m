function W=VB3_createPrior(runinput,N)
% W=V3_createPrior(runinput,N)
%
% Creates a model structure W with N states, and prior distributions
% according to runinput, either a runinputfile or an options struct. By
% default no aggregated states are generated.

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
if(ischar(runinput) && exist(runinput,'file'))
    runinputfile = runinput;
    opt=VB3_getOptions(runinputfile);
    % if an option struct, read in the runinputfilename
elseif(isstruct(runinput))
    opt=runinput;
else
    error(['Not a valid input, aborting VB3_createPrior']);
end
%% start of actual code

%% default prior type choices
if(~isfield(opt,'prior_type_D'))
    warning('VB3_createPrior: prior_type_D not specified. Using default : mean_strengt')
    prior_type_D='mean_strength';
end
if(~isfield(opt,'prior_type_Pi'))
    warning('VB3_createPrior: prior_type_Pi not specified. Using default : flat')
    prior_type_Pi='flat';
end
if(~isfield(opt,'prior_type_A'))
    warning('VB3_createPrior: prior_type_A not specified. Using default : dwell_Bflat')
    prior_type_A='dwell_Bflat';
end
%% diffusion constant prior
if(strcmp(opt.prior_type_D,'mean_strength'))
    [W.PM.n,W.PM.c]=priorD_mean_strength(opt,N);
else
    error(['VB3_createPrior: did not recognize prior_type_D : ' opt.prior_type_D])
end
%% initial state prior
if(strcmp(opt.prior_type_Pi,'flat'))
    W.PM.wPi=ones(1,N); % flat initial state distribution
elseif(strcmp(opt.prior_type_Pi,'natmet13'))
    W.PM.wPi=ones(1,N)*opt.prior_piStrength/N; % the choice used in the nat. meth. 2013 paper.
else
    error(['VB3_createPrior: did not recognize prior_type_Pi : ' opt.prior_type_Pi])
end
%% transition prior
if(strcmp(opt.prior_type_A,'dwell_Bflat'))
    [W.PM.wa,W.PM.wB]=priorA_dwell_Bflat(opt,N);
elseif(strcmp(opt.prior_type_A,'natmet13'))
    [W.PM.wa,W.PM.wB]=priorA_natmet13(opt,N);
else
    error(['VB3_createPrior: did not recognize prior_type_A : ' opt.prior_type_A])
end
%% misc
W.dim=opt.dim;
W.N=N;
W.M.SA=1:N; % aggregate assignments
end
%% more complicated prior choices
function [n,c]=priorD_mean_strength(opt,N)
% Default diffusion constant prior, specified using the mean value and
% total strength per state.

timestep=opt.timestep;              % sampling time step
D0=opt.prior_D;                     % prior diffusion constant
Dn=opt.prior_Dstrength;             % strength of diffusion constant prior

% each emission variable gets same strength independent of model size
n=Dn*ones(1,N);
c=4*D0*timestep*n; % match <gamma>=<1/D> rather than <D> to avoid Inf

end

function [wa,wB]=priorA_dwell_Bflat(opt,N)
% Default transition prior, which specifies mean and std of the mean
% dwell times, and uses a flat prior for each row of the jump matrix B.

timestep=opt.timestep;          % sampling time step
t0=opt.prior_tD/timestep;       % prior mean dwell time
if(t0<1)
    error('VB3_createPrior/priorA_dwell_Bflat: prior dwell time must exceed sampling time step')
end

t0Var=(opt.prior_tDstd/timestep)^2; % prior dwell time variance

% conditional jump probabilities are also flat
wB=ones(N,N)-eye(N);

% mean dwell times are Gamma-distributed
wa=ones(N,2)+t0*(t0-1)/t0Var;
wa(:,2)=wa(:,1)*(t0-1);

% in case there is only one hidden states, B and a are not defined in the
% model:
if(N==1)
    wB=0;
    %wa=[0 0]; % there really should not be an a-variable for N=1, but
    %this is taken care of elsewhere (mostly in VBEMiterator).
end
end

function [wa,wB]=priorA_natmet13(opt,N)
% transition prior used in the nat. meth. (2013) paper. This is a prior
% on the transition matrix directly, which means that
% wa(j,1)=sum(wB(j,:)) is enforced. The prior is specified in terms of
% a mean dwell time and an overall row strength (same for all states).

timestep=opt.timestep;              % sampling time step
At0=opt.prior_tD/timestep;          % prior mean dwell time
An =opt.prior_tDstrength;           % transition probability strength (per row)

% transition matrix
A0=1;
if(N>1)
    A0=(1-1/At0)*eye(N)+1/At0/(N-1)*(ones(N,N)-eye(N));
end
wA =An*A0;           % each row gets prior strength An

wB=wA-diag(diag(wA));
wa=[sum(wB,2) diag(wA)];

end

