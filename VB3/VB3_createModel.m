function M=VB3_createModel(D,A,p0,W0,dt,strength)
% function M=VB3_createModel(D,A,p0,W0,dt,strength)
% 
% Creates a VB3 model from parameters and strength. The function either
% returns an M-field, or (if an input model W0 is given), adds the M-field
% to an input model. The model parameters are chosen so that the
% vatiational mean values coincide with the input parameters. This is good
% when constructing initial guesses for synthetic data with known
% parameters.
%
% Input:
% D  : vector of diffusion constants * timestep
% A  : transition matrix, needs to be properly normalized etc
% p0 : initial state probability (default = approx. stationary state of A)
% strength  : strength parameter, corresponding to the number of counts that go
%      into each parameter distribution (default = 1e4). 
% W0 : optional starting model. Only the PM and dim fields of this
%      model is inherited, if present. If not given, only the M-field
%      corresponding to the other input parameters are returned.
% dt : data timestep (default=1, which means that one can omit dt and input
%      D*dt for the diffusion constants directly)

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_createModel, creates a VB3 model from parameters and strengths.
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

%% check input parameters and size compatibility
if(exist('dt','var') && ~isempty(dt) && length(dt)==1 && dt>0)
   Ddt=D*dt;
else
    if(length(dt)~=1 || dt<=0)
        error('VB3_createModel: dt must be a positive scalar')
    end
   Ddt=D;
end
clear dt D;

lambda=sort(eig(A));
if(abs(lambda(end)-1)>1e-15) % then there is no proper stationary state
    error('VB3_createModel: A does not seem to be a proper transition matrix')
end

if(~exist('p0','var') || ~isempty(p0)) % then use stationary state of A
    % compute an equilibration time from second largest eigenvalue
    % at this point, one knows that lambda(end-1)<1
    if(size(A,1)>1)
        if(abs(real(lambda(end-1)))<1)
            neq=-1/log(abs(real(lambda(end-1))));
            p0=A^ceil(100*neq);
            p0=p0(1,:);
        else
            error('VB3_createModel: cannot compute equilibration time for transition matrix')
            endq
        end
    else
        p0=1;
    end
end
%clear lambda

N_Ddt=length(Ddt);
N_A=size(A);
N_p0=length(p0);
if( N_Ddt==N_A(1) && N_A(1)==N_A(2) && N_A(2)==N_p0)
    N=N_Ddt;
else
    error('VB3_createModel: incompatible input parameter sizes!')
end
N=N_Ddt;
clear N_Ddt N_A N_p0

if(~exist('strength','var') || isempty(strength))
    strength=1e4;
end
if(length(strength)~=1 || strength<=1)
    error('VB3_createModel: strength must be a scalar > 1')
end

%% compute M field
M.n=ones(1,N)*strength;
M.c=4*Ddt*(strength-1);
wA=A*strength;

M.wa=[sum(wA,2)-diag(wA) diag(wA)];
M.wB=wA-diag(diag(wA));

M.wPi=p0*strength;

M.SA=1:N;
%clear strength Ddt D A p0

% assemble result structure
if(exist('W0','var')) % thgen we shall return an output model
    Wout=struct;
    Wout.M=M;
    if(isfield(W0,'PM'))
        Wout.PM=W0.PM;
        fn=fieldnames(W0.PM);
        for n=1:length(fn)
           Wout.M.(fn{n})=Wout.M.(fn{n})+Wout.PM.(fn{n});           
        end
    end
    if(isfield(W0,'dim'))
        Wout.dim=W0.dim;
    end
    if(isfield(W0,'N') && W0.N~=N)
        error('VB3_createModel: input model size incompatible with parameters')
    end    
    M=Wout;
end
