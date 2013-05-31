function [x,s]=VB3_switching_diffusion_ddim(A,g,T,d)
% [x,s]=VB3_switching_diffusion_ddim(A,g,T,d)
%
% simulate hidden states s(t) and position x(t) for a diffusion that
% switches between different diffusion constant according to a Markov
% process. Displaced so that x(1)=0;
%
% A(i,j)=p(s_t=j|s_{t-1}=i), so the evolution of the hidden states is
% given by p(t)=p(t-1)*A, where p(t) is a row vector
%
% g is a vector of precision parameters: if s(t)=i, then 
% 2*D(i)*dt = Var[x(t+1)-x(t)] = 1/2/g(i). 
%
% d is the dimension of the output data (default:2)
%
% T is the number of time points to simulate.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_switching_diffusion_ddim.m, simulate d-dimensional diffusion with 
% multiple diffusion constants, part of the vbSPT package
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
%% start of actual code


%% parameter handling
if(~exist('d','var') || isempty(d)); d=2;end
g=reshape(g,length(g),1);

N=length(g);
%% start of actual code
%% stationary distribution for initial state
A0=A^10000;
p0=A0(1,:);
if(abs(sum(p0)-1)>1e-10)
    error('stationary state might not be normalized')
end
p0=p0/sum(p0);
%% generate state trajectories
s=zeros(T-1,1);
s(1)=find(rand<cumsum(p0),1); % initial state close to the stationary distribution
cumA=cumsum(A,2);
for t=2:T-1
    ras=rand;
    s(t)=find(ras<cumA(s(t-1),:),1);
end
%% generate noise trajectory
dx=zeros(T-1,d);
for k=1:d
    dx(1:end,k)=randn(T-1,1).*sqrt(0.5./g(s));
end
x=[zeros(1,d); cumsum(dx,1)];
