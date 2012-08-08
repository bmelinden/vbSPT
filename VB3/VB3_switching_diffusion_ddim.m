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
 

%% change log
% M.L. 2012-01-31 : simulate d-dimensional free diffusion
% M.L. 2011-12-05 : started

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
