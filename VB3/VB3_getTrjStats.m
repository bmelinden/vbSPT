function [M,Ptot,Ntot]=VB3_getTrjStats(X,dim,stateColumn,trjLmin)
% [M,Ptot,Ntot]=VB3_getTrjStats(X,dim,stateColumn,trjLmin)
% assemble optimal model from data with known hidden state sequence.
%
% X     : cell vector with trajectory data. In each cell, the first
%         dim columns are position data, and column stateColumn are the
%         hidden state sequence. dim=0 means no position stats are
%         collected.
% M     : M-field of the model resulting from the data and hidden state
%         sequence in X.
% two extra variables are computed:
% Ntot(j) = total occupancy of state s=j, Ptot=Ntot/sum(Ntot);
%
% trjLmin: minimum trajectory length (default 2)
% The current code cannot handle more than 100 states.
%
% M.L. 2012-07-05   : started

warning('VB3_getTrjSttats has not been thouroughly tested. Use at your own risk...')

%% parameters
if(~exist('trjLmin','var') || isempty(trjLmin))
    trjLmin=2;
end
if(trjLmin<2)
    error('minimum trajectory length cannot be < 2!')
end

%% compute model
n=zeros(1,100);
c=zeros(1,100);
wA=zeros(100,100);
wPi=zeros(1,100);
Ptot=zeros(1,100);
N=0;
for t=1:length(X)
    if(size(X{t},1)>=trjLmin)
        % state statistics
        s=X{t}(1:end-1,stateColumn);
        N=max(N,max(s));
        wPi(s(1))=wPi(s(1))+1;
        
        for m=1:length(s)-1
            s0=s(m);
            s1=s(m+1);
            wA(s0,s1)=wA(s0,s1)+1;
            Ptot(s0)=Ptot(s0)+1;
        end
        Ptot(s(end))=Ptot(s(end))+1;
        % position statitics
        %E(m).c=zeros(1,N);             % sum_{t=2}^T P(s(t)=j,c(t)=1).*dx(t)^2
        %for j=1:N
        %    E(m).c(j)=sum(pst(:,j).*dx2{m});
        %end
     
        
        
        if(dim>0)
            x=X{t}(:,1:dim);
            dx2=sum(diff(x,1,1).^2,2);
            for m=1:length(s)
                n(s(m))=n(s(m))+dim/2;
                c(s(m))=c(s(m))+dx2(m);
            end
        end
    end
end

M.n=n(1:N);
M.c=c(1:N);
M.wA=wA(1:N,1:N);
M.wPi=wPi(1:N);
Ntot=Ptot(1:N);
Ptot=Ptot(1:N)/sum(Ptot(1:N));
end
