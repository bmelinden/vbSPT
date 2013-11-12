% [logL,W]=MC1_logL(p0,A,g,dat,estimates)
%
% Evaluate log likelihood for the parameter set p0,A,D for diffusive data
% structure dat (see SPT_preprocess).
%
% estimates: 0 - compute only logL (fastest)
%            1 - compute light-weight estimates (requires a backward sweep, default)
%            2 - compute some memory intensive estimates (requires a
%            Viterbi sweep)
% 
% ML 2013

light_estimates=true;
heavy_estimates=false;
if(exist('estimates','var'))
    if(estimates==0)
        light_estimates=false;
    else
        light_estimates=true;
    end
    if(estimates<2)
        heavy_estimates=false;
    else
        heavy_estimates=true;
    end
end
Ntrj=length(dat.T);
N=length(p0);
dim=dat.dim;
%% Actual computation
lnQ=log(A);
lnQmax=max(max(lnQ));
Q=exp(lnQ-lnQmax);

clear lnH H lnH0
lnZQ=zeros(1,Ntrj);lnZq=zeros(1,Ntrj);lnZz=zeros(1,Ntrj);
W=struct;
E=struct;
E.wPi=zeros(1,N);
E.wA=zeros(N,N);
E.n=zeros(1,N);
E.c=zeros(1,N);
%% assemble hidden state distribution
trjstarts=[1 dat.end(1:end-1)+1];
% we compute log(H) first, and then exponentiate
lnH1=log(p0);                   % initial state term
lnH0=dim/2*log(g/pi);           % data-independent term
lnH  =zeros(length(dat.dx2),N); % data-dependent terms
for j=1:N
    lnH(:,j)=lnH0(j)-g(j)*dat.dx2;
    lnH(trjstarts,j)=lnH(trjstarts,j)+lnH1(j);
end
lnHMax=max(lnH,[],2);
H=zeros(length(dat.dx2),N);
for j=1:N % now we compute the actual pointwise emission contribution
    H(:,j)=exp(lnH(:,j)-lnHMax);
end
%% forward sweep  (depends only on Q and H)
if(light_estimates)
    [lnZz,E.wA,pst]=HMM_multiForwardBackward(Q,H,dat.end);
else
    [lnZz]=HMM_multiForwardBackward(Q,H,dat.end);    
end

%[lnZz,E.wA,pst]=VB_multiForwardBackward(Q,H,dat.end,light_estimates);

% forward sweep normalization constant (same as VB3)
lnZQ=(sum(dat.T-2))*lnQmax;
lnZq=sum(lnHMax);
% assemble free energy
logL=lnZQ+lnZq+sum(lnZz);
if(~isfinite(logL))
    error('MC1_logL not finite (lnZ)')
end

if(light_estimates)
    % transition counts
    E.wPi=sum(pst(trjstarts,:),1);          % <s_1>=<\delta_{j,s_1}>_{q(s)}
    
    % for the emission models
    E.n = E.n+dim/2*sum(pst,1);   % sum_{t=2}^T p(s_t)
    for j=1:N                     % sum_{t=2}^T P(s(t)=j).*dx(t)^2
        E.c(j) = E.c(j)+sum(pst(:,j).*dat.dx2);
    end
    
    % trajectory specific light-weight estimates (always)
    % occupation probabilities
    WestTs=sum(pst,1);           % time spent in each state
    WestPs=rowNormalize(WestTs); % occupation probability
    
    % check for problems
    isNanInf=(sum(~isfinite([E.n E.c E.wPi E.wA(1:end)]))>1);
    if(isNanInf)
        error('VB2_VBEMiter:Efield_not_finite','Nan/Inf generated in E step')
    end
    W.E=E;
    W.Fterms.lnZQ=lnZQ;
    W.Fterms.lnZq=lnZq;
    W.Fterms.lnZz=lnZz;
    W.est.Ttot=sum(pst,1);    
    W.est.Ptot=WestTs/sum(WestTs); % total average
    W.est.Q=Q;
    W.est.lnQ=log(Q);
    if(heavy_estimates)
        W.est2.pst   =pst;
        W.est2.H     =H;
        W.est2.lnH   =lnH;
        [~,W.est2.sMaxP] = max(pst,[],2);
    end
    W.F=logL;
    W.Fterms.Fterms=[ W.Fterms.lnZQ+W.Fterms.lnZq+W.Fterms.lnZz];
    W.Fterms.FtermsNames='[lnZQ+lnZq+lnZz]';
    if(~isfinite(W.F))
        disp(W.Fterms)
        error('ML1_logL:F_not_finite','Nan/Inf generated in likelihhod')
    end
end
%% some large estimates that can be done on the last iteration.
if(heavy_estimates)
    %W.est2.viterbiM=HMM_multiViterbi_log0(log(Q),log(H),dat.end);
    W.est2.viterbi=HMM_multiViterbi_log(log(Q),log(H),dat.end);
end

