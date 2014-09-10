function [W,C,F]=VB3_VBEMiterator(W,dat,varargin)
%% [W,C,F]=VB3_VBEMiterator(W,dat,varargin)
%
% Perform VBEM iterations, on the VB structure W, with data (structure)
% X, until convergence. This version accepts d-dimensional data (T by d
% matrices), and uses a 'short forward' model, i.e., hidden states are
% associated with steps, not positions, and the hidden state at time t is
% associated with the particle motion on the interval t -> t+dt.
%
% options:
% 'estimate'   : if given, add a field est2, which contains memory and
%                computer intensive estimates, such as the Viterbi path.
% 'slim'       : if given, remove some potentially bulky fields to decrease
%                storage footprint of the model.
%
% 'maxIter',n  : run at most n VBE iterations.
% 'relTolF',tf : convergence criterion for relative change in likelihood bound.
% 'tolPar' ,tp : convergence criterion for M-step parameters.
%                iterate until
%                |dF/F| <= tf and max|dPar(i)/Par(i)| <= tp,
%                where Par(i) are all parameters for the variational
%                distributions. Default values are
%                maxIter=1000, relTolF=1e-8, tolPar=1e-2.
% 'outputLevel', {0,1,2}
%                0: no output, 1: display convergence, not progress (default),
%                2: display convergence measures for each iteration
%
% Fields in the VB1 tructure W:
% W.M.wPi;         : initial state
% W.M.wa,W.M.wB;         : transition counts for s(t) transition matrix
% W.M.n; W.M.c;  : B-distribution for s(t)
% Prior parameters are stored in W.PM
% E-step parameters are under W.E
% W.F              : lower bound on the likelihood
% W.est            : various useful estimates that do not take up a lot of
%                    memory
% W.est2           : more memory intensensive estimates
%
% This function uses the mex file HMM_multiForwardBackward.mexXXX from
% HMMcore for the computer intensive nner loops, where XXX is a platform
% dependent extension. Please refer to HMMcore/compile_code.m to
% (re)compile for your system. 

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_VBEMiterator, variational EM iterations in the vbSPT package
% =========================================================================
% 
% Copyright (C) 2014 Martin Lindén and Fredrik Persson
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
% process options
C.maxIter=1000;
C.relTolF=1e-8;
C.tolPar=1e-2;
do_estimates=false;
do_slim=false;
displayProgress=false;
displayExit=true;

if(nargin>2)        % then parse options
    k=1;            % argument counter
    kmax=nargin-2;  % stopping criterion
    while(k<=kmax)
        option=varargin{k};
        if(strcmpi(option,'estimate'))
            do_estimates=true;
            k=k+1;
        elseif(strcmpi(option,'slim'))
            do_slim=true;            
            k=k+1;
        elseif(strcmpi(option,'maxIter'))
            if(~isempty(varargin{k+1}))
                C.maxIter=varargin{k+1};
                if(~isnumeric(C.maxIter) || C.maxIter~=round(C.maxIter) || C.maxIter<=0)
                    error('VB1_VBEMiter: maxIter option must be followed by a positive integer.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'relTolF'))
            if(~isempty(varargin{k+1}))
                C.relTolF=varargin{k+1};
                if(~isnumeric(C.relTolF) || C.relTolF<=0)
                    error('VB1_VBEMiter: relTolF option must be followed by a positive number.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'tolPar'))
            if(~isempty(varargin{k+1}))
                C.tolPar=varargin{k+1};
                if(~isnumeric(C.tolPar) || C.tolPar<=0)
                    error('VB1_VBEMiter: tolPar option must be followed by a positive number.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'outputLevel'))
            if(~isempty(varargin{k+1}))
                outputLevel=varargin{k+1};
                if(~isnumeric(outputLevel) || ~ismember(outputLevel,[0 1 2]))
                    error('VB1_VBEMiter: outpuLevel option must be followed by 0, 1, or 2.')
                end
                if(outputLevel==0) % adjust output settings
                    displayExit=false;
                elseif(outputLevel==2)
                    displayProgress=true;
                end
            end
            k=k+2;
        else
            error(['VB1_VBEMiter: option ' option ' not recognized.'])
        end
    end
end
%% remove old estimates
if(isfield(W,'est'))
    W=rmfield(W,'est');
end
if(isfield(W,'est2'))
    W=rmfield(W,'est2');
end
if(isfield(W,'Fterms'))
    W=rmfield(W,'Fterms');
end
%% finish convergence structure
C.iter    =0;
C.converged=false;
C.W0=W;
%% preprocess the data and insert default aggregation
%%%Ntrj=length(X);
%%%dx2=cell(1,Ntrj);
W.T=dat.T;
W.dim=dat.dim;
N=size(W.PM.wB,1);
W.N=N;
if(~isfield(W.M,'SA')) % add default state aggregation (no aggregation)
    W.M.SA=1:W.N;
end
%% initialize VBEM iterations
runMore=true;
C.exitStatus='';
Wm3=[];Wm2=[];Wm1=[];
%% iterate
E=struct;
E.wPi=zeros(1,N);
E.wA=zeros(N,N);
E.n=zeros(1,N);
E.c=zeros(1,N);
while(runMore)
    % keep a short history in case something goes wrong...
    Wm3=Wm2;Wm2=Wm1;Wm1=W;
    %try % try one round of VBEM iterations
    %% M-step: contributions from all trajectories are added
    if(isfield(W,'E')) % then start with M step
        
        % transitions and initial conditions
        W.M.wPi = W.PM.wPi + W.E.wPi;
        wB=W.E.wA.*(1-eye(W.N)).*(W.PM.wB>0); % only allowed transitions included
        W.M.wa =  W.PM.wa  + [sum(wB,2) diag(W.E.wA)];
        W.M.wB =  W.PM.wB  + wB;

        % emission model part, with aggregated states
        for a=1:max(W.M.SA)
            % all states in aggregate a gets emission statistics from
            % all states in the same aggregate
            ind=find(a==W.M.SA);
            W.M.n(ind)  = W.PM.n(ind)  + sum(W.E.n(ind));
            W.M.c(ind)  = W.PM.c(ind)  + sum(W.E.c(ind));
        end
    end
    clear wB
    % check for problems
    isNanInf=(sum(~isfinite([W.M.wPi W.M.wa(1:end) W.M.wB(1:end) W.M.n  W.M.c]))>1);
    if(isNanInf)
        error('VB1_VBEMiter:Mfield_not_finite','Nan/Inf generated in VBM step')
    end
    %% E-steps starts here
    % coupling matrix is the same for all trajectories
    %lnQ=psi(W.M.wA)-psi(sum(W.M.wA,2))*ones(1,W.N); % old version
    lnQ=zeros(W.N,W.N);
    wB0=sum(W.M.wB,2);
    wa0=sum(W.M.wa,2);
    for i=1:W.N
        lnQ(i,i)=psi(W.M.wa(i,2))-psi(wa0(i));
        for j=[1:(i-1) (i+1):W.N]
            lnQ(i,j)=psi(W.M.wa(i,1))-psi(wa0(i))...
                    +psi(W.M.wB(i,j))-psi(wB0(i));
        end
    end    
    
    lnQmax=max(max(lnQ));
    Q=exp(lnQ-lnQmax)+0*(eps);
    % for estimates of transition rates and dwell times: s(t)
    W.est.Q=Q;
    W.est.lnQ=lnQ;
    
    lnH1=psi(W.M.wPi)-psi(sum(W.M.wPi)); % initial state probability    
    T=W.T;N=W.N;dim=W.dim;
    WMn=W.M.n;
    WMc=W.M.c;

    %% trial emission priobability q(s,t)
    % variational pointwise contributions, notation as in ML1 notes
    % except perhaps an off-by-one difference in time index
        
    trjStarts=[1 dat.end(1:end-1)+1];
    trjEnds=dat.end;
    
    % we compute log(H) first, and then exponentiate
    lnH1=psi(W.M.wPi)-psi(sum(W.M.wPi)); % initial state probability
    lnH0=dim/2*(psi(WMn)-log(pi*WMc));   % data-independent term, same for all t
    lnH  =zeros(length(dat.dx2),N); % data-dependent terms
    for j=1:N
        lnH(:,j)=lnH0(j)-W.M.n(j)/W.M.c(j)*dat.dx2;        
        lnH(trjStarts,j)=lnH(trjStarts,j)+lnH1(j);
    end
    lnHMax=max(lnH,[],2);
    H=zeros(length(dat.dx2),N);
    for j=1:N % now we compute the actual pointwise emission contribution
        H(:,j)=exp(lnH(:,j)-lnHMax);
    end 
    %% forward sweep  (depends only on Q and H)
    [lnZz,E.wA,pst]=HMM_multiForwardBackward(Q,H,trjEnds);
    % forward sweep normalization constant (same as VB3)
    lnZQ=(sum(dat.T-2))*lnQmax;
    lnZq=sum(lnHMax);

    % compute quantities for next M-step
    %% transition counts
    E.wPi=sum(pst(trjStarts,:),1);          % <s_1>=<\delta_{j,s_1}>_{q(s)}    
    % E.wA: already done!
    %% for the emission models
    E.n=dim/2*sum(pst,1);   % sum_{t=2}^T p(s_t)
    E.c=zeros(1,N);             % sum_{t=2}^T P(s(t)=j,c(t)=1).*dx(t)^2
    for j=1:N
        E.c(j)=sum(pst(:,j).*dat.dx2);
    end
    W.E=E;
    %% check for problems
    %isNanInf=(sum(~isfinite([E(m).n E(m).c]))>1);
    %if(isNanInf)
    %    error('VB2_VBEMiter:Efield_not_finite','Nan/Inf generated in VBE step')
    %end

    % check for problems
    isNanInf=(sum(~isfinite([E.n E.c]))>1);
    if(isNanInf)
        error('VB2_VBEMiter:Efield_not_finite','Nan/Inf generated in VBE step')
    end
    %% assemble free energy
    F=lnZQ+lnZq+lnZz;
    if(~isfinite(F))
        error('VB3_VBEM: F not finite (lnZ)')
    end    
    % KL divergence of transition probabilities of s(t), new
    % parameterization
    KL_a=zeros(W.N,1);
    if(W.N>1) % a is only defined if N>1
        wa0=sum(W.M.wa,2);
        ua0=sum(W.PM.wa,2);
        KL_a=gammaln(wa0)-gammaln(ua0)...
            -(wa0-ua0).*psi(wa0)-(...
            gammaln(W.M.wa(:,1))-gammaln(W.PM.wa(:,1))...
            -(W.M.wa(:,1)-W.PM.wa(:,1)).*psi(W.M.wa(:,1))...
            +gammaln(W.M.wa(:,2))-gammaln(W.PM.wa(:,2))...
            -(W.M.wa(:,2)-W.PM.wa(:,2)).*psi(W.M.wa(:,2)));
    end
    W.Fterms.aterms=-KL_a;
    F=F-sum(KL_a);
    if(~isfinite(F))
        error('VB3_VBEM: F not finite (KL_a)')
    end
    clear wa0 ua0;    
    % jump probabilities
    KL_B=zeros(1,W.N);
    if(W.N>1) % B is only defined for N>1        
        for k=1:W.N
            %ind=setdiff(1:N,k); % only include non-diagonal elements
            ind=find(W.PM.wB(k,:)>0); % only include non-zero elements
            wB0=sum(W.M.wB(k,ind));
            uB0=sum(W.PM.wB(k,ind));
            KL_B(k)=gammaln(wB0)-gammaln(uB0)-(wB0-uB0)*psi(wB0)...
                +sum((W.M.wB(k,ind)-W.PM.wB(k,ind)).*psi(W.M.wB(k,ind))...
                -gammaln(W.M.wB(k,ind))+gammaln(W.PM.wB(k,ind)));
        end
    end
    W.Fterms.Bterms=-KL_B;
    F=F-sum(KL_B);
    if(~isfinite(F))
        error('VB3_VBEM: F not finite (KL_B)')
    end    
    clear wA0 uA0 ind;
    % KL divergence of initial state probability 
    u0Pi=sum(W.PM.wPi);
    w0Pi=sum(W.M.wPi);
    KL_pi=gammaln(w0Pi)-gammaln(u0Pi)...
          +sum((gammaln(W.PM.wPi)-gammaln(W.M.wPi))...
               +(W.M.wPi-W.PM.wPi).*(psi(W.M.wPi)-psi(w0Pi)));
    W.Fterms.piTerms=-KL_pi;
    F=F-KL_pi;    
    if(~isfinite(F))
        error('VB3_VBEM: F not finite (KL_pi)')
    end    
    % KL divergence of emission parameters
    KL_gj= W.PM.n.*log(W.M.c./W.PM.c)...
        -W.M.n.*(1-W.PM.c./W.M.c)...
        -gammaln(W.M.n)+gammaln(W.PM.n)...
        +(W.M.n-W.PM.n).*psi(W.M.n);
    % remove duplicate terms in each aggregate
    for a=1:max(W.M.SA)
       ind=find(a==W.M.SA);
       KL_gj(ind(2:end))=0;
    end
    W.Fterms.gTerms=-KL_gj;
    F=F-sum(KL_gj);
    if(~isfinite(F))
        error('VB3_VBEM: F not finite (KL_gj)')
    end
    %% assembly of the free energy
    W.F=F;    
    W.Fterms.lnZQ=lnZQ;
    W.Fterms.lnZq=lnZq;
    W.Fterms.lnZz=lnZz;
            
    if(~isfinite(W.F))
        error('VB_VBEMiter:F_not_finite','Nan/Inf generated in lower bound')
    end
    %catch me
    %% catch potential errors
    %me.getReport
    %crashfile=sprintf('VB1_VBEMiterator_error_%f.mat',now);
    %save(crashfile)
    %runMore=false;
    %C.exitStatus=[C.exitStatus 'VB1_VBEMiterator generated error'];
    %error(['VB1_VBEMiterator: VB iterations returned error. Saved state to ' crashfile])
    %end
    C.iter=C.iter+1;
    %% do convergence check!
    if(isfield(Wm1,'F')) % then we can check convergence
        %% converge criterion in terms of relative changes in F and parameter values
        if(~isfinite(W.F)) % check for problem
            crashfile=sprintf('VB1_VBEMiterator_NaNInf_%f.mat',now);
            save(crashfile);
            runMore=false;
            error(['VB1_VBEMiterator found W.F=NaN or Inf. Saving state to ' crashfile])
            C.exitStatus=[C.exitStatus 'W.F=NaN or Inf. '];
            pause(0.5)
        end
        C.dFrel=(W.F-Wm1.F)/abs(W.F); % convergence statistic for F
        % convergence statistics for parameters
        fM=fields(W.M);
        C.dPrel=-Inf;
        C.limitPar='';
        for k=1:length(fM)
            Pnew=W.M.(fM{k})(1:end);
            Pold=Wm1.M.(fM{k})(1:end);
            ind=find(Pnew~=0);
            dPrel=max(abs((Pnew(ind)-Pold(ind))./Pnew(ind)));
            if(dPrel>C.dPrel)
                C.limitPar=['M.' fM{k}];
                C.dPrel=dPrel;
            end
        end
        
        Fconverged=(abs(C.dFrel)<C.relTolF);
        Pconverged=C.dPrel<C.tolPar;
        if(Fconverged && Pconverged)
            C.exitStatus=['Converged normally after ' int2str(C.iter) ' iterations, relTolF = ' num2str(C.relTolF) ', tolPar = ' num2str(C.tolPar) ', (' C.limitPar ' limiting).'];
            C.converged=true;
            runMore=false;
        end
        if(C.iter>=C.maxIter) % check for max iterations
            runMore=false;
            C.exitStatus=[C.exitStatus 'Maximum number of iterations reached. '];
        end
        if(displayProgress) % display convergence progress
            if(mod(C.iter,10)==1)
                displayHeader();
            end
            fm=displayConvergence();
            if(fm)
                disp('---------- lower bound decrease?!!! ----------')
            end
        end
        if(C.iter<2 && C.iter<C.maxIter) % run at least 2 iterations unless specifi
            runMore=true;
            C.exitStatus='';
        end
    end
    %% estimates (only last iteration)
    if(~runMore)
        
        % global light-weight estimates (always)
        wa0=sum(W.M.wa,2);
        W.est.aMean=W.M.wa(:,1)./wa0;
        W.est.aMode=(W.M.wa(:,1)-1)./(wa0-2);
        W.est.aVar=W.M.wa(:,1).*W.M.wa(:,2)./(wa0.^2.*(1+wa0));
        
        wB0=sum(W.M.wB,2)*ones(1,W.N);
        eyeB=1-eye(W.N);
        W.est.Bmean=W.M.wB./wB0;
        W.est.Bmode=(W.M.wB-1+eye(W.N))./(wB0-W.N+1);
        W.est.Bvar=W.M.wB.*(wB0.*eyeB-W.M.wB)./(wB0.^2.*(1+wB0));
        B2=W.M.wB.*(eyeB+W.M.wB)./(wB0.*(1+wB0));  % <Bjk^2>
        
        W.est.Amean=diag(W.M.wa(:,2)./sum(W.M.wa,2))...
            +(W.M.wa(:,1)./sum(W.M.wa,2)./sum(W.M.wB,2))*ones(1,W.N).*W.M.wB;
        %W.est.Amode : have not figured that one out yet (ML 2014-05.02)
        W.est.Astd=diag(W.est.aVar)...
            +(W.est.aVar*ones(1,W.N)).*B2...
            +(W.est.aMean.^2*ones(1,W.N)).*W.est.Bvar;
        
        W.est.dwellMean=1./W.est.aMean;
        W.est.dwellMode=wa0./(1+W.M.wa(:,1));
        clear wB0 eyeB B2 wa0
        
        % emission parameters
        W.est.gMean=W.M.n./W.M.c;
        W.est.gMode=(W.M.n-1)./W.M.c;
        W.est.gStd=sqrt(W.M.n./W.M.c.^2); % sqrt(Var(g))
        W.est.DdtMean=W.M.c/4./(W.M.n-1);
        W.est.DdtMode=W.M.c/4./(W.M.n+1);
        W.est.Ddtstd=W.M.c/4./(W.M.n-1)./sqrt(W.M.n-2);

        % occupation        
        W.est.Ttot=sum(pst,1);
        W.est.Ptot=W.est.Ttot/sum(W.est.Ttot);

        W.Fterms.Fterms=[ W.Fterms.lnZQ+W.Fterms.lnZq+W.Fterms.lnZz -sum(KL_a) -sum(KL_B) -sum(KL_pi) -sum(KL_gj)];
        W.Fterms.FtermsNames='[lnZQ+lnZq+lnZz -sum(KL_a) -sum(KL_B) -sum(KL_pi) -sum(KL_gj)]';

        
        
        W.est.Ts=sum(pst,1);
        W.est.Ps=W.est.Ts/sum(W.est.Ts);
        %% potentially demanding estimates (only if asked)
        if(do_estimates)
            % extract trajectory estimates
            Wviterbi=uint8(HMM_multiViterbi_log(lnQ,lnH,trjEnds)); % Viterbi path
            [~,WsMaxP]=max(pst,[],2);
            
            for kk=1:length(trjStarts)
                W.est2.pst{kk}    =           pst(trjStarts(kk):trjEnds(kk),:);
                W.est2.H{kk}      =             H(trjStarts(kk):trjEnds(kk),:);
                W.est2.lnH{kk}    =           lnH(trjStarts(kk):trjEnds(kk),:);
                W.est2.lnHMax{kk} =        lnHMax(trjStarts(kk):trjEnds(kk),:);
                W.est2.sMaxP{kk}  =uint8(  WsMaxP(trjStarts(kk):trjEnds(kk)));
                W.est2.viterbi{kk}=uint8(Wviterbi(trjStarts(kk):trjEnds(kk)));
                %W.est2.start=trjstarts;
                %W.est2.end=trjEnds;
            end
            clear Wviterbi WsMaxP
            
            W.est2.Ts=zeros(length(trjStarts),N);
            W.est2.Ps=zeros(length(trjStarts),N);
            for m=1:length(trjStarts)
                W.est2.Ts(m,:)=sum(pst(trjStarts(m):trjEnds(m),:),1); % time spent in each state
                W.est2.Ps(m,:)=W.est2.Ts(m,:)/sum(W.est2.Ts(m,:));
            end
            
            try
                W.est.lnAmean=logm(W.est.Amean);
                W.est.lnA_error='none';
            catch me
                W.est.lnAmean=0*W.est.Amean;
                W.est.lnA_error=me;
            end
        end
    end
end
%% exit message
if(displayExit) % display exit message
    displayHeader();
    displayConvergence();
    fprintf('%s \n',C.exitStatus)
end
%% slim down the model, on request, by deleting the bulky E-field
if(do_slim)
    W=rmfield(W,'E');
    W.est=rmfield(W.est,{'Ts','Ps'});
end
%% auxiliary functions
    function dFminus=displayConvergence()
        fprintf('%02d % 5d % 0.2e ',[W.N C.iter C.dFrel]);
        fprintf('%0.2e  %s ',C.dPrel,C.limitPar)
        %fprintf(', dFterms :')
        %fprintf('% 0.2e ',[W.Fterms.Fterms-Wm1.Fterms.Fterms]);
        fprintf(', F : %0.3e \n',W.F)
        dFminus=(C.dFrel<-1e-11);
        if(dFminus)
            %keyboard
        end
    end
    function displayHeader()
        disp([' N     it   dF/|F|   d<p>/p   p '])
    end
end

