function [W,C,F]=VB3_VBEMiterator(W,X,varargin)
%% [W,C,F]=VB3_VBEMiterator(W,X,varargin)
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
% W.M.wA;         : transition counts for s(t) transition matrix
% W.M.n; W.M.c;  : B-distribution for s(t)
% Prior parameters are stored in W.PM
% E-step parameters are under W.E
% W.F              : lower bound on the likelihood
% W.est            : various useful estimates that do not take up a lot of
%                    memory
% W.est2           : more memory intensensive estimates
%
% This function uses mex files from HMMcore for the computer intensive
% inner loops: VB_forwardbackward.mexXXX VB_wAcount.mexXX, where XXX is a
% platform dependent extension. Please refer to HMMcore/compile_code.m to
% (re)compile for your system.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_VBEMiterator, variational EM iterations in the vbSPT package
% =========================================================================
% 
% Copyright (C) 2012 Martin Lind??n and Fredrik Persson
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
%% preprocess the data
Ntrj=length(X);
dx2=cell(1,Ntrj);
W.T=[];
dim0=size(X{1},2);
for m=1:Ntrj
    dx2{m}=sum(diff(X{m},1,1).^2,2);
    [W.T(m),dim]=size(X{m});
    W.T(m)=W.T(m)-1;
    if(dim~=dim0)
        error('VB3_ddim: all trajectories must have the same dimensionality (=same number of columns)')
    end
    dim0=dim;
end
W.dim=dim;
W.N=size(W.PM.wB,1);
%% initialize VBEM iterations
runMore=true;
C.exitStatus='';
Wm3=[];Wm2=[];Wm1=[];
%% iterate
E(1:Ntrj)=struct;
while(runMore)
    % keep a short history in case something goes wrong...
    Wm3=Wm2;Wm2=Wm1;Wm1=W;
    %try % try one round of VBEM iterations
    %% M-step: contributions from all trajectories are added
    if(isfield(W,'E')) % then start with M step
        W.M.wPi = W.PM.wPi;
        W.M.wa =  W.PM.wa;
        W.M.wB =  W.PM.wB;
        W.M.n  = W.PM.n;
        W.M.c  = W.PM.c;
        for m=1:Ntrj
            % transitions and initial conditions
            W.M.wPi = W.M.wPi + W.E(m).wPi;
            wB=W.E(m).wA.*(1-eye(W.N));
            W.M.wa =  W.M.wa  + [sum(wB,2) diag(W.E(m).wA)];
            W.M.wB =  W.M.wB  + wB;
            % emission model
            W.M.n  = W.M.n  + W.E(m).n;
            W.M.c  = W.M.c  + W.E(m).c;
        end
        clear wB
        % check for problems
        isNanInf=(sum(~isfinite([W.M.wPi W.M.wa(1:end) W.M.wB(1:end) W.M.n  W.M.c]))>1);
        if(isNanInf)
            error('VB1_VBEMiter:Mfield_not_finite','Nan/Inf generated in VBM step')
        end
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

    % here starts a weird optimization scheme that makes it possible to run
    % the loop m=1:Ntrj in parallell. However, it turned out that most of
    % the speedup remained even in serial mode, perhaps it is expensive to
    % use nested struct vector constructs?
    lnH1=psi(W.M.wPi)-psi(sum(W.M.wPi)); % initial state probability
    clear lnH H lnH0
    T=W.T;N=W.N;dim=W.dim;
    WMn=W.M.n;
    WMc=W.M.c;
    WestTs=zeros(Ntrj,W.N);WestPs=zeros(Ntrj,W.N);
    lnZQ=zeros(1,Ntrj);lnZq=zeros(1,Ntrj);lnZz=zeros(1,Ntrj);
    West2lnHMax=cell(1,Ntrj);
    West2pst=cell(1,Ntrj); 
    West2H=cell(1,Ntrj); 
    West2lnH=cell(1,Ntrj);
    West2sMaxP=cell(1,Ntrj);
    for m=1:Ntrj
        %% trial emission priobability q(s,t)
        % variational pointwise contributions, notation as in ML1 notes
        % except perhaps an off-by-one difference in time index

        % we compute log(H) first, and then exponentiate
        lnH0=dim/2*(psi(WMn)-log(pi*WMc)); % same for all t
        lnH  =zeros(T(m),N);
        for j=1:N
            lnH(:,j)=lnH0(j)-WMn(j)/WMc(j)*dx2{m};
        end
        lnH(1,:)=lnH(1,:)+lnH1;
        lnHMax=max(lnH,[],2);
        H=zeros(T(m),N);
        for j=1:N % now we compute the actual pointwise emission contribution
            H(:,j)=exp(lnH(:,j)-lnHMax);
        end
        %% forward sweep  (depends only on Q and H)
        [Za,alpha,~,beta,pst]=VB_forwardbackward(Q,H);
        % compute quantities for next M-steps, same as VB5      
        %% transition counts
        E(m).wPi=pst(1,:);          % <s_1>=<\delta_{j,s_1}>_{q(s)}
        E(m).wA =VB_wAcount(alpha,beta,H,Q);
        %% for the emission models
        E(m).n=dim/2*sum(pst,1);   % sum_{t=2}^T p(s_t)
        E(m).c=zeros(1,N);             % sum_{t=2}^T P(s(t)=j,c(t)=1).*dx(t)^2
        for j=1:N
            E(m).c(j)=sum(pst(:,j).*dx2{m});
        end
        %% check for problems
        %isNanInf=(sum(~isfinite([E(m).n E(m).c]))>1);
        %if(isNanInf)
        %    error('VB2_VBEMiter:Efield_not_finite','Nan/Inf generated in VBE step')
        %end
        % forward sweep normalization constant (same as VB4)
        lnZQ(m)=(T(m)-1)*lnQmax;
        lnZq(m)=sum(lnHMax);
        lnZz(m)=sum(log(Za));
        % trajectory specific light-weight estimates (always)
        % occupation probabilities
        WestTs(m,:)=sum(pst,1); % time spent in each state 
        WestPs(m,:)=WestTs(m,:)/sum(WestTs(m,:));
        %% potentially demanding estimates (only if asked, and only on last iteration)
        if(do_estimates)
            West2pst{m}=pst;
            West2H{m}=H;
            West2lnH{m}=lnH;
            [~,West2sMaxP{m}]=max(pst,[],2);
            West2sMaxP{m}=uint8(West2sMaxP{m});
        end                
    end
    for m=1:Ntrj
        % check for problems
        isNanInf=(sum(~isfinite([E(m).n E(m).c]))>1);
        if(isNanInf)
            error('VB2_VBEMiter:Efield_not_finite','Nan/Inf generated in VBE step')
        end
    end
    W.E=E;
    W.Fterms.lnZQ=lnZQ;
    W.Fterms.lnZq=lnZq;
    W.Fterms.lnZz=lnZz;
    W.est.Ts=WestTs;
    W.est.Ps=WestPs;
    W.est.Ptot=sum(WestTs)/sum(sum(WestTs)); % total average
    if(do_estimates)
        W.est2.lnHMax=West2lnHMax;
        W.est2.pst   =West2pst;
        W.est2.H     =West2H;
        W.est2.lnH   =West2lnH;
        W.est2.sMaxP =West2sMaxP;
        clear West2lnHMax West2pst West2H West2lnH
    end
    % global light-weight estimates (always)
    error('ML: code conversion continues here!')
    W.est.Amean=zeros(size(W.M.wA));
    W.est.Astd=zeros(size(W.M.wA));
    W.est.Amode=zeros(size(W.M.wA));
    a0=sum(W.M.wA,2);
    for k=1:W.N
        W.est.Amean(k,:)=W.M.wA(k,:)/a0(k);
        W.est.Astd(k,:)=sqrt((W.M.wA(k,:).*(a0(k)-W.M.wA(k,:)))./(a0(k)^2*(1+a0(k))));
        W.est.Amode(k,:)=(W.M.wA(k,:)-1)/(a0(k)-W.N);
    end
    
    W.est.dwellMean=1./(1-diag(W.est.Amean));
    W.est.dwellMode=1./(1-diag(W.est.Amode));
    % emission parameters
    W.est.gMean=W.M.n./W.M.c;
    W.est.gMode=(W.M.n-1)./W.M.c;
    W.est.gStd=sqrt(W.M.n./W.M.c.^2); % sqrt(Var(g))
    W.est.DdtMean=W.M.c/4./(W.M.n-1);
    W.est.DdtMode=W.M.c/4./(W.M.n+1);
    W.est.Ddtstd=W.M.c/4./(W.M.n-1)./sqrt(W.M.n-2);
    %% assemble free energy
    F=sum(W.Fterms.lnZQ+W.Fterms.lnZq+W.Fterms.lnZz);
    % KL divergence of transition probabilities of s(t) (same as VB4)
    KL_A=zeros(1,W.N);
    for k=1:W.N
        wA0=sum(W.M.wA(k,:));
        uA0=sum(W.PM.wA(k,:));
        KL_A(k)=gammaln(wA0)-gammaln(uA0)-(wA0-uA0)*psi(wA0)...
            +sum((W.M.wA(k,:)-W.PM.wA(k,:)).*psi(W.M.wA(k,:))...
            -gammaln(W.M.wA(k,:))+gammaln(W.PM.wA(k,:)));
        %-sum(gammaln(W.M.wA(k,:))-gammaln(W.PM.wA(k,:))...
        %-(W.M.wA(k,:)-W.PM.wA(k,:)).*psi(W.M.wA(k,:)));
    end
    W.Fterms.Aterms=-KL_A;
    F=F-sum(KL_A);
    clear wA0 uA0;
    % KL divergence of initial state probability 
    u0Pi=sum(W.PM.wPi);
    w0Pi=sum(W.M.wPi);
    KL_pi=gammaln(w0Pi)+gammaln(u0Pi)...
          +sum((gammaln(W.PM.wPi)-gammaln(W.M.wPi))...
               +(W.M.wPi-W.PM.wPi).*(psi(W.M.wPi)-psi(w0Pi)));
    W.Fterms.piTerms=-KL_pi;
    F=F-sum(KL_pi);
    % KL divergence of emission parameters
    KL_gj= W.PM.n.*log(W.M.c./W.PM.c)...
        -W.M.n.*(1-W.PM.c./W.M.c)...
        -gammaln(W.M.n)+gammaln(W.PM.n)...
        +(W.M.n-W.PM.n).*psi(W.M.n);
    W.Fterms.gTerms=-KL_gj;
    F=F-sum(KL_gj);
    if(~isfinite(F))
        error('VB3_VBEM: F not finite (KL_gj)')
    end
    %% assembly of the free energy
    W.F=F;
    %W.Fterms.Fterms=[ W.Fterms.lnZQ W.Fterms.lnZq W.Fterms.lnZz -KL_A -KL_pi -KL_gj];
    W.Fterms.Fterms=[ W.Fterms.lnZQ+W.Fterms.lnZq+W.Fterms.lnZz -sum(KL_A) -sum(KL_pi) -sum(KL_gj)];
    W.Fterms.FtermsNames='[lnZQ+lnZq+lnZz -sum(KL_A) -sum(KL_pi) -sum(KL_gj)]';
    if(~isfinite(W.F))
        disp(W.Fterms)
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
end
%% exit message
if(displayExit) % display exit message
    displayHeader();
    displayConvergence();
    fprintf('%s \n',C.exitStatus)
end
%% some large estimates that can be done on the last iteration.
if(do_estimates)
    try
        W.est.lnAmean=logm(W.est.Amean);
        W.est.lnAmode=logm(W.est.Amode);
        W.est.lnA_error='none';
    catch me
        W.est.lnAmean=0*W.est.Amean;
        W.est.lnAmode=0*W.est.Amean;
        W.est.lnA_error=me;
    end
    for m=1:Ntrj
        W.est2.viterbi{m}=uint8(VBviterbi_log(W.est.lnQ,W.est2.lnH{m})); % Viterbi path
    end    
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

