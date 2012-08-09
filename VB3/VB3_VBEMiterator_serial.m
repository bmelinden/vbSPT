function [W,C,F]=VB3_VBEMiterator_serial(W,X,varargin)
%% [W,C,F]=VB3_VBEMiterator(W,X,varargin)
%
% Perform VBEM iterations, on the VB structure W, with data (structure)
% X, until convergence. This version accepts d-dimensional data (T by d
% matrices), and uses a 'short' model, i.e., hidden states are associated
% with steps, not positions. 
%
% options:
% 'estimate'   : if given, add a field est2, which contains memory and
%                computer intensive estimates, such as the Viterbi path.
%
% 'maxIter',n  : run at most n VBE iterations.
% 'relTolF',tf : convergence criterion for relative change in likelihood bound.
% 'tolPar' ,tp : convergence criterion for M-step parameters.
%                iterate until
%                |dF/F| <= tf and max|dPar(i)/Par(i)| <= tp,
%                where Par(i) are all parameters for the variational
%                distributions. Default values are
%                maxIter=1000, relTolF=1e-6, tolPar=1e-3.
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
% this function uses mex files for the computer intensive inner loops:
% VB_forwardbackward.mexXXX VB_wAcount.mexXX, where XXX is a platform
% deoendent extension.
% To (re)compile binaries for your system, issue
% mex -setup
% mex('VB_forwardbackward.c',['-' computer('arch')],'-O')
% mex('VB_wAcount.c',['-' computer('arch')],'-O')


%% change log and notes
% M.L. 2012-02-26   : changed to d-dimensional data (second attempt), and
%                     then to the shorter model representation. Put this as
%                     new VB3 iterator, from VB2_VBEMiterator_ddim_short.
% M.L. 2012-01-27   : cleaning up notations in the code in order to prepare
%                     for jump to d-dimensional data. 
% M.L. 2011-12-15   : tested forward-backward against vbFRET
%                     implementation: same results. Corrected sign error in
%                     KL-divergence for g. F seems to increase always now.
%                     Added exit message.
% M.L. 2011-12-08   : first test: runs, but F does not always increase,
%                     especially for short trajectories. seem lnZ is the
%                     culprit?
% M.L. 2011-12-07   : started this diffusion version, with convergence
%                     check included. No testing done, some variable
%                     renaming remains, and the iteration part should
%                     perhaps be put in a separate subfunction to make the
%                     code easier to read.
% M.L. 2011-02-02   : started from VB5_VBEMiter.m
% M.L. 2011-02-21   : ~20% optimization by skipping some steps in the
%                     transitions count loop when the model has no spurious
%                     states.
% M.L. 2011-02-25   : Further optimizations by rewriting large parts of the
%                     VBE step in c.
% M.L. 2011-03-01   : successfully renamed E.C -> E.U, E.S -> E.C E.M ->
%                     E.V, E.P -> E.M, and made corresponding changes in
%                     the Ec field names.
% M.L. 2011-07-25   : removed the internal function minimalStorage from being saved to the VB7 structure, since
%                     it turns out to be very bulky to store. Created an
%                     external VB7 function instead.
% M.L. 2011-08-15   : reformulated VBE step slightly to avoid division by
%                     zero in case of points x=(0,0).

%% process options
C.maxIter=1000;
C.relTolF=1e-6;
C.tolPar=1e-3;
do_estimates=false;
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
        elseif(strcmpi(option,'maxIter'))
            C.maxIter=varargin{k+1};
            if(~isnumeric(C.maxIter) || C.maxIter~=round(C.maxIter) || C.maxIter<=0)
                error('VB1_VBEMiter: maxIter option must be followed by a positive integer.')
            end
            k=k+2;
        elseif(strcmpi(option,'relTolF'))
            C.relTolF=varargin{k+1};
            if(~isnumeric(C.relTolF) || C.relTolF<=0)
                error('VB1_VBEMiter: relTolF option must be followed by a positive number.')
            end
            k=k+2;
        elseif(strcmpi(option,'tolPar'))
            C.tolPar=varargin{k+1};
            if(~isnumeric(C.tolPar) || C.tolPar<=0)
                error('VB1_VBEMiter: tolPar option must be followed by a positive number.')
            end
            k=k+2;
        elseif(strcmpi(option,'outputLevel'))
            outputLevel=varargin{k+1};
            if(~isnumeric(outputLevel) || ~ismember(outputLevel,[0 1 2]))
                error('VB1_VBEMiter: outpuLevel option must be followed by 0, 1, or 2.')
            end
            if(outputLevel==0) % adjust output settings
                displayExit=false;
            elseif(outputLevel==2)
                displayProgress=true;
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
W.N=size(W.PM.wA,1);
%% initialize VBEM iterations
runMore=true;
C.exitStatus='';
Wm5=[];Wm4=[];Wm3=[];Wm2=[];Wm1=[];
%% iterate
while(runMore)
    % keep a short history in case something goes wrong...
    Wm5=Wm4;Wm4=Wm3;Wm3=Wm2;Wm2=Wm1;Wm1=W;
    %try % try one round of VBEM iterations
    %% M-step: contributions from all trajectories are added
    if(isfield(W,'E')) % then start with M step
        W.M.wPi = W.PM.wPi;
        W.M.wA =  W.PM.wA;
        W.M.n  = W.PM.n;
        W.M.c  = W.PM.c;
        for m=1:Ntrj
            % transitions and initial conditions
            W.M.wPi = W.M.wPi + W.E(m).wPi;
            W.M.wA =  W.M.wA  + W.E(m).wA;
            % emission model
            W.M.n  = W.M.n  + W.E(m).n;
            W.M.c  = W.M.c  + W.E(m).c;
        end
        % check for problems
        isNanInf=(sum(~isfinite([W.M.wPi W.M.wA(1:end) W.M.n  W.M.c]))>1);
        if(isNanInf)
            error('VB1_VBEMiter:Mfield_not_finite','Nan/Inf generated in VBM step')
        end
    end
    %% E-steps starts here
    % coupling matrix is the same for all trajectories
    lnQ=psi(W.M.wA)-psi(sum(W.M.wA,2))*ones(1,W.N);
    lnQmax=max(max(lnQ));
    Q=exp(lnQ-lnQmax);
    for m=1:Ntrj
        %% trial emission priobability q(s,t)
        % variational pointwise contributions, notation as in ML1 notes
        % except for index rul for which step 'belongs' to which hidden state
        % initial state t=1
        lnH1=psi(W.M.wPi)-psi(sum(W.M.wPi)); % initial state probability
        % t>1: we compute log(H) first, and then exponentiate
        lnH   =zeros(W.T(m),W.N);
        lnH0=W.dim/2*(psi(W.M.n)-log(pi*W.M.c)); % same for all t>1
        for j=1:W.N
            lnH(:,j)=lnH0(j)-W.M.n(j)/W.M.c(j)*dx2{m};
        end
        lnH(1,:)=lnH(1,:)+lnH1;
        lnHMax=max(lnH,[],2);
        H=zeros(W.T(m),W.N);
        for j=1:W.N % now we compute the actual pointwise emission contribution
            H(:,j)=exp(lnH(:,j)-lnHMax);
        end
        %% forward sweep  (depends only on Q and H)
        [Za,alpha,~,beta,pst]=VB_forwardbackward(Q,H);
        % compute quantities for next M-steps, same as VB5
       
        %% transition counts
        W.E(m).wPi=pst(1,:);          % <s_1>=<\delta_{j,s_1}>_{q(s)}
        W.E(m).wA =VB_wAcount(alpha,beta,H,Q);
        %% for the emission models
        W.E(m).n=W.dim/2*sum(pst,1);   % sum_{t=2}^T p(s_t)
        W.E(m).c=zeros(1,W.N);             % sum_{t=2}^T P(s(t)=j,c(t)=1).*dx(t)^2
        for j=1:W.N
            W.E(m).c(j)=sum(pst(:,j).*dx2{m});
        end
        % check for problems
        isNanInf=(sum(~isfinite([W.E(m).n W.E(m).c]))>1);
        if(isNanInf)
            error('VB2_VBEMiter:Efield_not_finite','Nan/Inf generated in VBE step')
        end
        % forward sweep normalization constant (same as VB4)
        W.Fterms.lnZQ(m)=(W.T(m)-1)*lnQmax;
        W.Fterms.lnZq(m)=sum(lnHMax);
        W.Fterms.lnZz(m)=sum(log(Za));
        % trajectory specific light-weight estimates (always)
        % occupation probabilities
        W.est.Ts(m,:)=sum(pst,1); % time spent in each state 
        W.est.Ps(m,:)=W.est.Ts(m,:)/sum(W.est.Ts(m,:));
        % transition rates and dwell times: s(t)
        W.est.Q=Q;
        W.est.lnQ=lnQ;                
        %% potentially demanding estimates (only if asked, and only on last iteration)
        if(do_estimates)
            W.est2.pst{m}=pst;
            W.est2.lnHMax{m}=lnHMax;
            W.est2.alpha{m}=alpha;
            W.est2.beta{m}=beta;
            W.est2.H{m}=H;
            W.est2.lnH{m}=lnH;
            [~,W.est2.sMaxP{m}]=max(pst,[],2);
            W.est2.sMaxP{m}=uint8(W.est2.sMaxP{m});
        end                
    end
    % global light-weight estimates (always)
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
    % KL divergence of initial state probability for s(t),c(t) (same as VB4)
    u0Pi=sum(W.PM.wPi);
    w0Pi=sum(W.M.wPi);
    KL_pi=gammaln(w0Pi)+gammaln(u0Pi)...
          +sum((gammaln(W.PM.wPi)-gammaln(W.M.wPi))...
               +(W.M.wPi-W.PM.wPi).*(psi(W.M.wPi)-psi(w0Pi)));
    % this expression assumes sum(W.M.wPi)=1+sum(W.PM.wPi), not true for
    % multiple data sets
    %KL_pi = log(u0Pi)-psi(u0Pi)-1/u0Pi...
    %    +sum((W.M.wPi-W.PM.wPi).*psi(W.M.wPi)...
    %    -gammaln(W.M.wPi)+gammaln(W.PM.wPi));
    W.Fterms.piTerms=-KL_pi;
    F=F-sum(KL_pi);
    % KL divergence of emission parameters for s(t)
    KL_gj= W.PM.n.*log(W.M.c./W.PM.c)...
        -W.M.n.*(1-W.PM.c./W.M.c)...
        -gammaln(W.M.n)+gammaln(W.PM.n)...
        +(W.M.n-W.PM.n).*psi(W.M.n);
    W.Fterms.gTerms=-KL_gj;
    F=F-sum(KL_gj);
    %% assembly of the free energy (same as VB5)
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
%% potentially demanding estimates (only if asked, and only on last iteration)
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

