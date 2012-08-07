function res=VB3_HMManalysis(runinputfile)
% res=VB3_HMManalysis(runinputfile)
%
% Run the variational HMM analysis specified in the runinputfile, which
% should be in the current directory.  
% it is also possible to use an options structure, e.g., from 
% opt=VB3_getOptions(runinputfile) instead, which might be useful for
% parameter sweeps.
%
% M.L. 2012-03-30

% change-log
% F.P. 2012-04-16 : fixed a bug concerning the bootstrapping arguments
% M.L. 2012-04-04 : added the possibility to pass an options structure
%                   instead of a runinput filename.
% M.L. 2012-03-31 : added some parameter corrections, and forwarded
%                   convergence criteria to the VBEM iterations.

% to do:
% - 

%% read analysis parameters
% if an existing file, generate options structure
if(isstr(runinputfile) && exist(runinputfile)==2)
    opt=VB3_getOptions(runinputfile);
    disp(['Read runinput file ' runinputfile])
elseif(isstruct(runinputfile))
    opt=runinputfile;
    runinputfile=opt.runinputfile;
    disp(['Read options structure based on runinput file ' runinputfile ])
end
disp(['jobID: ' opt.jobID])
disp('Starting greedy optimization to find best model...')
disp('----------')

%% start of actual analysis code
%% set up model structure: prior distributions
%% converge models
% load data
X=VB3_readData(opt);

nTot=0; % total number of time steps in data
for n=1:length(X) % count total number of time steps
    nTot=nTot+length(X{n})-1;                        
end
Ntrj=length(X{1}); % number of trajectories
        
if(0) %% start with size 1 not necessary for greedy search
w=VB3_createPrior(opt,1);
w.M.wPi=Ntrj; % initial state distribution: ~flat
w.M.wA=nTot;  % 'transition' matrix
w.M.n=prior_n+nTot/2;            % for 1D data
w.M.c=w.M.n*4*sqrt(prod(opt.D)); % a value somewhere in the middle of the 
                                 % diffusion constant range.         
WbestN=cell(1);
INF=cell(1);

WbestN{1}=VB3_VBEMiterator(w,X);
for k=2:Nmax
    WbestN{k}.F=-inf;
end
INF(1,:)=[0 1 WbestN{1}.F]; % iteration modelsize lowerbound        
clear w
end

Witer  =cell(1,opt.runs);
INFiter=cell(1,opt.runs);
% setup distributed computation toolbox
if(opt.parallelize_config)
    if(matlabpool('size')>0) % disable existing setting
        matlabpool close
    end
    eval(opt.parallel_start)
end
parfor iter=1:opt.runs
    % Greedy search strategy is probably more efficient than to start over
    % at each model size. We simply start with a large model, and systematically remove the least occupied statate until things start to get worse.
    tic
    
    Nmax=opt.Nmax;
    dt=opt.dt;
    
    w=VB3_createPrior(opt,Nmax);
    %% initial guess
    w.M.wPi=w.PM.wPi+Ntrj*ones(1,Nmax)/Nmax; % uniform initial state distr.
    
    ntD=opt.init_tD/dt;                               % init dwell time interval
    ntD=exp(log(ntD(1))+diff(log(ntD))*rand(1,Nmax)); % log(ntD) has flat distribution
    A0=diag(1-1./ntD)*eye(Nmax)+diag(1./ntD/(Nmax-1))*(ones(Nmax,Nmax)-eye(Nmax));
    w.M.wA=w.PM.wA+nTot*A0;
    %clear ntD A0;
        
    w.M.n=w.PM.n+nTot*ones(1,Nmax);
    Ddt=opt.init_D;
    Ddt=exp(log(Ddt(1))+diff(log(Ddt))*rand(1,Nmax)); % log(Ddt) has flat distribution
    w.M.c=w.M.n*4.*Ddt;
    %clear Ddt

    % converge largest model
    w=VB3_VBEMiterator(w,X,'outputLevel',0,'maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar);
    %% greedy search
    INFiter{iter}=[iter Nmax w.F];  % keep track of all results
    Witer{iter}=w;                  % best model so far
    for k=Nmax-1:-1:1
        w0=w;
        [~,sMin]=min(w0.est.Ptot);
        w=VB3_createPrior(opt,k);        
        w.M=VB3_removeState(w0,sMin);
        try
            w=VB3_VBEMiterator(w,X,'outputLevel',0,'maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar);
        catch me
            disp('VB3_HMManalysis encountered an error:')
            disp(me.message)
            disp(me.stack)
            errFile=['errorlog_VB3_HMManalysis' num2str(rand) '.mat'];
            save(errFile);
            rethrow(me);
        end
        INFiter{iter}(end+1,:)=[iter w.N w.F];
        if(w.F>Witer{iter}.F) % then this helped, and we should go on
            Witer{iter}=w;
            continue
        end
        % if this did not help, and some significant occupation probability
        % was removed, try adding some extra transitions and reconverge.
        if(w0.est.Ptot(sMin)>0.1/nTot)
            w.M.wA=w.M.wA+2;
            try
                w=VB3_VBEMiterator(w,X,'outputLevel',0,'maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar);
            catch me
                disp('VB3_HMManalysis encountered an error:')
                disp(me.message)
                disp(me.stack)
                errFile=['errorlog_VB3_HMManalysis' num2str(rand) '.mat'];
                save(errFile);
                rethrow(me);
            end

            INFiter{iter}(end+1,:)=[iter w.N w.F];
            if(w.F>Witer{iter}.F) % then this helped, and we should go on
                Witer{iter}=w;
                continue
            end
        end
        % if both these attempts failed, then it is time to give in and
        % hope for luck
        break
    end
    disp(['Finished greedy search, iter ' int2str(iter) ', with model size ' int2str(Witer{iter}.N) ...
        ' : ' num2str(toc) ' s.'] )
end
%% compare best results from the different iterations
Wbest.F=-Inf;
bestIter=0;
INF=[];
for iter=1:opt.runs
    INF=[INF;INFiter{iter}];
    if(Witer{iter}.F>Wbest.F)
        Wbest=Witer{iter};
        bestIter=iter;
    end
end  
Wbest=VB3_sortModel(Wbest);     % sort in order of increasing diffusion constants
% regenerate unsorted fields
disp('regenerating best model and estimates')
if(opt.pathestimate)
    Wbest=VB3_VBEMiterator(Wbest,X,'estimate','maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar);
else
    Wbest=VB3_VBEMiterator(Wbest,X,'maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar);
end
% 
disp(['Best model size: ' int2str(Wbest.N) ', from attempt ' int2str(bestIter) '.'])
%% write results to savefile
res=struct;
res.options=opt;
res.Wbest=Wbest;
res.INF=INF;
%% bootstrap?
if(opt.BSnum>0)
    [wbs,Wmean,Wstd]=VB3_bootstrap(Wbest,X,opt,opt.BSnum);
    res.bootstrap.help='run "help VB3_bootstrap" for an explanation of these fields';
    res.bootstrap.wbs=wbs;     % this field will take up most of the storage.
    res.bootstrap.Wmean=Wmean;
    res.bootstrap.Wstd=Wstd;
end
save(opt.savefile,'-struct','res');
if(opt.parallelize_config)
    eval(opt.parallel_end)
end

end
