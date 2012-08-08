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
% M.L. 2012-06-18 : corrected unit error in initial guess for diffusion
%                   constant
% M.L. 2012-06-13 : decreased scope of full bottstrap to optimal size +-3,
%                   and introduced the 'slim' option to decrease computer
%                   time and storage footprint of the models. 
% M.L. 2012-06-11 : Save the best model for each size. Merged updated
%                   variable names. Changed viterbiEstimate ->
%                   stateEstimate
% F.P. 2012-06-07 : changed variable names Nmax -> maxHidden, dt ->
%                   timestep, pathestimate -> viterbiEstimate, BSnum ->
%                   bootstrapNum and savefile -> outputfile
% M.L. 2012-05-16 : changed to a deeper greedy serach that tries to remove
%                   all states before giving up. Also made a more ambitious
%                   bootstrap, that bootstraps all models sizes, and
%                   estimates the probability that a certain model is
%                   optimal. 
% F.P. 2012-04-16 : fixed a bug concerning the bootstrapping arguments
% M.L. 2012-04-04 : added the possibility to pass an options structure
%                   instead of a runinput filename.
% M.L. 2012-03-31 : added some parameter corrections, and forwarded
%                   convergence criteria to the VBEM iterations.


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
maxHidden=opt.maxHidden;
timestep=opt.timestep;

nTot=0; % total number of time steps in data
for n=1:length(X) % count total number of time steps
    nTot=nTot+length(X{n})-1;                        
end
Ntrj=length(X); % number of trajectories
        
Witer  =cell(1,opt.runs); % save all models generated in each run

% setup distributed computation toolbox
if(opt.parallelize_config)
    if(matlabpool('size')>0) % disable existing setting
        matlabpool close
    end
    eval(opt.parallel_start)
end
parfor iter=1:opt.runs
    % Greedy search strategy is probably more efficient than to start over
    % at each model size. We simply start with a large model, and
    % systematically remove the least occupied statate until things start
    % to get worse.  
    tic
    w=VB3_createPrior(opt,maxHidden);
    w0=w;
    %% initial guess
    w.M.wPi=w.PM.wPi+Ntrj*ones(1,maxHidden)/maxHidden; % uniform initial state distr.
    
    ntD=opt.init_tD/timestep;                               % init dwell time interval
    ntD=exp(log(ntD(1))+diff(log(ntD))*rand(1,maxHidden)); % log(ntD) has flat distribution
    A0=diag(1-1./ntD)*eye(maxHidden)+diag(1./ntD/(maxHidden-1))*(ones(maxHidden,maxHidden)-eye(maxHidden));
    w.M.wA=w.PM.wA+nTot*A0;
    %clear ntD A0;
        
    w.M.n=w.PM.n+nTot*ones(1,maxHidden);
    Ddt=opt.init_D*timestep;
    Ddt=exp(log(Ddt(1))+diff(log(Ddt))*rand(1,maxHidden)); % log(Ddt) has flat distribution
    w.M.c=w.M.n*4.*Ddt;
    %clear Ddt
    
    % converge largest model
    w=VB3_VBEMiterator(w,X,'outputLevel',0,'maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar,'slim');
    %% re-written more comprehensive greedy search
    
    Witer{iter}{1}=w;
    oneMoreTry=true;
    move=0;
    while(oneMoreTry) % try successive removal of low-occupancy states
        foundImprovement=false;
        move=move+1;
        w0=w; % reference state
        [~,h]=sort(w0.est.Ptot); % order in increasing occupancy
        for k=1:length(h)
            if(w0.N>1)
                %% try to remove one state
                w=VB3_createPrior(opt,w0.N-1);
                w.M=VB3_removeState(w0,h(k));
                try
                    w=VB3_VBEMiterator(w,X,'outputLevel',0,'maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar,'slim');
                catch me
                    disp('VB3_HMManalysis encountered an error:')
                    disp(me.message)
                    disp(me.stack)
                    errFile=['errorlog_VB3_HMManalysis' num2str(rand) '.mat'];
                    save(errFile);
                    rethrow(me);
                end
                Witer{iter}{end+1}=w; % save this attempt
                if(w.F>w0.F) % then this helped, and we should go on
                    w0=w;
                    foundImprovement=true;
                    %disp(['Iter ' int2str(iter) ': removing state ' int2str(h(k)) ' helped, new size N = ' int2str(w0.N)])
                    break % do not make further attempts at this model size
                end
                % if this did not help, and some significant occupation probability
                % was removed, try adding some extra transitions and reconverge.
                if(w0.est.Ptot(h(k))>0.1/nTot)
                    %disp(['Iter ' int2str(iter) ': simple removal did not help. Trying to add some extra transitions'])
                    w=VB3_createPrior(opt,w0.N-1);
                    w.M=VB3_removeState(w0,h(k));
                    w.M.wA=w.M.wA+2;
                    try
                        w=VB3_VBEMiterator(w,X,'outputLevel',0,'maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar,'slim');
                    catch me
                        disp('VB3_HMManalysis encountered an error:')
                        disp(me.message)
                        disp(me.stack)
                        errFile=['errorlog_VB3_HMManalysis' num2str(rand) '.mat'];
                        save(errFile);
                        rethrow(me);
                    end
                    Witer{iter}{end+1}=w; % save this attempt too
                end
                if(w.F>w0.F)
                    w0=w;
                    foundImprovement=true;
                    %disp(['Iter ' int2str(iter) ': removing state ' int2str(h(k)) ' + adding extra transitions helped, new size N = ' int2str(w0.N)'])
                    break % do not make further attempts at this model size
                end
            end
            % give up if W0 could not be improved upon
            if(~foundImprovement)
                oneMoreTry=false;
            end
        end
    end
 
    disp(['Finished greedy search, iter ' int2str(iter) ', with model size ' int2str(w0.N) ...
        ' : ' num2str(toc) ' s.'] )
end
%% collect best models for all sizes
INF=[];
Wbest.F=-inf;
WbestN=cell(1,maxHidden);
bestIter=0;
for k=1:maxHidden
    WbestN{k}.F=-inf;
end
for iter=1:opt.runs
    for k=1:length(Witer{iter})
        w=VB3_sortModel(Witer{iter}{k});
        INF(end+1,1:3)=[iter w.N w.F];
        if(w.F>Wbest.F)
            Wbest=w;
        end
        if(w.F>WbestN{w.N}.F)
            WbestN{w.N}=w;
            bestIter=iter;
        end
    end
end
% regenerate unsorted fields
disp('regenerating best model and estimates')
if(opt.stateEstimate)
    Wbest=VB3_VBEMiterator(Wbest,X,'estimate','maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar,'outputLevel',0);
else
    Wbest=VB3_VBEMiterator(Wbest,X,'maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar,'outputLevel',0);
end
parfor k=1:length(WbestN)
    if(WbestN{k}.F>-inf);
        WbestN{k}=VB3_VBEMiterator(WbestN{k},X,'maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar,'slim');
    end
end
% 
disp(['Best model size: ' int2str(Wbest.N) ', from attempt ' int2str(bestIter) '.'])
%% write results to outputfile
res=struct;
res.options=opt;
res.Wbest=Wbest;
res.WbestN=WbestN;
res.INF=INF;
for k=1:length(WbestN)
   res.dF(k)=WbestN{k}.F-Wbest.F; 
end
%% bootstrap?
if(opt.bootstrapNum>0)
    % bootstrap best model
    [wbs,Wmean,Wstd]=VB3_bootstrap(Wbest,X,opt,opt.bootstrapNum);
    res.bootstrap.wbs=wbs;     % this field will take up most of the storage.
    res.bootstrap.Wmean=Wmean;
    res.bootstrap.Wstd=Wstd;
    
    if(~isfield(opt,'fullBootstrap'))
        warning(['the flag fullBootstrap is missing in ' opt.runinputfile '. default=false']);
        opt.fullBootstrap=false;
    end
    if(opt.fullBootstrap)
    % bootstrap all other models with the same data selection
    clear WmeanN WstdN dFmean dFstd
    Fbs=zeros(opt.bootstrapNum,maxHidden);
    for k=1:maxHidden
        if(WbestN{k}.F>-inf && abs(WbestN{k}.N-Wbest.N)<=3) 
            [W,Wm,Ws]=VB3_bootstrap(WbestN{k},X,opt,opt.bootstrapNum,wbs); %VB3_bootstrap is parallelized
            WmeanN(k)=Wm;
            WstdN(k) =Ws;
            dF=[W.F]'-[wbs.F]';
            Fbs(:,k)=[W.F]';
            dFmean(k)=mean(dF);
            dFstd(k)=std(dF);
        else
            dFmean(k)=-inf;
            dFstd(k)=0;
            pBest(k)=0;
            Fbs(:,k)=-inf;
        end
    end
    % best model size in each bootstrap 
    nBest=zeros(1,maxHidden);
    for bs=1:opt.bootstrapNum
        [fb,nb]=max(Fbs(bs,:));
        nBest(nb)=nBest(nb)+1;
    end
    
    res.bootstrap.WmeanN=WmeanN;
    res.bootstrap.Fbootstrap=Fbs;
    res.bootstrap.WstdN=WstdN;
    res.bootstrap.dFmean=dFmean;
    res.bootstrap.dFstd=dFstd;
    res.bootstrap.pBest=nBest/sum(nBest);
    end
end
save(opt.outputfile,'-struct','res');
if(opt.parallelize_config)
    eval(opt.parallel_end)
end

end
