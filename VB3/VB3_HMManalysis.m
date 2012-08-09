function res=VB3_HMManalysis(runinputfile)
% res=VB3_HMManalysis(runinputfile)
%
% Run the variational HMM analysis specified in the runinputfile, which
% should be in the current directory. 
% It is also possible to use an options structure, e.g., from
% opt=VB3_getOptions(runinputfile) instead, which might be useful for
% scripting things like parameter sweeps from a single runinput file.
%
% res : structure containing the result of the analysis. The same
% information is also written to the outputfile specified in the runinput
% file or options structure. 

%% change-log
% M.L. 2012-08-09 : Write command line output to a log file next to the
%                   output file (but with a .log extension instead), and
%                   made minor changes to the log output.
%                   logfile stored as opt.logfile
%                   Explicitely added .mat extension to outputfile   
% F.P. 2012-07-03 : fixed bug in by expanding the if (w.N>1) statement for
%                   adding transitions to include the reports if it helped 
%                   (timer tx0 could be undefined in the reports before)
% F.P. 2012-07-03 : added saving step prior to bootstrapping and put
%                   bootstrapping in an external function with no redundancy
% M.L. 2012-06-29 : added timing of extra transition count convergences
% M.L. 2012-06-28 : added another search step with extra transition counts
% F.P. 2012-06-27 : fixed bug when printing dF/abs(F)
% M.L. 2012-06-27 : added a search step that fills in missing transitions and
%                   reconverges after no more states can be removed.
% M.L. 2012-06-18 : created a version with simpler search: only attempts to
%                   remove the least occupied state
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

%% start of actual code
tstart=tic;
%% read analysis parameters
% if an existing file, generate options structure
if(ischar(runinputfile) && exist(runinputfile)==2)
    opt=VB3_getOptions(runinputfile);
    disp(['Read runinput file ' runinputfile])
elseif(isstruct(runinputfile))
    opt=runinputfile;
    runinputfile=opt.runinputfile;
    disp(['Read options structure based on runinput file ' runinputfile ])
end

% add .mat extension to output file if not present
[outpath,outfile]=fileparts(opt.outputfile);
opt.outputfile=fullfile(outpath,[outfile '.mat']);
clear outfile outpath;

% construct log file name from outputfile, with extension .log
[logpath,logfile]=fileparts(opt.outputfile);
opt.logfile=fullfile(logpath,[logfile '.log']);
clear logfile logpath;

% start the diary, and clear old entries (the output is overwritten too). 
if(exist(opt.logfile,'file'))
    delete(opt.logfile)
end
diary(opt.logfile);
diary on
disp('Copyright notice should be printed here!')

disp('----------')
disp([ datestr(now) ' : Starting greedy optimization to find best model.'])
disp(['jobID        : ' opt.jobID])
disp(['runinput file: ' opt.runinputfile])
disp(['output file  : ' opt.outputfile])
disp(['log file     : ' opt.logfile])
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
    od=[]; tx0=[];
    % Greedy search strategy is probably more efficient than to start over
    % at each model size. We simply start with a large model, and
    % systematically remove the least occupied statate until things start
    % to get worse.
    titer=tic;
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
    %% greedy search
    
    Witer{iter}{1}=w;
    oneMoreTry=true;
    move=0;
    while(oneMoreTry) % try successive removal of low-occupancy states
        foundImprovement=false;
        move=move+1;
        w0=w; % reference state
        [~,h]=sort(w0.est.Ptot); % order in increasing occupancy
        for k=1:1 % length(h) %%% k=1 only is where one decides to only try the least occupied state
            if(w0.N>1)
                %% try to remove a looping state
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
                % if this did not help, try adding more transition counts and reconverge.
                
                if(w.N>1)
                    tx0=tic;
                    %disp(['Iter ' int2str(iter) ': simple removal did not help. Trying to add some extra transitions'])
                    w=VB3_createPrior(opt,w0.N-1);
                    w.M=VB3_removeState(w0,h(k));
                    od=max(max((w.M.wA-w.PM.wA).*(1-eye(w.N)))); % largest off-diagonal element
                    w.M.wA= diag(diag(w.M.wA))+od*(1-eye(w.N));
                    if(isfield(w,'E')) % make sure that the new M field is used
                        w=rmfield(w,'E');
                    end
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
                    
                    if(w.F>w0.F)
                        disp(['Iter ' int2str(iter) '. Removing state ' int2str(h(k)) ' (of ' int2str(w0.N) ...
                            ' ) + adding  ' int2str(od) ' extra transitions helped, dF/|F| = ' ...
                            num2str((w.F-w0.F)/abs(w.F)) ', time : ' num2str(toc(tx0)) ' s.']);
                        w0=w;
                        foundImprovement=true;
                        break % do not make further attempts at this model size
                    else
                        disp(['Iter ' int2str(iter) '. Removing state ' int2str(h(k)) ' (of ' int2str(w0.N) ...
                            ' ) + adding ' int2str(od) ' extra transitions did not help, dF/|F| = ' ...
                            num2str((w.F-w0.F)/abs(w.F)) ', time : ' num2str(toc(tx0)) ' s.']);
                    end
                end
            end
            % give up if W0 could not be improved upon
            if(~foundImprovement)
                oneMoreTry=false;
            end
            %% last attempt: add transition pseudocounts to w0 and reconverge.
            if(w0.N>1)
                tx0=tic;
                w=w0;
                od=max(max((w.M.wA-w.PM.wA).*(1-eye(w.N)))); % largest off-diagonal element
                w.M.wA= diag(diag(w.M.wA))+od*(1-eye(w.N));
                if(isfield(w,'E')) % make sure that the new M field is used
                    w=rmfield(w,'E');
                end
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
                if(w.F>w0.F)
                    disp(['Iter ' int2str(iter) '. Adding ' int2str(od) ...
                        ' extra transitions to final model helped : dF/|F| = ' ...
                        num2str((w.F-w0.F)/abs(w.F)) ', time : ' num2str(toc(tx0)) ' s.']);
                else
                    disp(['Iter ' int2str(iter) '. Adding ' int2str(od) ...
                        ' extra transitions to final model did not help : dF/|F| = ' ...
                        num2str((w.F-w0.F)/abs(w.F)) ', time : ' num2str(toc(tx0)) ' s.']);
                end
            end
        end
    end
    
    disp(['Iter ' int2str(iter) '. Finished greedy search in '  num2str(toc(titer)) ' s.'] )
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
disp('Sorting and regenerating best model and estimates')
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
disp(['Best model size: ' int2str(Wbest.N) ', from iteration ' int2str(bestIter) '.'])
%% write results to outputfile
res=struct;
res.options=opt;
res.Wbest=Wbest;
res.WbestN=WbestN;
res.INF=INF;
for k=1:length(WbestN)
    res.dF(k)=WbestN{k}.F-Wbest.F;
end

% saving the models prior to bootstrapping them
disp(['Saving ' opt.outputfile ' after ' num2str(toc(tstart)/60) ' min.']);
save(opt.outputfile,'-struct','res');


%% bootstrapping

if(opt.bootstrapNum>0)
bootstrap = VB3_bsResult(opt, 'HMM_analysis');
res.bootstrap=bootstrap;

% save again after bootstrapping
save(opt.outputfile,'-struct','res');
end

% End parallel computing
if(opt.parallelize_config)
    eval(opt.parallel_end)
end

disp([datestr(now) ' : Finished ' opt.runinputfile '. Total run time ' num2str(toc(tstart)/60) ' min.'])
diary off
end


