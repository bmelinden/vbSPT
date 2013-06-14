function bootstrap=VB3_bsResult(runinput, varargin)
% res=VB3_bsResult(runinputfile,option1,option2,...)
%
% Run bootstrapping on the HMM analysis result specified in the runinputfile 
% which should be in the current directory.
% It is also possible to use an options structure, e.g., from
% opt=VB3_getOptions(runinputfile) instead. Note that bootstrapping 
% parameters have to be set in the options struct or runinputfile before 
% running this function.
% If called from VB3_HMManalysis as part of the initial
% analysis 'HMM_analysis' should exist as an option
%
% options:
% 'save'        : if given, the results will be saved in a new .mat file with 
%                 '_bootstrapped' added on the old outputfilenameontains memory and
% 'overwrite'   : if given, the 'save' option will overwrite the old outputfile 
%                 defined in the runinputfile/option struct.
% 'HMM_analysis' : if given, ignore the local treatment of parallel computing.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_HMManalysis, runs data analysis in the vbSPT package
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

tbootstrap=tic;
%% read bootstrapping input

% if an existing file, generate options structure
if(ischar(runinput) && exist(runinput, 'file')==2)
    runinputfile = runinput;
    opt=VB3_getOptions(runinputfile);
    disp(['Read runinput file ' runinputfile])
    % if an option struct, read in the runinputfilename
elseif(isstruct(runinput))
    opt=runinput;
    runinputfile=opt.runinputfile;
    disp(['Read options structure based on runinput file ' runinputfile ])
else
    error(['Not a valid input, aborting bootstrapping']);
end
%% load data
X=VB3_readData(opt);
res=load(opt.outputfile);

Wbest = res.Wbest;
WbestN = res.WbestN;

%% read options
do_save = false;
do_overwrite = false;
do_HMM_analysis = false;

if(nargin>1)        % then parse options
    k=1;            % argument counter
    kmax=nargin-1;  % stopping criterion
    while(k<=kmax)
        option=varargin{k};
        if(strcmpi(option,'save'))
            do_save=true;
            k=k+1;
        elseif(strcmpi(option,'overwrite'))
            do_overwrite=true;            
            k=k+1;
        elseif(strcmpi(option,'HMM_analysis'))
            do_HMM_analysis=true;            
            k=k+1;
        else
            error(['VB1_VBEMiter: option ' option ' not recognized.'])
        end
    end
end



%% initiate distributed computation
% setup distributed computation toolbox
if(opt.parallelize_config && ~do_HMM_analysis)
    if(matlabpool('size')>0) % disable existing setting
        matlabpool close
    end
    eval(opt.parallel_start)
end



%% bootstrap
disp(['jobID: ' opt.jobID])
disp('Starting bootstrapping...')
disp('----------')

if(opt.bootstrapNum>0)
    
    if(~isfield(opt,'fullBootstrap'))
        warning(['the flag fullBootstrap is missing in ' opt.runinputfile '. default=false']);
        opt.fullBootstrap=false;
    end
    if(opt.fullBootstrap)
        % bootstrap the best global model
        [wbs,Wmean,Wstd]=VB3_bootstrap(Wbest,X,opt,opt.bootstrapNum);
        res.bootstrap.wbs=wbs;     % this field will take up most of the storage.
        res.bootstrap.Wmean=Wmean;
        res.bootstrap.Wstd=Wstd;
        
        % bootstrap all other models with the same data selection
        clear WmeanN WstdN dFmean dFstd
        Fbs=zeros(opt.bootstrapNum,opt.maxHidden);
        for k=1:opt.maxHidden
            disp(['Bootstrapping the ' num2str(k) ' state model of ' num2str(opt.maxHidden) '.'])
            if k == Wbest.N
                % only fill in the parameters
                WmeanN(k)=Wmean;
                WstdN(k)=Wstd;
                Fbs(:,k)=[wbs.F]';
                dFmean(k)=0;
                dFstd(k)=0;
            elseif(WbestN{k}.F>-inf && abs(WbestN{k}.N-Wbest.N)<=2)
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
        nBest=zeros(1,opt.maxHidden);
        for bs=1:opt.bootstrapNum
            [fb,nb]=max(Fbs(bs,:));
            nBest(nb)=nBest(nb)+1;
        end
        
        
        % put the parameters in the res struct
        res.bootstrap.WmeanN=WmeanN;
        res.bootstrap.Fbootstrap=Fbs;
        res.bootstrap.WstdN=WstdN;
        res.bootstrap.dFmean=dFmean;
        res.bootstrap.dFstd=dFstd;
        res.bootstrap.pBest=nBest/sum(nBest);
        
    else
        % only bootstrap best global model
        [wbs,Wmean,Wstd]=VB3_bootstrap(Wbest,X,opt,opt.bootstrapNum);
        res.bootstrap.wbs=wbs;     % this field will take up most of the storage.
        res.bootstrap.Wmean=Wmean;
        res.bootstrap.Wstd=Wstd;
    end
    
    %% Save result
    res.options = opt;
    if do_save && ~do_overwrite
        save([opt.outputfile '_Bootstrapped'],'-struct','res');
    elseif do_overwrite
        save(opt.outputfile,'-struct','res');
    end

else
    warning(['No bootstrapping requested, check ' runinput]);
    return;
end


%% end distributed computing
if(opt.parallelize_config && ~do_HMM_analysis)
    eval(opt.parallel_end)
end

bootstrap = res.bootstrap;

disp(['Finished ' opt.runinputfile ' bootstrapping in ' num2str(toc(tbootstrap)/60) ' min.'])
end
