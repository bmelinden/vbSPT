function [wbs,Wmean,Wstd]=VB3_bootstrap(W,X,opt,NB,varargin)
% [wbs,Wmean,Wstd]=VB3_bootstrap(W,dat,opt,NB,wbs0)
%
% Run NB bootstrap fits starting from model W and data X. Each bootstrap
% fit resamples the data with replacement, and then converges the VB3 model
% using W as a starting point. No other model sizes are attempted, so if
% the data set is small, or if some states are only present in very few
% trajectories, strange things may happen...
%
% No path estimates are done, even if opt.pathestimates=true
% This function uses a parfor loop, but does not open or close the
% matlabpool by itself. 
%
% W     : Initial guess for bootstrap fit. For good results, use the best
%         model for the inriginal data set X. 
% X     : diffusion data 
% opt   : VB3 parameter object ( obj=VB3_getOptions(runinputfile) ).
% NB    : optional number of bootstrap iterations to run. Default: 100.
% wbs0  : wbs struct array, containing bootstrap indices in wBS.ind. This
%         is used to duplicate a bootstrap analysis on a different starting
%         model.
%
% wbs   : vector of bootstrapped models. To save some space, only three fields are stored. 
%         wbs(k).M, wbm(k).est (from the bootstrapped models), and
%         wbs(k).ind, the resampling indices which can be used to recreate
%         the resampled data, as Y={X{wbs(k).ind}};, and from there the
%         rest of the model.
% Wmean, Wstd: bootstrap average and standard deviation models, taken
%              elements-wise for some chosen fields. Only M and est fields
%              are analyzed.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_bootstrap.m, runs a bootstrap analysis in the vbSPT package
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


%% parameters
hasindices=false;
if(nargin>4)
    wbs0=varargin{1};
    NB=length(wbs0);
    hasindices=true;
else
    wbs0=struct('ind',cell(1,NB)); % trick to avoid complaints from parfor
end
if(~exist('NB','var') || isempty(NB)); NB=100;end
%% bootstrap fits
L=length(X);
tBootStrapFit=tic;
parfor k=1:NB
    %tic
    if(hasindices)
        ind=wbs0(k).ind;
    else
        ind=sort(ceil(L*rand(1,L))); % sample with replacement
    end
    %disp(['bootstrap iter ' int2str(k) ' 1'])
    % produce resampled data set
    Y=X(ind);
    dat=VB3_preprocess(Y,opt.dim);
    
    %disp(['bootstrap iter ' int2str(k) ' 2'])
    ww=VB3_VBEMiterator(W,dat,'outputLevel',0,'maxIter',opt.maxIter,...
        'relTolF',opt.relTolF,'tolPar',opt.tolPar,'slim');
    %disp(['bootstrap iter ' int2str(k) ' 3'])    
    wbs(k).M=ww.M;
    wbs(k).est=ww.est;
    wbs(k).F=ww.F;    
    wbs(k).ind=ind;
    %disp(['bootstrap iter ' int2str(k) ' finished'])
end
disp(['bootstrap : ' int2str(NB) ' resampling in ' num2str(toc(tBootStrapFit)) ' s.'])
%% default bootstrap analysis
Wmean=struct;
Wstd=struct;
f={'M' 'est'}; % fields on which to compute bootstrap statistics

for m=1:length(f)
   g=fieldnames(wbs(1).(f{m}));
   for j=1:length(g)
       [a,b]=size(wbs(1).(f{m}).(g{j}));
       bsnum=zeros(NB,a*b);
       
       % construct table
       for k=1:NB
            bsnum(k,:)=wbs(k).(f{m}).(g{j})(1:end);
       end
       % compute mean and std
       Wmean.(f{m}).(g{j})       =zeros(a,b);
       Wmean.(f{m}).(g{j})(1:end)=mean(bsnum,1);       
       Wstd.(f{m}).(g{j})       =zeros(a,b);
       Wstd.(f{m}).(g{j})(1:end)=std(bsnum,[],1);
   end
end
Wmean.F=mean([wbs(1:end).F]);
Wstd.F=std([wbs(1:end).F]);


