function dat=SPT_preprocess(runinput,dim)
% dat=SPT_preprocess(runinput,dim)
%
% assemble single particle diffusion data in a form for fast EM iterations
%
% runinput  : a VB3/VBSPT runinput file, or a VB3/VBSPT runinput structure,
%             or a cell vector of diffusion trajectories
% optional parameters with cell vector input (ignored otherwise):
% dim       : optional specification of data dimension when input is cell
%             vector; only the first dim columns of each trajectory will be
%             analysed. Default: use all columns.

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
%% parse input
% if an existing file, generate options structure
if(ischar(runinput) && exist(runinput, 'file')==2)
    runinputfile = runinput;
    opt=VB3_getOptions(runinputfile);
    disp(['Read runinput file ' runinputfile])
    X=VB3_readData(opt);
    dim=opt.dim;
elseif(isstruct(runinput)) % if already an options structure
    opt=runinput;
    runinputfile=opt.runinputfile;
    X=VB3_readData(opt);
    dim=opt.dim;
elseif(iscell(runinput))
    X=runinput;
    if(~exist('dim','var'))
       dim=size(X{1},2);
    end
else
    error(['Not a valid input, aborting SPT_preprocess']);
end
clear opt runinput 

%% assemble output structure
dat=struct;
dat.dim=dim;

% count output size
T=zeros(size(X));
for k=1:length(X)
    T(k)=size(X{k},1);
end
if(~isempty(find(T<2,1)))
   warning('SPT_preprocess: data contains traces with no steps.')
end
dat.T=T;
dat.dx2=zeros(sum(T-1),1);
trjstarts=zeros(1,length(X),'double');
dat.end  =zeros(1,length(X),'double');

ind=1;
for k=1:length(X)
    dx2= sum(diff(X{k}(:,1:dim),1,1).^2,2);
    trjstarts(k)=ind;    
    dat.end(k)  =ind+length(dx2)-1;
    ind=ind+length(dx2);
    dat.dx2(trjstarts(k):dat.end(k))=dx2;
end




