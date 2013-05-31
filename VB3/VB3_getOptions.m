function opt=VB3_getOptions(runinputfile)
% opt=VB3_getOptions(runinputfile)
%
% convert HMM runinput parameters from a runinput file into an options
% structure opt. In fact, all variables created by the command
% eval(runinputfile) are stored in the opt structure. 
% An (incomplete) sanity check of some parameter values is also performed.
%
% M.L. and F.P. 2012

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_getOptions.m, read runinout parameters in the vbSPT package
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


% Split the runinput filename
[path, name, ext] = fileparts(runinputfile);
if(isempty(path))
    path='.';
else
    warning('VB3_getOptions warning: runinput file not in the current folder,')
end
if(ismember('.',name))
   error(['runinputfile : ' name '. It is currently not possible to handle runinputfile names containing a .'])
end

% read the raw options from the runinput file
oldFolder = cd(path);
eval(name)
cd(oldFolder);
clear oldFolder; % forget what folder the options file happend to be called from
vv=whos;
opt=struct;

for m=1:length(vv)
    opt.(vv(m).name)=eval(vv(m).name);    
end
opt.localroot=pwd;

% temporary translations from old (pre 2012-06-11) parameter names
old={'Nmax','dt','pathestimate','BSnum','savefile','sourcefile'};
new={'maxHidden','timestep','stateEstimate','bootstrapNum','outputfile','inputfile'};
for kk=1:length(old)
    if(isfield(opt,old{kk}))
        opt.(new{kk})=opt.(old{kk});
        opt=rmfield(opt,old{kk});
        warning(['This runinput file uses an outdated parameter name: ' old{kk} '. Use ' new{kk} ' instead.'])
    end
end

% make some sanity checks
% 
if(opt.prior_tD<=opt.timestep)
    error('VB3: prior mean dwell time must be greater than time step')
end
if(length(opt.init_D)~=2)
    error('VB3: needs init_D to be an interval')
end
% none discovered so far
