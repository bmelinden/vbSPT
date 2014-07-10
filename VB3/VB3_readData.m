function X=VB3_readData(runinput)
% X=VB3_readData(runinput)
% 
% Read diffusion data set as specified in a runinputfile. It is also possible
% to use an options structure, e.g., from opt=VB3_getOptions(runinputfile) instead.
% Only trajectories longer than trjLmin (specified in opt) are returned. 
%

% ML 2014-07-10 : handle runinput files without the .m suffix

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_readData.m, reads position data fo use with the vbSPT package
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
if(ischar(runinput))% && exist(runinput, 'file')==2)
    runinputfile = runinput;
    opt=VB3_getOptions(runinputfile);
    disp(['Read runinput file ' runinputfile])
    % if an option struct, read in the runinputfilename
elseif(isstruct(runinput))
    opt=runinput;
    runinputfile=opt.runinputfile;
    disp(['Read options structure based on runinput file ' runinputfile ])
else
    error(['Not a valid input, aborting VB3_readData']);
end

%% start of actual code
foo=load(opt.inputfile,opt.trajectoryfield);

if(isfield(opt,'trjLmin'))
    Lmin=opt.trjLmin;
else
    warning('VB3_readData: cannot find minimum trajectory length in opt structure. Using minimum value: trjLmin=2');
    Lmin=2;
end
% extract the relevant columns
X=cell(1,length(foo.(opt.trajectoryfield)));
k=0;

for m=1:length(foo.(opt.trajectoryfield))
    if(size(foo.(opt.trajectoryfield){m}(:,1:opt.dim),1)>=Lmin)
        k=k+1;
        X{k}=foo.(opt.trajectoryfield){m}(:,1:opt.dim);
   end
end
% remove empty elements
X={X{1:k}};
