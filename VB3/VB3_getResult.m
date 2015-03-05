function res=VB3_getResult(runinput)
% res=VB3_getResult(runinput)
%
% Find saved outputfile from a runinput file or options struct, load the analysis
% results, and print a short description to the command line.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_getResult.m, loads analysis results for the vbSPT package
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

%% Parse input

% if an existing file, generate options structure
if(isstr(runinput) && exist(runinput)==2)
    runinputfile = runinput;
    opt=VB3_getOptions(runinputfile);
    disp(['Read runinput file ' runinputfile])
    % if an option struct, read in the runinputfilename
elseif(isstruct(runinput))
    opt=runinput;
    runinputfile=opt.runinputfile;
    disp(['Read options structure based on runinput file ' runinputfile ])
else
    error(['Not a valid input, aborting']);
end

%% Start of actual code
try
    res=load(opt.outputfile);
catch me
    if(~exist(opt.outputfile)) % then try to build the full path
        res=load(fullfile(fileparts(fullfile(opt.localroot,opt.runinputfile)),opt.outputfile));
    else
        me.rethrow;
    end
end

%% Present results in the Matlab prompt

disp(['The best global model for ' opt.runinputfile ':']);
disp(sprintf('\n'));
disp(['Number of states: ' num2str(res.Wbest.N)]);
disp(sprintf('\n'));
disp(['Diffusion rate constants: ']);
disp(num2str(res.Wbest.est.DdtMean/opt.timestep, 3));
disp(sprintf('\n'));
disp(['Occupancy: ']);
disp(num2str(res.Wbest.est.Ptot, 3));
disp(sprintf('\n'));
disp(['Transition matrix [per timestep]: ']);
disp(num2str(res.Wbest.est.Amean, 3));


end



