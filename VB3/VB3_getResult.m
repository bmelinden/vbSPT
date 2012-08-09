function res=VB3_getResult(runinputfile)
% res=VB3_getResult(runinputfile)
%
% Find saved outputfile from a runinput file, load the analysis
% results, and print a short description to the command line.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_getResult.m, loads analysis results for the vbSPT package
% =========================================================================
% 
% Copyright (C) 2012 Martin Lindén and Fredrik Persson
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
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%% start of actual code

if(isstruct(runinput))
    opt=runinput;
else    
    opt=VB3_getOptions(runinput);
end

res=load(opt.outputfile);

%% Present results in the Matlab prompt

disp(['The best global model for ' opt.runinputfile ':']);
disp(sprintf('\n'));
disp(['Number of states: ' num2str(res.Wbest.N)]);
disp(sprintf('\n'));
disp(['Diffusion rate constants: ']);
disp(num2str(res.Wbest.est.DdtMean/opt.timestep/1e6, 3));
disp(sprintf('\n'));
disp(['Occupancy: ']);
disp(num2str(res.Wbest.est.Ptot, 3));
disp(sprintf('\n'));
disp(['Transition matrix [per timestep]: ']);
disp(num2str(res.Wbest.est.Amean, 3));


end



