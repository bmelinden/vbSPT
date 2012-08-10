function VB3_license(funcName)
% res=VB3_license(funcName)
%?
% Prints a short license text to the command line.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_license.m, prints a short license text for the vbSPT package
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

fprintf(...
['\nvbSPT, %s, Copyright (C) 2012 Martin Lindén and Fredrik Persson \n' ...
 'This program comes with ABSOLUTELY NO WARRANTY. \n' ...
 'This is free software, and you are welcome to redistribute it \n' ...
 'under certain conditions. See license.txt for details. \n\n'],funcName);
end



