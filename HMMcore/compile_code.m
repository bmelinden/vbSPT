% Matrlab script to compile all c files in this folder.
% M.L. 2011

mex -setup

ff=dir('*.c');

for k=1:length(ff)
    disp(ff(k).name)
    mex(ff(k).name,['-' computer('arch')],'-O')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                 %
%       rowNormalize for VB HMM rate estimates    %
% =============================================== %
%                                                 %
% Copyright (C) 2012 Martin Lindén                %
%                                                 %
%        E-mail: bmelinden@gmail.com              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
