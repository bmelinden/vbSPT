% A=rowNormalize(Q)
% normalize the rows of a transition matrix to enforce sum_j Q(i,j) = 1
% If each row is the parameters of a Dirichlet distribution, then each row
% in the output is the mode (most likely vlues) of that distribution.
% M.L. 2011

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rowNormalize.m, part of HMMcore/
% =========================================================================
% 
% Copyright (C) 2013 Martin Lindén, E-mail: bmelinden@gmail.com
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
% Additional permission under GNU GPL version 3 section 7
% 
% If you modify this Program, or any covered work, by linking or combining it
% with Matlab or any Matlab toolbox, the licensors of this Program grant you 
% additional permission to convey the resulting work.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%% start of actual code
function A=rowNormalize(Q)

A=Q;
for k=1:size(Q,1)
    A(k,:)=Q(k,:)/sum(Q(k,:));
end
