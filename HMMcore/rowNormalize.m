function A=rowNormalize(Q)
% A=rowNormalize(Q)
% normalize the rows of a transition matrix to enforce sum_j Q(i,j) = 1

% M.L. 2011

A=Q;
for k=1:size(Q,1)
    A(k,:)=Q(k,:)/sum(Q(k,:));
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
