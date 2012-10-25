function w=VB3_reorder(w0,ind)
% function w=VB3_reorder(w0,ind)
%
% Permute the state indices in model in the order specified in ind, and 
% return PM and M fields.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_reorder.m, permute state ordering in the vbSPT package
% =========================================================================
% 
% Copyright (C) 2012 Martin Lind??n and Fredrik Persson
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


w.dim=w0.dim;
w.N=w0.N;

fn=fieldnames(w0.M);
for m=1:length(fn)
    if(strcmp(fn{m},'wA'))
        w.M.wA=w0.M.wA(ind,ind);
        w.PM.wA=w0.PM.wA(ind,ind);        
    else
        w.M.(fn{m}) = w0.M.(fn{m})(ind);
        w.PM.(fn{m})=w0.PM.(fn{m})(ind);
    end
end

