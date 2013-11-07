function W1=VB3_sortModel(W,ind)
% W1=VB3_sortModel(W,ind)
%
% sort the states of the VB3 model W in order ind, where in is a
% permutation of 1:W.N. If ind is not given, the model is sorter in order
% of increasing most likely diffusion constant, 
% W.est.DdtMode=W.M.c/4./(W.M.n+1)
%
% Only the M,PM, and est fields are sorted, other fields must be recreated
% by running further VB iterations.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_sortModel.m, reorders state indices of models in the vbSPT package.
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


% M.L. 2012-04-14

%% check parameters
if(~exist('ind','var')|| isempty(ind))
    [~,ind]=sort(W.M.c/4./(W.M.n+1));
end
ind0=sort(union(ind(1),ind));
if(length(ind0)~=W.N || sum(ind0==1:W.N)~=W.N)
   error('VB3_sortModel error: incorrect order.')
end
clear ind0;
%% actual code


% fields that do not need reordering
W1.dim= W.dim;
W1.N  = W.N;
W1.T  = W.T;
W1.F  =W.F;
% reorder some fields
f={'PM','M','est'};
%g={'wPi','n','c'};
for a=1:length(f)    
    F=W.(f{a});
    g=fieldnames(F);
    for b=1:length(g)
        [R,C]=size(F.(g{b}));
        
        % only sort dimensions that have length W.N
        ri=1:R;
        ci=1:C;
        if(R==W.N)
            ri=ind;            
        end        
        if(C==W.N)
            ci=ind;
        end
        if(strcmp(g{b},'wa')) % special case: only rows of wa should be sorted
            ri=ind;
            ci=1:2;
        end
        W1.(f{a}).(g{b})=F.(g{b})(ri,ci);        
    end
end

