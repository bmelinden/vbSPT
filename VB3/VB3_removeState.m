function M=VB3_removeState(w,s)
% M=VB3_removeState(w,s)
% create a new model parameter field M by removing state s from model w.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VB3_removeState.m, removes state from models in the vbSPT package
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


M=struct;

sk=[1:s-1 s+1:w.N]; % states to keep
M.wPi=w.M.wPi(sk);
M.n  =  w.M.n(sk);
M.c  =  w.M.c(sk);

% transfer observed transitions
wA =w.M.wA -w.PM.wA;
% try to compensate for transitions that went via the removed state
toS=wA(sk,s);
frS=wA(s,sk);
M.wA=wA(sk,sk)+toS*frS+w.PM.wA(sk,sk);
