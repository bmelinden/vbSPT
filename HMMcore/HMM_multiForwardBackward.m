% [lnZ,wA,pst]=HMM_multiForwardBackard_m(Q,H,iEnd) 
%
% performs forward-backward sweeps on HMM-type time series stacked on top
% of each other. 
%
% Q    : transition matrix (not necessarily normalized) transition matrix
% H    : emission likelihood, including initial state probabilities on
%        appropriate rows 
% iEnd : indices to ends of individual time series, e.g., individual
%        trajectories run from  1 -> iEnd(1), iEnd(1)+2 -> iEnd(2), etc.
%        This is so that the algorithm know at what lines the transition
%        matrix should be omitted from the calculations.
%        
% lnZ  : log of normalization constant for the forward sweep
% wA   : transition count matrix (requires backward sweep)
% pst  : pst(t,j)= P(s(t)==j), averaged over all state sequences (some
%        extra computational cost in addition to backward sweep).
%
% The algorithm performs only those calculations needed for the number of
% outputs. E.g., lnZ=HMM_multiForwardBackard_m(Q,H,iEnd) is about twice as
% fast as a call with all three output variables.
%
% M.L. 2013-05-21
