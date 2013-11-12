function w=VB3_convertOld(w0)
% w1=VB3_convertOld(w0) 
% Convert a model in the 1.0 format (Nat. Meth. 2013) to the new
% parameterization, choosing prior parameters to make the models
% equivalent.

% ML 2013-05-02 : got numerically identical results with old and new
% models on some test data after conversion with this function

w=struct;

w.dim=w0.dim;
w.N=w0.N;
w.M.SA=1:w0.N;

w.PM.n=w0.PM.n;
w.PM.c=w0.PM.c;
w.PM.wPi=w0.PM.wPi;
w.PM.wB=w0.PM.wA.*(1-eye(w.N));% state changes
w.PM.wa(:,1)=sum(w.PM.wB,2);   % state changes
w.PM.wa(:,2)=diag(w0.PM.wA);   % state non-changes

w.M.n=w0.M.n;
w.M.c=w0.M.c;
w.M.wPi=w0.M.wPi;
w.M.wB=w0.M.wA.*(1-eye(w.N)); % state changes
w.M.wa(:,1)=sum(w.M.wB,2);    % state changes
w.M.wa(:,2)=diag(w0.M.wA);    % state non-changes

