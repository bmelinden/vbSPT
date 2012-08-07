function M=VB3_removeState(w,s)
% M=VB3_removeState(w,s)
% create a new model parameter field M by removing state s from model w.

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
