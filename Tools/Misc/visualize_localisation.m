% a script to plot histogram of accumulated transverse positions for the
% different states. The idea is to try to distinguish cytoplasmic states
% from membrane bound ones: there are more cytoplasm per projected xy-area
% in the middle of the cell, but more membrane area per projected area near
% the edges. 
%


if(1) % read data
    opt=VB3_getOptions('runinput.m');
    opt2=opt;opt2.dim=2;
    X=VB3_readData(opt2);
    hmm=load(opt.savefile);
    W=hmm.Wbest;
end

%% total positions with state weights
xw=[];%zeros(sum(W.T),2+W.N);

for t=1:length(W.T)
    xw=[xw;X{t}(1:end-1,:) W.est2.pst{t} double(W.est2.sMaxP{t})];
end

%% y-histograms for different states
ind=cell(1,W.N);
for s=1:W.N
    ind{s}=find(xw(:,2+W.N+1)==s);    
end

figure(1)
clf

for s=1:W.N
   subplot(W.N,1,s)
   hist(xw(ind{s},2))   
end