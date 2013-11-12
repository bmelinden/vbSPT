% A little script to illustrate the meaning of the Dirichlet prior strength
% for the Bij matrix, assuming equal weight to all components




dim=[ 3 6 9 12 100 100]
%dim=[ 3:10]
str=logspace(-1,1,3);
col='krbgmkrbgmkrbgm';
mar='.+x*.+x*.+x*.+x*';


figure(1)
clf
hold on
for d=1:length(dim) % dimension of Dirichlet variable
    for s=1:length(str)      % total strength
        b=dirrnd(str(s)*ones(1,dim(d))/dim(d)); % a Dirichlet vector with d components of total strength s
        b=-sort(-b); % sort in descending order
        
        plot(1:dim(d),b,'color',col(s),'marker','none','linew',2)
       
       
       
   end
end

legend(num2str(str'))
set(gca,'yscale','log','ylim',[1e-5 1],'xscale','log','xlim',[1 100])
box on
xlabel('rank')
ylabel('p_j')