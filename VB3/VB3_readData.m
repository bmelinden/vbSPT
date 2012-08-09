function X=VB3_readData(opt)
% X=VB3_readData(opt)
% 
% read diffusion data set as specified in the options structure opt, from
% e.g., opt=VB3_getOptions(runinputfile). Only trajectories longer than
% opt.trjLmin are returned. 
%
% M.L. 2012-04-17

foo=load(opt.inputfile,opt.trajectoryfield);

if(isfield(opt,'trjLmin'))
    Lmin=opt.trjLmin;
else
    warning('VB3_readData: cannot find minimum trajectory length in opt structure. Using minimum value: trjLmin=2');
    Lmin=2;
end
% extract the relevant columns
X=cell(1,length(foo.(opt.trajectoryfield)));
k=0;
for m=1:length(foo.(opt.trajectoryfield))
    if(size(foo.(opt.trajectoryfield){m}(:,1:opt.dim),1)>=Lmin)
        k=k+1;
        X{k}=foo.(opt.trajectoryfield){m}(:,1:opt.dim);
   end
end
% remove empty elements
X={X{1:k}};
