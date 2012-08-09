function opt=VB3_getOptions(runinputfile)
% opt=VB3_getOptions(runinputfile)
%
% convert HMM runinput parameters from a runinput file into an options
% structure opt. In fact, all variables created by the command
% eval(runinputfile) are stored in the opt structure. 
% An (incomplete) sanity check of some parameter values is also performed.
%
% M.L. and F.P. 2012

% Split the runinput filename
[path, name, ext] = fileparts(runinputfile);
if(isempty(path))
    path='.';
else
    warning('VB3_getOptions warning: runinput file not in the current folder,')
end
    
% read the raw options from the runinput file
oldFolder = cd(path);
eval(name)
cd(oldFolder);
clear oldFolder; % forget what folder the options file happend to be called from
vv=whos;
opt=struct;

for m=1:length(vv)
    opt.(vv(m).name)=eval(vv(m).name);    
end
opt.localroot=pwd;

% temporary translations from old (pre 2012-06-11) parameter names
old={'Nmax','dt','pathestimate','BSnum','savefile','sourcefile'};
new={'maxHidden','timestep','stateEstimate','bootstrapNum','outputfile','inputfile'};
for kk=1:length(old)
    if(isfield(opt,old{kk}))
        opt.(new{kk})=opt.(old{kk});
        opt=rmfield(opt,old{kk});
        warning(['This runinput file uses an outdated parameter name: ' old{kk} '. Use ' new{kk} ' instead.'])
    end
end

% make some sanity checks
% 
if(opt.prior_tD<=opt.timestep)
    error('VB3: prior mean dwell time must be greater than time step')
end
if(length(opt.init_D)~=2)
    error('VB3: needs init_D to be an interval')
end
% none discovered so far
