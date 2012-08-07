function opt=VB3_getOptions(runinputfile)
% opt=VB3_getOptions(runinputfile)
%
% read HMM parameters into an options structure opt.
% M.L. 2012-03-30
% F.P. 2012-06-05 Added functionaltity of running a runinput file not in
% the current path.
% M.L. 2012-06-07 Tried to tweak path recognition to deal with paths in current directory 
% ML 2012-06-11: added translations to handle transitions to new variable
%                names
% M.L. 2012-06-12 : updated sourcefile -> inputfile


% Split the filename
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
