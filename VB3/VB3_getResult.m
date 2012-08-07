function res=VB3_getResult(runinput)
% res=VB3_getResult(runinputfile)
%
% Find saved outputfile from a runinput file, and load the analysis
% results.

% M.L. 2012-06-12
% 
% Change log:
% F.P. 2012-07-04 : Prints the important results to the Matlab prompt
% M.L. 2012-06-14 : handle opt structure as input
%

if(isstruct(runinput))
    opt=runinput;
else    
    opt=VB3_getOptions(runinput);
end

res=load(opt.outputfile);

%% Present results in the Matlab prompt

disp(['The best global model for ' opt.runinputfile ':']);
disp(sprintf('\n'));
disp(['Number of states: ' num2str(res.Wbest.N)]);
disp(sprintf('\n'));
disp(['Diffusion rate constants: ']);
disp(num2str(res.Wbest.est.DdtMean/opt.timestep/1e6, 3));
disp(sprintf('\n'));
disp(['Occupancy: ']);
disp(num2str(res.Wbest.est.Ptot, 3));
disp(sprintf('\n'));
disp(['Transition matrix [per timestep]: ']);
disp(num2str(res.Wbest.est.Amean, 3));


end



