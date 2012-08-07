function w=VB3_reorder(w0,ind)
% permute the state indices in mode w0, and return PM and fields.

w.dim=w0.dim;
w.N=w0.N;

fn=fieldnames(w0.M);
for m=1:length(fn)
    if(strcmp(fn{m},'wA'))
        w.M.wA=w0.M.wA(ind,ind);
        w.PM.wA=w0.PM.wA(ind,ind);        
    else
        w.M.(fn{m}) = w0.M.(fn{m})(ind);
        w.PM.(fn{m})=w0.PM.(fn{m})(ind);
    end
end
