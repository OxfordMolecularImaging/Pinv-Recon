function [F_slice, pinv_time] = getF_slice(E_slice,thresh,mode,useGPU)
%GETF Calculates reconstruction matrix slicewise from 
%      E_slice encoding matrix (cell array)
% Kylie Yeung 06/2025

nslices=length(E_slice);
F_slice=cell(1,nslices);

for ls=1:nslices
    clear progressBar
    progressBar(ls,nslices,'Reconstruction matrix for slice')
    [F_slice{ls},pinv_time]=getF(E_slice{ls},thresh,mode,useGPU,0);
end

F_slice=cat(3,F_slice{:});

end