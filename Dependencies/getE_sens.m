function [E_sens] = getE_sens(E,sens)
%GETE Calculate encoding matrix from 
%            E encoding matrix
%         sens coil sensitivity map
% Kylie Yeung 06/2025

[~,~,nslices,ncoils]=size(sens);

if ~iscell(E);E=repmat({E},1,nslices);end

E_sens=cell(1,nslices);
for ls=1:nslices    
progressBar(ls,nslices,'Calculating encoding matrix with coil sensitivity encoding for slice')
E_sens_slice=cell(1,ncoils);
    for lc=1:ncoils
        coil_sense_c=sens(:,:,ls,lc);
        E_sens_slice{lc}=E{ls}.*coil_sense_c(:)';
    end
E_sens{ls}=cat(1,E_sens_slice{:});
end

end