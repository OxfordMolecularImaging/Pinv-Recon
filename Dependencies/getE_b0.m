function [E_b0] = getE_b0(E0,b0,t)
%GETE Calculate encoding matrix from 
%            E encoding matrix
%           b0 b0 off-resonance
% Kylie Yeung 06/2025

[~,~,nslices]=size(b0);

E_b0=cell(1,nslices);
for ls=1:nslices    
    progressBar(ls,nslices,'Calculating encoding matrix with b0 encoding for slice')
    b0_s=b0(:,:,ls);
    D=exp(-1i*2*pi*t.'*(b0_s(:)).'); %% check if the - sign should be here
    E_b0{ls}=E0.*D;
end

end