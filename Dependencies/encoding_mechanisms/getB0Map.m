function [fieldmap] = getB0Map(iField,TE, mask)
%getB0Map gets off-resonance fieldmap from known echo times of acquisition
% fit for b0 map
if ~exist("mask","var");mask=[];end
if isempty(mask);    mask=ones(size(iField,1),size(iField,2),size(iField,3));end
ncoils=size(iField,5);
for lc=1:ncoils
    fieldmap_rad(:,:,:,lc) = Fit_ppm_complex_TE(iField(:,:,:,:,lc).*repmat(mask,1,1,1,size(iField,4)),TE);
end
deltaTE=TE(2)-TE(1);
fieldmap=fieldmap_rad./(2*pi*(deltaTE));
fieldmap=mean(fieldmap,4); % b0 map can be assumed to be same for all coils
fieldmap=fiex3d(fieldmap,2,4);
end