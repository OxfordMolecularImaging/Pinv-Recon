function [rel_coil_sense] = getSENSmap(data,sum_pts,plt)
%GETSENSMAP get coil sensitvity map from image
% Inputs
%            data   complex data from which to extract map (x,y,z,t,coil)           
%         sum_pts   number of time points to take weighted sum of     ([1])
%             plt   plotting
% Output
%  rel_coil_sense   coil sensitivity map
% 
% 08/2024   Kylie Yeung

if ~exist('sum_pts','var');sum_pts=[];end
if isempty(sum_pts);sum_pts=1;end
if ~exist('plt','var');sum_pts=[];end
if isempty(plt);plt=0;end

ncoils=size(data,5);
max_sig_mag=squeeze(sum(abs(data),[1,2,3,5])); % signal over time
[~,max_lt]=maxk(max_sig_mag,sum_pts); % take the maximum k time points

% get sensitivity from each time point
for lt=1:length(max_lt)
    image_t=squeeze(mean(data(:,:,:,max_lt(lt),:),4)); % x,y,z,coil

    coil_rss=sqrt(sum((image_t.*conj(image_t)),4));
    rel_coil_sense_lt=conj(image_t)./repmat(coil_rss,1,1,1,ncoils);

    rel_coil_sense(:,:,:,:,lt)=rel_coil_sense_lt*max_sig_mag(max_lt(lt));
end

% take weighted sum
rel_coil_sense=sum(rel_coil_sense,5)./sum(max_sig_mag(max_lt));

% polynomial smooting
for lc=1:ncoils;rel_coil_sense(:,:,:,lc)=fiex3d(rel_coil_sense(:,:,:,lc),1,4);end

% normalize
rel_coil_sense=bsxfun(@rdivide,rel_coil_sense,sqrt(mean(rel_coil_sense.*conj(rel_coil_sense),4)));

% plotting
if plt
    figure;
    for lc=1:ncoils
        nexttile;mat2montage(flipud(abs(rel_coil_sense(:,:,:,lc))));colormap gray
        title(strcat('coil ',num2str(lc)))
    end
end

end